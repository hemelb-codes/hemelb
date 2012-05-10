#include <cmath>
#include <list>
#include <map>
#include <algorithm>
#include <zlib.h>

#include "debug/Debugger.h"
#include "io/formats/geometry.h"
#include "io/writers/xdr/XdrMemReader.h"
#include "geometry/decomposition/BasicDecomposition.h"
#include "geometry/GeometryReader.h"
#include "lb/lattices/D3Q27.h"
#include "net/net.h"
#include "topology/NetworkTopology.h"
#include "log/Logger.h"
#include "util/utilityFunctions.h"
#include "constants.h"

namespace hemelb
{
  namespace geometry
  {

    GeometryReader::GeometryReader(const bool reserveSteeringCore,
                                   const lb::lattices::LatticeInfo& latticeInfo,
                                   reporting::Timers &atimings) :
        latticeInfo(latticeInfo), timings(atimings)
    {
      // Get the group of all procs.
      MPI_Group worldGroup;
      MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

      // This rank should participate in the domain decomposition if
      //  - there's no steering core (then all ranks are involved)
      //  - we're not on core 0 (the only core that might ever not participate)
      //  - there's only one processor (so core 0 has to participate)
      participateInTopology = !reserveSteeringCore || topology::NetworkTopology::Instance()->GetLocalRank() != 0
          || topology::NetworkTopology::Instance()->GetProcessorCount() == 1;

      // Create our own group, without the root node if we're not running with it.
      if (reserveSteeringCore && topology::NetworkTopology::Instance()->GetProcessorCount() > 1)
      {
        int lExclusions[1] = { 0 };
        MPI_Group_excl(worldGroup, 1, lExclusions, &topologyGroup);
      }
      else
      {
        topologyGroup = worldGroup;
      }

      // Create a communicator just for the domain decomposition.
      MPI_Comm_create(MPI_COMM_WORLD, topologyGroup, &topologyCommunicator);

      // Each rank needs to know its rank wrt the domain
      // decomposition.
      if (participateInTopology)
      {
        topologyComms = topology::Communicator(topologyCommunicator);
      }
    }

    GeometryReader::~GeometryReader()
    {
      MPI_Group_free(&topologyGroup);

      // Note that on rank 0, this is the same as MPI_COMM_WORLD.
      if (participateInTopology)
      {
        MPI_Comm_free(&topologyCommunicator);
      }
    }

    Geometry GeometryReader::LoadAndDecompose(const std::string& dataFilePath)
    {
      log::Logger::Log<log::Info, log::OnePerCore>("Starting file read timer");
      timings[hemelb::reporting::Timers::fileRead].Start();

      // Open the file using the MPI parallel I/O interface at the path
      // given, in read-only mode.
      MPI_Info_create(&fileInfo);

      // Create hints about how we'll read the file. See Chapter 13, page 400 of the MPI 2.2 spec.
      std::string accessStyle = "access_style";
      std::string accessStyleValue = "sequential";
      std::string buffering = "collective_buffering";
      std::string bufferingValue = "true";

      MPI_Info_set(fileInfo, const_cast<char*>(accessStyle.c_str()), const_cast<char*>(accessStyleValue.c_str()));
      MPI_Info_set(fileInfo, const_cast<char*>(buffering.c_str()), const_cast<char*>(bufferingValue.c_str()));

      // Open the file.
      // Stupid C-MPI lack of const-correctness
      int error = MPI_File_open(MPI_COMM_WORLD,
                                const_cast<char *>(dataFilePath.c_str()),
                                MPI_MODE_RDONLY,
                                fileInfo,
                                &file);

      currentComms = topology::Communicator(MPI_COMM_WORLD);

      if (error != MPI_SUCCESS)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Unable to open file %s, exiting", dataFilePath.c_str());
        fflush(0x0);
        exit(0x0);
      }
      else
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Opened config file %s", dataFilePath.c_str());
      }
      fflush(NULL);

      // Set the view to the file.
      std::string mode = "native";
      MPI_File_set_view(file, 0, MPI_CHAR, MPI_CHAR, const_cast<char*>(mode.c_str()), fileInfo);

      log::Logger::Log<log::Info, log::OnePerCore>("Reading file preamble");
      Geometry geometry = ReadPreamble();

      // Read the file header.
      log::Logger::Log<log::Info, log::OnePerCore>("Reading file header");

      principalProcForEachBlock = std::vector<proc_t>(geometry.GetBlockCount());

      ReadHeader(geometry.GetBlockCount());
      timings[hemelb::reporting::Timers::initialDecomposition].Start();
      // Perform an initial decomposition, of which processor should read each block.
      log::Logger::Log<log::Info, log::OnePerCore>("Beginning initial decomposition");

      if (!participateInTopology)
      {
        for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
        {
          principalProcForEachBlock[block] = -1;
        }
      }
      else
      {
        // Get an initial base-level decomposition of the domain macro-blocks over processors.
        // This will later be improved upon by ParMetis.
        decomposition::BasicDecomposition basicDecomposer(geometry, latticeInfo, topologyComms, fluidSitesOnEachBlock);
        basicDecomposer.Decompose(principalProcForEachBlock);

        if (ShouldValidate())
        {
          basicDecomposer.Validate(principalProcForEachBlock);
        }
      }
      timings[hemelb::reporting::Timers::initialDecomposition].Stop();
      // Perform the initial read-in.
      log::Logger::Log<log::Info, log::OnePerCore>("Reading in my blocks");

      // Close the file - only the ranks participating in the topology need to read it again.
      MPI_File_close(&file);

      if (participateInTopology)
      {
        // Reopen in the file just between the nodes in the topology decomposition. Read in blocks
        // local to this node.
        MPI_File_open(topologyCommunicator, const_cast<char*>(dataFilePath.c_str()), MPI_MODE_RDONLY, fileInfo, &file);

        currentComms = topologyComms;

        ReadInBlocksWithHalo(geometry, principalProcForEachBlock, topologyComms.GetRank());

        if (ShouldValidate())
        {
          ValidateGeometry(geometry);
        }
      }

      timings[hemelb::reporting::Timers::fileRead].Stop();

      hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Begin optimising the domain decomposition.");
      timings[hemelb::reporting::Timers::domainDecomposition].Start();

      // Having done an initial decomposition of the geometry, and read in the data, we optimise the
      // domain decomposition.
      if (participateInTopology)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Beginning domain decomposition optimisation");
        OptimiseDomainDecomposition(geometry, principalProcForEachBlock);
        log::Logger::Log<log::Info, log::OnePerCore>("Ending domain decomposition optimisation");

        if (ShouldValidate())
        {
          ValidateGeometry(geometry);
        }
        MPI_File_close(&file);
      }

      // Finish up - close the file, set the timings, deallocate memory.
      MPI_Info_free(&fileInfo);

      timings[hemelb::reporting::Timers::domainDecomposition].Stop();

      return geometry;
    }

    /**
     * Read in the section at the beginning of the config file.
     */
    Geometry GeometryReader::ReadPreamble()
    {
      const unsigned preambleBytes = io::formats::geometry::PreambleLength;

      // Read in the file preamble into a buffer. We read this on a single core then broadcast it.
      // This has proven to be more efficient than reading in on every core (even using a collective
      // read).
      char preambleBuffer[preambleBytes];

      if (currentComms.GetRank() == HEADER_READING_RANK)
      {
        MPI_File_read(file, preambleBuffer, preambleBytes, MpiDataType(preambleBuffer[0]), MPI_STATUS_IGNORE);
      }

      MPI_Bcast(preambleBuffer,
                preambleBytes,
                MpiDataType<char>(),
                HEADER_READING_RANK,
                currentComms.GetCommunicator());

      // Create an Xdr translator based on the read-in data.
      hemelb::io::writers::xdr::XdrReader preambleReader = hemelb::io::writers::xdr::XdrMemReader(preambleBuffer,
                                                                                                  preambleBytes);

      unsigned hlbMagicNumber, gmyMagicNumber, version;
      // Read in housekeeping values
      preambleReader.readUnsignedInt(hlbMagicNumber);
      preambleReader.readUnsignedInt(gmyMagicNumber);
      preambleReader.readUnsignedInt(version);

      // Check the value of the HemeLB magic number.
      if (hlbMagicNumber != io::formats::HemeLbMagicNumber)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("This file starts with %d, not the HemeLB magic number %d.",
                                                     hlbMagicNumber,
                                                     io::formats::HemeLbMagicNumber);
        exit(1);
      }

      // Check the value of the geometry file magic number.
      if (gmyMagicNumber != io::formats::geometry::MagicNumber)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("This file is not a geometry file: had %d, not the geometry magic number %d.",
                                                     gmyMagicNumber,
                                                     io::formats::geometry::MagicNumber);
        exit(1);
      }

      log::Logger::Log<log::Warning, log::OnePerCore>("Geometry file version: %d", version);

      // Variables we'll read.
      // We use temporary vars here, as they must be the same size as the type in the file
      // regardless of the internal type used.
      unsigned int blocksX, blocksY, blocksZ, blockSize;
      double voxelSize;
      util::Vector3D<double> origin;

      // Read in the values.
      preambleReader.readUnsignedInt(blocksX);
      preambleReader.readUnsignedInt(blocksY);
      preambleReader.readUnsignedInt(blocksZ);
      preambleReader.readUnsignedInt(blockSize);
      preambleReader.readDouble(voxelSize);
      for (unsigned int i = 0; i < 3; ++i)
      {
        preambleReader.readDouble(origin[i]);
      }

      // Read the padding unsigned int.
      unsigned paddingValue;
      preambleReader.readUnsignedInt(paddingValue);

      return Geometry(util::Vector3D<site_t>(blocksX, blocksY, blocksZ), blockSize, voxelSize, origin);
    }

    /**
     * Read the header section, with minimal information about each block.
     *
     * Results are placed in the member arrays fluidSitesPerBlock,
     * bytesPerCompressedBlock and bytesPerUncompressedBlock.
     */
    void GeometryReader::ReadHeader(site_t blockCount)
    {
      site_t headerByteCount = GetHeaderLength(blockCount);
      // Allocate a buffer to read into, then do the reading.
      char* headerBuffer = new char[headerByteCount];

      if (currentComms.GetRank() == HEADER_READING_RANK)
      {
        MPI_File_read(file, headerBuffer, (int) headerByteCount, MpiDataType(headerBuffer[0]), MPI_STATUS_IGNORE);
      }

      MPI_Bcast(headerBuffer,
                (int) headerByteCount,
                MpiDataType<char>(),
                HEADER_READING_RANK,
                currentComms.GetCommunicator());

      // Create a Xdr translation object to translate from binary
      hemelb::io::writers::xdr::XdrReader preambleReader =
          hemelb::io::writers::xdr::XdrMemReader(headerBuffer, (unsigned int) headerByteCount);

      // Read in all the data.
      for (site_t block = 0; block < blockCount; block++)
      {
        unsigned int sites, bytes, uncompressedBytes;
        preambleReader.readUnsignedInt(sites);
        preambleReader.readUnsignedInt(bytes);
        preambleReader.readUnsignedInt(uncompressedBytes);

        fluidSitesOnEachBlock.push_back(sites);
        bytesPerCompressedBlock.push_back(bytes);
        bytesPerUncompressedBlock.push_back(uncompressedBytes);
      }

      delete[] headerBuffer;
    }

    /**
     * Read in the necessary blocks from the file.
     */
    void GeometryReader::ReadInBlocksWithHalo(Geometry& geometry,
                                              const std::vector<proc_t>& unitForEachBlock,
                                              const proc_t localRank)
    {
      // Create a list of which blocks to read in.
      timings[hemelb::reporting::Timers::readBlocksPrelim].Start();

      // Populate the list of blocks to read (including a halo one block wide around all
      // local blocks).
      log::Logger::Log<log::Info, log::OnePerCore>("Determining blocks to read");
      std::vector<bool> readBlock = DecideWhichBlocksToReadIncludingHalo(geometry, unitForEachBlock, localRank);

      if (ShouldValidate())
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating block sizes");

        // Validate the uncompressed length of the block on disk fits out expectations.
        for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
        {
          if (bytesPerUncompressedBlock[block]
              > io::formats::geometry::GetMaxBlockRecordLength(geometry.GetBlockSize(), fluidSitesOnEachBlock[block]))
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Block %i is %i bytes when the longest possible block should be %i bytes",
                                                          block,
                                                          bytesPerUncompressedBlock[block],
                                                          io::formats::geometry::GetMaxBlockRecordLength(geometry.GetBlockSize(),
                                                                                                         fluidSitesOnEachBlock[block]));
          }
        }
      }

      // Next we spread round the lists of which blocks each core needs access to.
      log::Logger::Log<log::Info, log::OnePerCore>("Informing reading cores of block needs");
      net::Net net = net::Net(currentComms.GetCommunicator());
      Decomposition decomposition(geometry.GetBlockCount(),
                                  readBlock,
                                  util::NumericalFunctions::min(READING_GROUP_SIZE, currentComms.GetSize()),
                                  net,
                                  currentComms.GetCommunicator(),
                                  currentComms.GetRank(),
                                  currentComms.GetSize());

      timings[hemelb::reporting::Timers::readBlocksPrelim].Stop();
      log::Logger::Log<log::Info, log::OnePerCore>("Reading blocks");
      timings[hemelb::reporting::Timers::readBlocksAll].Start();

      // Set the initial offset to the first block, which will be updated as we progress
      // through the blocks.
      MPI_Offset offset = io::formats::geometry::PreambleLength + GetHeaderLength(geometry.GetBlockCount());

      // Iterate over each block.
      for (site_t nextBlockToRead = 0; nextBlockToRead < geometry.GetBlockCount(); ++nextBlockToRead)
      {
        // Read in the block on all cores (nothing will be done if this core doesn't need the block).
        ReadInBlock(offset,
                    geometry,
                    decomposition.ProcessorsNeedingBlock(nextBlockToRead),
                    nextBlockToRead,
                    readBlock[nextBlockToRead]);

        // Update the offset to be ready for the next block.
        offset += bytesPerCompressedBlock[nextBlockToRead];
      }

      timings[hemelb::reporting::Timers::readBlocksAll].Stop();
    }

    void GeometryReader::ReadInBlock(MPI_Offset offsetSoFar,
                                     Geometry& geometry,
                                     const std::vector<proc_t>& procsWantingThisBlock,
                                     const site_t blockNumber,
                                     const bool neededOnThisRank)
    {
      // Easy case if there are no sites on the block.
      if (fluidSitesOnEachBlock[blockNumber] <= 0)
      {
        return;
      }
      std::vector<char> compressedBlockData;
      proc_t readingCore = GetReadingCoreForBlock(blockNumber);

      net::Net net = net::Net(currentComms.GetCommunicator());

      if (readingCore == currentComms.GetRank())
      {
        timings[hemelb::reporting::Timers::readBlock].Start();
        // Read the data.
        compressedBlockData.resize(bytesPerCompressedBlock[blockNumber]);
        MPI_File_read_at(file,
                         offsetSoFar,
                         &compressedBlockData.front(),
                         bytesPerCompressedBlock[blockNumber],
                         MPI_CHAR,
                         MPI_STATUS_IGNORE);

        // Spread it.
        for (std::vector<proc_t>::const_iterator receiver = procsWantingThisBlock.begin();
            receiver != procsWantingThisBlock.end(); receiver++)
        {
          if (*receiver != currentComms.GetRank())
          {
            net.RequestSend(&compressedBlockData.front(), bytesPerCompressedBlock[blockNumber], *receiver);
          }
        }
        timings[hemelb::reporting::Timers::readBlock].Stop();
      }
      else if (neededOnThisRank)
      {
        compressedBlockData.resize(bytesPerCompressedBlock[blockNumber]);
        net.RequestReceive(&compressedBlockData.front(), bytesPerCompressedBlock[blockNumber], readingCore);
      }
      else
      {
        return;
      }
      timings[hemelb::reporting::Timers::readNet].Start();
      net.Send();
      net.Receive();
      net.Wait();
      timings[hemelb::reporting::Timers::readNet].Stop();
      timings[hemelb::reporting::Timers::readParse].Start();
      if (neededOnThisRank)
      {
        // Create an Xdr interpreter.
        std::vector<char> blockData = DecompressBlockData(compressedBlockData, bytesPerUncompressedBlock[blockNumber]);
        io::writers::xdr::XdrMemReader lReader(&blockData.front(), blockData.size());

        ParseBlock(geometry, blockNumber, lReader);

        // If debug-level logging, check that we've read in as many sites as anticipated.
        if (ShouldValidate())
        {
          // Count the sites read,
          site_t numSitesRead = 0;
          for (site_t site = 0; site < geometry.GetSitesPerBlock(); ++site)
          {
            if (geometry.Blocks[blockNumber].Sites[site].targetProcessor != BIG_NUMBER2)
            {
              ++numSitesRead;
            }
          }
          // Compare with the sites we expected to read.
          if (numSitesRead != fluidSitesOnEachBlock[blockNumber])
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Was expecting %i fluid sites on block %i but actually read %i",
                                                          fluidSitesOnEachBlock[blockNumber],
                                                          blockNumber,
                                                          numSitesRead);
          }
        }
      }
      else if (!geometry.Blocks[blockNumber].Sites.empty())
      {
        geometry.Blocks[blockNumber].Sites = std::vector<GeometrySite>(0, GeometrySite(false));
      }
      timings[hemelb::reporting::Timers::readParse].Stop();
    }

    std::vector<char> GeometryReader::DecompressBlockData(const std::vector<char>& compressed,
                                                          const unsigned int uncompressedBytes)
    {
      timings[hemelb::reporting::Timers::unzip].Start();
      // For zlib return codes.
      int ret;

      // Set up the buffer for decompressed data. We know how long the the data is
      std::vector<char> uncompressed(uncompressedBytes);

      // Set up the inflator
      z_stream stream;
      stream.zalloc = Z_NULL;
      stream.zfree = Z_NULL;
      stream.opaque = Z_NULL;
      stream.avail_in = compressed.size();
      stream.next_in = reinterpret_cast<unsigned char*>(const_cast<char*>(&compressed.front()));

      ret = inflateInit(&stream);
      if (ret != Z_OK)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Decompression error for block");
        std::exit(1);
      }
      stream.avail_out = uncompressed.size();
      stream.next_out = reinterpret_cast<unsigned char*>(&uncompressed.front());

      ret = inflate(&stream, Z_FINISH);
      if (ret != Z_STREAM_END)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Decompression error for block");
        std::exit(1);
      }

      uncompressed.resize(uncompressed.size() - stream.avail_out);
      ret = inflateEnd(&stream);
      if (ret != Z_OK)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Decompression error for block");
        std::exit(1);
      }
      timings[hemelb::reporting::Timers::unzip].Stop();
      return uncompressed;
    }

    void GeometryReader::ParseBlock(Geometry& geometry, const site_t block, io::writers::xdr::XdrReader& reader)
    {
      // We start by clearing the sites on the block. We read the blocks twice (once before
      // optimisation and once after), so there can be sites on the block from the previous read.
      geometry.Blocks[block].Sites.clear();

      for (site_t localSiteIndex = 0; localSiteIndex < geometry.GetSitesPerBlock(); ++localSiteIndex)
      {
        geometry.Blocks[block].Sites.push_back(ParseSite(reader));
      }
    }

    GeometrySite GeometryReader::ParseSite(io::writers::xdr::XdrReader& reader)
    {
      // Read the fluid property.
      unsigned isFluid;
      bool success = reader.readUnsignedInt(isFluid);

      if (!success)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Error reading site type");
      }

      GeometrySite readInSite(isFluid != 0);

      // If solid, there's nothing more to do.
      if (!readInSite.isFluid)
      {
        return readInSite;
      }

      const io::formats::geometry::DisplacementVector& neighbourhood = io::formats::geometry::Get().GetNeighbourhood();
      // Prepare the links array to have enough space.
      readInSite.links.resize(latticeInfo.GetNumVectors() - 1);

      // For each link direction...
      for (Direction readDirection = 0; readDirection < neighbourhood.size(); readDirection++)
      {
        // read the type of the intersection and create a link...
        unsigned intersectionType;
        reader.readUnsignedInt(intersectionType);

        GeometrySiteLink link;
        link.type = (GeometrySiteLink::IntersectionType) intersectionType;

        // walls have a floating-point distance to the wall...
        if (link.type == GeometrySiteLink::WALL_INTERSECTION)
        {
          float distance;
          reader.readFloat(distance);
          link.distanceToIntersection = distance;
        }
        // inlets and outlets (which together with none make up the other intersection types)
        // have an iolet id and a distance float...
        else if (link.type != GeometrySiteLink::NO_INTERSECTION)
        {
          float distance;
          unsigned ioletId;
          reader.readUnsignedInt(ioletId);
          reader.readFloat(distance);

          link.ioletId = ioletId;
          link.distanceToIntersection = distance;
        }

        // Now, attempt to match the direction read from the local neighbourhood to one in the
        // lattice being used for simulation. If a match is found, assign the link to the read
        // site.
        for (Direction usedLatticeDirection = 1; usedLatticeDirection < latticeInfo.GetNumVectors();
            usedLatticeDirection++)
        {
          if (latticeInfo.GetVector(usedLatticeDirection) == neighbourhood[readDirection])
          {
            // If this link direction is necessary to the lattice in use, keep the link data.
            readInSite.links[usedLatticeDirection - 1] = link;
            break;
          }
        }
      }

      return readInSite;
    }

    proc_t GeometryReader::GetReadingCoreForBlock(site_t blockNumber)
    {
      return proc_t(blockNumber % util::NumericalFunctions::min(READING_GROUP_SIZE, currentComms.GetSize()));
    }

    /**
     * This function is only called if in Debug mode.
     * @param geometry
     */
    void GeometryReader::ValidateGeometry(const Geometry& geometry)
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Validating the GlobalLatticeData");

      // We check the isFluid property and the link type for each direction
      site_t blockSiteDataLength = geometry.GetSitesPerBlock() * (1 + latticeInfo.GetNumVectors() - 1);

      std::vector<proc_t> myProcForSite;
      std::vector<unsigned> dummySiteData;

      std::vector<proc_t> procForSiteRecv(geometry.GetSitesPerBlock());
      std::vector<unsigned> siteDataRecv(blockSiteDataLength);

      // We also validate that each processor has the same beliefs about each site.
      for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
      {
        // Clear vectors
        myProcForSite.clear();
        dummySiteData.clear();

        if (geometry.Blocks[block].Sites.size() == 0)
        {
          for (site_t localSite = 0; localSite < geometry.GetSitesPerBlock(); ++localSite)
          {
            myProcForSite.push_back(BIG_NUMBER2);
            dummySiteData.push_back(std::numeric_limits<unsigned>::max());
            for (Direction direction = 1; direction < latticeInfo.GetNumVectors(); ++direction)
            {
              dummySiteData.push_back(std::numeric_limits<unsigned>::max());
            }
          }
        }
        else
        {
          for (site_t localSite = 0; localSite < geometry.GetSitesPerBlock(); ++localSite)
          {
            myProcForSite.push_back(geometry.Blocks[block].Sites[localSite].targetProcessor);

            dummySiteData.push_back(geometry.Blocks[block].Sites[localSite].isFluid);

            for (Direction direction = 1; direction < latticeInfo.GetNumVectors(); ++direction)
            {
              if (geometry.Blocks[block].Sites[localSite].isFluid)
              {
                dummySiteData.push_back(geometry.Blocks[block].Sites[localSite].links[direction - 1].type);
              }
              else
              {
                dummySiteData.push_back(std::numeric_limits<unsigned>::max());
              }
            }
          }
        }

        // Reduce using a minimum to find the actual processor for each site (ignoring the
        // BIG_NUMBER2 entries).
        MPI_Allreduce(&myProcForSite[0],
                      &procForSiteRecv[0],
                      (int) geometry.GetSitesPerBlock(),
                      MpiDataType(procForSiteRecv[0]),
                      MPI_MIN,
                      topologyCommunicator);

        MPI_Allreduce(&dummySiteData[0],
                      &siteDataRecv[0],
                      (int) blockSiteDataLength,
                      MpiDataType(dummySiteData[0]),
                      MPI_MIN,
                      topologyCommunicator);

        for (site_t site = 0; site < geometry.GetSitesPerBlock(); ++site)
        {
          if (procForSiteRecv[site] == ConvertTopologyRankToGlobalRank(topologyComms.GetRank())
              && (myProcForSite[site] != ConvertTopologyRankToGlobalRank(topologyComms.GetRank())))
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Other cores think this core has site %li on block %li but it disagrees.",
                                                          site,
                                                          block);
          }
          else if (myProcForSite[site] != BIG_NUMBER2 && procForSiteRecv[site] != myProcForSite[site])
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("This core thought that core %li has site %li on block %li but others think it's on core %li.",
                                                          myProcForSite[site],
                                                          site,
                                                          block,
                                                          procForSiteRecv[site]);
          }

          if (geometry.Blocks[block].Sites.size() > 0)
          {
            if (dummySiteData[site * latticeInfo.GetNumVectors()] != siteDataRecv[site * latticeInfo.GetNumVectors()])
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Different fluid state was found for site %li on block %li. One: %li, Two: %li .",
                                                            site,
                                                            block,
                                                            dummySiteData[site * latticeInfo.GetNumVectors()],
                                                            siteDataRecv[site * latticeInfo.GetNumVectors()]);
            }

            for (Direction direction = 1; direction < latticeInfo.GetNumVectors(); ++direction)
            {
              if (dummySiteData[site * latticeInfo.GetNumVectors() + direction]
                  != siteDataRecv[site * latticeInfo.GetNumVectors() + direction])
              {
                log::Logger::Log<log::Debug, log::OnePerCore>("Different link type was found for site %li, link %i on block %li. One: %li, Two: %li .",
                                                              site,
                                                              direction,
                                                              block,
                                                              dummySiteData[site * latticeInfo.GetNumVectors()
                                                                  + direction],
                                                              siteDataRecv[site * latticeInfo.GetNumVectors()
                                                                  + direction]);
              }
            }

          }
        }
      }
    }

    std::vector<bool> GeometryReader::DecideWhichBlocksToReadIncludingHalo(const Geometry& geometry,
                                                                           const std::vector<proc_t>& unitForEachBlock,
                                                                           proc_t localRank)
    {
      std::vector<bool> shouldReadBlock(geometry.GetBlockCount(), false);

      // Read a block in if it has fluid sites and is to live on the current processor. Also read
      // in any neighbours with fluid sites.
      for (site_t blockI = 0; blockI < geometry.GetBlockDimensions().x; ++blockI)
      {
        for (site_t blockJ = 0; blockJ < geometry.GetBlockDimensions().y; ++blockJ)
        {
          for (site_t blockK = 0; blockK < geometry.GetBlockDimensions().z; ++blockK)
          {
            site_t lBlockId = geometry.GetBlockIdFromBlockCoordinates(blockI, blockJ, blockK);

            if (unitForEachBlock[lBlockId] != localRank)
            {
              continue;
            }

            // Read in all neighbouring blocks.
            for (site_t neighI = util::NumericalFunctions::max<site_t>(0, blockI - 1);
                (neighI <= (blockI + 1)) && (neighI < geometry.GetBlockDimensions().x); ++neighI)
            {
              for (site_t neighJ = util::NumericalFunctions::max<site_t>(0, blockJ - 1);
                  (neighJ <= (blockJ + 1)) && (neighJ < geometry.GetBlockDimensions().y); ++neighJ)
              {
                for (site_t neighK = util::NumericalFunctions::max<site_t>(0, blockK - 1);
                    (neighK <= (blockK + 1)) && (neighK < geometry.GetBlockDimensions().z); ++neighK)
                {
                  site_t lNeighId = geometry.GetBlockIdFromBlockCoordinates(neighI, neighJ, neighK);

                  shouldReadBlock[lNeighId] = true;
                }
              }
            }
          }
        }
      }

      return shouldReadBlock;
    }

    void GeometryReader::OptimiseDomainDecomposition(Geometry& geometry, const std::vector<proc_t>& procForEachBlock)
    {
      timings[hemelb::reporting::Timers::dbg5].Start(); //overall dbg timing
      // Get some arrays that ParMetis needs.
      std::vector<idx_t> vtxDistribn = GetSiteDistributionArray(geometry.GetBlockCount(),
                                                                procForEachBlock,
                                                                fluidSitesOnEachBlock);

      std::vector<idx_t> firstSiteIndexPerBlock(geometry.GetBlockCount());

      GetFirstSiteIndexOnEachBlock(firstSiteIndexPerBlock,
                                   geometry.GetBlockCount(),
                                   vtxDistribn,
                                   procForEachBlock,
                                   fluidSitesOnEachBlock);

      idx_t localVertexCount = vtxDistribn[topologyComms.GetRank() + 1] - vtxDistribn[topologyComms.GetRank()];

      std::vector<idx_t> adjacenciesPerVertex;
      std::vector<idx_t> localAdjacencies;

      GetAdjacencyData(adjacenciesPerVertex,
                       localAdjacencies,
                       geometry,
                       localVertexCount,
                       procForEachBlock,
                       firstSiteIndexPerBlock);

      log::Logger::Log<log::Debug, log::OnePerCore>("Adj length %i", localAdjacencies.size());

      if (ShouldValidate())
      {
        ValidateGraphData(vtxDistribn, localVertexCount, adjacenciesPerVertex, localAdjacencies);
      }

      timings[hemelb::reporting::Timers::dbg5].Stop();

      // Call parmetis.
      std::vector<idx_t> partitionVector(localVertexCount);
      timings[hemelb::reporting::Timers::parmetis].Start();
      log::Logger::Log<log::Info, log::OnePerCore>("Making the call to Parmetis");
      // Reserve on the vectors to be certain they're at least 1 in capacity (so &vector[0] works)
      partitionVector.reserve(1);
      vtxDistribn.reserve(1);
      adjacenciesPerVertex.reserve(1);
      localAdjacencies.reserve(1);
      CallParmetis(partitionVector, localVertexCount, vtxDistribn, adjacenciesPerVertex, localAdjacencies);
      timings[hemelb::reporting::Timers::parmetis].Stop();
      log::Logger::Log<log::Info, log::OnePerCore>("Parmetis has finished.");

      // Convert the ParMetis results into a nice format.
      std::vector<idx_t> allMoves(topologyComms.GetSize());
      timings[hemelb::reporting::Timers::dbg4].Start();
      log::Logger::Log<log::Info, log::OnePerCore>("Getting moves lists for this core.");
      idx_t* movesList = GetMovesList(allMoves,
                                      geometry,
                                      firstSiteIndexPerBlock,
                                      procForEachBlock,
                                      fluidSitesOnEachBlock,
                                      vtxDistribn,
                                      partitionVector);
      log::Logger::Log<log::Info, log::OnePerCore>("Done getting moves lists for this core");
      timings[hemelb::reporting::Timers::dbg4].Stop();

      timings[hemelb::reporting::Timers::reRead].Start();
      log::Logger::Log<log::Info, log::OnePerCore>("Rereading blocks");
      // Reread the blocks based on the ParMetis decomposition.
      RereadBlocks(geometry, allMoves, movesList, procForEachBlock);
      timings[hemelb::reporting::Timers::reRead].Stop();

      timings[hemelb::reporting::Timers::moves].Start();
      // Implement the decomposition now that we have read the necessary data.
      log::Logger::Log<log::Info, log::OnePerCore>("Implementing moves");
      ImplementMoves(geometry, procForEachBlock, allMoves, movesList);
      timings[hemelb::reporting::Timers::moves].Stop();
      delete[] movesList;
    }

    // The header section of the config file contains a number of records.
    site_t GeometryReader::GetHeaderLength(site_t blockCount) const
    {
      return io::formats::geometry::HeaderRecordLength * blockCount;
    }

    /**
     * Get the cumulative count of sites on each processor.
     *
     * @param blockCount
     * @param procForEachBlock
     * @param fluidSitesPerBlock
     * @return Array of the number of vertices on each processor, suitable for passing
     * to ParMetis
     */
    std::vector<idx_t> GeometryReader::GetSiteDistributionArray(site_t blockCount,
                                                                const std::vector<proc_t>& procForEachBlock,
                                                                const std::vector<site_t>& fluidSitesPerBlock) const
    {
      std::vector<idx_t> vertexDistribn(topologyComms.GetSize() + 1, 0);

      // Firstly, count the sites per processor. Do this off-by-one
      // to be compatible with ParMetis.
      for (site_t block = 0; block < blockCount; ++block)
      {
        if (procForEachBlock[block] >= 0)
        {
          vertexDistribn[1 + procForEachBlock[block]] += (idx_t) fluidSitesPerBlock[block];
        }
      }

      // Now make the count cumulative, again off-by-one.
      for (proc_t rank = 0; rank < topologyComms.GetSize(); ++rank)
      {
        vertexDistribn[rank + 1] += vertexDistribn[rank];
      }

      // Validate if we're logging.
      if (ShouldValidate())
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the vertex distribution.");

        // vtxDistribn should be the same on all cores.
        std::vector<idx_t> vtxDistribnRecv(topologyComms.GetSize() + 1);

        MPI_Allreduce(&vertexDistribn[0],
                      &vtxDistribnRecv[0],
                      topologyComms.GetSize() + 1,
                      MpiDataType(vtxDistribnRecv[0]),
                      MPI_MIN,
                      topologyCommunicator);

        for (proc_t rank = 0; rank < topologyComms.GetSize() + 1; ++rank)
        {
          if (vertexDistribn[rank] != vtxDistribnRecv[rank])
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("vertexDistribn[%i] was %li but at least one other core had it as %li.",
                                                          rank,
                                                          vertexDistribn[rank],
                                                          vtxDistribnRecv[rank]);
          }
        }
      }

      return vertexDistribn;
    }

    void GeometryReader::GetFirstSiteIndexOnEachBlock(std::vector<idx_t>& firstSiteIndexPerBlock,
                                                      site_t blockCount,
                                                      const std::vector<idx_t>& vertexDistribution,
                                                      const std::vector<proc_t>& procForEachBlock,
                                                      const std::vector<site_t>& fluidSitesPerBlock) const
    {
      // First calculate the lowest site index on each proc - relatively easy.
      std::vector<idx_t> firstSiteOnProc(vertexDistribution);

      // Now for each block (in ascending order), the smallest site index is the smallest site
      // index on its processor, incremented by the number of sites observed from that processor
      // so far.
      for (site_t block = 0; block < blockCount; ++block)
      {
        proc_t proc = procForEachBlock[block];
        if (proc < 0)
        {
          firstSiteIndexPerBlock[block] = -1;
        }
        else
        {
          firstSiteIndexPerBlock[block] = firstSiteOnProc[proc];
          firstSiteOnProc[proc] += (idx_t) fluidSitesPerBlock[block];
        }
      }

      // Now, if logging debug info, we validate firstSiteIndexPerBlock
      if (ShouldValidate())
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the firstSiteIndexPerBlock values.");

        std::vector<idx_t> firstSiteIndexPerBlockRecv(blockCount);

        // Reduce finding the maximum across all nodes. Note that we have to use the maximum
        // because some cores will have -1 for a block (indicating that it has no neighbours on
        // that block.
        MPI_Allreduce(&firstSiteIndexPerBlock[0],
                      &firstSiteIndexPerBlockRecv[0],
                      (int) blockCount,
                      MpiDataType(firstSiteIndexPerBlock[0]),
                      MPI_MAX,
                      topologyCommunicator);

        for (site_t block = 0; block < blockCount; ++block)
        {
          if (firstSiteIndexPerBlock[block] >= 0 && firstSiteIndexPerBlock[block] != firstSiteIndexPerBlockRecv[block])
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("This core had the first site index on block %li as %li but at least one other core had it as %li.",
                                                          block,
                                                          firstSiteIndexPerBlock[block],
                                                          firstSiteIndexPerBlockRecv[block]);
          }
        }
      }
    }

    void GeometryReader::GetAdjacencyData(std::vector<idx_t>& adjacenciesPerVertex,
                                          std::vector<idx_t>& localAdjacencies,
                                          const Geometry& geometry,
                                          const idx_t localVertexCount,
                                          const std::vector<proc_t>& procForEachBlock,
                                          const std::vector<idx_t>& firstSiteIndexPerBlock) const
    {
      adjacenciesPerVertex.push_back(0);

      // For each block (counting up by lowest site id)...
      for (site_t blockI = 0; blockI < geometry.GetBlockDimensions().x; blockI++)
      {
        for (site_t blockJ = 0; blockJ < geometry.GetBlockDimensions().y; blockJ++)
        {
          for (site_t blockK = 0; blockK < geometry.GetBlockDimensions().z; blockK++)
          {
            const site_t blockNumber = geometry.GetBlockIdFromBlockCoordinates(blockI, blockJ, blockK);

            // ... considering only the ones which live on this proc...
            if (procForEachBlock[blockNumber] != topologyComms.GetRank())
            {
              continue;
            }

            const BlockReadResult& blockReadResult = geometry.Blocks[blockNumber];

            site_t m = -1;

            // ... iterate over sites within the block...
            for (site_t localSiteI = 0; localSiteI < geometry.GetBlockSize(); localSiteI++)
            {
              for (site_t localSiteJ = 0; localSiteJ < geometry.GetBlockSize(); localSiteJ++)
              {
                for (site_t localSiteK = 0; localSiteK < geometry.GetBlockSize(); localSiteK++)
                {
                  ++m;

                  // ... only looking at non-solid sites...
                  if (blockReadResult.Sites[m].targetProcessor == BIG_NUMBER2)
                  {
                    continue;
                  }

                  // ... for each lattice direction...
                  for (unsigned int l = 1; l < latticeInfo.GetNumVectors(); l++)
                  {
                    // ... which leads to a valid neighbouring site...
                    site_t neighbourI = blockI * geometry.GetBlockSize() + localSiteI + latticeInfo.GetVector(l).x;
                    site_t neighbourJ = blockJ * geometry.GetBlockSize() + localSiteJ + latticeInfo.GetVector(l).y;
                    site_t neighbourK = blockK * geometry.GetBlockSize() + localSiteK + latticeInfo.GetVector(l).z;

                    if (neighbourI < 0 || neighbourJ < 0 || neighbourK < 0
                        || neighbourI >= (geometry.GetBlockSize() * geometry.GetBlockDimensions().x)
                        || neighbourJ >= (geometry.GetBlockSize() * geometry.GetBlockDimensions().y)
                        || neighbourK >= (geometry.GetBlockSize() * geometry.GetBlockDimensions().z))
                    {
                      continue;
                    }

                    // ... (that is actually being simulated and not a solid)...
                    site_t neighbourBlockI = neighbourI / geometry.GetBlockSize();
                    site_t neighbourBlockJ = neighbourJ / geometry.GetBlockSize();
                    site_t neighbourBlockK = neighbourK / geometry.GetBlockSize();

                    site_t neighbourSiteI = neighbourI % geometry.GetBlockSize();
                    site_t neighbourSiteJ = neighbourJ % geometry.GetBlockSize();
                    site_t neighbourSiteK = neighbourK % geometry.GetBlockSize();

                    site_t neighbourBlockId = geometry.GetBlockIdFromBlockCoordinates(neighbourBlockI,
                                                                                      neighbourBlockJ,
                                                                                      neighbourBlockK);

                    const BlockReadResult& neighbourBlock = geometry.Blocks[neighbourBlockId];

                    site_t neighbourSiteId = geometry.GetSiteIdFromSiteCoordinates(neighbourSiteI,
                                                                                   neighbourSiteJ,
                                                                                   neighbourSiteK);

                    if (neighbourBlock.Sites.size() == 0
                        || neighbourBlock.Sites[neighbourSiteId].targetProcessor == BIG_NUMBER2)
                    {
                      continue;
                    }

                    // Calculate the site's id over the whole geometry,
                    site_t neighGlobalSiteId = firstSiteIndexPerBlock[neighbourBlockId];

                    for (site_t neighSite = 0; neighSite < geometry.GetSitesPerBlock(); ++neighSite)
                    {
                      if (neighSite == neighbourSiteId)
                      {
                        break;
                      }
                      else if (neighbourBlock.Sites[neighSite].targetProcessor != BIG_NUMBER2)
                      {
                        ++neighGlobalSiteId;
                      }
                    }

                    // then add this to the list of adjacencies.
                    localAdjacencies.push_back((idx_t) neighGlobalSiteId);
                  }

                  // The cumulative count of adjacencies for this vertex is equal to the total
                  // number of adjacencies we've entered.
                  // NOTE: The prefix operator is correct here because
                  // the array has a leading 0 not relating to any site.
                  adjacenciesPerVertex.push_back(localAdjacencies.size());
                }
              }
            }
          }
        }
      }
    }

    void GeometryReader::CallParmetis(std::vector<idx_t>& partitionVector,
                                      idx_t localVertexCount,
                                      std::vector<idx_t>& vtxDistribn,
                                      std::vector<idx_t>& adjacenciesPerVertex,
                                      std::vector<idx_t>& adjacencies)
    {
      // From the ParMETIS documentation:
      // --------------------------------
      // Processor Pi holds ni consecutive vertices and mi corresponding edges
      //
      // xadj[ni+1] has the cumulative number of adjacencies per vertex (with a leading 0 on each processor)
      // vwgt[ni] has vertex weight coefficients and can be NULL
      // adjncy[mi] has the adjacent vertices for each edge (using a global index, starting at 0)
      // adjwgt[mi] has edge weights and can be NULL
      // vtxdist[P+1] has an identical array of the number of the vertices on each processor, cumulatively.
      //           So Pi has vertices from vtxdist[i] to vtxdist[i+1]-1
      // wgtflag* is 0 with no weights (1 on edges, 2 on vertices, 3 on edges & vertices)
      // numflag* is 0 for C-style numbering (1 for Fortran-style)
      // ncon* is the number of weights on each vertex
      // nparts* is the number of sub-domains (partition domains) desired
      // tpwgts* is the fraction of vertex weight to apply to each sub-domain
      // ubvec* is an array of the imbalance tolerance for each vertex weight
      // options* is an int array of options
      //
      // edgecut[1] will contain the number of edges cut by the partitioning
      // part[ni] will contain the partition vector of the locally-stored vertices
      // comm* is a pointer to the MPI communicator of the processes involved

      // Initialise the partition vector.
      for (idx_t vertexIndex = 0; vertexIndex < localVertexCount; ++vertexIndex)
      {
        partitionVector[vertexIndex] = topologyComms.GetRank();
      }

      // Weight all vertices evenly.
      std::vector<idx_t> vertexWeight(localVertexCount, 1);

      // Set the weights of each partition to be even, and to sum to 1.
      idx_t desiredPartitionSize = topologyComms.GetSize();

      std::vector<real_t> domainWeights(desiredPartitionSize, (real_t)(1.0) / ((real_t) desiredPartitionSize));

      // A bunch of values ParMetis needs.
      idx_t noConstraints = 1;
      idx_t weightFlag = 2;
      idx_t numberingFlag = 0;
      idx_t edgesCut = 0;
      idx_t options[4] = { 0, 0, 0, 0 };

      if (ShouldValidate())
      {
        // Specify that some options are set and that we should
        // debug everything.

        // Specify that we have set some options
        options[0] = 1;
        // From parmetis.h
        // We get timing info (1)
        // more timing info (2)
        // details of the graph-coarsening process (4)
        // info during graph refinement (8)
        // NOT info on matching (16)
        // info on communication during matching (32)
        // info on remappining (64)
        options[1] = 1 | 2 | 4 | 8 | 32 | 64;
      }

      real_t tolerance = 1.001F;

      log::Logger::Log<log::Info, log::OnePerCore>("Calling ParMetis");

      // Reserve 1 on these vectors so that the reference to their first element
      // exists (even if it's unused).
      partitionVector.reserve(1);
      vertexWeight.reserve(1);

      ParMETIS_V3_PartKway(&vtxDistribn[0],
                           &adjacenciesPerVertex[0],
                           &adjacencies[0],
                           &vertexWeight[0],
                           NULL,
                           &weightFlag,
                           &numberingFlag,
                           &noConstraints,
                           &desiredPartitionSize,
                           &domainWeights[0],
                           &tolerance,
                           options,
                           &edgesCut,
                           &partitionVector[0],
                           &topologyCommunicator);
      log::Logger::Log<log::Info, log::OnePerCore>("ParMetis returned.");
    }

    void GeometryReader::ValidateGraphData(const std::vector<idx_t>& vtxDistribn,
                                           idx_t localVertexCount,
                                           const std::vector<idx_t>& adjacenciesPerVertex,
                                           const std::vector<idx_t>& adjacencies)
    {
      // If we're using debugging logs, check that the arguments are consistent across all cores.
      // To verify: vtxDistribn, adjacenciesPerVertex, adjacencies
      if (ShouldValidate())
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the graph adjacency structure");

        // Create an array of lists to store all of this node's adjacencies, arranged by the
        // proc the adjacent vertex is on.
        std::vector<std::multimap<idx_t, idx_t> > adjByNeighProc(topologyComms.GetSize(),
                                                                 std::multimap<idx_t, idx_t>());

        // The adjacency data should correspond across all cores.
        for (idx_t index = 0; index < localVertexCount; ++index)
        {
          idx_t vertex = vtxDistribn[topologyComms.GetRank()] + index;

          // Iterate over each adjacency (of each vertex).
          for (idx_t adjNumber = 0; adjNumber < (adjacenciesPerVertex[index + 1] - adjacenciesPerVertex[index]);
              ++adjNumber)
          {
            idx_t adjacentVertex = adjacencies[adjacenciesPerVertex[index] + adjNumber];
            proc_t adjacentProc = -1;

            // Calculate the proc of the neighbouring vertex.
            for (proc_t proc = 0; proc < topologyComms.GetSize(); ++proc)
            {
              if (vtxDistribn[proc] <= adjacentVertex && vtxDistribn[proc + 1] > adjacentVertex)
              {
                adjacentProc = proc;
                break;
              }
            }

            // If it doesn't appear to belong to any proc, something's wrong.
            if (adjacentProc == -1)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("The vertex %li has a neighbour %li which doesn\'t appear to live on any processor.",
                                                            vertex,
                                                            adjacentVertex);
              continue;
            }

            // Store the data if it does belong to a proc.
            adjByNeighProc[adjacentProc].insert(std::pair<idx_t, idx_t>(adjacentVertex, vertex));
          }
        }

        // Create variables for the neighbour data to go into.
        std::vector<idx_t> counts(topologyComms.GetSize());
        std::vector<std::vector<idx_t> > data(topologyComms.GetSize());
        std::vector<MPI_Request> requests(2 * topologyComms.GetSize());

        log::Logger::Log<log::Debug, log::OnePerCore>("Validating neighbour data");

        // Now spread and compare the adjacency information. Larger ranks send data to smaller
        // ranks which receive the data and compare it.
        for (proc_t neigh = 0; neigh < (proc_t) topologyComms.GetSize(); ++neigh)
        {
          requests[2 * neigh] = MPI_REQUEST_NULL;
          requests[2 * neigh + 1] = MPI_REQUEST_NULL;

          if (neigh < topologyComms.GetRank())
          {
            // Send the array length.
            counts[neigh] = 2 * adjByNeighProc[neigh].size();
            MPI_Isend(&counts[neigh], 1, MpiDataType(counts[0]), neigh, 42, topologyCommunicator, &requests[2 * neigh]);

            MPI_Wait(&requests[2 * neigh], MPI_STATUS_IGNORE);

            // Create a sendable array (std::lists aren't organised in a sendable format).
            data[neigh].resize(counts[neigh]);

            unsigned int adjacencyIndex = 0;

            for (std::multimap<idx_t, idx_t>::iterator it = adjByNeighProc[neigh].begin();
                it != adjByNeighProc[neigh].end(); ++it)

            {
              data[neigh][2 * adjacencyIndex] = it->first;
              data[neigh][2 * adjacencyIndex + 1] = it->second;
              ++adjacencyIndex;
            }

            // Send the data to the neighbour.
            MPI_Isend(&data[neigh][0],
                      (int) counts[neigh],
                      MpiDataType<idx_t>(),
                      neigh,
                      43,
                      topologyCommunicator,
                      &requests[2 * neigh + 1]);

            MPI_Wait(&requests[2 * neigh + 1], MPI_STATUS_IGNORE);

            // Sending arrays don't perform comparison.
            continue;
          }

          // If this is a greater rank number than the neighbour, receive the data.
          if (neigh > topologyComms.GetRank())
          {
            MPI_Irecv(&counts[neigh],
                      1,
                      MpiDataType(counts[neigh]),
                      neigh,
                      42,
                      topologyCommunicator,
                      &requests[2 * neigh]);

            MPI_Wait(&requests[2 * neigh], MPI_STATUS_IGNORE);

            data[neigh].resize(counts[neigh]);

            MPI_Irecv(&data[neigh][0],
                      (int) counts[neigh],
                      MpiDataType<idx_t>(),
                      neigh,
                      43,
                      topologyCommunicator,
                      &requests[2 * neigh + 1]);

            MPI_Wait(&requests[2 * neigh] + 1, MPI_STATUS_IGNORE);
          }
          // Neigh == mTopologyRank, i.e. neighbouring vertices on the same proc
          // Duplicate the data.

          else
          {
            counts[neigh] = 2 * adjByNeighProc[neigh].size();
            data[neigh].resize(counts[neigh]);

            int adjacencyIndex = 0;
            for (std::multimap<idx_t, idx_t>::iterator it = adjByNeighProc[neigh].begin();
                it != adjByNeighProc[neigh].end(); ++it)

            {
              data[neigh][2 * adjacencyIndex] = it->first;
              data[neigh][2 * adjacencyIndex + 1] = it->second;
              ++adjacencyIndex;
            }
          }

          // Now we compare. First go through the received data which is ordered as (adjacent
          // vertex, vertex) wrt the neighbouring proc.
          for (idx_t ii = 0; ii < counts[neigh]; ii += 2)
          {
            bool found = false;

            // Go through each neighbour we know about on this proc, and check whether it
            // matches the current received neighbour-data.
            for (std::multimap<idx_t, idx_t>::iterator it = adjByNeighProc[neigh].find(data[neigh][ii + 1]);
                it != adjByNeighProc[neigh].end(); ++it)
            {
              idx_t recvAdj = it->first;
              idx_t recvAdj2 = it->second;

              if (data[neigh][ii] == recvAdj2 && data[neigh][ii + 1] == recvAdj)
              {
                adjByNeighProc[neigh].erase(it);
                found = true;
                break;
              }
            }

            // No neighbour data on this proc matched the data received.
            if (!found)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Neighbour proc %i had adjacency (%li,%li) that wasn't present on this processor.",
                                                            neigh,
                                                            data[neigh][ii],
                                                            data[neigh][ii + 1]);
            }
          }

          // The local store of adjacencies should now be empty, if there was complete matching.
          std::multimap<idx_t, idx_t>::iterator it = adjByNeighProc[neigh].begin();
          while (it != adjByNeighProc[neigh].end())
          {
            idx_t adj1 = it->first;
            idx_t adj2 = it->second;
            ++it;

            // We had neighbour-data on this proc that didn't match that received.
            log::Logger::Log<log::Debug, log::OnePerCore>("The local processor has adjacency (%li,%li) that isn't present on neighbouring processor %i.",
                                                          adj1,
                                                          adj2,
                                                          neigh);
          }
        }
      }
    }

    /**
     * Returns a list of the fluid sites to be moved.
     *
     * NOTE: This function's return value is a dynamically-allocated array of all the moves to be
     * performed, ordered by (origin processor [with a count described by the content of the first
     * parameter], site id on the origin processor). The contents of the array are contiguous
     * triplets of ints: (block id, site id on block, destination rank).
     *
     * @param movesFromEachProc
     * @param vtxDistribn
     * @param partitionVector
     * @return
     */
    idx_t* GeometryReader::GetMovesList(std::vector<idx_t>& movesFromEachProc,
                                        const Geometry& geometry,
                                        const std::vector<idx_t>& firstSiteIndexPerBlock,
                                        const std::vector<proc_t>& procForEachBlock,
                                        const std::vector<site_t>& fluidSitesPerBlock,
                                        const std::vector<idx_t>& vtxDistribn,
                                        const std::vector<idx_t>& partitionVector)
    {
      timings[hemelb::reporting::Timers::dbg1].Start();

      // Right. Let's count how many sites we're going to have to move. Count the local number of
      // sites to be moved, and collect the site id and the destination processor.
      std::vector<idx_t> moveData;

      const idx_t myLowest = vtxDistribn[topologyComms.GetRank()];
      const idx_t myHighest = vtxDistribn[topologyComms.GetRank() + 1] - 1;

      // Create a map for looking up block Ids: the map is from the contiguous site index
      // of the first fluid site on the block, to the block id.
      std::map<site_t, site_t> blockIdLookupByLastSiteIndex;

      for (site_t blockId = 0; blockId < geometry.GetBlockCount(); ++blockId)
      {
        if (procForEachBlock[blockId] >= 0 && procForEachBlock[blockId] != BIG_NUMBER2)
        {
          site_t lastFluidSiteId = firstSiteIndexPerBlock[blockId] + fluidSitesPerBlock[blockId] - 1;
          blockIdLookupByLastSiteIndex[lastFluidSiteId] = blockId;
        }
      }

      // For each local fluid site...
      for (idx_t ii = 0; ii <= (myHighest - myLowest); ++ii)
      {
        // ... if it's going elsewhere...
        if (partitionVector[ii] != topologyComms.GetRank())
        {
          // ... get it's id on the local processor...
          idx_t localFluidSiteId = myLowest + ii;

          timings[hemelb::reporting::Timers::dbg3].Start();
          // ... find out which block it's on, using our lookup map...
          // A feature of std::map::equal_range is that if there's no equal key, both iterators
          // returned will point to the entry with the next greatest key. Since we store block
          // ids by last fluid site number, this immediately gives us the block id.
          std::pair<std::map<site_t, site_t>::iterator, std::map<site_t, site_t>::iterator> rangeMatch =
              blockIdLookupByLastSiteIndex.equal_range(localFluidSiteId);

          idx_t fluidSiteBlock = rangeMatch.first->second;

          // Check the block id is correct
          if (ShouldValidate())
          {
            if (procForEachBlock[fluidSiteBlock] < 0)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Found block %i for site %i but this block has a processor of %i assigned",
                                                            fluidSiteBlock,
                                                            localFluidSiteId,
                                                            procForEachBlock[fluidSiteBlock]);
            }
            if (firstSiteIndexPerBlock[fluidSiteBlock] > localFluidSiteId)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Found block %i for site %i but sites on this block start at number %i",
                                                            fluidSiteBlock,
                                                            localFluidSiteId,
                                                            firstSiteIndexPerBlock[fluidSiteBlock]);
            }
            if (firstSiteIndexPerBlock[fluidSiteBlock] + fluidSitesPerBlock[fluidSiteBlock] - 1 < localFluidSiteId)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Found block %i for site %i but there are %i sites on this block starting at %i",
                                                            fluidSiteBlock,
                                                            localFluidSiteId,
                                                            fluidSitesPerBlock[fluidSiteBlock],
                                                            firstSiteIndexPerBlock[fluidSiteBlock]);
            }
          }

          timings[hemelb::reporting::Timers::dbg3].Stop();

          // ... and find its site id within that block. Start by working out how many fluid sites
          // we have to pass before we arrive at the fluid site we're after...
          idx_t fluidSitesToPass = localFluidSiteId - firstSiteIndexPerBlock[fluidSiteBlock];
          idx_t siteIndex = 0;

          while (true)
          {
            // ... then keep going through the sites on the block until we've passed as many fluid
            // sites as we need to.
            if (geometry.Blocks[fluidSiteBlock].Sites[siteIndex].targetProcessor != BIG_NUMBER2)
            {
              fluidSitesToPass--;
            }
            if (fluidSitesToPass < 0)
            {
              break;
            }
            siteIndex++;
          }

          // The above code could go wrong, so in debug logging mode, we do some extra tests.
          if (ShouldValidate())
          {
            // If we've ended up on an impossible block, or one that doesn't live on this rank,
            // inform the user.
            if (fluidSiteBlock >= geometry.GetBlockCount()
                || procForEachBlock[fluidSiteBlock] != topologyComms.GetRank())
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Partition element %i wrongly assigned to block %u of %i (block on processor %i)",
                                                            ii,
                                                            fluidSiteBlock,
                                                            geometry.GetBlockCount(),
                                                            procForEachBlock[fluidSiteBlock]);
            }

            // Similarly, if we've ended up with an impossible site index, or a solid site,
            // print an error message.
            if (siteIndex >= geometry.GetSitesPerBlock()
                || geometry.Blocks[fluidSiteBlock].Sites[siteIndex].targetProcessor == BIG_NUMBER2)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Partition element %i wrongly assigned to site %u of %i (block %i%s)",
                                                            ii,
                                                            siteIndex,
                                                            fluidSitesPerBlock[fluidSiteBlock],
                                                            fluidSiteBlock,
                                                            geometry.Blocks[fluidSiteBlock].Sites[siteIndex].targetProcessor
                                                                == BIG_NUMBER2 ?
                                                              " and site is solid" :
                                                              "");
            }
          }

          // Add the block, site and destination rank to our move list.
          moveData.push_back(fluidSiteBlock);
          moveData.push_back(siteIndex);
          moveData.push_back(partitionVector[ii]);
        }
      }

      // Spread the move data around
      timings[hemelb::reporting::Timers::dbg1].Stop();
      log::Logger::Log<log::Info, log::OnePerCore>("Starting to spread move data");
      timings[hemelb::reporting::Timers::dbg2].Start();

      // First, for each core, gather a list of which blocks the current core wants to
      // know more data about.
      // Handily, the blocks we want to know about are exactly those for which we already
      // have some data.
      std::map<proc_t, std::vector<site_t> > blockIdsIRequireFromX;

      for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
      {
        if (!geometry.Blocks[block].Sites.empty())
        {
          proc_t residentProc = procForEachBlock[block];

          blockIdsIRequireFromX[residentProc].push_back(block);
          log::Logger::Log<log::Debug, log::OnePerCore>("I require block %i from proc %i (running total %i)",
                                                        block,
                                                        residentProc,
                                                        blockIdsIRequireFromX[residentProc].size());
        }
      }

      timings[hemelb::reporting::Timers::moveForcingNumbers].Start();

      // We also need to force some data upon blocks, i.e. when they're receiving data from a new
      // block they didn't previously want to know about.
      std::map<proc_t, std::vector<site_t> > blockForcedUponX;
      std::vector<proc_t> numberOfBlocksIForceUponX(topologyComms.GetSize(), 0);

      for (idx_t moveNumber = 0; moveNumber < (idx_t) moveData.size(); moveNumber += 3)
      {
        proc_t target_proc = moveData[moveNumber + 2];
        site_t blockId = moveData[moveNumber];

        if (std::count(blockForcedUponX[target_proc].begin(), blockForcedUponX[target_proc].end(), blockId) == 0)
        {
          blockForcedUponX[target_proc].push_back(blockId);
          ++numberOfBlocksIForceUponX[target_proc];

          log::Logger::Log<log::Debug, log::OnePerCore>("I'm ensuring proc %i takes data about block %i",
                                                        target_proc,
                                                        blockId);
        }
      }

      // Now find how many blocks are being forced upon us from every other core.
      net::Net netForMoveSending(topologyCommunicator);

      std::vector<proc_t> blocksForcedOnMe(topologyComms.GetSize(), 0);

      log::Logger::Log<log::Info, log::OnePerCore>("Moving forcing block numbers");

      MPI_Alltoall(&numberOfBlocksIForceUponX[0],
                   1,
                   MpiDataType<proc_t>(),
                   &blocksForcedOnMe[0],
                   1,
                   MpiDataType<proc_t>(),
                   topologyCommunicator);

      timings[hemelb::reporting::Timers::moveForcingNumbers].Stop();
      timings[hemelb::reporting::Timers::moveForcingData].Start();

      // Now get all the blocks being forced upon me.
      std::map<proc_t, std::vector<site_t> > blocksForcedOnMeByEachProc;
      for (proc_t otherProc = 0; otherProc < (proc_t) topologyComms.GetSize(); ++otherProc)
      {
        if (blocksForcedOnMe[otherProc] > 0)
        {
          blocksForcedOnMeByEachProc[otherProc] = std::vector<site_t>(blocksForcedOnMe[otherProc]);
          netForMoveSending.RequestReceive(&blocksForcedOnMeByEachProc[otherProc][0],
                                           blocksForcedOnMe[otherProc],
                                           otherProc);
        }

        if (numberOfBlocksIForceUponX[otherProc] > 0)
        {
          netForMoveSending.RequestSend(&blockForcedUponX[otherProc][0],
                                        numberOfBlocksIForceUponX[otherProc],
                                        otherProc);
        }

        log::Logger::Log<log::Debug, log::OnePerCore>("I'm forcing %i blocks on proc %i.",
                                                      numberOfBlocksIForceUponX[otherProc],
                                                      otherProc);
      }

      log::Logger::Log<log::Info, log::OnePerCore>("Moving forcing block ids");

      netForMoveSending.Receive();
      netForMoveSending.Send();
      netForMoveSending.Wait();

      // Now go through every block forced upon me and add it to the list of ones I want.
      for (proc_t otherProc = 0; otherProc < (proc_t) topologyComms.GetSize(); ++otherProc)
      {
        if (blocksForcedOnMe[otherProc] > 0)
        {
          for (std::vector<site_t>::iterator it = blocksForcedOnMeByEachProc[otherProc].begin();
              it != blocksForcedOnMeByEachProc[otherProc].end(); ++it)
          {
            if (std::count(blockIdsIRequireFromX[otherProc].begin(), blockIdsIRequireFromX[otherProc].end(), *it) == 0)
            {
              blockIdsIRequireFromX[otherProc].push_back(*it);

              log::Logger::Log<log::Debug, log::OnePerCore>("I'm being forced to take block %i from proc %i",
                                                            *it,
                                                            otherProc);
            }

            // We also need to take all neighbours of the forced block from their processors.
            BlockLocation blockCoords = geometry.GetBlockCoordinatesFromBlockId(*it);

            // Iterate over every direction we might need (except 0 as we obviously already have
            // that block in the list).
            for (Direction direction = 1; direction < lb::lattices::D3Q27::NUMVECTORS; ++direction)
            {
              // Calculate the putative neighbour's coordinates...
              BlockLocation neighbourCoords = blockCoords
                  + BlockLocation(lb::lattices::D3Q27::CX[direction],
                                  lb::lattices::D3Q27::CY[direction],
                                  lb::lattices::D3Q27::CZ[direction]);

              // If the neighbour is a real block...
              if (geometry.AreBlockCoordinatesValid(neighbourCoords))
              {
                // Get the block id, and check whether it has any fluid sites...
                site_t neighbourBlockId = geometry.GetBlockIdFromBlockCoordinates(neighbourCoords.x,
                                                                                  neighbourCoords.y,
                                                                                  neighbourCoords.z);
                proc_t neighbourBlockProc = procForEachBlock[neighbourBlockId];

                if (neighbourBlockProc >= 0)
                {
                  // Check whether this is a block we're already interested in from that neighbour.
                  if (std::count(blockIdsIRequireFromX[neighbourBlockProc].begin(),
                                 blockIdsIRequireFromX[neighbourBlockProc].end(),
                                 neighbourBlockId) == 0)
                  {
                    // Then add it to the list of blocks we're getting from that neighbour.
                    blockIdsIRequireFromX[neighbourBlockProc].push_back(neighbourBlockId);

                    log::Logger::Log<log::Debug, log::OnePerCore>("I need to also take block %i from proc %i",
                                                                  neighbourBlockId,
                                                                  neighbourBlockProc);

                  }
                }
              }
            }
          }
        }
      }

      timings[hemelb::reporting::Timers::moveForcingData].Stop();
      timings[hemelb::reporting::Timers::blockRequirements].Start();

      log::Logger::Log<log::Info, log::OnePerCore>("Calculating block requirements");

      // Now we want to spread this info around so that each core knows which blocks each other
      // requires from it.
      std::vector<site_t> numberOfBlocksRequiredFrom(topologyComms.GetSize(), 0);
      std::vector<site_t> numberOfBlocksXRequiresFromMe(topologyComms.GetSize(), 0);

      // Populate numberOfBlocksRequiredFrom
      for (proc_t otherProc = 0; otherProc < (proc_t) topologyComms.GetSize(); ++otherProc)
      {
        numberOfBlocksRequiredFrom[otherProc] = blockIdsIRequireFromX.count(otherProc) == 0 ?
          0 :
          blockIdsIRequireFromX[otherProc].size();

        log::Logger::Log<log::Debug, log::OnePerCore>("I require a total of %i blocks from proc %i",
                                                      numberOfBlocksRequiredFrom[otherProc],
                                                      otherProc);
      }

      // Now perform the exchange s.t. each core knows how many blocks are required of it from
      // each other core.
      MPI_Alltoall(&numberOfBlocksRequiredFrom[0],
                   1,
                   MpiDataType<site_t>(),
                   &numberOfBlocksXRequiresFromMe[0],
                   1,
                   MpiDataType<site_t>(),
                   topologyCommunicator);

      // Awesome. Now we need to get a list of all the blocks wanted from each core by each other
      // core.
      std::map<proc_t, std::vector<site_t> > blockIdsXRequiresFromMe;

      for (proc_t otherProc = 0; otherProc < (proc_t) topologyComms.GetSize(); ++otherProc)
      {
        blockIdsXRequiresFromMe[otherProc] = std::vector<site_t>(numberOfBlocksXRequiresFromMe[otherProc]);

        log::Logger::Log<log::Debug, log::OnePerCore>("Proc %i requires %i blocks from me",
                                                      otherProc,
                                                      blockIdsXRequiresFromMe[otherProc].size());

        netForMoveSending.RequestReceive(&blockIdsXRequiresFromMe[otherProc][0],
                                         numberOfBlocksXRequiresFromMe[otherProc],
                                         otherProc);

        netForMoveSending.RequestSend(&blockIdsIRequireFromX[otherProc][0],
                                      numberOfBlocksRequiredFrom[otherProc],
                                      otherProc);
      }

      netForMoveSending.Receive();
      netForMoveSending.Send();
      netForMoveSending.Wait();

      timings[hemelb::reporting::Timers::blockRequirements].Stop();

      timings[hemelb::reporting::Timers::moveCountsSending].Start();

      // OK, now to get on with the actual sending of the data...
      // Except first, it'll be helpful to organise things by blocks.
      std::map<site_t, std::vector<proc_t> > coresInterestedInEachBlock;
      std::map<site_t, std::vector<idx_t> > moveDataForEachBlock;
      std::map<site_t, idx_t> movesForEachLocalBlock;

      // Initialise the moves for each local block to 0. This handles an edge case where a local
      // block has no moves.
      for (site_t blockId = 0; blockId < geometry.GetBlockCount(); ++blockId)
      {
        if (procForEachBlock[blockId] == topologyComms.GetRank())
        {
          movesForEachLocalBlock[blockId] = 0;
        }
      }

      for (proc_t otherProc = 0; otherProc < (proc_t) topologyComms.GetSize(); ++otherProc)
      {
        for (site_t blockNum = 0; blockNum < (site_t) blockIdsXRequiresFromMe[otherProc].size(); ++blockNum)
        {
          site_t blockId = blockIdsXRequiresFromMe[otherProc][blockNum];

          log::Logger::Log<log::Debug, log::OnePerCore>("Proc %i requires block %i from me", otherProc, blockId);

          if (coresInterestedInEachBlock.count(blockId) == 0)
          {
            coresInterestedInEachBlock[blockId] = std::vector<proc_t>();
          }

          coresInterestedInEachBlock[blockId].push_back(otherProc);
        }
      }

      for (site_t moveNumber = 0; moveNumber < (site_t) moveData.size(); moveNumber += 3)
      {
        site_t blockId = moveData[moveNumber];

        if (moveDataForEachBlock.count(blockId) == 0)
        {
          moveDataForEachBlock[blockId] = std::vector<idx_t>();
        }

        moveDataForEachBlock[blockId].push_back(blockId);
        moveDataForEachBlock[blockId].push_back(moveData[moveNumber + 1]);
        moveDataForEachBlock[blockId].push_back(moveData[moveNumber + 2]);
        movesForEachLocalBlock[blockId]++;
      }

      // And it'll also be super-handy to know how many moves we're going to have locally.
      std::vector<idx_t> movesForEachBlockWeCareAbout(geometry.GetBlockCount(), 0);

      for (proc_t otherProc = 0; otherProc < (proc_t) topologyComms.GetSize(); ++otherProc)
      {
        for (std::vector<site_t>::iterator it = blockIdsIRequireFromX[otherProc].begin();
            it != blockIdsIRequireFromX[otherProc].end(); ++it)
        {
          netForMoveSending.RequestReceive(&movesForEachBlockWeCareAbout[*it], 1, otherProc);
          log::Logger::Log<log::Debug, log::OnePerCore>("I want the move count for block %i from proc %i",
                                                        *it,
                                                        otherProc);
        }

        for (std::vector<site_t>::iterator it = blockIdsXRequiresFromMe[otherProc].begin();
            it != blockIdsXRequiresFromMe[otherProc].end(); ++it)
        {
          netForMoveSending.RequestSend(&movesForEachLocalBlock[*it], 1, otherProc);
          log::Logger::Log<log::Debug, log::OnePerCore>("I'm sending move count for block %i to proc %i",
                                                        *it,
                                                        otherProc);
        }
      }

      log::Logger::Log<log::Info, log::OnePerCore>("Sending move counts");

      netForMoveSending.Receive();
      netForMoveSending.Send();
      netForMoveSending.Wait();

      timings[hemelb::reporting::Timers::moveCountsSending].Stop();
      timings[hemelb::reporting::Timers::moveDataSending].Start();

      idx_t totalMovesToReceive = 0;

      for (site_t blockId = 0; blockId < geometry.GetBlockCount(); ++blockId)
      {
        totalMovesToReceive += movesForEachBlockWeCareAbout[blockId];
      }

      log::Logger::Log<log::Debug, log::OnePerCore>("I'm expecting a total of %i moves", totalMovesToReceive);

      // Gather the moves to the places they need to go to.
      // Moves list has block, site id, destination proc
      idx_t* movesList = new idx_t[totalMovesToReceive * 3];

      idx_t localMoveId = 0;

      for (proc_t otherProc = 0; otherProc < (proc_t) topologyComms.GetSize(); ++otherProc)
      {
        movesFromEachProc[otherProc] = 0;

        for (std::vector<site_t>::iterator it = blockIdsIRequireFromX[otherProc].begin();
            it != blockIdsIRequireFromX[otherProc].end(); ++it)
        {
          if (movesForEachBlockWeCareAbout[*it] > 0)
          {
            netForMoveSending.RequestReceive(&movesList[localMoveId * 3],
                                             3 * movesForEachBlockWeCareAbout[*it],
                                             otherProc);
            localMoveId += movesForEachBlockWeCareAbout[*it];
            movesFromEachProc[otherProc] += movesForEachBlockWeCareAbout[*it];

            log::Logger::Log<log::Debug, log::OnePerCore>("Expect %i moves from from proc %i about block %i",
                                                          movesForEachBlockWeCareAbout[*it],
                                                          otherProc,
                                                          *it);
          }
        }

        for (std::vector<site_t>::iterator it = blockIdsXRequiresFromMe[otherProc].begin();
            it != blockIdsXRequiresFromMe[otherProc].end(); ++it)
        {
          if (moveDataForEachBlock[*it].size() > 0)
          {
            netForMoveSending.RequestSend(&moveDataForEachBlock[*it][0], moveDataForEachBlock[*it].size(), otherProc);

            log::Logger::Log<log::Debug, log::OnePerCore>("Sending %i moves from to proc %i about block %i",
                                                          moveDataForEachBlock[*it].size() / 3,
                                                          otherProc,
                                                          *it);
          }
        }

        log::Logger::Log<log::Debug, log::OnePerCore>("%i moves from proc %i", movesFromEachProc[otherProc], otherProc);
      }

      log::Logger::Log<log::Info, log::OnePerCore>("Sending move data");

      netForMoveSending.Receive();
      netForMoveSending.Send();
      netForMoveSending.Wait();

      timings[hemelb::reporting::Timers::moveDataSending].Stop();
      timings[hemelb::reporting::Timers::dbg2].Stop();

      // ... and return the list of moves.
      return movesList;
    }

    void GeometryReader::RereadBlocks(Geometry& geometry,
                                      const std::vector<idx_t>& movesPerProc,
                                      const idx_t* movesList,
                                      const std::vector<int>& procForEachBlock)
    {
      // Initialise the array (of which proc each block belongs to) to what it was before.
      std::vector<int> newProcForEachBlock(geometry.GetBlockCount());

      for (site_t blockNumber = 0; blockNumber < geometry.GetBlockCount(); ++blockNumber)
      {
        newProcForEachBlock[blockNumber] = procForEachBlock[blockNumber];
      }

      // Set the proc for each block to be the current proc whenever a site on that block is
      // going to be moved to the current proc.
      idx_t moveIndex = 0;

      for (proc_t fromProc = 0; fromProc < topologyComms.GetSize(); ++fromProc)
      {
        for (idx_t moveNumber = 0; moveNumber < movesPerProc[fromProc]; ++moveNumber)
        {
          idx_t block = movesList[3 * moveIndex];
          idx_t toProc = movesList[3 * moveIndex + 2];
          ++moveIndex;

          if (toProc == (idx_t) topologyComms.GetRank())
          {
            newProcForEachBlock[block] = topologyComms.GetRank();
          }
        }
      }

      // Reread the blocks into the GlobalLatticeData now.
      ReadInBlocksWithHalo(geometry, newProcForEachBlock, topologyComms.GetRank());
    }

    void GeometryReader::ImplementMoves(Geometry& geometry,
                                        const std::vector<proc_t>& procForEachBlock,
                                        const std::vector<idx_t>& movesFromEachProc,
                                        const idx_t* movesList) const
    {
      // First all, set the proc rank for each site to what it originally was before
      // domain decomposition optimisation. Go through each block...
      for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
      {
        // If this proc has owned a fluid site on this block either before or after optimisation,
        // the following will be non-null.
        if (geometry.Blocks[block].Sites.size() > 0)
        {
          // Get the original proc for that block.
          proc_t originalProc = procForEachBlock[block];

          // For each site on that block...
          for (site_t siteIndex = 0; siteIndex < geometry.GetSitesPerBlock(); ++siteIndex)
          {
            // ... if the site is non-solid...
            if (geometry.Blocks[block].Sites[siteIndex].targetProcessor != BIG_NUMBER2)
            {
              // ... set its rank to be the rank it had before optimisation.
              geometry.Blocks[block].Sites[siteIndex].targetProcessor = ConvertTopologyRankToGlobalRank(originalProc);
            }
          }
        }
      }

      // Now implement the moves suggested by parmetis.
      idx_t moveIndex = 0;

      // For each source proc, go through as many moves as it had.
      for (proc_t fromProc = 0; fromProc < topologyComms.GetSize(); ++fromProc)
      {
        for (idx_t moveNumber = 0; moveNumber < movesFromEachProc[fromProc]; ++moveNumber)
        {
          // For each move, get the block, site and destination proc.
          idx_t block = movesList[3 * moveIndex];
          idx_t site = movesList[3 * moveIndex + 1];
          idx_t toProc = movesList[3 * moveIndex + 2];

          // Only implement the move if we have read that block's data.
          if (geometry.Blocks[block].Sites.size() > 0)
          {
            // Some logging code - the unmodified rank for each move's site should equal
            // lFromProc.
            if (ShouldValidate())
            {
              if (geometry.Blocks[block].Sites[site].targetProcessor
                  != ConvertTopologyRankToGlobalRank((proc_t) fromProc))
              {
                log::Logger::Log<log::Debug, log::OnePerCore>("Block %ld, site %ld from move %u was originally on proc %i, not proc %u.",
                                                              block,
                                                              site,
                                                              moveIndex,
                                                              geometry.Blocks[block].Sites[site].targetProcessor,
                                                              fromProc);
              }
            }

            // Implement the move.
            geometry.Blocks[block].Sites[site].targetProcessor = ConvertTopologyRankToGlobalRank((proc_t) toProc);
          }

          ++moveIndex;
        }
      }
    }

    proc_t GeometryReader::ConvertTopologyRankToGlobalRank(proc_t topologyRankIn) const
    {
      // If the global rank is not equal to the topology rank, we are not using rank 0 for
      // LBM.
      return (topology::NetworkTopology::Instance()->GetLocalRank() == topologyComms.GetRank()) ?
        topologyRankIn :
        (topologyRankIn + 1);
    }

    bool GeometryReader::ShouldValidate() const
    {
      return log::Logger::ShouldDisplay<log::Debug>();
    }
  }
}
