#include "extraction/LocalDistributionInput.h"
#include "extraction/OutputField.h"
#include "extraction/LocalPropertyOutput.h"
#include "geometry/LatticeData.h"
#include "io/formats/formats.h"
#include "io/formats/extraction.h"
#include "io/formats/offset.h"
#include "io/writers/xdr/XdrFileReader.h"

namespace hemelb
{
  namespace extraction
  {
    LocalDistributionInput::LocalDistributionInput(const std::string dataFilePath,
						   const net::IOCommunicator& ioComm) :
      comms(ioComm), filePath(dataFilePath)
    {
    }

    LocalDistributionInput::~LocalDistributionInput()
    {
    }

    namespace {
      // Helper class to simplify pulling data out of the file and XDR
      // decoding it
      struct MpiFileReader : public io::writers::xdr::XdrReader {

	MpiFileReader(net::MpiFile& f, size_t len) : file(f), buffer(len) {
	  file.Read(buffer);
	  xdrmem_create(&mXdr, buffer.data(), buffer.size(), XDR_DECODE);
	}

      private:
	net::MpiFile& file;
	std::vector<char> buffer;
      };
    }

    void LocalDistributionInput::LoadDistribution(geometry::LatticeData* latDat)
    {
      // We could supply hints regarding how the file should be read
      // but we are not doing so yet.

      // Open the file as read-only.
      // TODO: raise an exception if the file does not exist.
      auto inputFile = net::MpiFile::Open(comms, filePath, MPI_MODE_RDONLY);
      // Set the view to the file.
      inputFile.SetView(0, MPI_CHAR, MPI_CHAR, "native");
      ReadExtractionHeaders(inputFile);

      // Now open the offset file.
      ReadOffsets(io::formats::offset::ExtractionToOffset(filePath));

      // Read the local part of the checkpoint
      unsigned readLength = localStop - localStart;

      std::vector<char> dataBuffer(readLength);
      inputFile.ReadAt(localStart, dataBuffer);
      io::writers::xdr::XdrMemReader dataReader(&dataBuffer[0], readLength);

      // Read the timestep
      if (comms.OnIORank()) {
	dataReader.read(timestep);
      }
      comms.Broadcast(timestep, comms.GetIORank());

      site_t iSite = 0;
      // while (dataReader.GetPosition() < readLength) {
      for (; dataReader.GetPosition() < readLength; iSite++) {
	// Read the grid coord and check it's consistent with latDat.
	[&](){
	  // Stored as 32 b unsigned
	  util::Vector3D<uint32_t> tmp;
	  dataReader.read(tmp.x);
	  dataReader.read(tmp.y);
	  dataReader.read(tmp.z);

	  // Convert to canonical type
	  util::Vector3D<site_t> grid = tmp;
	  // Look up the site ID and rank, as decomposed by this run
	  // of HemeLB, for the grid coordinate read from the
	  // checkpoint file.
	  proc_t rank; site_t index;
	  if (!latDat->GetContiguousSiteId(grid, rank, index)) {
	    // function returns a 'valid' flag
	    throw Exception() << "Cannot get valid site from extracted site coordinate";
	  }
	  if (rank != comms.Rank())
	    throw Exception() << "Site read on rank " << comms.Rank()
			      << " but should be read on " << rank;
	  if (index != iSite)
	    throw Exception() << "Site read at index " << iSite
			      << " but should be read at " << index;
	}();

	distribn_t* f_old_p = latDat->GetFOld(iSite * LatticeType::NUMVECTORS);
	distribn_t* f_new_p = latDat->GetFNew(iSite * LatticeType::NUMVECTORS);
	// distField.numberOfFloats is read on IO rank and checked to
	// be equal to LatticeType::NUMVECTORS so we use that instead
	// of broadcasting and storing.
	for (int i = 0; i < LatticeType::NUMVECTORS; i++) {
	  float field_val;
	  dataReader.read(field_val);
	  field_val += distField.offset;
	  f_new_p[i] = f_old_p[i] = field_val;
	}
      }

      if (iSite != latDat->GetLocalFluidSiteCount())
	throw Exception() << "Read " << iSite
			  << " sites but expected " << latDat->GetLocalFluidSiteCount();
    }

    void LocalDistributionInput::ReadExtractionHeaders(net::MpiFile& inputFile) {
      namespace fmt = hemelb::io::formats;

      if (comms.OnIORank()) {
	MpiFileReader preambleReader(inputFile, fmt::extraction::MainHeaderLength);

	// Read the magic numbers.
	uint32_t hlbMagicNumber, extMagicNumber, version;
	preambleReader.read(hlbMagicNumber);
	preambleReader.read(extMagicNumber);
	preambleReader.read(version);

	// Check the value of the HemeLB magic number.
	if (hlbMagicNumber != io::formats::HemeLbMagicNumber)
	{
	  throw Exception() << "This file does not start with the HemeLB magic number."
			    << " Expected: " << unsigned(io::formats::HemeLbMagicNumber)
			    << " Actual: " << hlbMagicNumber;
	}

	// Check the value of the extraction file magic number.
	if (extMagicNumber != io::formats::extraction::MagicNumber)
        {
	  throw Exception() << "This file does not have the extraction magic number."
			    << " Expected: " << unsigned(io::formats::extraction::MagicNumber)
			    << " Actual: " << extMagicNumber;
	}

	// Check the version number.
	if (version != io::formats::extraction::VersionNumber)
	{
	  throw Exception() << "Version number incorrect."
			    << " Supported: " << unsigned(io::formats::extraction::VersionNumber)
			    << " Input: " << version;
	}

	{
	  // Obtain the size of voxel in metres.
	  double voxelSize;
	  preambleReader.read(voxelSize);

	  // Obtain the origin.
	  double origin[3];
	  preambleReader.read(origin[0]);
	  preambleReader.read(origin[1]);
	  preambleReader.read(origin[2]);
	}
	// Obtain the total number of sites, fields & header len
	uint64_t numberOfSites;
	uint32_t numberOfFields, lengthOfFieldHeader;
	preambleReader.read(numberOfSites);
	preambleReader.read(numberOfFields);
	preambleReader.read(lengthOfFieldHeader);

	if (numberOfFields != 1 )
	  throw Exception() << "Checkpoint file must contain exactly one field, the distributions, but has "
			    << numberOfFields;
	if (lengthOfFieldHeader != 32 )
	  throw Exception() << "Checkpoint file's field header must be 32 B long, but is "
			    <<  lengthOfFieldHeader;

	MpiFileReader fieldHeaderReader(inputFile, lengthOfFieldHeader);
	fieldHeaderReader.read(distField.name);
	fieldHeaderReader.read(distField.numberOfFloats);
	fieldHeaderReader.read(distField.offset);
	if (distField.name != "distributions")
	  throw Exception() << "Checkpoint file must contain field named 'distributions', but has '"
			    << distField.name << "'";

	if (distField.numberOfFloats != LatticeType::NUMVECTORS)
	  throw Exception() << "Checkpoint field distributions contains " << distField.numberOfFloats
			    << " distributions but this build of HemeLB requires " << LatticeType::NUMVECTORS;

      }
      comms.Broadcast(distField.offset, comms.GetIORank());
    }

    void LocalDistributionInput::ReadOffsets(const std::string& offsetFileName) {
      std::vector<uint64_t> offsets;

      // Only actually read on IO rank
      if (comms.OnIORank()) {
	namespace fmt = hemelb::io::formats;
	io::writers::xdr::XdrFileReader offsetReader(offsetFileName);
	uint32_t hlbMagicNumber, offMagicNumber, version, nRanks;
	offsetReader.read(hlbMagicNumber);
	offsetReader.read(offMagicNumber);
	offsetReader.read(version);
	offsetReader.read(nRanks);

	if (hlbMagicNumber != io::formats::HemeLbMagicNumber)
	  throw Exception() << "This file does not start with the HemeLB magic number."
			    << " Expected: " << unsigned(io::formats::HemeLbMagicNumber)
			    << " Actual: " << hlbMagicNumber;

	if (offMagicNumber != io::formats::offset::MagicNumber)
	  throw Exception() << "This file does not have the offset magic number."
			    << " Expected: " << unsigned(io::formats::offset::MagicNumber)
			    << " Actual: " << offMagicNumber;

	if (version != io::formats::offset::VersionNumber)
	  throw Exception() << "Version number incorrect."
			    << " Supported: " << unsigned(io::formats::offset::VersionNumber)
			    << " Input: " << version;

	if (nRanks != comms.Size())
	  throw Exception() << "Offset file has wrong number of MPI ranks."
			    << " Running with: " << comms.Size()
			    << " Input: " << nRanks;

	// Now read the encoded nProcs+1 values
	// We are going to duplicate these into a flattened array of shape (nRanks, 2)
	// [start0, end0, start1, end1, ...]
	// where end_i == start_i+1 (except for the start finish obvs)
	offsets.resize(2*nRanks);
	offsetReader.read(offsets[0]);
	for (unsigned i = 1; i < nRanks; ++i) {
	  offsetReader.read(offsets[2*i]);
	  offsets[2*i - 1] = offsets[2*i];
	}
	offsetReader.read(offsets[2*nRanks-1]);
	// Compute the total length of a record
	allCoresWriteLength = offsets[2*nRanks-1] - offsets[0];
      }
      // Now bcast/scatter from IO rank to all
      comms.Broadcast(allCoresWriteLength, comms.GetIORank());
      auto start_finish = comms.Scatter(offsets, 2, comms.GetIORank());
      localStart = start_finish[0];
      localStop = start_finish[1];
    }
  }
}
