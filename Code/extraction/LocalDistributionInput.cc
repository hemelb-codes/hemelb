#include "extraction/LocalDistributionInput.h"
#include "extraction/OutputField.h"
#include "extraction/LocalPropertyOutput.h"
#include "geometry/LatticeData.h"
#include "io/formats/formats.h"
#include "io/formats/extraction.h"

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

    void LocalDistributionInput::LoadDistribution(geometry::LatticeData* latDat)
    {
      // We could supply hints regarding how the file should be read but
      // We are not doing so yet.
      MPI_Info fileInfo;
      HEMELB_MPI_CALL(MPI_Info_create, (&fileInfo));

      // Open the file as read-only.
      // TO DO: raise an exception if the file does not exist.
      // TO DO: there seems to be a missing argument in the call to Open in LocalPropertyOutput.cc.
      inputFile = net::MpiFile::Open(comms, filePath, MPI_MODE_RDONLY, fileInfo);
      // Set the view to the file.
      inputFile.SetView(0, MPI_CHAR, MPI_CHAR, "native", fileInfo);

      // Now open the offset file.
      std::string offsetFileName;
      int32_t pathLength = filePath.length();
      offsetFileName = filePath.substr(0, pathLength-3) + "off";
      offsetFile = net::MpiFile::Open(comms, offsetFileName, MPI_MODE_RDONLY, fileInfo);
      offsetFile.SetView(0, MPI_CHAR, MPI_CHAR, "native", fileInfo);

      CheckPreamble();

      ReadHeaderInfo();

      // Read the data section.
      // Read just one timestep.
      if (comms.Rank() == HEADER_READING_RANK)
      {
        const unsigned timestepBytes = sizeof(uint64_t);

	std::vector<char> timestepBuffer(timestepBytes);
	inputFile.Read(timestepBuffer);

	// Create an XDR translator based on the read buffer.
	io::writers::xdr::XdrReader timestepReader = io::writers::xdr::XdrMemReader(&timestepBuffer[0],
										    timestepBytes);

	// Obtain the timestep.
	uint64_t timestep;
	timestepReader.readUnsignedLong(timestep);
      }

      if (!comms.OnIORank())
      {
        // Read the offset for this rank and the subsequent rank.
        const unsigned offsetBytes = 2*sizeof(uint64_t);
        std::vector<char> offsetBuffer(offsetBytes);
        offsetFile.ReadAt((comms.Rank()+1)*sizeof(uint64_t), offsetBuffer);
        io::writers::xdr::XdrReader offsetReader = io::writers::xdr::XdrMemReader(&offsetBuffer[0],
                                                                                  offsetBytes);
        uint64_t thisOffset, nextOffset;
        offsetReader.readUnsignedLong(thisOffset);
        offsetReader.readUnsignedLong(nextOffset);

        // Read the grid and distribution data.
        unsigned readLength = nextOffset - thisOffset;
        std::vector<char> dataBuffer(readLength);
        inputFile.ReadAt(thisOffset, dataBuffer);
        io::writers::xdr::XdrReader dataReader = io::writers::xdr::XdrMemReader(&dataBuffer[0],
                                                                                readLength);
	// TO DO: is this the best way to do this?
        uint32_t numberOfFloats = LocalPropertyOutput::GetFieldLength(hemelb::extraction::OutputField::Distributions);
        uint32_t lengthOfSegment = 3*sizeof(uint32_t) + numberOfFloats*sizeof(float);

        uint32_t numberOfLocalSites = 0;
        unsigned position = thisOffset;
        while (position < nextOffset)
        {
          uint32_t x, y, z;
          dataReader.readUnsignedInt(x);
          dataReader.readUnsignedInt(y);
          dataReader.readUnsignedInt(z);

	  distribn_t* f_old_p = latDat->GetFOld(numberOfLocalSites * LatticeType::NUMVECTORS);
	  distribn_t* f_new_p = latDat->GetFNew(numberOfLocalSites * LatticeType::NUMVECTORS);

	  float offset = LocalPropertyOutput::GetOffset(hemelb::extraction::OutputField::Distributions);
          for (int i = 0; i < numberOfFloats; i++)
          {
            float field_val;
            dataReader.readFloat(field_val);
            field_val += offset;
	    f_new_p[i] = f_old_p[i] = field_val;
          }

          position += lengthOfSegment;
          numberOfLocalSites++;
        }
      }
    }

    void LocalDistributionInput::CheckPreamble()
    {
      if (comms.Rank() == HEADER_READING_RANK)
      {
	const unsigned preambleBytes = 3*sizeof(uint32_t) + 4*sizeof(double) + sizeof(uint64_t);

	std::vector<char> preambleBuffer(preambleBytes);
	inputFile.Read(preambleBuffer);

	// Create an XDR translator based on the read buffer.
	io::writers::xdr::XdrReader preambleReader = io::writers::xdr::XdrMemReader(&preambleBuffer[0],
										    preambleBytes);

	// Read the magic numbers.
	uint32_t hlbMagicNumber, extMagicNumber, version;
	preambleReader.readUnsignedInt(hlbMagicNumber);
	preambleReader.readUnsignedInt(extMagicNumber);
	preambleReader.readUnsignedInt(version);

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

	// Obtain the size of voxel in metres.
	double voxelSize;
	preambleReader.readDouble(voxelSize);

	// Obtain the origin.
	double origin[3];
	preambleReader.readDouble(origin[0]);
	preambleReader.readDouble(origin[1]);
	preambleReader.readDouble(origin[2]);

	// Obtain the total number of sites.
	preambleReader.readUnsignedLong(numberOfSites);
      }
    }

    void LocalDistributionInput::ReadHeaderInfo()
    {
      if (comms.Rank() == HEADER_READING_RANK)
      {
	const unsigned infoBytes = 2*sizeof(uint32_t);

	std::vector<char> infoBuffer(infoBytes);
	inputFile.Read(infoBuffer);

	// Create an XDR translator based on the read buffer.
	io::writers::xdr::XdrReader infoReader = io::writers::xdr::XdrMemReader(&infoBuffer[0],
										infoBytes);
	uint32_t numberOfFields;
	infoReader.readUnsignedInt(numberOfFields);

	uint32_t lengthOfFieldHeader;
	infoReader.readUnsignedInt(lengthOfFieldHeader);

	std::vector<char> fieldHeaderBuffer(lengthOfFieldHeader);
	inputFile.Read(fieldHeaderBuffer);

	// Create an XDR translator based on the read buffer.
	io::writers::xdr::XdrReader fieldHeaderReader =
	  io::writers::xdr::XdrMemReader(&fieldHeaderBuffer[0], lengthOfFieldHeader);

	for (int i = 0; i < numberOfFields; i++)
	{
	  // When encoding a string XDR places an unsigned int at its head which gives the
	  // length of the string.
	  uint32_t lengthOfFieldName;
	  fieldHeaderReader.readUnsignedInt(lengthOfFieldName);
          uint32_t position = fieldHeaderReader.GetPosition();
	  std::string name;
	  for (int j = position; j < lengthOfFieldName + position; j++)
	  {
	    name += fieldHeaderBuffer[j];
	  }

	  // TO DO: This does not work as expected so change it.
	  if ((i == 1) & (name != "distributions"))
	  {
	    throw Exception() << "The first fields must be 'distributions'."
			      << " Actual: " << name;
	  }

	  uint32_t lengthOfPaddedFieldName = ((lengthOfFieldName +3)/4)*4;
	  fieldHeaderReader.SetPosition(fieldHeaderReader.GetPosition() + lengthOfPaddedFieldName);
	  uint32_t numberOfFloats;
	  fieldHeaderReader.readUnsignedInt(numberOfFloats);
	  double offset;
	  fieldHeaderReader.readDouble(offset);

	  extraction::InputField field;
	  field.name = name;
	  field.numberOfFloats = numberOfFloats;
	  field.offset = offset;

	  fields.push_back(field);
	}
      }
    }
  }
}
