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
      // We could supply hints regarding how the file should be read
      // but we are not doing so yet.

      // Open the file as read-only.
      // TO DO: raise an exception if the file does not exist.
      inputFile = net::MpiFile::Open(comms, filePath, MPI_MODE_RDONLY);
      // Set the view to the file.
      inputFile.SetView(0, MPI_CHAR, MPI_CHAR, "native");

      // Now open the offset file.
      std::string offsetFileName;
      int32_t pathLength = filePath.length();
      offsetFileName = filePath.substr(0, pathLength-3) + "off";
      offsetFile = net::MpiFile::Open(comms, offsetFileName, MPI_MODE_RDONLY);
      offsetFile.SetView(0, MPI_CHAR, MPI_CHAR, "native");

      CheckPreamble();

      ReadHeaderInfo();

      // Read the data section.
      // Read just one timestep.
      if (comms.OnIORank())
      {
        const unsigned timestepBytes = sizeof(uint64_t);

	std::vector<char> timestepBuffer(timestepBytes);
	inputFile.Read(timestepBuffer);

	// Create an XDR translator based on the read buffer.
	io::writers::xdr::XdrMemReader timestepReader(&timestepBuffer[0], timestepBytes);

	// Obtain the timestep.
	uint64_t timestep;
	timestepReader.read(timestep);
      }
      else
      {
        // Read the offset for this rank and the subsequent rank.
        const unsigned offsetBytes = 2*sizeof(uint64_t);
        std::vector<char> offsetBuffer(offsetBytes);
        offsetFile.ReadAt((comms.Rank()+1)*sizeof(uint64_t), offsetBuffer);
	io::writers::xdr::XdrMemReader offsetReader(&offsetBuffer[0], offsetBytes);
        uint64_t thisOffset, nextOffset;
        offsetReader.read(thisOffset);
        offsetReader.read(nextOffset);

        // Read the grid and distribution data.
        unsigned readLength = nextOffset - thisOffset;
        std::vector<char> dataBuffer(readLength);
        inputFile.ReadAt(thisOffset, dataBuffer);
	io::writers::xdr::XdrMemReader dataReader(&dataBuffer[0], readLength);
	// TO DO: is this the best way to do this?
        uint32_t numberOfFloats = LocalPropertyOutput::GetFieldLength(hemelb::extraction::OutputField::Distributions);
        uint32_t lengthOfSegment = 3*sizeof(uint32_t) + numberOfFloats*sizeof(float);

        uint32_t numberOfLocalSites = 0;
        unsigned position = thisOffset;
        while (position < nextOffset)
        {
          uint32_t x, y, z;
          dataReader.read(x);
          dataReader.read(y);
          dataReader.read(z);

	  distribn_t* f_old_p = latDat->GetFOld(numberOfLocalSites * LatticeType::NUMVECTORS);
	  distribn_t* f_new_p = latDat->GetFNew(numberOfLocalSites * LatticeType::NUMVECTORS);

	  float offset = LocalPropertyOutput::GetOffset(hemelb::extraction::OutputField::Distributions);
          for (int i = 0; i < numberOfFloats; i++)
          {
            float field_val;
            dataReader.read(field_val);
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
      if (comms.OnIORank())
      {
	const unsigned preambleBytes = 3*sizeof(uint32_t) + 4*sizeof(double) + sizeof(uint64_t);

	std::vector<char> preambleBuffer(preambleBytes);
	inputFile.Read(preambleBuffer);

	// Create an XDR translator based on the read buffer.
	io::writers::xdr::XdrMemReader preambleReader(&preambleBuffer[0], preambleBytes);

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

	// Obtain the size of voxel in metres.
	double voxelSize;
	preambleReader.read(voxelSize);

	// Obtain the origin.
	double origin[3];
	preambleReader.read(origin[0]);
	preambleReader.read(origin[1]);
	preambleReader.read(origin[2]);

	// Obtain the total number of sites.
	preambleReader.read(numberOfSites);
      }
    }

    void LocalDistributionInput::ReadHeaderInfo()
    {
      if (comms.OnIORank())
      {
	const unsigned infoBytes = 2*sizeof(uint32_t);

	std::vector<char> infoBuffer(infoBytes);
	inputFile.Read(infoBuffer);

	// Create an XDR translator based on the read buffer.
	io::writers::xdr::XdrMemReader infoReader(&infoBuffer[0], infoBytes);
	uint32_t numberOfFields;
	infoReader.read(numberOfFields);

	uint32_t lengthOfFieldHeader;
	infoReader.read(lengthOfFieldHeader);

	std::vector<char> fieldHeaderBuffer(lengthOfFieldHeader);
	inputFile.Read(fieldHeaderBuffer);

	// Create an XDR translator based on the read buffer.
	io::writers::xdr::XdrMemReader fieldHeaderReader(&fieldHeaderBuffer[0], lengthOfFieldHeader);

	for (int i = 0; i < numberOfFields; i++)
	{
	  // When encoding a string XDR places an unsigned int at its head which gives the
	  // length of the string.
	  uint32_t lengthOfFieldName;
	  fieldHeaderReader.read(lengthOfFieldName);
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
	  fieldHeaderReader.read(numberOfFloats);
	  double offset;
	  fieldHeaderReader.read(offset);

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
