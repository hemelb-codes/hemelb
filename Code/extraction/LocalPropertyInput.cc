#include "extraction/LocalPropertyInput.h"
#include "io/formats/formats.h"
#include "io/formats/extraction.h"
#include "io/writers/xdr/XdrMemReader.h"

namespace hemelb
{
  namespace extraction
  {
    LocalPropertyInput::LocalPropertyInput(const std::string dataFilePath,
                                           const net::IOCommunicator& ioComm) :
      comms(ioComm), filePath(dataFilePath)
    {
    }

    LocalPropertyInput::~LocalPropertyInput()
    {
    }

    void LocalPropertyInput::LoadDistribution()
    {
      // We could supply hints regarding how the file should be read but
      // We are not doing so yet.
      MPI_Info fileInfo;
      HEMELB_MPI_CALL(MPI_Info_create, (&fileInfo));

      // Open the file as read-only.
      // TO DO: raise an exception if the file does not exist.
      // TO DO: there seems to be a missing argument in the call to Open in LocalPropertyOutput.cc.
      inputFile = net::MpiFile::Open(comms, filePath, MPI_MODE_RDONLY, fileInfo);
      //inputFile = net::MpiFile::Open(comms, inputSpec->filename, MPI_MODE_RDONLY, fileInfo);
      //net::MpiFile inputFile;
      // Set the view to the file.
      inputFile.SetView(0, MPI_CHAR, MPI_CHAR, "native", fileInfo);

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
	std::cout << "timestep = " << timestep << std::endl;

	//for (uint64_t s = 0; s < numberOfSites; s++)
	for (uint64_t s = 0; s < 3; s++)
        {
	  // Deal with the grid position first.

	  const unsigned gridBytes = 3*sizeof(uint32_t);

	  std::vector<char> gridBuffer(gridBytes);
	  inputFile.Read(gridBuffer);

	  // Obtain the grid position.
	  uint32_t x, y, z;
	  io::writers::xdr::XdrReader dataReader = io::writers::xdr::XdrMemReader(&gridBuffer[0],
										  gridBytes);
	  dataReader.readUnsignedInt(x);
	  dataReader.readUnsignedInt(y);
	  dataReader.readUnsignedInt(z);

	  std::cout << "gridPosition = (" << x << ", " << y << ", " <<  z << ")" << std::endl;

	  // Obtain the fields.
	  for (auto &f : fields)
	  {
	    std::cout << "Reading " << f.name << std::endl;

	    uint32_t fieldBytes = f.numberOfFloats * sizeof(float);
	    //std::cout << fieldBytes << std::endl;

	    std::vector<char> fieldBuffer(fieldBytes);
	    inputFile.Read(fieldBuffer);

	    io::writers::xdr::XdrReader fieldReader = io::writers::xdr::XdrMemReader(&fieldBuffer[0],
										     fieldBytes);

	    for (int i = 0; i < f.numberOfFloats; i++)
	    {
	      float field_val;
	      fieldReader.readFloat(field_val);
	      //field_val -= f.offset;
	      std::cout << i << ": " << field_val << std::endl;
	    }
	  }
	}
      }
    }

    void LocalPropertyInput::CheckPreamble()
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

    void LocalPropertyInput::ReadHeaderInfo()
    {
      if (comms.Rank() == HEADER_READING_RANK)
      {
	const unsigned infoBytes = 2*sizeof(uint32_t);

	std::vector<char> infoBuffer(infoBytes);
	inputFile.Read(infoBuffer);

	// Create an XDR translator based on the read buffer.
	io::writers::xdr::XdrReader infoReader = io::writers::xdr::XdrMemReader(&infoBuffer[0],
										infoBytes);
	uint32_t numberOfFields, lengthOfFieldHeader;
	infoReader.readUnsignedInt(numberOfFields);
	std::cout << "nFields " << numberOfFields << std::endl;
	infoReader.readUnsignedInt(lengthOfFieldHeader);
	std::cout << "lHeader " << lengthOfFieldHeader << std::endl;

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
	  std::cout << name << std::endl;
	  uint32_t lengthOfPaddedFieldName = ((lengthOfFieldName +3)/4)*4;
	  fieldHeaderReader.SetPosition(fieldHeaderReader.GetPosition() + lengthOfPaddedFieldName);
	  uint32_t numberOfFloats;
	  fieldHeaderReader.readUnsignedInt(numberOfFloats);
	  std::cout << "number of floats: " << numberOfFloats << std::endl;
	  double offset;
	  fieldHeaderReader.readDouble(offset);
	  std::cout << "offset: " << offset << std::endl;

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
