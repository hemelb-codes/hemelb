// In extraction/XdrPropertyOutput.h

#pragma once

#include "extraction/LocalPropertyOutput.h"
#include "net/MpiFile.h"

namespace hemelb::extraction
{
  class XdrPropertyOutput : public LocalPropertyOutput
  {
  public:
    XdrPropertyOutput(IterableDataSource& dataSource,
                      const PropertyOutputFile& outputSpec,
                      const net::IOCommunicator& ioComms);
    ~XdrPropertyOutput() = default;

    void Write(unsigned long timestepNumber, unsigned long totalSteps) override;

  private:
  // from LocalPropertyOutput.h private

      // How many bytes are written for a single site?
      std::uint64_t CalcSiteWriteLen(std::vector<OutputField> const& fields) const;

      // Make the XTR header
      std::vector<char> PrepareHeader() const;

      // Open the file specified and write the header. Collective.
      void StartFile(std::string const& fn);

      // Write the offset file. Collective on the communicator.
      // for XDR private
      void WriteOffsetFile();

      // For single-timestep-per-file mode, hold the pattern we'll
      // pass to printf.
      std::string output_file_pattern;
      // The MPI file to write into.
      net::MpiFile outputFile;

      // The length, in bytes, of the whole header
      std::uint64_t header_length;
      // The data that makes up the header (only used on rank 0)
      std::vector<char> header_data;

      // The length, in bytes, of the local/global data write for one timestep
      std::uint64_t local_data_write_length;
      std::uint64_t global_data_write_length;

      // Where, in bytes, to begin writing into the file.
      std::uint64_t local_write_start;

      // Buffer to serialise into before writing to disk.
      std::vector<char> buffer;

      // The MPI file to write the offsets into.
      std::string offset_file_name;
  };
} // namespace hemelb::extraction
