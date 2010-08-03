#include <fstream>

#include "asciiFileWriter.h"

// Implement a constructor that opens the file and creates the Xdr
// object to write to it.
AsciiFileWriter::AsciiFileWriter(char* fileName) {
  outStream = new std::ofstream(fileName);
}

// Destructor closes the file and cleans up the stream object.
AsciiFileWriter::~AsciiFileWriter() {
  
  std::ofstream* outFileStream = 
    static_cast<std::ofstream *>(outStream);
  
  outFileStream->close();
  delete outFileStream;
}
