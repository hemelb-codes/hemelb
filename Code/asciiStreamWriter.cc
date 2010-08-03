#include <iostream>

#include "asciiStreamWriter.h"

using namespace std;

// Method to get the current position in the stream.
unsigned int AsciiStreamWriter::getCurrentStreamPosition() const {
  return outStream->tellp();
}

// No field/record separators in XDR files
void AsciiStreamWriter::writeFieldSeparator() {
  *outStream << ' ';
}
void AsciiStreamWriter::writeRecordSeparator() {
  *outStream << endl;
}

template <typename T>
void AsciiStreamWriter::_write(T const & value) {
  *outStream << value;
}
