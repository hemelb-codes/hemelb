#ifndef __xdrReader_h_
#define __xdrReader_h_

#include <stdio.h>

#include <rpc/types.h>
#include <rpc/xdr.h>

// Class to read an Xdr-style file from disk.
class XdrReader {
  public:
    // Constructor and destructor.
    XdrReader(FILE* xdrFile);
    ~XdrReader();

    // Functions for reading the next bit of the file.
    void readDouble(double& outDouble);
    void readInt(int& outInt);
    void readUnsignedInt(unsigned int& outUInt);
    
  private:
    XDR  myXdr;

};

#endif// __xdrReader_h_
