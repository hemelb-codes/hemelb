#include <rpc/xdr.h>
#include <string>

using namespace std;

// Class to read an Xdr-style file from disk.
class XdrReader
{
  private:
    XDR  myXdr;

  public:
    // Constructor and destructor.
    XdrReader(FILE* xdrFile);
    ~XdrReader();

    // Functions for reading the next bit of the file.
    void readDouble(double& outDouble);
    void readInt(int& outInt);
    void readUnsignedInt(unsigned int& outUInt);
};
