//#include <rpc/types.h>
#include <rpc/xdr.h>
#include <string>

using namespace std;

class XdrReader
{
  private:
    XDR  myXdr;

  public:
    XdrReader(FILE* xdrFile);
    ~XdrReader();

    void readDouble(double& outDouble);
    void readInt(int& outInt);
    void readUnsignedInt(unsigned int& outUInt);
};