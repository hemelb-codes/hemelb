#include "xdrWriter.h"

// Class to write Xdr to a file. The actual write functions are implemented in the base class, XdrWriter.
class XdrFileWriter : public XdrWriter
{
  private:
    FILE *myFile;

  // Implement the constructor and destructor to deal with the FILE and Xdr objects.
  public:
    XdrFileWriter(char* fileName);
    ~XdrFileWriter();
};
