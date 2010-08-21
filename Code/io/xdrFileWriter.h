#ifndef HEME_IO_XDRFILEWRITER_H
#define HEME_IO_XDRFILEWRITER_H

#include <stdio.h>

#include "io/XdrWriter.h"

namespace heme
{
  namespace io 
  {
    // Class to write Xdr to a file. The actual write functions are implemented in the base class, XdrWriter.
    class XdrFileWriter : public XdrWriter {
  
      // Implement the constructor and destructor to deal with the FILE
      // and XDR objects.
    public:
      XdrFileWriter(char* fileName);
      ~XdrFileWriter();
  
    private:
      FILE *myFile;
  
    };
  }
}
#endif // HEME_IO_XDRFILEWRITER_H
