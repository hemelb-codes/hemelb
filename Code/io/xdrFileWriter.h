#ifndef __xdrFileWriter_h_
#define __xdrFileWriter_h_

#include <stdio.h>

#include "xdrWriter.h"

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

#endif //__xdrFileWriter_h_
