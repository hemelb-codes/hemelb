#ifndef __xdrFileWriter_h_
#define __xdrFileWriter_h_

#include <stdio.h>

#include "xdrWriter.h"

// Class to write Xdr to a file. The actual write functions are implemented in the base class, XdrWriter.
class XdrFileWriter : public XdrWriter {
 private:
  FILE *myFile;
  
  // Implement the constructor and destructor to deal with the FILE and Xdr objects.
 public:
  XdrFileWriter(char* fileName);
  ~XdrFileWriter();
};

#endif //__xdrFileWriter_h_
