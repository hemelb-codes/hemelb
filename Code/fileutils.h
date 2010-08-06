#ifndef __fileutils_h_
#define __fileutils_h_

#include <dirent.h>

// Define a suitable type for the system we're on
// (scandir has slightly different definitions on different
// architectures)
namespace util {
#ifdef DARWIN
  typedef struct direct direct_t;
#else
  typedef const struct direct direct_t;
#endif // DARWIN
  
  // Exits if the named file doesn't exist or can't be opened for
  // reading.
  void check_file(const char * filename);
  
  // Delete all files within a directory.
  int DeleteDirContents (char *pathname);
  
  // Function to create the directory of given path, which user group
  // and anyone can read write and execute.
  void MakeDirAllRXW(char* dirPath);
  
}

#endif // __fileutils_h_
