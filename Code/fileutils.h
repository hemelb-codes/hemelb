#ifndef HEME_FILEUTILS_H
#define HEME_FILEUTILS_H

#include <dirent.h>

// Define a suitable type for the system we're on
// (scandir has slightly different definitions on different
// architectures)
namespace util {
#ifdef HEME_CFG_ON_BSD
  typedef struct direct direct_t;
#else // HEME_CFG_ON_BSD
  typedef const struct direct direct_t;
#endif // HEME_CFG_ON_BSD
  
  // Exits if the named file doesn't exist or can't be opened for
  // reading.
  void check_file(const char * filename);
  
  // Delete all files within a directory.
  int DeleteDirContents (char *pathname);
  
  // Function to create the directory of given path, which user group
  // and anyone can read write and execute.
  void MakeDirAllRXW(char* dirPath);
  
}

#endif // HEME_FILEUTILS_H
