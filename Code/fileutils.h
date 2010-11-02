#ifndef HEMELB_FILEUTILS_H
#define HEMELB_FILEUTILS_H

#include <dirent.h>

// Define a suitable type for the system we're on
// (scandir has slightly different definitions on different
// architectures)
namespace hemelb
{
  namespace util
  {
#ifdef HEMELB_CFG_ON_BSD
    typedef struct direct direct_t;
#else // HEMELB_CFG_ON_BSD
    typedef const struct direct direct_t;
#endif // HEMELB_CFG_ON_BSD
    // Copy a file from the source path destination path to the destination
    // path.
    bool FileCopy(const char* iOriginalPath, const char* iNewPath);

    // Exits if the named file doesn't exist or can't be opened for
    // reading.
    void check_file(const char * filename);

    // Delete all files within a directory.
    int DeleteDirContents(char *pathname);

    // Check whether a directory exists.
    bool DoesDirectoryExist(const char *pathname);

    // Function to create the directory of given path, which user group
    // and anyone can read write and execute.
    void MakeDirAllRXW(const char* dirPath);

  }
}
#endif // HEMELB_FILEUTILS_H
