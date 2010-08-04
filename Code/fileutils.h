#ifndef __fileutils_h_
#define __fileutils_h_

#include <dirent.h>

// Define a suitable type for the system we're on
// (scandir has slightly different definitions on different
// architectures)
#ifdef DARWIN
typedef struct direct direct_t;
#else
typedef const struct direct direct_t;
#endif // DARWIN

class FileUtils
{
  private:
    // Returns true if the file with the given name exists for reading, false otherwise.
    static bool file_exists(const char * filename);

    // Function to select directory contents that are not "." or ".."
    static int SelectOnlyContents (direct_t *entry);
    
  public:
    // Exits if the named file doesn't exist or can't be opened for reading.
    static void check_file(const char * filename);

    // Delete all files within a directory.
    static int DeleteDirContents (char *pathname);

    // Function to create the directory of given path, which user group and anyone
    // can read write and execute.
    static void MakeDirAllRXW(char* dirPath);
};

#endif // __fileutils_h_
