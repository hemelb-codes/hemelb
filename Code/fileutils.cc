#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <sys/dir.h>
#include <sys/stat.h>

#include "fileutils.h"

// Return true if file exists for reading, false if not.
bool FileUtils::file_exists(const char * filename) {
  if (FILE * file = fopen(filename, "r")) {
    fclose(file);
    return true;
  }
  return false;
}

// Check the existence of a critical file - exit if it's not there
void FileUtils::check_file(const char * filename) {
        if(!file_exists(filename)) {
                fprintf(stderr,"Cannot open file %s\nExiting.\n", filename);
                exit(0);
        }
}

// Function to select directory contents that are not "." or ".."
// The hack is necessary because of different versions of 'scandir'
// used by different compilers.
#ifdef DARWIN
int FileUtils::SelectOnlyContents (struct direct *entry)
#else
int FileUtils::SelectOnlyContents (const struct direct *entry)
#endif
{
  if ((strcmp(entry->d_name, ".") == 0) || (strcmp(entry->d_name, "..") == 0))
    {
      return 0;
    }
  else
    {
      return 1;
    }
}

// Delete all files within a directory.
int FileUtils::DeleteDirContents (char *pathname)
{
  struct direct **files;
  
  int file_count = scandir(pathname, &files, SelectOnlyContents, alphasort);
  
  char filename[1024];
  
  for (int i = 0; i < file_count; i++)
    {
      snprintf (filename, 1024, "%s/%s", pathname, files[i]->d_name);	
      unlink (filename);
    }
  return 0;
}

// Function to create the directory of given path, which user group and anyone
// can read write and execute.
void FileUtils::MakeDirAllRXW(char* dirPath)
{
  mkdir(dirPath, 0777);
}

