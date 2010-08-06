#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>

#include <sys/dir.h>
#include <sys/stat.h>

#include "fileutils.h"

namespace util {
  namespace {
    // Returns true if the file with the given name exists for reading,
    // false otherwise.
    //bool file_exists(const char * filename);
    
    // Function to select directory contents that are not "." or ".."
    // int selectOnlyContents (direct_t *entry);
    
    
    // Return true if file exists for reading, false if not.
    bool file_exists(const char * filename) {
      
      if (access(filename, R_OK) == -1) {
	return false;
      }
      return true;
    }
    
  }
  
  // Check the existence of a critical file - exit if it's not there
  void check_file(const char * filename) {
    if(!file_exists(filename)) {
      fprintf(stderr,"Cannot open file %s\nExiting.\n", filename);
      exit(0);
    }
  }

  // Function to select directory contents that are not "." or ".."
  int selectOnlyContents (direct_t *entry)
  {
    if ((strcmp(entry->d_name, ".") == 0) ||
	(strcmp(entry->d_name, "..") == 0)) {
      return 0;
      
    } else {
      return 1;
    }
    
  }

  // Delete all files within a directory.
  int DeleteDirContents (char *pathname)
  {
    struct direct **files;
  
    int file_count = scandir(pathname, &files,
			     selectOnlyContents, alphasort);
  
    char filename[1024];
    
    for (int i = 0; i < file_count; i++) {
      snprintf(filename, 1024,
	       "%s/%s", pathname, files[i]->d_name);
      unlink (filename);
    }
    return 0;
  }

  // Function to create the directory of given path, which user group and anyone
  // can read write and execute.
  void MakeDirAllRXW(char* dirPath)
  {
    mkdir(dirPath, 0777);
  }

}
