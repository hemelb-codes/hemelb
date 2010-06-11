#include <stdlib.h>
#include <stdio.h>

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

