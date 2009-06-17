#include <stdlib.h>
#include <stdio.h>

#include "fileutils.h"

int file_exists(const char * filename) {
	if (FILE * file = fopen(filename, "r")) {
                fclose(file);
                return 0;
    }
    return -1;
}

void check_file(const char * filename) {
        if(file_exists(filename) < 0 ) {
                fprintf(stderr,"Cannot open file %s\nExiting.\n", filename);
                exit(0);
        } else {
               // fprintf(stderr,"Located file %s\n", filename);
        }
}

