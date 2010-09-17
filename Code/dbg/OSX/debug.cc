#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <unistd.h>

#include "mpiInclude.h"

#include "dbg/debug.h"
#include "dbg/OSX/debug.h"

void hemelb::dbg::attach(char *executable) {
  int amWaiting = 1;
  int rank; MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  int nProcs; MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  int *pIds = new int[nProcs];
  int pId = getpid();
  
  MPI_Gather((void*)&pId, 1, MPI_INT,
	     (void*)pIds, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (rank==0) {
    std::ostringstream oStream;
    std::string include (__FILE__);
    std::string osxDebugDir = include.substr(0, include.rfind('/'));
    char *cwd = getcwd(NULL, 0);
    std::string binaryPath (cwd);
    binaryPath += "/";
    binaryPath += executable;
    std::free(cwd);
    
    //("/Users/rupert/working/hemelb-clean/Code/build/hemelb");
    
    oStream << "osascript "
	    << osxDebugDir << "/MPIdebug.scpt "
	    << osxDebugDir << "/resume.gdb "
	    << binaryPath;
    
    for (rank=0; rank<nProcs; rank++) {
      oStream << ' ' << pIds[rank];
    }
    std::string cmd = oStream.str();
    
    std::cerr << cmd << std::endl;
    std::cerr.flush();
    
    std::system(cmd.c_str());
    
  }
  while (amWaiting) sleep(5);
  delete pIds;
}

