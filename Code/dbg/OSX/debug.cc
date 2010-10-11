#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <unistd.h>

#include "mpiInclude.h"

#include "dbg/debug.h"
#include "dbg/OSX/debug.h"

std::string intToString(int i) {
  // Convert an int to a string.
  // Remember you will have to delete it!
  std::stringstream ss;
  ss << i;
  return ss.str();//return a string with the contents of the stream
}

void hemelb::dbg::attach(char *executable) {
  int amWaiting = 1;
  int rank; MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  int nProcs; MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  int *pIds = new int[nProcs];
  int pId = getpid();
  pid_t childPid = 0;
  
  MPI_Gather((void*)&pId, 1, MPI_INT,
	     (void*)pIds, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if (rank==0) {    
    //std::system(cmd.c_str());
    childPid = fork();
    
    if (childPid == 0) {
      // We're the child, so we are DOOMED.
      // (Either to exec() or if that fails exit().)
      
      std::string include (__FILE__);
      std::string osxDebugDir = include.substr(0, include.rfind('/'));
      std::string binaryPath (getcwd(NULL, 0)); // This will leak, but don't care
      binaryPath += "/";
      binaryPath += executable;
      
      typedef std::vector<std::string> VoS; // Vector of Strings
      VoS args;
      
      args.push_back(std::string("osascript"));
      args.push_back(osxDebugDir + "/MPIdebug.scpt");
      std::string gdbScript;
      
      {
	char* gdbScriptEnv = std::getenv("HEMELB_DEBUG_SCRIPT");
	if (gdbScriptEnv == NULL) {
	  gdbScript = osxDebugDir + "/resume.gdb";
	} else {
	  gdbScript = std::string(gdbScriptEnv);
	}
      }
      
      args.push_back(gdbScript);
      args.push_back(binaryPath);
      
      for (int i=0; i<nProcs; i++) {
	// This leaks memory
	args.push_back(intToString(pIds[i]));
      }
      
      args.push_back("");
      
      char **argv = new char *[args.size()];
      int i = 0;
      
      for(VoS::iterator it = args.begin(); it < args.end(); it++) {
	argv[i]  = new char[it->length()]; // for terminating null
	std::strcpy(argv[i], it->c_str());
	++i;
      }
      
      // Exec to replace hemelb with osascript
      execvp(argv[0], argv);
      // OK- that didn't work if we get here, better die, but print
      // command first
      i = 0;
      for(VoS::iterator it = args.begin(); it < args.end(); it++) {
	std::cerr << *it << " ";
      }
      std::cerr.flush();

      exit(1);
    } // if child
    
  } // if rank 0
  
  while (amWaiting) sleep(5);
  delete pIds;
  
  if (rank == 0) {
    int deadPid = waitpid(childPid, NULL, 0);
  }
  
}

