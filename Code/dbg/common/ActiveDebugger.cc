#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <unistd.h>
#include <cerrno>
#include "mpiInclude.h"

#include "dbg/common/ActiveDebugger.h"

namespace hemelb
{
  
  namespace dbg
  {
    ActiveDebugger::ActiveDebugger(char* executable) : Debugger(executable),
						       mAmAttached(false),
						       mPIds(NULL) {}
    
    std::string ActiveDebugger::ConvertIntToString(int i) {
      // Convert an int to a string.
      // Remember you will have to delete it!
      std::stringstream ss;
      ss << i;
      return ss.str();
    }
    
    void ActiveDebugger::BreakHere(void) {
      return;
    }
    
    void ActiveDebugger::Attach(void) {
      // Start up a the debuggers, tell them to attach to the
      // processes and wait for them to attach. This function forks
      // another process on the rank 0 task and waits for it.
      if (mAmAttached) return;
      
      int amWaiting = 1;
      int rank; MPI_Comm_rank (MPI_COMM_WORLD, &rank);
      int nProcs; MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
      mPIds = new VoI(nProcs);
      int pId = getpid();
      pid_t childPid = 0;
      
      MPI_Gather((void*)&pId, 1, MPI_INT,
		 (void*)&(mPIds->front()), 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      if (rank==0) {    
	childPid = fork();
	
	if (childPid == 0) {
	  // This function won't return.
	  SpawnDebuggers();
	}
      } // if rank 0
      
      // Sit here waiting for GDB to attach
      while (amWaiting) sleep(5);
      
      if (rank == 0) {
	// Reap the spawner
	int deadPid = waitpid(childPid, NULL, 0);
      }
      
      mAmAttached = true;
    }
    
    void ActiveDebugger::GatherProcessIds(){
      if (mAmAttached) {
	// Return a vector of the process ids
	int rank; MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	int nProcs; MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	
	mPIds = new VoI(nProcs);
	
	int pId = getpid();
	
	MPI_Gather((void*)&pId, 1, MPI_INT,
		   (void*)&(mPIds->front()), 1, MPI_INT, 0, MPI_COMM_WORLD);
	
      }
    }
    
    void ActiveDebugger::SpawnDebuggers(void) {
      // Run the script to Tell the OS to start appropriate
      // terminal(s), launch the debuggers and attach them to our
      // processes.  We're the child, so we are DOOMED.  (Either to
      // exec() or if that fails exit().)
      std::string srcFile (__FILE__);
      std::string dbgCommonDir = srcFile.substr(0, srcFile.rfind('/'));
      
      std::string binaryPath (getcwd(NULL, 0)); // This will leak, 
                                                // but don't care
      binaryPath += "/";
      binaryPath += mExecutable;
      
      
      VoS args;
      
      args.push_back(GetPlatformInterpreter());
      
      args.push_back(GetPlatformScript());
      
      // Get the GDB script to exec
      // This will either be the value of the environment variable
      // HEMELB_DEBUG_SCRIPT or resume.gdb
      std::string gdbScript;
      {
	char* gdbScriptEnv = std::getenv("HEMELB_DEBUG_SCRIPT");
	if (gdbScriptEnv == NULL) {
	  gdbScript = dbgCommonDir + "/resume.gdb";
	} else {
	  gdbScript = std::string(gdbScriptEnv);
	}
      }
      
      args.push_back(gdbScript);
      args.push_back(binaryPath);
      
      for (VoI::iterator i = mPIds->begin(); i<mPIds->end(); ++i) {
	// This leaks memory
	args.push_back(ConvertIntToString(*i));
      }
      
      // +1 to include required NULL pointer for execvp()
      char **argv = new char *[args.size()+1];
      int i = 0;
      
      // convert to C array of char arrays.
      for(VoS::iterator it = args.begin(); it < args.end(); it++) {
	argv[i]  = new char[it->length()]; // for terminating null
	std::strcpy(argv[i], it->c_str());
	++i;
      }
      // Terminating NULL
      argv[i] = NULL;
      
      // Exec to replace hemelb with osascript
      int code = execvp(argv[0], argv);
      
      // OK- that didn't work if we get here, better die (since we're
      // the extra process). Print the error code too.
      std::cerr << "Couldn't exec() script to launch debuggers, error code "
		<< errno << std::endl;
      // Now print the command we wanted to exec()
      for(VoS::iterator it = args.begin(); it < args.end(); it++) {
	std::cerr << *it << " ";
      }
      std::cerr.flush();
      // Die
      exit(1);
    }

  }
}
