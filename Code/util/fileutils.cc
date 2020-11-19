// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <unistd.h>
#include <dirent.h>
#include <sys/dir.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>

#include "Exception.h"
#include "log/Logger.h"
#include "util/fileutils.h"

namespace hemelb
{
  namespace util
  {

    // Returns true if the file with the given name exists for reading,
    // false otherwise.
    //bool file_exists(const char * filename);

    // Function to select directory contents that are not "." or ".."
    // int selectOnlyContents (direct_t *entry);

    // Return true if file exists for reading, false if not.
    bool file_exists(const char * filename)
    {

      if (access(filename, R_OK) == -1)
      {
        return false;
      }
      return true;
    }

    // Check the existence of a critical file - exit if it's not there
    void check_file(const char * filename)
    {
      if (!file_exists(filename))
      {
        log::Logger::Log<log::Critical, log::OnePerCore>("Cannot open file %s\nExiting.", filename);
        std::exit(0);
      }
    }

    // Function to select directory contents that are not "." or ".."
    int selectOnlyContents(direct_t *entry)
    {
      if ( (std::strcmp(entry->d_name, ".") == 0) || (std::strcmp(entry->d_name, "..") == 0))
      {
        return 0;

      }
      else
      {
        return 1;
      }

    }

    void ChangeDirectory(const char * target)
    {
      chdir(target);
    }

    void ChangeDirectory(const std::string & target)
    {
      chdir(target.c_str());
    }

    void GetCurrentDir(char * result, int bufflength)
    {
      getcwd(result, bufflength);
    }

    std::string GetCurrentDir()
    {
      char buff[1000];
      GetCurrentDir(buff, 1000);
      return std::string(buff); // return by copy.
    }

    // This copied from BOOST. TODO: Use boost
    std::string GetTemporaryDir()
    {
      const char *dirname;
      dirname = std::getenv("HEME_TMP");
      if (nullptr == dirname)
        dirname = std::getenv("TMP");
      if (nullptr == dirname)
        dirname = std::getenv("TMPDIR");
      if (nullptr == dirname)
        dirname = std::getenv("TEMP");
      if (nullptr == dirname)
      {
        //assert(false); // no temp directory found
        return GetCurrentDir();
      }
      return std::string(dirname); // return by copy
    }

    // Delete all files within a directory.
    int DeleteDirContents(std::string pathname)
    {
      struct direct **files;

      int file_count = scandir(pathname.c_str(), &files, selectOnlyContents, alphasort);

      for (int i = 0; i < file_count; i++)
      {
        std::stringstream filename;
        filename << pathname.c_str() << "/" << files[i]->d_name << std::flush;
        unlink(filename.str().c_str());
      }
      return 0;
    }

    void DeleteDirTree(const std::string& pathname) {
      return DeleteDirTree(pathname.c_str());
    }

    namespace {
      // Helper for DeleteDirTree that wraps opendir/readdir/closedir in
      // RAII style.
      // 
      // Note that it also changes directory into its directory to make
      // the operations simpler, but changes back in d'tor.
      struct DirIter {
	char* startdir = nullptr;
	DIR* dstream = nullptr;

	DirIter() = default;
	DirIter(const char* name) {
	  dstream = opendir(name);
	  if (dstream == nullptr) {
	    throw Exception() << std::strerror(errno);
	  }
	  startdir = getcwd(NULL, 0); // this must be free()'d
	  chdir(name);
	}

	~DirIter() {
	  chdir(startdir);
	  std::free(startdir);
	  if (dstream)
	    closedir(dstream);
	}

	struct dirent* next() {
	  return readdir(dstream);
	}
      };
    }

    void DeleteDirTree(const char* pathname) {
      struct stat statbuf;
      if (stat(pathname, &statbuf) != 0)
	throw Exception() << std::strerror(errno);

      if (S_ISDIR(statbuf.st_mode)) {
	// is directory
	DirIter dir(pathname);
	for (struct dirent* entry = dir.next(); entry != nullptr; entry = dir.next()) {
	  if (entry->d_name[0] == '.' && (entry->d_name[1] == 0 || (entry->d_name[1] == '.' && entry->d_name[2] == 0)))
	    continue;
	  DeleteDirTree(entry->d_name);
	}
      }
      // this works for file or directory
      if (remove(pathname) != 0)
	throw Exception() << std::strerror(errno); 
    }

    // Detect whether a directory exists.
    bool DoesDirectoryExist(const char *pathname)
    {
      struct stat st;
      return stat(pathname, &st) == 0;
    }

    bool FileCopy(const char* iOriginalPath, const char* iNewPath)
    {
      std::ifstream lSource;
      std::ofstream lDestination;

      // open in binary to prevent jargon at the end of the buffer
      lSource.open(iOriginalPath, std::ios::binary);
      lDestination.open(iNewPath, std::ios::binary);

      if (!lSource.is_open() || !lDestination.is_open())
      {
        return false;
      }

      lDestination << lSource.rdbuf();

      lDestination.close();
      lSource.close();

      return true;
    }

    // Function to create the directory of given path, which user group and anyone
    // can read write and execute.
    bool MakeDirAllRXW(std::string const &dirPath)
    {
      int returnValue = mkdir(dirPath.c_str(), 0777);
      return (returnValue == 0);
    }

    std::string NormalizePathRelativeToPath(std::string inPath, std::string basePath)
    {
      // If it's an absolute path, just return it
      if (inPath[0] == '/')
      {
        return inPath;
      }

      // Going to check if it's a directory
      std::string baseDir;
      struct stat st;
      stat(basePath.c_str(), &st);
      // Assume it's a regular file in case it doesn't exist
      st.st_mode = S_IFREG;

      if (st.st_mode == S_IFDIR)
      {
        // It's a directory
        baseDir = basePath;
      }
      else
      {
        // Not a dir, find the last slash
        unsigned long lastSlash = basePath.rfind('/');
        if (lastSlash == basePath.npos)
        {
          // No slashes, so the baseDir is just the working dir
          baseDir = ".";
        }
        else
        {
          // Has slashes, return up to the last
          baseDir = basePath.substr(0, lastSlash);
        }
      }

      // Make sure it ends in a slash
      if (baseDir[baseDir.size() - 1] != '/')
      {
        baseDir += "/";
      }

      //Append the path of interest
      return baseDir + inPath;
    }

  }
}
