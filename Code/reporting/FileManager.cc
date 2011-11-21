#include "FileManager.h"
#include <sstream>
namespace hemelb
{
  namespace reporting
  {
    FileManager::FileManager(configuration::CommandLine &commandLine, const bool & io, const int & processorCount):
                                                            options(commandLine), ok(false),doIo(io)
    {

      inputFile=options.GetInputFile();
      outputDir=options.GetOutputDir();

      GuessOutputDir();

      imageDirectory = outputDir + "/Images/";
      snapshotDirectory = outputDir + "/Snapshots/";

      if (doIo)
      {
        if (hemelb::util::DoesDirectoryExist(outputDir.c_str()))
        {
          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("\nOutput directory \"%s\" already exists. Exiting.",
                                                                              outputDir.c_str());
          return;
        }


        hemelb::util::MakeDirAllRXW(outputDir);
        hemelb::util::MakeDirAllRXW(imageDirectory);
        hemelb::util::MakeDirAllRXW(snapshotDirectory);
        std::stringstream timings_name_stream;
        timings_name_stream<< outputDir << "/timings" << processorCount << ".asc" << std::flush;
        timings_name=timings_name_stream.str();
      }

      ok=true;
    }

    const std::string & FileManager::GetInputFile() const {
      return(inputFile);
    }
    const std::string & FileManager::GetSnapshotDirectory() const {
      return(snapshotDirectory);
    }
    const std::string & FileManager::GetImageDirectory() const {
      return(imageDirectory);
    }
    const std::string & FileManager::GetReportPath() const {
          return(timings_name);
    }

    void FileManager::EmptyOutputDirectories()
    {
      hemelb::util::DeleteDirContents(snapshotDirectory);
      hemelb::util::DeleteDirContents(imageDirectory);
    }



    hemelb::io::XdrFileWriter * FileManager::XdrImageWriter(const long int time){
      char filename[255];
      snprintf(filename, 255, "%08li.dat", time);
      return(new hemelb::io::XdrFileWriter(imageDirectory + std::string(filename)));
    }

    const std::string FileManager::SnapshotPath(unsigned long time) const {
      char snapshot_filename[255];
      snprintf(snapshot_filename, 255, "snapshot_%06li.dat", time);
      return(snapshotDirectory+std::string(snapshot_filename)); // by copy
    }

    void FileManager::SaveConfiguration(configuration::SimConfig * simConfig)
    {
      if (doIo) {
        simConfig->Save(outputDir + "/" + configLeafName);
      }
    }

    void FileManager::GuessOutputDir()
    {
      unsigned long lLastForwardSlash = inputFile.rfind('/');
      if (lLastForwardSlash == std::string::npos)
      {
        // input file supplied is in current folder
        configLeafName= inputFile;
        if (outputDir.length() == 0)
        {
          // no output dir given, defaulting to local.
          outputDir="./results";
        }
      } else {
        // input file supplied is a path to the input file
        configLeafName=  inputFile.substr(lLastForwardSlash);
        if (outputDir.length() == 0) {
          // no output dir given, defaulting to location of input file.
          outputDir=inputFile.substr(0, lLastForwardSlash)+"results";
        }
      }
    }
  }
}


