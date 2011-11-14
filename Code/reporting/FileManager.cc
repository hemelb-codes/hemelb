#include "FileManager.h"
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

        InitialiseReport(processorCount);
      }
      ok=true;
    }

    FileManager::~FileManager()
    {
      if (doIo)
      {
        fclose(mTimingsFile);
      }
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


    void FileManager::EmptyOutputDirectories()
    {
      hemelb::util::DeleteDirContents(snapshotDirectory);
      hemelb::util::DeleteDirContents(imageDirectory);
    }

    void FileManager::InitialiseReport(const int & processorCount)
    {
      char timings_name[256];
      char procs_string[256];

      sprintf(procs_string, "%i", processorCount);
      strcpy(timings_name, outputDir.c_str());
      strcat(timings_name, "/timings");
      strcat(timings_name, procs_string);
      strcat(timings_name, ".asc");

      mTimingsFile = fopen(timings_name, "w");
      ReportHeader();
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

    void FileManager::ReportHeader()
    {
      fprintf(mTimingsFile, "***********************************************************\n");
      fprintf(mTimingsFile, "Opening config file:\n %s\n", inputFile.c_str());
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

    void FileManager::ReportProcessorTimings(std::string *const names,double *const mins,double *const means,double *const maxes){
       fprintf(ReportFile(),

                  "\n\nPer-proc timing data (secs per [cycle,cycle,cycle,image,snapshot]): \n\n");
          fprintf(ReportFile(), "\t\tMin \tMean \tMax\n");
          for (int ii = 0; ii < 5; ii++)
          {
            fprintf(ReportFile(),
                    "%s\t\t%.3g\t%.3g\t%.3g\n",
                    names[ii].c_str(),
                    mins[ii],
                    means[ii],
                    maxes[ii]);
          }
     }


    void FileManager::ReportPhase1(long int site_count, int total_time_steps, long int cycle_id,
                                   double simulation_time, bool unstable, unsigned long time_steps_per_cycle,
                                   double decomposition_time, double initialise_time, double read_time, double creation_time){
      fprintf(ReportFile(), "\n");
      fprintf(ReportFile(),
              "threads: %i, machines checked: %i\n\n",
              hemelb::topology::NetworkTopology::Instance()->GetProcessorCount(),
              hemelb::topology::NetworkTopology::Instance()->GetMachineCount());
      fprintf(ReportFile(),
              "topology depths checked: %i\n\n",
              hemelb::topology::NetworkTopology::Instance()->GetDepths());
      fprintf(ReportFile(), "fluid sites: %li\n\n", site_count);
      fprintf(ReportFile(),
              "cycles and total time steps: %li, %i \n\n",
              (cycle_id - 1), // Note that the cycle-id is 1-indexed.
              total_time_steps);
      fprintf(ReportFile(), "time steps per second: %.3f\n\n", total_time_steps / simulation_time);

      if (unstable)
      {
        fprintf(ReportFile(),
                "Attention: simulation unstable with %li timesteps/cycle\n", time_steps_per_cycle
        );
        fprintf(ReportFile(), "Simulation terminated\n");
      }

      fprintf(ReportFile(),
                    "time steps per cycle: %li\n",
                    time_steps_per_cycle);
            fprintf(ReportFile(), "\n");

            fprintf(ReportFile(), "\n");

            fprintf(ReportFile(), "\n");
            fprintf(ReportFile(), "decomposition optimisation time (s):       %.3f\n", decomposition_time);
            fprintf(ReportFile(), "pre-processing buffer management time (s): %.3f\n", initialise_time);
            fprintf(ReportFile(), "input configuration reading time (s):      %.3f\n", read_time);

            fprintf(ReportFile(),
                    "total time (s):                            %.3f\n\n",
                    (hemelb::util::myClock() - creation_time));

            fprintf(ReportFile(), "Sub-domains info:\n\n");

            for (hemelb::proc_t n = 0; n
            < hemelb::topology::NetworkTopology::Instance()->GetProcessorCount(); n++)
            {
              fprintf(ReportFile(),
                      "rank: %lu, fluid sites: %lu\n",
                      (unsigned long) n,
                      (unsigned long) hemelb::topology::NetworkTopology::Instance()->FluidSitesOnEachProcessor[n]);
            }
    }
  }
}


