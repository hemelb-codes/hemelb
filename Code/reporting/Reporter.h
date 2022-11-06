// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REPORTING_REPORTER_H
#define HEMELB_REPORTING_REPORTER_H

#include <string>
#include <filesystem>
#include <vector>
#include "configuration/SimConfig.h"
#include "configuration/CommandLine.h"
#include "log/Logger.h"
#include "reporting/Policies.h"
#include "reporting/Reportable.h"
#include "resources/Resource.h"

namespace hemelb::reporting
{
    /**
     * Report generator class.
     * Class defining the creation of a report, intended for long-term archiving, describing what happened during a HemeLB run.
     * Accepts three policies as template arguments, defining:
     * @todo CommsPolicy and BroadcastPolicy should be unified
     */
    class Reporter
    {
      public:
        /**
         * Build a reporter capable of reporting on timings and other aspects of a HemeLB run.
         * @param path Path to write the report file to.
         * @param inputFile Input XML config file used to initialise the run.
         * @param aSiteCount Total count of sites used in the simulation.
         * @param timers Reference to list of timers used to measure performance.
         * @param aState Reference to state of ongoing simulation.
         */
        Reporter(const std::filesystem::path& report_dir, const std::string &inputFile);

        void AddReportable(Reportable* reportable);

        inline void WriteXML()
        {
          Write(resources::Resource("report.xml.ctp").Path(), "report.xml");
        }
        inline void WriteTxt()
        {
          Write(resources::Resource("report.txt.ctp").Path(), "report.txt");
        }
        void FillDictionary();
        inline void Write()
        {
          WriteXML();
          WriteTxt();
        }
        inline const Dict& GetDictionary() const
        {
          return dictionary;
        }
      private:
        std::filesystem::path report_dir;
        void Write(const std::filesystem::path& ctemplate, const std::string &as); //! Write the report to disk, (or wherever the WriterPolicy decides.)
        Dict dictionary;
        std::vector<Reportable*> reportableObjects;
    };

}

#endif // HEMELB_REPORTING_REPORTER_H
