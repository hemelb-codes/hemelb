// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_REPORTING_REPORTER_H
#define HEMELB_REPORTING_REPORTER_H

#include <string>
#include <vector>
#include "configuration/SimConfig.h"
#include "configuration/CommandLine.h"
#include "log/Logger.h"
#include "util/fileutils.h"
#include "reporting/Policies.h"
#include "ctemplate/template.h"
#include "reporting/Reportable.h"
#include "resources/Resource.h"

namespace hemelb
{
  namespace reporting
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
        Reporter(const std::string &path, const std::string &inputFile);
        void Image(); //! Inform the reporter that an image has been saved.
        void Snapshot(); //! Inform the reporter that a simulation snapshot has been taken.

        void AddReportable(Reportable* reportable);

        void WriteXML()
        {
          Write(resources::Resource("report.xml.ctp").Path(), "report.xml");
        }
        void WriteTxt()
        {
          Write(resources::Resource("report.txt.ctp").Path(), "report.txt");
        }
        void FillDictionary();
        void Write()
        {
          WriteXML();
          WriteTxt();
        }
        ctemplate::TemplateDictionary const & GetDictionary()
        {
          return dictionary;
        }
      private:
        const std::string &path;
        void Write(const std::string &ctemplate, const std::string &as); //! Write the report to disk, (or wherever the WriterPolicy decides.)
        bool doIo; //! Is this the processor which should write the report.
        unsigned int snapshotCount; //! Number of snapshots taken.
        unsigned int imageCount; //! Number of images written.
        bool stability; //! Stability of the simulation.
        ctemplate::TemplateDictionary dictionary;
        std::vector<Reportable*> reportableObjects;
    };
  }
}

#endif // HEMELB_REPORTING_REPORTER_H
