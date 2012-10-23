// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_RESOURCES_RESOURCE_H
#define HEMELB_RESOURCES_RESOURCE_H

#include "util/fileutils.h"
#include "resources/path_parameters.h"

namespace hemelb
{
  namespace resources
  {
    /***
     * Define how to find resources paths for tests which use resources
     */
    class Resource
    {
      public:
        Resource(const std::string &aResourceName) :
            resourceName(aResourceName)
        {
        }
        std::string Path() const
        {
          if (util::file_exists(BuildPath().c_str()))
          {
            return BuildPath();
          }
          if (util::file_exists(InstallPath().c_str()))
          {
            return InstallPath();
          }
          std::cerr << "Resource " << resourceName  << " not found either at: " << BuildPath() << " or: " << InstallPath() << std::endl;
          return "";
        }
        std::string BuildPath() const
        {
          return build_resource_path + "/" + resourceName;
        }
        std::string InstallPath() const
        {
          return install_resource_path + "/" + resourceName;
        }
      private:
        std::string resourceName;
    };
  }
}
#endif // ONCE
