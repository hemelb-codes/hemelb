// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_RESOURCES_RESOURCE_H
#define HEMELB_RESOURCES_RESOURCE_H

#include <iostream>
#include <filesystem>

#include "resources/path_parameters.h"

namespace hemelb::resources
{
    /***
     * Define how to find resources paths for tests which use resources
     */
    class Resource
    {
    public:
        inline explicit Resource(const std::string &aResourceName) :
                resourceName(aResourceName)
        {
        }
        [[nodiscard]] inline std::filesystem::path Path() const
        {
            namespace fs = std::filesystem;
            if (fs::exists(BuildPath()))
            {
                return BuildPath();
            }
            if (fs::exists(InstallPath()))
            {
                return InstallPath();
            }
            std::cerr << "Resource " << resourceName << " not found either at: " << BuildPath()
                      << " or: " << InstallPath() << std::endl;
            return "";
        }
        [[nodiscard]] inline std::string BuildPath() const
        {
            return build_resource_path / resourceName;
        }
        [[nodiscard]] inline std::string InstallPath() const
        {
            return install_resource_path / resourceName;
        }
    private:
        std::string resourceName;
    };
}
#endif
