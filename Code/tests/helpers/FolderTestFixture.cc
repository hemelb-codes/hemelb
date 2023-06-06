// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/FolderTestFixture.h"

#include <sstream>
#include <fstream>
#include <iomanip>
#include <filesystem>

#include <unistd.h>

#include <catch2/catch.hpp>
#include <tinyxml.h>

#include "resources/Resource.h"
#include "util/clock.h"

namespace fs = std::filesystem;

namespace hemelb::tests::helpers
{
    void ModifyXMLInput(TiXmlDocument &document, std::vector<std::string> const& elements,
                        std::string const& _value)
    {
        std::string const& attribute = elements.back();
        // Point to the *actual last elem*
        auto end = --elements.end();
        auto child = document.FirstChildElement("hemelbsettings");
        for (auto iter = elements.begin(); iter != end; ++iter) {
            auto& name = *iter;
            auto next_child = child->FirstChildElement(name);
            if(next_child  == nullptr) {
                next_child = new TiXmlElement(name);
                child->LinkEndChild(next_child);
            }
            child = next_child;
        }
        child->SetAttribute(attribute, _value);
    }

    //! \brief Modify XML document by deleting an element if it exists
    //! \details HemeLB parameters cannot be modified programmatically, so we have to jump
    //! these hoops to test it.
    //! \param[in] document: Document to modify
    //! \param[in] elements: hierarchy of elements, last item will be removed.
    //!   Should not include "hemelbsettings"
    //! \param[in] value: Value to set the attribute to
    void DeleteXMLInput(TiXmlDocument &document, std::vector<std::string> const& elements)
    {
        auto element = document.FirstChildElement("hemelbsettings");
        for (std::string const &name : elements)
        {
            auto next_child = element->FirstChildElement(name);
            if(next_child  == nullptr)
            {
                return;
            }
            element = next_child;
        }
        element->Parent()->RemoveChild(element);
    }

    FolderTestFixture::FolderTestFixture()
    {
        auto comm = Comms();
        // We will do the MPI comms using std::string instead of std::filesystem::path
        std::string temp_string;
        if (comm.OnIORank()) {
            tempPath = ConstructTempPath();
            // create a folder to work in
            fs::create_directory(tempPath);
            ownsTempPath = true;
            temp_string = tempPath;
        }
        comm.Broadcast(temp_string, comm.GetIORank());
        if (!comm.OnIORank()) {
            tempPath = temp_string;
        }

        // store current location
        origin = fs::current_path();
        MoveToTempdir();
    }

    FolderTestFixture::~FolderTestFixture() {
        ReturnToOrigin();
        if (ownsTempPath) {
            fs::remove_all(tempPath);
        }
    }

    namespace {
        std::int64_t HackyUUID() {
            static std::int32_t uuid[2] = {getpid(), std::int32_t(util::clock() * 100000)};
            ++uuid[1];
            return *reinterpret_cast<std::int64_t*>(uuid);
        }
    }

    auto FolderTestFixture::ConstructTempPath() const -> path {
        std::stringstream tmp_name;
        tmp_name << "HemeLBTest-" << std::fixed << HackyUUID();
        // TODO: find a portable uuid solution. BOOST?
        return fs::temp_directory_path() / std::move(tmp_name).str();
    }

    void FolderTestFixture::ReturnToOrigin()
    {
        // return to origin
        fs::current_path(origin);
    }

    void FolderTestFixture::CopyResourceToTempdir(const std::string & resource) const
    {
        bool ok = fs::copy_file(resources::Resource(resource).Path(),
                                tempPath / resource);
        REQUIRE(ok);
    }

    void FolderTestFixture::ModifyXMLInput(std::string const &resource, std::vector<std::string> const& elements,
                                           std::string const &_value) const
    {
        auto const filename = tempPath / resource;
        TiXmlDocument document(filename.c_str());
        document.LoadFile();
        helpers::ModifyXMLInput(document, elements, _value);
        std::ofstream output(filename);
        [&](std::ostream& o, TiXmlDocument const& doc) {
            o << doc;
        } (output, document);
    }

    void FolderTestFixture::DeleteXMLInput(std::string const &resource, std::vector<std::string> const& elements) const
    {
        std::string const filename = tempPath / resource;
        TiXmlDocument document(filename.c_str());
        document.LoadFile();
        helpers::DeleteXMLInput(document, elements);
        std::ofstream output(filename);
        output << document;
    }

    void FolderTestFixture::MoveToTempdir()
    {
        fs::current_path(tempPath);
    }

    void FolderTestFixture::AssertPresent(const std::string &fname)
    {
        // "does directory exist" actually works for files too.
        REQUIRE(fs::exists(fname));
    }
    const fs::path& FolderTestFixture::GetTempdir()
    {
        return tempPath;
    }
}

