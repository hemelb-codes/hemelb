#include "io/h5md/H5MD.h"
#include "io/h5md/H5MDException.h"

namespace hemelb
{
  namespace io
  {
    namespace h5md
    {

      void AttachAttribute(const H5::H5Location & location, const std::string & name,
                           const std::string & value)
      {
        const H5::DataSpace space;
        const H5::Attribute attribute = location.createAttribute(name.c_str(),
                                                                 H5::PredType::C_S1,
                                                                 space);
        attribute.write(H5::PredType::C_S1, value.c_str());
      }

      H5MD H5MD::Create(const H5::CommonFG & location, const std::string & creator_name,
                        const std::string & creator_version, const std::string & author_name,
                        const std::string & author_email)
      {
        const H5::Group h5md = location.createGroup("h5md");

        const hsize_t dims[] = { 2 };
        const H5::DataSpace space(1, dims);
        const H5::Attribute attribute = h5md.createAttribute("version",
                                                             H5::PredType::NATIVE_INT,
                                                             space);
        attribute.write(H5::PredType::NATIVE_INT, VERSION);

        const H5::Group creator = h5md.createGroup("creator");
        AttachAttribute(creator, "name", creator_name);
        AttachAttribute(creator, "version", creator_version);

        const H5::Group author = h5md.createGroup("author");
        AttachAttribute(author, "name", author_name);
        if (author_email != "")
          AttachAttribute(author, "email", author_email);

        return H5MD(location, creator_name, creator_version, author_name, author_email);
      }

      H5MD H5MD::Open(const H5::CommonFG & location)
      {
        const H5::Group h5md = location.openGroup("h5md");

        int version[2];
        h5md.openAttribute("version").read(H5::PredType::C_S1, version);
        if (version[0] != VERSION[0] || version[1] > VERSION[1])
          throw H5MDException("H5MD file is incompatible version");

        std::string creator_name, creator_version;
        const H5::Group creator = h5md.openGroup("creator");
        creator.openAttribute("name").read(H5::PredType::C_S1, creator_name);
        creator.openAttribute("version").read(H5::PredType::C_S1, creator_version);

        std::string author_name, author_email;
        const H5::Group author = h5md.openGroup("author");
        author.openAttribute("name").read(H5::PredType::C_S1, author_name);
        if (author.attrExists("email"))
          author.openAttribute("email").read(H5::PredType::C_S1, author_email);

        return H5MD(location, creator_name, creator_version, author_name, author_email);
      }

      bool H5MD::IsH5MD(const H5::CommonFG & location)
      {
        try
        {
          const H5::Group h5md = location.openGroup("h5md");
          if (h5md.attrExists("version"))
          {
            const H5::Attribute version = h5md.openAttribute("version");
            if (version.getTypeClass() != H5::PredType::NATIVE_INT ||
                version.getSpace().getSimpleExtentNdims() != 1)
            {
              return false;
            }
            hsize_t dims[1];
            version.getSpace().getSimpleExtentDims(dims);
            if (dims[0] != 2)
              return false;
          }

          const H5::Group creator = h5md.openGroup("creator");
          if (!creator.attrExists("name") || !creator.openAttribute("name").getSpace().isSimple()
              || creator.openAttribute("name").getTypeClass() != H5::PredType::C_S1)
          {
            return false;
          }
          if (!creator.attrExists("version")
              || !creator.openAttribute("version").getSpace().isSimple()
              || creator.openAttribute("version").getTypeClass() != H5::PredType::C_S1)
          {
            return false;
          }

          const H5::Group author = h5md.openGroup("author");
          if (!author.attrExists("name") || !author.openAttribute("name").getSpace().isSimple()
              || author.openAttribute("name").getTypeClass() != H5::PredType::C_S1)
          {
            return false;
          }
          if (author.attrExists("email"))
          {
            if (!author.openAttribute("email").getSpace().isSimple()
                || author.openAttribute("email").getTypeClass() != H5::PredType::C_S1)
            {
              return false;
            }
          }
        }
        catch (H5::Exception & e)
        {
          return false;
        }
        return true;
      }

      void H5MD::SetCreator(const std::string & c)
      {
        if (creator != c)
        {
          const H5::Group h5md = location->openGroup("h5md");
          const H5::Group group = h5md.openGroup("creator");
          H5::Attribute attribute = group.openAttribute("name");
          attribute.write(H5::PredType::C_S1, c.c_str());
          creator = c;
        }
      }

      void H5MD::SetCreatorVersion(const std::string & v)
      {
        if (version != v)
        {
          const H5::Group h5md = location->openGroup("h5md");
          const H5::Group group = h5md.openGroup("creator");
          H5::Attribute attribute = group.openAttribute("version");
          attribute.write(H5::PredType::C_S1, v.c_str());
          version = v;
        }
      }

      void H5MD::SetAuthor(const std::string & a)
      {
        if (author != a)
        {
          const H5::Group h5md = location->openGroup("h5md");
          const H5::Group group = h5md.openGroup("author");
          H5::Attribute attribute = group.openAttribute("name");
          attribute.write(H5::PredType::C_S1, a.c_str());
          author = a;
        }
      }

      void H5MD::SetAuthorEmail(const std::string & e)
      {
        if (email != e)
        {
          const H5::Group h5md = location->openGroup("h5md");
          const H5::Group group = h5md.openGroup("author");
          if (group.attrExists("email"))
          {
            if (e == "")
            {
              group.removeAttr("email");
            }
            else
            {
              H5::Attribute attribute = group.openAttribute("email");
              attribute.write(H5::PredType::C_S1, e.c_str());
            }
          }
          else
          {
            const H5::DataSpace space;
            H5::Attribute attribute = group.createAttribute("email", H5::PredType::C_S1, space);
            attribute.write(H5::PredType::C_S1, e.c_str());
          }
          email = e;
        }
      }

    }
  }
}
