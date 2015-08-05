#ifndef HEMELB_IO_H5MD_H5MD_H
#define HEMELB_IO_H5MD_H5MD_H

namespace hemelb
{
  namespace io
  {
    namespace h5md
    {

      class ParticlesGroup;
      class TimeData;
      class FixedData;

      class H5MD
      {

        public:

          static constexpr int VERSION[] = { 1, 1 };

          static H5MD Create(const H5::CommonFG &, const std::string &, const std::string &,
                             const std::string &, const std::string & = "");

          static H5MD Open(const H5::CommonFG &);

          static bool IsH5MD(const H5::CommonFG &);

          const std::string & GetCreator() const
          {
            return creator;
          }

          const std::string & GetCreatorVersion() const
          {
            return version;
          }

          const std::string & GetAuthor() const
          {
            return author;
          }

          const std::string & GetAuthorEmail() const
          {
            return email;
          }

          void SetCreator(const std::string &);
          void SetCreatorVersion(const std::string &);
          void SetAuthor(const std::string &);
          void SetAuthorEmail(const std::string &);

          ParticlesGroup AddParticlesGroup(const std::string & name);
          void RemoveParticlesGroup(const std::string & name);
          bool HasParticlesGroup(const std::string & name) const;
          ParticlesGroup GetParticlesGroup(const std::string & name) const;

          H5MDData AddObservable(const std::string & name, int rank, const hsize_t * dims,
                                 const H5::DataType & type, const std::string & unit,
                                 bool timeDependent = true);
          void RemoveObservable(const std::string &);
          bool HasObservable(const std::string &) const;
          H5MDData GetObservable(const std::string &) const;

        private:
          H5MD(const H5::CommonFG &, const std::string &, const std::string &, const std::string &,
               const std::string &);
          std::shared_ptr<H5::CommonFG> location;
          std::string creator, version, author, email;

      };

    }
  }
}

#endif  // HEMELB_IO_H5MD_H5MD_H
