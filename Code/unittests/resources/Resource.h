#ifndef HEMELB_UNITTESTS_RESOURCES_RESOURCE_H
#define HEMELB_UNITTESTS_RESOURCES_RESOURCE_H
namespace hemelb
{
  namespace unittests
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
          std::string Path()
          {
            return ResourceFolder() + resourceName;
          }
        private:
          std::string FileName() const
          {
            return std::string(__FILE__);
          }
          std::string ResourceFolder() const
          {
            size_t lastForwardSlash = FileName().rfind('/');
            return FileName().substr(0, lastForwardSlash + 1);
          }
          std::string resourceName;
      };
    }
  }
}
#endif // ONCE
