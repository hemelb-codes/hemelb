#ifndef HEMELB_UNITTESTS_RESOURCES_H
#define HEMELB_UNITTESTS_RESOURCES_H
namespace hemelb{
  namespace unittests{
    /***
     * Define how to find resources paths for tests which use resources
     */
    class Resource {
      public:
        Resource(const std::string &aResourceName):resourceName(aResourceName){}
        std::string Path() {
          return ResourceFolder()+resourceName;
        }
      private:
        std::string FileName() const {
          return std::string(__FILE__);
        }
        std::string ResourceFolder() const {
          size_t lLastForwardSlash = FileName().rfind('/');
          return FileName().substr(0,lLastForwardSlash+1);
        }
        std::string resourceName;
    };
  }
}
#endif // ONCE
