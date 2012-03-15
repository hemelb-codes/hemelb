#ifndef HEMELB_MULTISCALE_INTERCOMMUNICAND_TYPE_H
#define HEMELB_MULTISCALE_INTERCOMMUNICAND_TYPE_H
#include <string>
#include <vector>

namespace hemelb
{
  namespace multiscale
  {
    // Different implementations of intercommunicator will use different systems to specify type at runtime
    // template on this choice.
    // Example: MPI Datatype
    template<class RuntimeTypeImplementation> class IntercommunicandType
    {
      public:
       
        IntercommunicandType(const std::string & alabel): fields(), label(alabel){}
        
        std::vector<std::pair<std::string, typename RuntimeTypeImplementation::RuntimeType> > &Fields(){return fields;}
        
        void RegisterSharedValue(const std::string &label,typename RuntimeTypeImplementation::RuntimeType type)
         {
            fields.push_back(std::make_pair(label,type));
         }
         template<class T> void RegisterSharedValue(const std::string &label){
            fields.push_back(std::make_pair(label,RuntimeTypeImplementation::template GetType<T>()));
         } 
      private:
        std::vector<std::pair<std::string, typename RuntimeTypeImplementation::RuntimeType> >fields;
        std::string label;
    };
  }
}

#endif // HEMELB_MULTISCALE_INTERCOMMUNICAND_H
