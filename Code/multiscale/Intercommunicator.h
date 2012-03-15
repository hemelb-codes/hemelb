#ifndef HEMELB_MULTISCALE_INTERCOMMUNICATOR_H
#define HEMELB_MULTISCALE_INTERCOMMUNICATOR_H

#include "IntercommunicandType.h"
#include "Intercommunicand.h"
#include <map>

namespace hemelb
{
  namespace multiscale
  {
    
    template<class RuntimeTypeImplementation> class Intercommunicator
    {
      public:
        typedef typename RuntimeTypeImplementation::RuntimeType RuntimeType;
        typedef IntercommunicandType<RuntimeTypeImplementation> IntercommunicandTypeT;
        typedef std::map<Intercommunicand *,std::pair<IntercommunicandTypeT *,std::string> > ContentsType;
        
        void AdvanceTime(double newtime);
        bool ShouldAdvance(); // return false if this simulation is ahead in time, and should wait.
        
        void RegisterIntercommunicand(IntercommunicandTypeT & resolver,Intercommunicand & intercommunicand, const std::string &label){
          registered_objects.insert(std::make_pair(&intercommunicand,std::make_pair(&resolver,label)));
        }
        void GetFromMultiscale();
        void SendToMultiscale();

       protected:
        ContentsType registered_objects;
    };
  }
}

#endif // HEMELB_MULTISCALE_INTERCOMMUNICATOR_H