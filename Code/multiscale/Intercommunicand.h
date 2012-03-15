#ifndef HEMELB_MULTISCALE_INTERCOMMUNICAND_H
#define HEMELB_MULTISCALE_INTERCOMMUNICAND_H
#include <vector>


namespace hemelb
{
  namespace multiscale
  {
    class BaseSharedValue;
    class Intercommunicand
    {
      public:
        void RegisterSharedValue(BaseSharedValue *value){
          values.push_back(value);
        }
        Intercommunicand(): values(){}
        std::vector<BaseSharedValue *> & Values(){return values;}
      private:
        std::vector<BaseSharedValue *> values;
    };
  }
}

#endif // HEMELB_MULTISCALE_INTERCOMMUNICAND_H
