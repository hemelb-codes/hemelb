#ifndef HEMELB_NET_PHASED_CONCERN_H
#define HEMELB_NET_PHASED_CONCERN_H
namespace hemelb{
  namespace net{
    namespace phased{
      class Concern{
        public:
          Concern(){}
          virtual ~Concern(){}
          virtual bool CallAction()=0;
      };
    }
  }
}
#endif //ONCE
