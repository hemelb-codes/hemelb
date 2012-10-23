// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_PHASED_CONCERN_H
#define HEMELB_NET_PHASED_CONCERN_H
namespace hemelb{
  namespace net{
    namespace phased{
      class Concern{
        public:
          Concern(){}
          virtual ~Concern(){}
          virtual bool CallAction(int action)=0;
      };
    }
  }
}
#endif //ONCE
