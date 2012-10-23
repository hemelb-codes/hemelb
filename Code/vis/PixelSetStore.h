// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_PIXELSETSTORE_H
#define HEMELB_VIS_PIXELSETSTORE_H

#include <queue>
#include "log/Logger.h"

namespace hemelb
{
  namespace vis
  {
    template<typename PixelSetType>
    class PixelSetStore
    {
      public:
        PixelSetType* GetUnusedPixelSet()
        {
          if (!store.empty())
          {
            // Take the first pixel set (and put it on the back of the queue for reuse later).
            PixelSetType* candidate = store.front();

            store.pop();
            store.push(candidate);

            // If we can use this one to write, return it.
            if (!candidate->IsInUse())
            {
              candidate->SetInUse();
              candidate->Clear();
              return candidate;
            }
          }

          // We've not got any inUse pixel sets, so create a new one, store it and return it.
          PixelSetType* ret = new PixelSetType();

          store.push(ret);

          ret->SetInUse();
          return ret;
        }

        ~PixelSetStore()
        {
          while (!store.empty())
          {
            if (store.front()->IsInUse())
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("This could be a problem: we've just "
                "cleared out a pixel set which appears to still be in use. If pixel sets are being "
                "managed properly there shouldn't be more than a few of these messages per core.");
            }

            delete store.front();
            store.pop();
          }
        }

      private:
        std::queue<PixelSetType*> store;
    };
  }
}

#endif // HEMELB_VIS_PIXELSETSTORE_H
