// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_STEERING_COMMON_SINGLETON_H
#define HEMELB_STEERING_COMMON_SINGLETON_H

namespace hemelb
{
  namespace steering
  {
    namespace _private
    {
      /**
       * Base class for singleton behaviour.
       * Uses CRTP to do static polymorphism.
       */
      template<typename Derived>
      class Singleton
      {
        protected:
          static Singleton* singleton;
          Singleton()
          {
          }

        public:
          static Derived* Init()
          {
            if (Singleton<Derived>::singleton == nullptr)
            {
              Singleton<Derived>::singleton = new Derived();
            }
            return static_cast<Derived*>(Singleton<Derived>::singleton);
          }
      };

      template<typename Derived>
      Singleton<Derived>* Singleton<Derived>::singleton = nullptr;

    }
  }
}

#endif // HEMELB_STEERING_COMMON_SINGLETON_H
