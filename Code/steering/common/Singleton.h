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
            if (Singleton<Derived>::singleton == NULL)
            {
              Singleton<Derived>::singleton = new Derived();
            }
            return static_cast<Derived*> (Singleton<Derived>::singleton);
          }
      };

      template<typename Derived>
      Singleton<Derived>* Singleton<Derived>::singleton = NULL;

    }
  }
}

#endif // HEMELB_STEERING_COMMON_SINGLETON_H
