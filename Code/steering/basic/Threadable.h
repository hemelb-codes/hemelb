#ifndef HEMELB_STEERING_BASIC_THREADABLE
#define HEMELB_STEERING_BASIC_THREADABLE

#include <pthread.h>

namespace hemelb
{
  namespace steering
  {

    class Threadable
    /**
     * An abstract class wrapping a pthread.
     * You must implement the pure virtual method
     * void DoWork(void) which is executed in the
     * new thread. It's up to you how to exchange
     * data with the thread.
     */
    {
      public:
        // Construct the wrapping object, but don't run it.
        Threadable();
        // D'tor
        virtual ~Threadable();

        // Start the new thread
        int Run(void);

      protected:
        // Pure virtual worker method, implement this in
        // your subclass.
        virtual void DoWork(void) = 0;

        // Override this to change the thread attributes
        virtual pthread_attr_t* GetPthreadAttributes(void);

        void Exit(void);

      private:
        // Helper method given to pthread_create
        static void* PthreadLaunch(void *);

        pthread_t mPthread;
        pthread_attr_t* mPthreadAttributes;
    };

  } // namespace steering
} //namespace hemelb

#endif //HEMELB_STEERING_BASIC_THREADABLE
