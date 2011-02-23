#include <pthread.h>
#include "steering/basic/Threadable.h"

namespace hemelb
{
  namespace steering
  {

    Threadable::Threadable()
    {
    }

    Threadable::~Threadable()
    {
    }

    int Threadable::Run(void)
    {
      return pthread_create(&mPthread, NULL, Threadable::PthreadLaunch, static_cast<void*> (this));
    }

    void Threadable::Exit(void)
    {
      pthread_exit(&mPthread);
    }

    void* Threadable::PthreadLaunch(void *ptr)
    {
      Threadable *instance = static_cast<Threadable*> (ptr);
      instance->DoWork();
      return NULL;
    }

  } // namespace steering
} //namespace hemelb
