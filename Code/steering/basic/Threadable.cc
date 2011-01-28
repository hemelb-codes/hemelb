#include <pthread.h>

#include "steering/basic/Threadable.h"

namespace hemelb
{
  namespace steering
  {

    Threadable::Threadable()
    {
      mPthreadAttributes = GetPthreadAttributes();
    }

    Threadable::~Threadable()
    {
      delete mPthreadAttributes;
    }

    pthread_attr_t* Threadable::GetPthreadAttributes(void)
    {
      return NULL;
    }

    int Threadable::Run(void)
    {
      return pthread_create(&mPthread, mPthreadAttributes,
                            Threadable::PthreadLaunch,
                            static_cast<void*> (this));
    }

    void Threadable::Exit(void) {
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
