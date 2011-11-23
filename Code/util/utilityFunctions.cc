#include "mpiInclude.h"
#include <uuid/uuid.h>
namespace hemelb
{
  namespace util
  {
    // Returns the number of seconds to 6dp elapsed since the Epoch
    double myClock()
    {
      return MPI_Wtime();
    }

    std::string GetUUID(){
      uuid_t uuid;
      uuid_generate_random ( uuid );
      char s[37];
      uuid_unparse ( uuid, s );
      return std::string(s);
    }

  }

}
