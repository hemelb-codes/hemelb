#include "vis/BasicPixel.h"

namespace hemelb
{
  namespace vis
  {
    BasicPixel::BasicPixel()
    {

    }
    BasicPixel::BasicPixel(int iIn, int jIn) :
      i(iIn), j(jIn)
    {

    }

    int BasicPixel::GetI() const
    {
      return i;
    }

    int BasicPixel::GetJ() const
    {
      return j;
    }

    bool BasicPixel::operator <(const BasicPixel& right) const
    {
      return (i < right.i || (i == right.i && j < right.j));
    }

    bool BasicPixel::operator ==(const BasicPixel& right) const
    {
      return (i == right.i && j == right.j);
    }

  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::BasicPixel>::RegisterMpiDataType()
  {
    MPI_Datatype ret = vis::BasicPixel::GetMPIType();
    MPI_Type_commit(&ret);
    return ret;
  }
}
