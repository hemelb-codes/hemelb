// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_UNITTESTS_REDBLOOD_INTEGRATIONS_H
#define HEMELB_UNITTESTS_REDBLOOD_INTEGRATIONS_H
namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
    }
  }
}
;
namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      using namespace hemelb::redblood;
    }
  }
}
namespace hemelb
{
  namespace unittests
  {
    namespace redblood_parallel
    {
      using namespace hemelb::redblood;
      using namespace hemelb::redblood::parallel;
    }
  }
}
#define HEMELB_DOING_UNITTESTS
#undef HEMELB_DOING_UNITTESTS
#endif

