// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "redblood/VTKError.h"

namespace hemelb {
  namespace redblood {

    ErrLogger* ErrLogger::New() {
      return new ErrLogger;
    }

    void ErrLogger::Clear() {
      *this = ErrLogger{};
    }

    void ErrLogger::Execute(vtkObject *vtkNotUsed(caller),
			    unsigned long event,
			    void *calldata)
    {
      if (event == vtkCommand::ErrorEvent) {
	err_occurred = true;
	err_message = static_cast<const char *>(calldata);
      }
    }

  }
}
