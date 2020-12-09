// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_VTKERROR_H
#define HEMELB_REDBLOOD_VTKERROR_H

#include <string>

#include <vtkCommand.h>
#include <vtkSmartPointer.h>

#include "Exception.h"

namespace hemelb {
  namespace redblood {
    
    // This class implements VTK's observer interface to be notified
    // in case of error.
    // Only holds the last error - can reset with Clear()
    class ErrLogger : public vtkCommand {
    public:
      bool err_occurred = false;
      std::string err_message;

      static ErrLogger* New();

      void Clear();

      void Execute(vtkObject *vtkNotUsed(caller),
		   unsigned long event,
		   void *calldata) override;
    };

    // Wrapper for VTK objects that installs an ErrLogger observer and
    // checks for errors on use, throwing if one occurred.
    //
    // Mostly treat just like a vtkSmartPointer<T> but you can always
    // get the members with the real pointer (val) and the ErrPointer
    // (log).
    template <typename T>
    struct ErrThrower {
      using ValType = T;
      using ValPointer = vtkSmartPointer<T>;
      using ErrPointer = vtkSmartPointer<ErrLogger>;

      ValPointer val = nullptr;
      ErrPointer log = nullptr;

      // Factory
      static ErrThrower New() {
	ErrThrower ans{ValPointer::New(), ErrPointer::New()};
	ans.val->AddObserver(vtkCommand::ErrorEvent, ans.log);
	return ans;
      }

      // Use the wrapped object
      ValPointer operator->() {
	if (log->err_occurred)
	  throw Exception() << "VTK error occurred: " << log->err_message;
	return val;
      }
    };

  }
}
#endif
