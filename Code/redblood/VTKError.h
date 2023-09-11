// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_VTKERROR_H
#define HEMELB_REDBLOOD_VTKERROR_H

#include <string>

#include <vtkOutputWindow.h>
#include <vtkSmartPointer.h>

#include "Exception.h"

namespace hemelb::redblood {

    // Subclass of the default output window that turns error
    // messages (e.g. from vtkErrorMacro) into HemeLB exceptions.
    class ErrThrowOutputWindow : public vtkOutputWindow
    {
    public:
    vtkTypeMacro(ErrThrowOutputWindow, vtkOutputWindow);
        static ErrThrowOutputWindow* New();

        void DisplayErrorText(const char* text) override;
    protected:
        ErrThrowOutputWindow();
        virtual ~ErrThrowOutputWindow();
    };

    // RAII type that replaces the VTK global output window with an instance of the above.
    class VtkErrorsThrow {
        vtkSmartPointer<vtkOutputWindow> originalWindow;
    public:
        VtkErrorsThrow();
        ~VtkErrorsThrow();
        VtkErrorsThrow(VtkErrorsThrow const&) = delete;
        VtkErrorsThrow& operator=(VtkErrorsThrow const&) = delete;
    };

}
#endif
