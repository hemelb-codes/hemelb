// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "redblood/VTKError.h"
#include <vtkObjectFactory.h>

namespace hemelb::redblood {

    vtkStandardNewMacro(ErrThrowOutputWindow);

    ErrThrowOutputWindow::ErrThrowOutputWindow() = default;

    ErrThrowOutputWindow::~ErrThrowOutputWindow() = default;

    void ErrThrowOutputWindow::DisplayErrorText(char const* txt) {
        throw Exception() << txt;
    }

    VtkErrorsThrow::VtkErrorsThrow() : originalWindow(vtkOutputWindow::GetInstance()) {
        vtkOutputWindow::SetInstance(ErrThrowOutputWindow::New());
    }

    VtkErrorsThrow::~VtkErrorsThrow() {
        vtkOutputWindow::SetInstance(originalWindow);
    }
}
