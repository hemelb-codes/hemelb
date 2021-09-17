# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import wx

# from wx.lib.masked import NumCtrl, EVT_NUM

from ..Bindings.Translators import FloatTranslator, NoneToValueTranslator
from ..Bindings.WxMappers import WxWidgetMapper, Mapper
from .Layout import H


def ForwardGet(func):
    def Get(self, val):
        return tuple(
            getattr(getattr(self, coord), func.__name__)() for coord in ("x", "y", "z")
        )

    Get.__name__ = func.__name__
    return Get


def ForwardSet(func):
    def Set(self, val):
        for coord in ("x", "y", "z"):
            setter = getattr(getattr(self, coord), func.__name__)
            setter(val)
            continue
        return

    Set.__name__ = func.__name__
    return Set


class VectorCtrl(wx.Panel):
    """Simple container of three TextCtrl's for a vector quantity."""

    def __init__(self, parent):
        wx.Panel.__init__(self, parent)

        # self.x = NumCtrl(parent, style=wx.TE_PROCESS_ENTER, integerWidth=3, fractionWidth=3)
        # self.y = NumCtrl(parent, style=wx.TE_PROCESS_ENTER, integerWidth=3, fractionWidth=3)
        # self.z = NumCtrl(parent, style=wx.TE_PROCESS_ENTER, integerWidth=3, fractionWidth=3)
        self.x = wx.TextCtrl(self, size=(50, 22))
        self.y = wx.TextCtrl(self, size=(50, 22))
        self.z = wx.TextCtrl(self, size=(50, 22))
        sizer = H(
            (self.x, 1, wx.EXPAND), (self.y, 1, wx.EXPAND), (self.z, 1, wx.EXPAND)
        ).create()
        self.SetSizer(sizer)

        return

    @ForwardSet
    def SetBackgroundColour():
        return

    @ForwardGet
    def GetBackgroundColour():
        return

    @ForwardSet
    def SetEditable():
        return

    pass


class VectorCtrlMapper(WxWidgetMapper):
    """Widget mapper for VectorCtrls."""

    def __init__(
        self,
        widget,
        key,
        event,
        translator=NoneToValueTranslator(float("nan"), inner=FloatTranslator()),
    ):
        # We want to skip the WxWidgetMapper's init for now as the
        # VectorCtrl typically won't have the required getters and
        # setters. On binding, this one mapper is turned into three
        # standard mappers anyway.
        Mapper.__init__(self, translator=translator)

        self.widget = widget
        self.key = key
        self.event = event

        return

    def CreateSubMapper(self, component):
        return WxWidgetMapper(
            getattr(self.widget, component),
            self.key,
            self.event,
            translator=self.translator,
        )

    pass
