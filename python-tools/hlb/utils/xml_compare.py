# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from fnmatch import fnmatch
from functools import partial
import xml.etree.ElementTree as ET
import numpy as np

QuantityKeys = {"units", "value"}


class QuantityCheck:
    def __init__(self, tol=1e-8):
        self.tol = tol

    def __call__(self, checker, ref, test):
        assert (
            set(ref.keys()) == QuantityKeys
        ), f"Invalid reference attributes for quantity element '{checker.path}'"
        assert (
            set(test.keys()) == QuantityKeys
        ), f"Invalid test attributes for quantity element '{checker.path}'"
        assert ref["units"] == test["units"], f"Units differ for '{checker.path}'"


class ScalarQuantityCheck(QuantityCheck):
    def __call__(self, checker, ref, test):
        super().__call__(checker, ref, test)
        rval = float(ref["value"])
        tval = float(test["value"])
        assert np.abs(rval - tval) < self.tol, f"Values differ for '{checker.path}'"


class VectorQuantityCheck(QuantityCheck):
    @staticmethod
    def vec3(txt):
        assert txt.startswith("(") and txt.endswith(")")
        return np.array(txt[1:-1].split(","), dtype=float)

    def __call__(self, checker, ref, test):
        rval = self.vec3(ref["value"])
        tval = self.vec3(test["value"])
        assert np.all(
            np.abs(rval - tval) < self.tol
        ), f"Values differ for '{checker.path}'"


class XmlChecker:
    """Test support class to compare XML files/streams.

    Decidedly not thread safe.

    You can add checks beyond string equality by adding key-value pairs
    to the attr_checks mapping. Keys are element paths (compared via
    fnmatch) and values are Callable[XmlChecker, ref_attr_dict,
    test_attr_dict].
    """

    def __init__(self, ref):
        """Construct with an ET.Element with the reference document."""
        self.ref_root = ref
        self.elem_path = []
        self.attr_checks = {}

    @classmethod
    def from_path(cls, ref_path):
        """Factory method to read from file."""
        return cls(ET.parse(ref_path).getroot())

    @classmethod
    def from_string(cls, ref_data):
        """Factory method to parse from string."""
        return cls(ET.fromstring(ref_data))

    @property
    def path(self):
        """Path to currently considered element."""
        return "/".join(self.elem_path)

    def check_attrib(self, rAttr, tAttr):
        """Assert that the attribute dicts are equivalent."""
        check = None
        for p, chk in self.attr_checks.items():
            if fnmatch(self.path, p):
                check = partial(chk, self)
                break

        if check is None:
            check = self.default_check_attrib
        check(rAttr, tAttr)

    def default_check_attrib(self, rAttr, tAttr):
        """Assert that the attribute dicts are equivalent."""
        for k in rAttr:
            assert k in tAttr, f"Missing attribute '{self.path}:{k}'"
            assert (
                rAttr[k] == tAttr[k]
            ), f"Attribute values differ for '{self.path}:{k}'"
            tAttr.pop(k)

        assert (
            len(tAttr) == 0
        ), f"Unexpected attribute(s) of '{self.path}': " + ", ".join(
            f"'{k}'" for k in tAttr.keys()
        )

    def check_elem(self, rEl, tEl):
        """Assert that the elements are equivalent, recursively, and
        ignoring text.
        """
        assert rEl.tag == tEl.tag
        self.elem_path.append(rEl.tag)
        try:
            self.check_attrib(rEl.attrib, tEl.attrib)

            for rChild in rEl:
                tChild = tEl.find(rChild.tag)
                assert tChild is not None, f"Missing element '{self.path}/{rChild.tag}'"
                self.check_elem(rChild, tChild)
                tEl.remove(tChild)
            assert (
                len(tEl) == 0
            ), f"Unexpected child element(s) of '{self.path}': " + ", ".join(
                f"'{el.tag}'" for el in tEl
            )
        finally:
            self.elem_path.pop()

    def check_path(self, test_path):
        """Entry point to check the path against the reference."""
        test_tree = ET.parse(test_path)
        test_root = test_tree.getroot()
        self.check_elem(self.ref_root, test_root)

    def check_string(self, test_str):
        """Entry point to check string against the reference."""
        test_root = ET.fromstring(test_str)
        self.check_elem(self.ref_root, test_root)
