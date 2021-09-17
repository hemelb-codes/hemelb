import xml.etree.ElementTree as ET


class XmlChecker:
    """Test support class to compare XML files/streams.

    Decidedly not thread safe.
    """

    def __init__(self, ref):
        """Construct with an ET.Element with the reference document."""
        self.ref_root = ref
        self.elem_path = []

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
