# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import pytest

from hlb.utils.xml_compare import XmlChecker

data = """<?xml version='1.0' encoding='US-ASCII'?>
<hemelbsettings version="3">
  <simulation>
    <step_length units="s" value="0.0008571428571" />
    <steps units="lattice" value="3000" />
    <stresstype value="1" />
    <voxel_size units="m" value="0.0001" />
    <origin units="m" value="(-0.001506801748275757,-0.0013002707004547119,-0.0015515305995941163)" />
  </simulation>
  <geometry>
    <datafile path="test.gmy" />
  </geometry>
  <inlets>
    <inlet>
      <condition subtype="cosine" type="pressure">
        <mean units="mmHg" value="80.0" />
        <amplitude units="mmHg" value="0.0" />
        <phase units="rad" value="0.0" />
        <period units="s" value="1" />
      </condition>
      <normal units="dimensionless" value="(0.5714906833845466,0.24691163832925012,-0.7825810128429067)" />
      <position units="m" value="(-0.0007175087835232063,-0.00035176118745805306,0.0008487027773518294)" />
    </inlet>
  </inlets>
  <outlets>
    <outlet>
      <condition subtype="cosine" type="pressure">
        <mean units="mmHg" value="80.0" />
        <amplitude units="mmHg" value="0.0" />
        <phase units="rad" value="0.0" />
        <period units="s" value="1" />
      </condition>
      <normal units="dimensionless" value="(-0.5838005926470043,-0.24178811710621645,0.7750583039044462)" />
      <position units="m" value="(0.0006813976915146474,0.00029353931278624915,-0.0008835672185368557)" />
    </outlet>
  </outlets>
  <visualisation>
    <centre units="m" value="(0.0,0.0,0.0)" />
    <orientation>
      <longitude units="deg" value="45.0" />
      <latitude units="deg" value="45.0" />
    </orientation>
    <display brightness="0.03" zoom="1.0" />
    <range>
      <maxvelocity units="m/s" value="0.1" />
      <maxstress units="Pa" value="0.1" />
    </range>
  </visualisation>
  <initialconditions>
    <pressure>
      <uniform units="mmHg" value="0.0" />
    </pressure>
  </initialconditions>
</hemelbsettings>
"""


def test_xmlchecker():
    xmlChecker = XmlChecker.from_string(data)

    # Self test
    xmlChecker.check_string(data)
    assert xmlChecker.elem_path == []

    # Totally different
    with pytest.raises(AssertionError):
        xmlChecker.check_string("<foo />")
    assert xmlChecker.elem_path == []

    with pytest.raises(AssertionError) as ex:
        xmlChecker.check_string("<hemelbsettings />")
    assert xmlChecker.elem_path == []
    assert ex.value.args[0].startswith("Missing attribute 'hemelbsettings:version'")

    with pytest.raises(AssertionError) as ex:
        xmlChecker.check_string('<hemelbsettings version="not 3" />')
    assert xmlChecker.elem_path == []
    assert ex.value.args[0].startswith(
        "Attribute values differ for 'hemelbsettings:version'"
    )

    with pytest.raises(AssertionError) as ex:
        xmlChecker.check_string('<hemelbsettings version="3" extra="whatever"/>')
    assert xmlChecker.elem_path == []
    assert ex.value.args[0].startswith(
        "Unexpected attribute(s) of 'hemelbsettings': 'extra'"
    )

    txt = "<root><child><granddaughter /><grandson /></child></root>"
    xmlChecker = XmlChecker.from_string(txt)
    # Self test
    xmlChecker.check_string(txt)
    # Equivalent
    xmlChecker.check_string(
        """
<root>
  <child>
    <grandson />
    <granddaughter></granddaughter>
  </child>
</root>"""
    )

    with pytest.raises(AssertionError) as ex:
        xmlChecker.check_string("<root><child><granddaughter /></child></root>")
    assert xmlChecker.elem_path == []
    assert ex.value.args[0].startswith("Missing element 'root/child/grandson'")
    with pytest.raises(AssertionError) as ex:
        xmlChecker.check_string(
            "<root><child><granddaughter /><grandson /></child><child /><branch /></root>"
        )
    assert xmlChecker.elem_path == []
    assert ex.value.args[0].startswith(
        "Unexpected child element(s) of 'root': 'child', 'branch'"
    )
