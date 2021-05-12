#!/usr/bin/env python3
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import argparse
import sys
import xml.etree.ElementTree as ET

expected_fail_names = {
    "hemelb::unittests::redblood::LoadingTimmMeshTests::testSimulationOutputFile",
    "hemelb::unittests::redblood::NodeIntegrationTests::testNodeWall<hemelb::redblood::stencil::FourPoint>",
    "hemelb::unittests::redblood::NodeIntegrationTests::testNodeWall<hemelb::redblood::stencil::CosineApprox>",
    "hemelb::unittests::redblood::NodeIntegrationTests::testNodeWall<hemelb::redblood::stencil::ThreePoint>",
}

p = argparse.ArgumentParser()
p.add_argument("testxml")
args = p.parse_args()

tree = ET.parse(args.testxml)
root = tree.getroot()

fails = root.findall('FailedTests/FailedTest')

expected = ET.SubElement(root, "ExpectedFailedTests")

for f in fails:
    name = f.find("Name").text
    if name in expected_fail_names:
        root.find("FailedTests").remove(f)
        expected.append(f)
        expected_fail_names.remove(name)

n_unexpected_fails = len(root.find("FailedTests").getchildren())
n_unexpected_passes = len(expected_fail_names)
if n_unexpected_passes:
    sys.stderr.write("Some tests unexpectedly passing: " + ", ".join(expected_fail_names) + "\n")

ET.dump(root)

sys.exit(n_unexpected_fails + n_unexpected_passes)
