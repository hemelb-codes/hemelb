# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

"""This tool upgrades an XML file and the referenced GMY file from GMY
version 3 to 4. Upgraded versions of both will be written to a new
directory "converted" next to the input XML file.

This is very much a rough-and-ready script, but perhaps worth keeping
around.
"""

import argparse
from pathlib import Path
import shutil
import xdrlib
from xml.etree import ElementTree as ET

from ..parsers.geometry import HemeLbMagicNumber, GeometryMagicNumber

pulsatile_period_s = 60.0 / 70.0
reference_pressure_mmHg = 80.0


def xyz2val(el, units=None):
    """Turn an element like
    <foo x="a" y="b" z="c" />
    into
    <foo value="(a,b,c)" />
    """
    x = el.attrib.pop("x")
    y = el.attrib.pop("y")
    z = el.attrib.pop("z")
    el.attrib["value"] = f"({x},{y},{z})"
    if units:
        el.attrib["units"] = units


def attrib2subel(el, name, units):
    """Transform from
    <foo attrib="val" />
    to
    <foo>
      <attrib value="val", units="units" />
    </foo>
    """
    val = el.attrib.pop(name)
    ET.SubElement(el, name, value=val, units=units)


def update_iolet_elem(ioletEl):
    """Update an inlet or outlet to newer XML"""
    pEl = ioletEl.find("pressure")
    amp = pEl.attrib.pop("amplitude")
    mean = pEl.attrib.pop("mean")
    phase = pEl.attrib.pop("phase")
    ioletEl.remove(pEl)
    condEl = ET.SubElement(ioletEl, "condition", type="pressure", subtype="cosine")
    ET.SubElement(condEl, "amplitude", value=amp, units="mmHg")
    ET.SubElement(condEl, "mean", value=mean, units="mmHg")
    ET.SubElement(condEl, "phase", value=phase, units="rad")
    ET.SubElement(condEl, "period", value=str(pulsatile_period_s), units="s")

    xyz2val(ioletEl.find("position"), "m")
    xyz2val(ioletEl.find("normal"), "dimensionless")


def convert_gmy(source_xml):
    """Do the upgrading"""
    tree = ET.ElementTree()
    root = tree.parse(source_xml)

    gmyEl = root.find("geometry/datafile")
    assert gmyEl is not None, "XML does not have element 'geometry/datafile'!"
    relpath = gmyEl.get("path")
    assert (
        relpath is not None
    ), "Element 'geometry/datafile' does not have attribute 'path'"

    source_dir = source_xml.parent
    source_gmy = source_dir / relpath

    dest_dir = source_dir / "converted"
    if dest_dir.exists():
        raise ValueError("output dir already exists!")

    dest_xml = dest_dir / source_xml.name
    dest_gmy = dest_dir / relpath

    # Read the old header
    v3_header_size = 64
    with open(source_gmy, "rb") as sgf:
        # Note we keep the file open until the copy of the body is done
        source_header = sgf.read(v3_header_size)
        sunp = xdrlib.Unpacker(source_header)

        heme_magic = sunp.unpack_uint()
        assert heme_magic == HemeLbMagicNumber

        gmy_magic = sunp.unpack_uint()
        assert gmy_magic == GeometryMagicNumber

        version = sunp.unpack_uint()
        assert version == 3

        blocks = [sunp.unpack_uint() for i in range(3)]
        block_size = sunp.unpack_uint()

        voxel_size = sunp.unpack_double()
        origin = tuple(sunp.unpack_double() for i in range(3))

        zero = sunp.unpack_uint()
        assert zero == 0

        # Create the new header
        dp = xdrlib.Packer()
        dp.pack_uint(HemeLbMagicNumber)
        dp.pack_uint(GeometryMagicNumber)
        dp.pack_uint(4)
        for b in blocks:
            dp.pack_uint(b)
        dp.pack_uint(block_size)
        dp.pack_uint(0)
        dest_header = dp.get_buffer()

        # Add voxel and origin to XML, and other updates
        root.attrib["version"] = "5"

        simEl = root.find("simulation")

        st = simEl.attrib.pop("stresstype")
        ET.SubElement(simEl, "stresstype", value=st)

        cycles = int(simEl.attrib.pop("cycles"))
        cyclesteps = int(simEl.attrib.pop("cyclesteps"))
        nsteps = cycles * cyclesteps
        timestep = pulsatile_period_s / cyclesteps
        ET.SubElement(simEl, "step_length", value=str(timestep), units="s")
        ET.SubElement(simEl, "steps", value=str(nsteps), units="lattice")

        ET.SubElement(simEl, "voxel_size", value=str(voxel_size), units="m")
        ET.SubElement(simEl, "origin", value="(%f,%f,%f)" % origin, units="m")

        ic = ET.SubElement(root, "initialconditions")
        p = ET.SubElement(ic, "pressure")
        ET.SubElement(p, "uniform", value="0.0", units="mmHg")

        for ioletEl in root.iterfind("inlets/inlet"):
            update_iolet_elem(ioletEl)
        for ioletEl in root.iterfind("outlets/outlet"):
            update_iolet_elem(ioletEl)

        xyz2val(root.find("visualisation/centre"), "m")
        orientEl = root.find("visualisation/orientation")
        attrib2subel(orientEl, "latitude", "deg")
        attrib2subel(orientEl, "longitude", "deg")

        rangeEl = root.find("visualisation/range")
        attrib2subel(rangeEl, "maxstress", "Pa")
        attrib2subel(rangeEl, "maxvelocity", "m/s")

        try:
            ET.indent(root)
        except AttributeError:
            # ET.indend only present in python >= 3.9
            pass

        dest_dir.mkdir()
        # Write new XML
        with open(dest_xml, "wb") as xmlFile:
            xmlFile.write(b'<?xml version="1.0" ?>\n')
            tree.write(xmlFile)

        # Copy body of GMY over
        with open(dest_gmy, "wb") as dgf:
            dgf.write(dest_header)
            shutil.copyfileobj(sgf, dgf)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("source_xml")
    args = p.parse_args()

    source_xml = Path(args.source_xml)
    assert source_xml.exists()
    convert_gmy(source_xml)


if __name__ == "__main__":
    main()
