# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

"""This module contains a VTKPythonAlgorithmBase subclass that will
read a HemeLB Geometry file (.gmy) and create a VTK data structure
with the points of the wall intersections.

This module can be run as a script on the command line to convert a
.gmy file to the corresponding .vtp (VTk PolyData XML file), for
usage, run the script with no arguments.
"""

import vtk
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
import numpy as np
from vtk.util import numpy_support

VTK_MAJOR_VERSION = vtk.vtkVersion().GetVTKMajorVersion()

from ..parsers.geometry import MooreNeighbourhoodDirections
from ..parsers.geometry.generic import Site
from ..parsers.geometry.simple import ConfigLoader


class GmyWallPointsReader(VTKPythonAlgorithmBase):
    """VTK-style reader for HemeLB Geometry files (.gmy). When run, it will
    create a VTK data structure with the wall intersection points from the file.

    The output units are metres and the object has no cell data or point data.
    """

    def __init__(self, FileName=""):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkPolyData"
        )
        self.FileName = FileName

    def SetFileName(self, name):
        """The file to read."""
        self.FileName = name
        self.Modified()

    def GetFileName(self):
        """The file to read."""
        return self.FileName

    def RequestData(self, request, inInfo, outInfo):
        loader = GeneratingLoader(self.FileName)
        loader.Load()
        output = vtk.vtkPolyData.GetData(outInfo)
        output.ShallowCopy(loader.Output)
        return 1


class ChunkedPointArray:
    def __init__(self, chunk_size=10000):
        self.chunk_size = chunk_size
        self.chunks = []
        self.next_pos = 0
        self._new_chunk()

    def _mk_n_pts(self, n):
        return np.empty((n, 3), dtype=float)

    def _new_chunk(self):
        self.chunks.append(self._mk_n_pts(self.chunk_size))
        self.next_pos = 0

    def append(self, pt):
        self.chunks[-1][self.next_pos] = pt
        self.next_pos += 1
        if self.next_pos == self.chunk_size:
            self._new_chunk()

    def __len__(self):
        if len(self.chunks):
            rest = self.chunks[:-1]
            len_rest = len(rest) * self.chunk_size
            len_last = self.next_pos
            return len_rest + len_last
        else:
            return 0

    def join(self):
        N = self.chunk_size
        n = len(self)

        ans = self._mk_n_pts(n)
        # Copy the complete chunks
        for i, ch in enumerate(self.chunks[:-1]):
            ans[i * N : (i + 1) * N] = ch
        # Copy the last chunk if needed
        if self.next_pos > 0:
            i = len(self.chunks) - 1
            lo = i * N
            hi = lo + self.next_pos
            ans[lo:hi] = self.chunks[-1][: self.next_pos]
        return ans


class GeneratingLoader(ConfigLoader):
    """This subclass will create a vtkPolyData of all the wall
    intersection points.
    """

    def OnEndHeader(self):
        """Get ready to add points"""
        self.Points = ChunkedPointArray()

    def OnEndSite(self, block, site):
        """Called once a site has been read from the file."""
        if site.IsSolid:
            return
        if not site.IsEdge:
            return

        # Site is edge, look at neighbours
        for iNeigh, delta in enumerate(MooreNeighbourhoodDirections):
            itype = site.IntersectionType[iNeigh]
            # Find the intersection and add (in voxel units - can transform later)
            if itype == Site.WALL_INTERSECTION:
                pt_grid = site.Index + site.IntersectionDistance[iNeigh] * delta
                self.Points.append(pt_grid)

    def OnEndBlock(self, bIdx, bIjk):
        """Since we're done with the block, we can free its sites."""
        self.Domain.GetBlock(bIdx).DeleteSites()
        return

    def OnEndBody(self):
        """We are done, so free all we can and make the VTK dataset."""
        p = self.Points.join()
        N = len(p)

        # XmlFileName is None => lattice units, else physical
        phys_units = False if self.XmlFileName is None else True
        if phys_units:
            p *= self.Domain.VoxelSize
            p += self.Domain.Origin

        # Free
        self.__dict__.clear()
        # To VTK
        # First, points
        p = numpy_support.numpy_to_vtk(p)
        pts = vtk.vtkPoints()
        pts.SetData(p)

        # Next cells (i.e. vertices)
        vert_ca = vtk.vtkCellArray()
        if VTK_MAJOR_VERSION >= 9:
            # VTK 9 uses new internal representation of cells in two arrays
            # offsets[i] = index of first connectivity ID for cell i
            # offsets[i+1] = index of past then end connectivity ID for cell i
            # (hence the N+1)
            # connectivity = concatenation of all the point IDs - use offsets to jump in with random access.
            offsets = np.arange(N + 1)
            connectivity = offsets[:-1]
            vert_ca.SetData(
                numpy_support.numpy_to_vtk(offsets),
                numpy_support.numpy_to_vtk(connectivity),
            )
        else:
            # Use old (npts, p0, p1, ...) representation of cells
            cells = np.empty((N, 2), dtype=int)
            cells[:, 0] = 1
            cells[:, 1] = np.arange(N)
            vert_ca.SetCells(N, numpy_support.numpy_to_vtkIdTypeArray(cells))

        # Finally build the polydata
        pd = vtk.vtkPolyData()
        pd.SetPoints(pts)
        pd.SetVerts(vert_ca)

        self.Output = pd
        return


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create a VTK polydata with the wall intersection points."
    )
    parser.add_argument(
        "input", help="the XML or GMY file - GMY files will be left in lattice units"
    )
    parser.add_argument(
        "output",
        nargs=argparse.OPTIONAL,
        help="file to write VTK data to, if omitted, use " 'input basename + ".vtp"',
        default=None,
    )

    args = parser.parse_args()

    import os.path

    input = args.input

    # If output wasn't given, use the input with '.vtp' extension.
    if args.output is None:
        base, ext = os.path.splitext(input)
        output = base + ".vtp"
    else:
        output = args.output
        pass

    reader = GmyWallPointsReader()
    reader.SetFileName(input)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output)

    writer.SetInputConnection(reader.GetOutputPort())

    writer.Write()
