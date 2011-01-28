The HemeLB segmentation tool is a graphical-editing software
application which effectively reconstructs vascular-related data to be
used by HemeLB and edit and manipulate boundary conditions which are
visualised with triangles and discs. It can operate in two different
modalities:

(1) The first execution mode is well-suited for the interactive
threshold-based segmentation of DICOM-based volumetric datasets and
the manipulation of boundaries. It can be employed by compiling the
code without defining any constant. The command line arguments next to
the executable name should be in order

a) the path where are located the DICOM images (or slices);

b) the name of the output file where the HemeLB-based configuration
will be written to; this is the input configuration of HemeLB;

c)the name of the output file where the pressure boundary
parameters are written to;

d) the name of the checkpoint file;

e) the distance between two subsequent slices (in mm);

f) the dimension of each pixel (in mm);

g) the flag indicating if the current checkpoint should be employed to
retrieve boundary data and a few parameters (like the HemeLB-grid
resolution enhancement), --and in this case the flag should be set to
1--, or not (flag = 0). If the aforementioned flag is one, the DICOM
images are still read and used: after files reading the segmentation
is carried out on the basis of the volumetric dataset and the data
supplied by the checkpoint.

Visualisation and editing tasks:

The visualisation begins by displaying the first slice or the one
specified in the checkpoint file (if the checkpoint flag was one); the
menu attached to the right mouse-button indicates the editing
facilities.  One can segment the system at 6 different resolutions by
taking into account the current mouse cursor position, adjust the
threshold interactively if the reconstructed system is not large),
change the slice, perform the 3-dimensional rendering, save the output
files or quit.

The various editing modalities are activated or adjusted by scrolling
the left mouse-button or clicking on it.

If one activates the 3-dimensional rendering he/she can operate as
before (apart from being able to change 2-dimensional image) plus
perform editing tasks related to the boundaries, like create them at
the lattice site pointed by the mouse (if any and if the estimation of
local parameters is not unsuccessful), rotate and scale them, and
change their pressure parameters. After creating a boundary, one can
select it by mouse-pointing one of its vertices. Finally, one should
revert the direction of the inlet normal if, with respect to the
segmented system, it points inward.


(2) The second execution mode of the tool bases its processing on a
triangle-mesh which is the result of the segmentation process carried
out by another preprocessing software. Specifically, that tool is
exploited to reconstruct the HemeLB-grid within the mesh. Furthermore,
the mesh is used to calculate surface-related data, like the normal
and the distance between any wall lattice site and the surface. This
second execution mode is enabled by defining the constant "MESH" at
compile time. The command line arguments next to the executable name
should be in order

a) the name of the file which comprises the triangulation;

b) the name of the output file where the HemeLB-based configuration
will be written to; this is the input configuration of HemeLB;

c)the name of the output file where the pressure boundary
parameters are written to;

d) the name of the checkpoint file;

e) the flag indicating if the current checkpoint should be employed to
retrieve boundary data and a few parameters (like the HemeLB-grid
resolution enhancement), --and in this case the flag should be set to
1--, or not (flag = 0). If the aforementioned flag is one, the mesh is
still read and used: after files reading the HemeLB-grid is
reconstructed on the basis of the triangulation and the data
supplied by the checkpoint.

Visualisation and editing tasks:

The visualisation begins by displaying the triangulation. Then, the
only difference with the editing-visualisation modality approached by
the first segmentation mode (achieved by not defining "MESH" at
compile time) whilst performing 3-dimensional rendering is the fact
that the threshold cannot be adjusted. The HemeLB-grid is
reconstructed by applying a clustering algorithm (at the menu-based
resolution selected by the user) to a seed point located at the center
of the artery location pointed by the mouse cursor coordinates.

Even if the software is very efficient, it has been optimised for
accuracy purposes. However, it might happen that, during the
reconstruction process, the clustering algorithm fails to contain the
HemeLB-grid within the triangulation.
