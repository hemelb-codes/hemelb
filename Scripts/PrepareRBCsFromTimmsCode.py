#!/usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
from argparse import ArgumentParser
import vtk
import os
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup


def write(file, thing_to_write):
    file.write(thing_to_write)
    return


def write_header(file, number_of_points, number_of_polys):
    file.write('<?xml version="1.0"?>\n')
    file.write('<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n')
    file.write('  <PolyData>\n')
    file.write('    <Piece NumberOfPoints="' + str(number_of_points) + '" NumberOfPolys="' + str(number_of_polys) + '">\n')
    file.write('      <Points>\n')
    file.write('        <DataArray NumberOfComponents="3" type="Float32">\n')


def write_points_footer(file):
    file.write('\n')
    file.write('        </DataArray>\n')
    file.write('      </Points>\n')


def write_connectivity_header(file):
    file.write('      <Polys>\n')
    file.write('        <DataArray type="Int32" Name="connectivity">\n')


def write_connectivity_footer(file):
    file.write('\n')
    file.write('        </DataArray>\n')
    file.write('        <DataArray type="Int32" Name="offsets">\n')


def write_footer(file):
    file.write('\n')
    file.write('        </DataArray>\n')
    file.write('      </Polys>\n')
    file.write('    </Piece>\n')
    file.write('  </PolyData>\n')
    file.write('</VTKFile>')


def write_list(file, list):
    for item in list:
        write(file, str(item))
        write(file, " ")


def get_node_array(polydata, centre_of_mass):
    cell_position_list = []
    points = polydata.GetPoints()
    number_of_points = polydata.GetNumberOfPoints()
    np_node_positions = np.array([points.GetPoint(i_1) for i_1 in range(points.GetNumberOfPoints())])
    for node_triplet in np_node_positions:
        node_triplet[0] = node_triplet[0] - centre_of_mass[0]
        node_triplet[1] = node_triplet[1] - centre_of_mass[1]
        node_triplet[2] = node_triplet[2] - centre_of_mass[2]
        # triplet = (-node_triplet[2], node_triplet[1], node_triplet[0])  # this rotates by 90 degrees to face z instead of x axis
        triplet = (-node_triplet[2], node_triplet[1], -node_triplet[0])  # need to mirro image around xy plane, test 1
        for node in triplet:
            cell_position_list.append(round(node, 9))
    return (cell_position_list, number_of_points)


def get_connectivity_and_offset_array(polydata):
    connectivity_list = []
    offset_list = []
    polys = polydata.GetPolys()
    number_of_polys = polydata.GetNumberOfPolys()
    poly_points = polys.GetData()
    poly_array = vtk_to_numpy(poly_points)
    array = np.resize(poly_array, (number_of_polys, 4))
    offset = 3
    for row in array:
        offset_list.append(offset)
        offset += 3
        connectivity_list.append(row[1]+1)
        connectivity_list.append(row[2]+1)
        connectivity_list.append(row[3]+1)
    return (connectivity_list, offset_list, number_of_polys)


def write_list_msh(file, list):
    resized_list = np.resize(np.array(list), (len(list)/3, 3))
    for count, item in enumerate(resized_list):
        counter = count + 1
        write(file, str(counter))
        write(file, " ")
        write(file, str(item[0]))
        write(file, " ")
        write(file, str(item[1]))
        write(file, " ")
        write(file, str(item[2]))
        write(file, "\n")


def write_connectivity_msh(file, list):
    resized_list = np.resize(np.array(list), (len(list)/3, 3))
    for count, item in enumerate(resized_list):
        counter = count + 1
        write(file, str(counter))
        write(file, " ")
        dummy = "2 3 0 1 0 "
        write(file, dummy)
        write(file, str(item[0]))
        write(file, " ")
        write(file, str(item[1]))
        write(file, " ")
        write(file, str(item[2]))
        write(file, "\n")


def write_vtp_file(file_name, output_directory, connectivity_list, offset_list, cell_position_list, number_of_polys, number_of_points):
    new_file_name = file_name.replace('.vtk', '.vtp')
    new_file_name = os.path.join(cell_output_directory, new_file_name)
    file = open(new_file_name, "w")
    write_header(file, number_of_points, number_of_polys)
    write_list(file, cell_position_list)
    write_points_footer(file)
    write_connectivity_header(file)
    write_list(file, connectivity_list)
    write_connectivity_footer(file)
    write_list(file, offset_list)
    write_footer(file)
    file.close()
    return


def write_msh_format(file_name, output_directory, connectivity_list, cell_position_list, number_of_polys, number_of_points, undeformed_msh_file):
    new_file_name = file_name.replace('.vtk', '.msh')
    new_file_name = os.path.join(cell_output_directory, new_file_name)
    file = open(new_file_name, "w")
    file.write('$MeshFormat\n')
    file.write('2 0 8\n')
    file.write('$EndMeshFormat\n')
    file.write('$Nodes\n')
    file.write(str(number_of_points) + "\n")
    write_list_msh(file, cell_position_list)
    file.write('$EndNodes\n')
    file.write('$Elements\n')
    file.write(str(number_of_polys) + "\n")
    write_connectivity_msh(file, connectivity_list)
    file.write('$EndElements')
    # with open(undeformed_msh_file) as read_msh: #This section here reads the ico_msh file and puts the connectivity from there. Unsure of how it all links together and wether this is maybe the source of error?
    #     start_point = 6 + number_of_points
    #     for counter, row in enumerate(read_msh):
    #         if counter >= start_point:
    #             file.write(row)
    file.close()
    return ()


def get_com(polydata):
    centreOfMassCalculator = vtk.vtkCenterOfMass()
    centreOfMassCalculator.SetInputData(polydata)
    centreOfMassCalculator.Update()
    centreOfMassCalculator.SetUseScalarsAsWeights(False)
    centreOfMass = centreOfMassCalculator.GetCenter()
    return centreOfMass


def xml_insert_cell_block(template, offset, x_value, y_value):
    new = ET.Element("insertcell", template=template.replace(".vtk", ""))
    ET.SubElement(new, "offset", units="LB", value=str(offset))
    ET.SubElement(new, "every", units="s", value="100")
    ET.SubElement(new, "x", units="LB", value=str(x_value))
    ET.SubElement(new, "y", units="LB", value=str(y_value))
    return new


def xml_cell_property_block(template, file_path, bending_modulus, strain_modulus):
    new = ET.Element("cell", name=template)
    moduli_sub = ET.SubElement(new, "moduli")
    ET.SubElement(moduli_sub, "bending", units="Nm", value=str(bending_modulus))
    ET.SubElement(moduli_sub, "surface", units="LB", value="1.0")
    ET.SubElement(moduli_sub, "volume", units="LB", value="1.0")
    ET.SubElement(moduli_sub, "dilation", units="LB", value="0.5")
    ET.SubElement(moduli_sub, "strain", units="N/m", value=str(strain_modulus))
    ET.SubElement(new, "shape", mesh_path=file_path)
    ET.SubElement(new, "scale", units="LB", value="1")
    return new


def get_cell_timestep(cell_file):
    intermediate_step = cell_file.replace('.vtk', '')
    timestep = intermediate_step.split("t")
    return int(timestep[-1])


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('xml_file', type=str, help='HemeLB XML file for simulation that will load RBCs')
    parser.add_argument('cell_directory', type=str, help='Directory containing RBCs generated with Timm\'s code')
    parser.add_argument('cell_output_directory', type=str, help='Directory containing RBCs generated from ths script in vtp format')
    parser.add_argument('bending_modulus', type=float, help='float value for the bending modulus of the RBCs')
    parser.add_argument('strain_modulus', type=float, help='float value for the strain modulus of the RBCs')
    parser.add_argument('undeformed_msh_file', type=str, help='the ico_msh_#ofFace.msh file to write down the connectivity')

    args = parser.parse_args()

    cell_directory = args.cell_directory
    xml_file = args.xml_file
    cell_output_directory = args.cell_output_directory  # need to drag and drop folder into terminal for this one to work, not too sure what the general solution is
    bending_modulus = args.bending_modulus
    strain_modulus = args.strain_modulus
    undeformed_msh_file = args.undeformed_msh_file

    tree = ET.parse(xml_file)
    root = tree.getroot()

    for cell in os.listdir(cell_directory):
        cell_path = os.path.join(cell_directory, cell)

        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(cell_path)
        reader.Update()
        polydata = reader.GetOutput()

        cell_timestep = get_cell_timestep(cell)

        centre_of_mass = get_com(polydata)
        cell_position_list, number_of_points = get_node_array(polydata, centre_of_mass)
        connectivity_list, offset_list, number_of_polys = get_connectivity_and_offset_array(polydata)

        write_msh_format(cell, cell_output_directory, connectivity_list, cell_position_list, number_of_polys, number_of_points, undeformed_msh_file)
        # write_vtp_file(cell, cell_output_directory, connectivity_list, offset_list, cell_position_list, number_of_polys, number_of_points) #cell needs to face in z direction

        cell_block = xml_insert_cell_block(cell, cell_timestep, (centre_of_mass[2]-25), (centre_of_mass[1]-25))
        inlets = root.find("inlets")
        inlet = inlets.find("inlet")
        inlet.append(cell_block)

        cell_path = os.path.join(os.path.basename(cell_output_directory), cell.replace(".vtk", ".msh"))
        cell_property_block = xml_cell_property_block(cell.replace(".vtk", ""), cell_path, bending_modulus, strain_modulus)
        redbloodcells = root.find("redbloodcells")
        cells = redbloodcells.find("cells")
        cells.append(cell_property_block)

    xml_file_name = os.path.join(cell_output_directory, "config.xml")
    tree.write(xml_file_name, xml_declaration=True, encoding="utf-8", method="xml")
