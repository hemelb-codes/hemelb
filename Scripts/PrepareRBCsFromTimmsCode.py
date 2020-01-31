#!/usr/bin/env python

from argparse import ArgumentParser


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('xml_file', type=str, help='HemeLB XML file for simulation that will load RBCs')
    parser.add_argument('folder_name', type=str, help='Directory containing RBCs generated with Timm\'s code')
    args = parser.parse_args()