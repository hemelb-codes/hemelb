# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from argparse import ArgumentParser

def main():
    # Parse command line arguments
    parser = ArgumentParser(description='Process an input STL file'
                            'into suitable input for HemeLB.')
    
    parser.add_argument('--profile', default=None, help='Load the '
                        'profile to use from an existing file. Other '
                        'options given override those in the profile '
                        'file.', metavar='PATH')
    
    parser.add_argument('--stl', default=None, dest='StlFile',
                        help='The STL file to use as input',
                        metavar='PATH')
    parser.add_argument('--geometry', default=None,
                        dest='OutputGeometryFile',
                        help='Config output file', metavar='PATH')
    parser.add_argument('--xml', default=None, dest='OutputXmlFile',
                        help='XML output file', metavar='PATH')

    args = parser.parse_args()
    
    # Separate the profile argument.
    profile = args.profile
    del args.profile

    # Import our module late to give erroneous args a chance to be
    # caught quickly
    from HemeLbSetupTool.App import SetupTool
    
    # Pass in the Namespace's __dict__, where it actually stores the
    # arguments
    app = SetupTool(args=args.__dict__, profile=profile, redirect=False)
    app.MainLoop()
