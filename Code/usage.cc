#include <stdio.h>
#include "usage.h"

// Prints the correct command-line usage to the terminal. Argument is the programme name.
void Usage::printUsage(char *progname)
{
  printf("-!-!-!-!-!-!-!-!-!-!-!-!");
  printf("Correct usage: %s [-<Parameter Name> <Parameter Value>]* \n", progname);
  printf("Parameter name and significance:\n");
  printf("-in \t Path to the configuration xml file (default is config.xml)\n");
  printf(
         "-out \t Path to the output folder (default is based on input file, e.g. config_xml_results)\n");
  printf("-s \t Number of snapshots to take per cycle (default 10)\n");
  printf("-i \t Number of images to create per cycle (default is 10)\n");
  printf("-ss \t Steering session identifier (default is 1)\n");
  printf("-!-!-!-!-!-!-!-!-!-!-!-!");
}

