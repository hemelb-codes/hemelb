#!/usr/bin/env python

# For all folders below that specified as first argument, ensure that all files with .h .hpp extensions have
# once guards matching the path and extension, underscore separated and capitalised
import sys
import glob
import re
root_path=sys.argv[1]
print 'Checking for correctness of once guards in sources below %s'% root_path
filenames=list(glob.glob('%s/**/*.h'%root_path))+list(glob.glob('%s/**/*.hpp'%root_path))
ok=True
for filename in filenames:
	with open(filename) as file:
		firstline=file.readline()
		secondline=file.readline()
		firstmatch=re.match('\s*#ifndef\s*HEMELB_([^\s]*)',firstline)
		secondmatch=re.match('\s*#define\s*HEMELB_([^\s]*)',secondline)
		if not firstmatch:
			print "First line of %s doesn't have a compliant #ifndef: %s" % (filename,firstline)
			ok=False
			continue
		if not secondmatch:
			print "First line of %s doesn't have a compliant #define: %s" % (filename,secondline)
			ok=False
			continue
		firstguard=firstmatch.group(1)
		secondguard=secondmatch.group(1)
		if not firstguard==secondguard:
			print "Once guard in %s doesn't have same value in #ifndef and #define: %s vs %s" % (filename, firstguard, secondguard)
			ok=False
		objective=filename.replace(root_path+'/','').replace('/','_').replace('.','_').upper()
		if not firstguard==objective:
			print "#ifndef in %s doesn't have file-name-based value: %s vs %s" % (filename, firstguard, objective)
			ok=False
		if not secondguard==objective:
			print "#define in %s doesn't have file-name-based value: %s vs %s" % (filename, secondguard, objective)
			ok=False
sys.exit(not ok)