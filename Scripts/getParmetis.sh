# Get the package.
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/ParMetis-3.1.1.tar.gz

# Unzip it
gunzip ParMetis-3.1.1.tar.gz

# Make a directory to untar into.
mkdir ../Code/parmetis

# Untar it
tar -xvf ParMetis-3.1.1.tar -C ../Code/parmetis

# Edit the makefile
emacs ../Code/parmetis/ParMetis-3.1.1/Makefile.in

# Make.
make -C ../Code/parmetis/ParMetis-3.1.1

