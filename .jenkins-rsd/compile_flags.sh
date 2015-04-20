case $os in
    "Legion")
        module load rsd-modules
        module load hemelb-dev/$compiler
        export CC=mpicc
        export CXX=mpiCC
        export MAKEFLAGS=-j4
    ;;
    "OSX")
        case $compiler in 
            "gnu")
                # inc=/Applications/Xcode.app/Contents/Developer//Toolchains/XcodeDefault.xctoolchain
                # inc=$(dirname $(find $inc -name type_traits))
                # export CC=gcc-4.9
                # # put everything in CXX because ctemplate + autotools
                # export CXX="c++-4.9 -nodefaultlibs -nostdinc++ -isystem$inc -std=c++11"
                # export CXX="$CXX -lc++ -lc -std=c++11"
                export CXX=g++-4.9
                export CC=gcc-4.9
            ;;
            "clang")
                export CC=/usr/bin/cc
                export CXX=/usr/bin/c++
                export CXXFLAGS="-std=c++11 -stdlib=libc++"
                export LDFLAGS="-stdlib=libc++"
            ;;
        esac
        export MAKEFLAGS=-j2
    ;;
esac
