# Building on DARWIN

[DARWIN, University of Delaware HPC](https://docs.hpc.udel.edu/abstract/darwin/darwin)

Author: [Dr. Ulf Schiller](https://github.com/uschille)

```
SW_DIR=${HOME}/sw
HEMELB_SRC_DIR=${SW_DIR}/hemelb/src
HEMELB_BUILD_DIR=${SW_DIR}/hemelb/build
HEMELB_INSTALL_DIR=${SW_DIR}/hemelb
WORKGROUP=your-workgroup-name-here
WORKDIR=$(workgroup -g ${WORKGROUP} --command 'echo ${WORKDIR}')
workgroup -g ${WORKGROUP} --command @ -- srun --partition=idle sh -s << EOF
vpkg_require cmake/3.28.3 gcc/14.2.0 openmpi/5.0.2 git/2.34.1
git clone git@github.com:hemelb-codes/hemelb.git ${HEMELB_SRC_DIR}
mkdir -p ${HEMELB_BUILD_DIR}
pushd ${HEMELB_BUILD_DIR}
cmake -D CMAKE_C_COMPILER=/opt/shared/gcc/14.2.0/bin/gcc -D CMAKE_CXX_COMPILER=/opt/shared/gcc/14.2.0/bin/g++ -D CMAKE_INSTALL_PREFIX=${HEMELB_INSTALL_DIR} -D HEMELB_DEPENDENCIES_INSTALL_PREFIX=${WORKDIR}/sw ${HEMELB_SRC_DIR}
make
popd
${HEMELB_INSTALL_DIR}/bin/hemelb-tests
EOF
```