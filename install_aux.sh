#!/usr/bin/env bash

curdir=$PWD
cd $(dirname $0)
srcdir=$PWD
cd ${curdir}

if [[ "$#" -eq "1" ]]; then
    LIBHINT=$1
fi

mkdir -p build_aux/local
cd build_aux
if [[ "$#" -eq "1" ]]; then
    cmake ${srcdir}/source/aux -DCMAKE_INSTALL_PREFIX=${PWD}/local -DCMAKE_BUILD_TYPE=Release -DBIOINFO_ROOT_DIR=$1 -DLOONLIB_ROOT_DIR=$1
else
    cmake ${srcdir}/source/aux -DCMAKE_INSTALL_PREFIX=${PWD}/local -DCMAKE_BUILD_TYPE=Release
fi
make install

echo "The aux binaries have been installed into \"${PWD}/local/bin\". Please add it the PATH before using it"
