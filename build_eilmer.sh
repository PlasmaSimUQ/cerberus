#!/usr/bin/env bash

if cd ${TOP}/.gdtk_git; then git pull; else git clone ${EILMER_URL} ${TOP}/.gdtk_git; fi

cd ${TOP}/.gdtk_git/src/eilmer
make prep-chem prep-gas MAKEFLAGS=${MF}

cd ${TOP}/.gdtk_git/src/gas
make build-libgas MAKEFLAGS=${MF}

mkdir -p ${EILMER_HOME}

cp -r ${TOP}/.gdtk_git/build/* ${EILMER_HOME}/

cp ${TOP}/.gdtk_git/extern/lua-5.4.3/install/bin/* ${EILMER_HOME}/bin
