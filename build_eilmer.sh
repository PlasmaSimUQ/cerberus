#!/usr/bin/env bash

git submodule update --init --recursive -- ${EILMER_SRC}

cd ${EILMER_SRC}/src/eilmer
make prep-chem prep-gas MAKEFLAGS=${MF}

cd ${EILMER_SRC}/src/gas
make build-libgas MAKEFLAGS=${MF}

mkdir -p ${EILMER_HOME}

cp -r ${EILMER_SRC}/build/* ${EILMER_HOME}/

cp -r ${EILMER_SRC}/extern/lua-5.4.3/install/* ${EILMER_HOME}/
