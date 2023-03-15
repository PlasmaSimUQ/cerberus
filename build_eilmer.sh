#!/usr/bin/env bash
if [ -d "${EILMER_SRC}" ]; then
  git pull
else
  git clone ${EILMER_URL} ${EILMER_SRC}
fi

(cd ${EILMER_SRC}/src/gas && make build-libgas FLAVOUR=${EILMER_FLAVOUR} MAKEFLAGS=${MF})

(cd ${EILMER_SRC}/src/eilmer && make prep-gas prep-chem FLAVOUR=${EILMER_FLAVOUR} MAKEFLAGS=${MF})

mkdir -p ${EILMER_HOME}

cp -r ${EILMER_SRC}/build/* ${EILMER_HOME}/

cp -r ${EILMER_SRC}/extern/lua-5.4.3/install/* ${EILMER_HOME}/

echo "Writing Eilmer environment variables to ./env"
echo "export DGD=${EILMER_HOME}" > env
echo "export DGD_REPO=${EILMER_SRC}" >> env
echo "export PATH=$PATH:${EILMER_HOME}/bin" >> env
echo "export DGD_LUA_PATH=${EILMER_HOME}/lib/?.lua" >> env
echo "export DGD_LUA_CPATH=${EILMER_HOME}/lib/?.so" >> env
