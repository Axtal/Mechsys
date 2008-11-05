#!/bin/bash

# SETTING PBUILDER
# ================
#   sudo apt-get install pbuilder fakeroot devscripts dh-make lintian
#   sudo pbuilder create
#   gvim ~/.pbuilderrc
#      COMPONENTS="main universe multiverse"
#      OTHERMIRROR="deb http://cneurocvs.rmki.kfki.hu /packages/binary/"
#   sudo pbuilder update --override-config
#   gpg --gen-key

if [ "$#" -ne 2 ]; then
	echo
	echo "Usage:"
	echo "        [1;34msh $0[0m [1;31mVERSION REV[0m"
	echo
	echo " where [1;31mVERSION=1.0[0m and [1;31mREV=1[0m, for example"
	echo
	exit 1
fi

VERSION=$1
REV=$2

set -e

echo
echo "[1;34m########################################## Creating source directory[0m"
echo
rm -rf mechsys-$VERSION
cp -r ~/mechsys/ ./mechsys-$VERSION/
cd mechsys-$VERSION
sh Reset_cmake.sh
rm -rf .hg/
rm .hgignore
find . -iname "*.swp" -exec rm {} \;
find . -iname "*.vtu" -exec rm {} \;
find . -iname "*.vtk" -exec rm {} \;
cd ..
tar czvf mechsys_$VERSION.orig.tar.gz mechsys-$VERSION/

echo
echo "[1;34m########################################## Generating source package[0m"
echo
cd mechsys-$VERSION
debuild -sa -S

echo
echo "[1;34m########################################### Verifying source package[0m"
echo
cd ..
lintian -i mechsys_$VERSION-$REV.dsc

echo
echo "[1;34m############################################# Generating deb package[0m"
echo
sudo pbuilder build mechsys_$VERSION-$REV.dsc
echo

echo
echo "[1;34m############################################# Cleaning up temp files[0m"
echo
rm -f *.build
rm -rf ./mechsys-$VERSION/

exit 0
