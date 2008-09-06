#!/bin/bash

# SETTING PBUILDER
# ================
#   sudo apt-get install pbuilder fakeroot devscripts dh-make lintian
#   sudo pbuilder create
#   gvim ~/.pbuilderrc
#      COMPONENTS="main universe"
#      OTHERMIRROR="deb http://cneurocvs.rmki.kfki.hu /packages/binary/"
#   sudo pbuilder update --override-config

if [ "$#" -ne 1 ]; then
	echo
	echo "Usage:"
	echo "        [1;34msh $0[0m [1;31mVERSION[0m"
	echo
	echo " where VERSION = 0.1, for example"
	echo
	exit 1
fi

VERSION=$1

set -e

echo
echo "[1;34m########################################## Creating source directory[0m"
echo
rm -rf mechsys_$VERSION
cp -r ~/mechsys/ ./mechsys_$VERSION/
cd mechsys_$VERSION
sh Reset_cmake.sh
rm -rf .hg/
rm .hgignore
find . -iname "*.swp" -exec rm {} \;
find . -iname "*.vtu" -exec rm {} \;
find . -iname "*.vtk" -exec rm {} \;

echo
echo "[1;34m########################################## Generating source package[0m"
echo
debuild -S

echo
echo "[1;34m########################################### Verifying source package[0m"
echo
cd ..
lintian -i mechsys_$VERSION.dsc

echo
echo "[1;34m############################################# Generating deb package[0m"
echo
sudo pbuilder build mechsys_$VERSION.dsc
cp /var/cache/pbuilder/result/mechsys_"$VERSION"_i386.deb .
echo
