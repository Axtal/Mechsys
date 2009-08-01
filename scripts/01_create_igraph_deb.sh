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

if [ "$#" -ne 1 ]; then
	echo
	echo "Usage:"
	echo "        $0  ARCH"
	echo
	echo "Where:"
	echo "	  ARCH =  i386  or  amd64"
	echo
	exit 1
fi

VERSION=0.5.2
REV=1
ARCH=$1

set -e

echo
echo "[1;34m########################################## Creating source directory[0m"
echo
rm -rf ./igraph-$VERSION
cp -r ~/pkg/igraph-$VERSION/ ./igraph-$VERSION/

echo "igraph ($VERSION-$REV) unstable; urgency=low"                                     >  ./igraph-$VERSION/debian/changelog
echo ""                                                                                 >> ./igraph-$VERSION/debian/changelog
echo "  * First Debian packaging by Dorival."                                           >> ./igraph-$VERSION/debian/changelog
echo ""                                                                                 >> ./igraph-$VERSION/debian/changelog
echo " -- Dorival Pedroso <dorival.pedroso@gmail.com>  Sat, 01 Aug 2009 22:05:28 +1000" >> ./igraph-$VERSION/debian/changelog

echo "Source: igraph"                                                                  >  ./igraph-$VERSION/debian/control
echo "Section: libs"                                                                   >> ./igraph-$VERSION/debian/control
echo "Priority: optional"                                                              >> ./igraph-$VERSION/debian/control
echo "Maintainer: Dorival Pedroso <dorival.pedroso@gmail.com>"                         >> ./igraph-$VERSION/debian/control
echo "Build-Depends: debhelper (>= 7.0)"                                               >> ./igraph-$VERSION/debian/control
echo "Standards-Version: 3.8.0"                                                        >> ./igraph-$VERSION/debian/control
echo ""                                                                                >> ./igraph-$VERSION/debian/control
echo "Package: libigraph"                                                              >> ./igraph-$VERSION/debian/control
echo "Architecture: any"                                                               >> ./igraph-$VERSION/debian/control
echo "Section: libs"                                                                   >> ./igraph-$VERSION/debian/control
echo "Priority: optional"                                                              >> ./igraph-$VERSION/debian/control
echo "Depends: \${shlibs:Depends}"                                                     >> ./igraph-$VERSION/debian/control
echo "Description: A library for creating and manipulating graphs"                     >> ./igraph-$VERSION/debian/control
echo " igraph is a library for creating and manipulating graphs."                      >> ./igraph-$VERSION/debian/control
echo " It is intended to be as powerful (ie. fast) as possible to enable the"          >> ./igraph-$VERSION/debian/control
echo " analysis of large graphs."                                                      >> ./igraph-$VERSION/debian/control
echo " ."                                                                              >> ./igraph-$VERSION/debian/control
echo " Homepage: http://cneurocvs.rmki.kfki.hu/igraph"                                 >> ./igraph-$VERSION/debian/control
echo ""                                                                                >> ./igraph-$VERSION/debian/control
echo "Package: libigraph-dev"                                                          >> ./igraph-$VERSION/debian/control
echo "Architecture: any"                                                               >> ./igraph-$VERSION/debian/control
echo "Section: libdevel"                                                               >> ./igraph-$VERSION/debian/control
echo "Priority: optional"                                                              >> ./igraph-$VERSION/debian/control
echo "Depends: libc6-dev, libxml2-dev, libigraph (= \${Source-Version})"               >> ./igraph-$VERSION/debian/control
echo "Description: A library for creating and manipulating graphs - development files" >> ./igraph-$VERSION/debian/control
echo " igraph is a library for creating and manipulating graphs."                      >> ./igraph-$VERSION/debian/control
echo " It is intended to be as powerful (ie. fast) as possible to enable the"          >> ./igraph-$VERSION/debian/control
echo " analysis of large graphs."                                                      >> ./igraph-$VERSION/debian/control
echo " ."                                                                              >> ./igraph-$VERSION/debian/control
echo " This package contains the include files and static library for igraph."         >> ./igraph-$VERSION/debian/control
echo " ."                                                                              >> ./igraph-$VERSION/debian/control
echo " Homepage: http://cneurocvs.rmki.kfki.hu/igraph"                                 >> ./igraph-$VERSION/debian/control

tar czvf igraph_$VERSION.orig.tar.gz igraph-$VERSION/


echo
echo "[1;34m########################################## Generating source package[0m"
echo
cd igraph-$VERSION
debuild -sa -S

echo
echo "[1;34m########################################### Verifying source package[0m"
echo
cd ..
lintian -i igraph_$VERSION-$REV.dsc

echo
echo "[1;34m############################################# Generating deb package[0m"
echo
sudo pbuilder build igraph_$VERSION-$REV.dsc
echo

echo
echo "[1;34m############################################# Cleaning up temp files[0m"
echo
rm -f *.build
rm -rf ./igraph-$VERSION/


exit 0
