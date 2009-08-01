#!/bin/bash

set -e

if [ "$#" -ne 1 ]; then
	echo
	echo "Usage:"
	echo "        $0  FLAG"
	echo
	echo "FLAG:"
	echo " 0 => only IGraph"
	echo " 1 => only MechSys"
	echo
	exit 1
fi

test -d binary  || mkdir binary
test -d sources || mkdir sources

if [ "$1" -eq 0 ]; then
	echo "[1;34m############################################## IGraph ####[0m"
	cp /var/cache/pbuilder/result/libigraph*.deb binary/
	mv igraph_*.diff.gz     sources/
	mv igraph_*.dsc         sources/
	mv igraph_*.changes     sources/
	mv igraph_*.orig.tar.gz sources/
fi

if [ "$1" -eq 1 ]; then
	echo "[1;34m############################################## MechSys ####[0m"
	cp /var/cache/pbuilder/result/mechsys*.deb binary/
	mv mechsys_*.diff.gz     sources/
	mv mechsys_*.dsc         sources/
	mv mechsys_*.changes     sources/
	mv mechsys_*.orig.tar.gz sources/
fi

# binary
dpkg-scanpackages binary /dev/null \
 | sed 's@Filename: binary/@Filename: /software/binary/@' \
 | gzip -9c > binary/Packages.gz

# sources
dpkg-scansources sources /dev/null \
 | sed 's@Directory: sources@Directory: /software/sources@' \
 | gzip -9c > sources/Sources.gz

exit 0
