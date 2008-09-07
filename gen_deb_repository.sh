#!/bin/bash

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

if [ -f mechsys_$VERSION.dsc ]; then
	echo
	echo "[1;34m############################################## Setting up repository[0m"
	echo
	test -d binary || mkdir binary
	test -d sources || mkdir sources
	if [ -f /var/cache/pbuilder/result/mechsys_"$VERSION"_i386.deb ]; then
		cp /var/cache/pbuilder/result/mechsys_"$VERSION"_i386.deb binary/
		mv mechsys_$VERSION.dsc sources/
		mv *.tar.gz sources/
		# binary
		dpkg-scanpackages binary \
		 | sed 's@Filename: binary/@Filename: releases/mechsys/binary/@' \
		 | sed 's@libboost-python1.34.1 (>= 1.34.1-2.1), @@' \
		 | sed 's@liblapack3gf | liblapack.so.3gf | libatlas3gf-base, @@' \
		 | sed 's@python (<< 2.6), @@' \
		 | gzip -9c > binary/Packages.gz
		# sources
		dpkg-scansources sources \
		 | sed 's@Directory: sources@Directory: releases/mechsys/sources@' \
		 | gzip -9c > sources/Packages.gz
		# upload
		echo
		echo "[1;34m################################################ Uploading to server[0m"
		echo
		echo "rsync -avz --delete ~/pkg/debs/ USERNAME@dl.sv.nongnu.org:/releases/mechsys"
		echo
	else
		echo
		echo "ERROR: file [1;31m/var/cache/pbuilder/result/mechsys_"$VERSION"_i386.deb[0m was not found"
		echo
		exit 1
	fi
else
	echo
	echo "ERROR: file [1;31mmechsys_$VERSION.dsc[0m was not found"
	echo
	exit 1
fi
