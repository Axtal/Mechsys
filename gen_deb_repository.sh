#!/bin/bash

if [ "$#" -ne 3 ]; then
	echo
	echo "Usage:"
	echo "        [1;34msh $0[0m [1;31mVERSION REV USER[0m"
	echo
	echo " where [1;31mVERSION=1.0[0m and [1;31mREV=1[0m, for example"
	echo
	exit 1
fi

VERSION=$1
REV=$2
USER=$3

set -e

if [ -f mechsys_$VERSION-$REV.dsc ]; then
	echo
	echo "[1;34m############################################## Setting up repository[0m"
	echo
	test -d binary || mkdir binary
	test -d sources || mkdir sources
	if [ -f /var/cache/pbuilder/result/mechsys_"$VERSION-$REV"_i386.deb ]; then
		cp /var/cache/pbuilder/result/mechsys_"$VERSION-$REV"_i386.deb binary/
		mv mechsys_$VERSION-$REV.dsc sources/
		mv mechsys_$VERSION-$REV.diff.gz sources/
		mv mechsys_$VERSION.orig.tar.gz sources/
		mv mechsys_"$VERSION-$REV"_source.changes sources/
		# binary
		dpkg-scanpackages binary /dev/null \
		 | sed 's@Filename: binary/@Filename: /releases/mechsys/binary/@' \
		 | gzip -9c > binary/Packages.gz
		# sources
		dpkg-scansources sources /dev/null \
		 | sed 's@Directory: sources@Directory: /releases/mechsys/sources@' \
		 | gzip -9c > sources/Sources.gz
		# upload
		echo
		echo "[1;34m################################################ Uploading to server[0m"
		echo
		rsync -avz --delete ~/pkg/debs/binary/ $USER@download.savannah.nongnu.org:/releases/mechsys/binary/
		rsync -avz --delete ~/pkg/debs/sources/ $USER@download.savannah.nongnu.org:/releases/mechsys/sources/
		echo
	else
		echo
		echo "ERROR: file [1;31m/var/cache/pbuilder/result/mechsys_"$VERSION-$REV"_i386.deb[0m was not found"
		echo
		exit 1
	fi
else
	echo
	echo "ERROR: file [1;31mmechsys_$VERSION-$REV.dsc[0m was not found"
	echo
	exit 1
fi
