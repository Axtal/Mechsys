#! /bin/bash

if [ "$#" -ne 2 ]; then
	echo
	echo "Usage:"
    echo "        $0 FileKey EraseFiles?"
	echo
	exit 1
fi

files=`find . -name "$1*.pov"`

for i in $files;
do
    povray -V -D +H640 +W640 +FN +I$i -V 2> NULL.txt
done
mencoder "mf://*.png" -mf fps=25 -ovc lavc -lavcopts vcodec=msmpeg4:vbitrate=14400:vhq -o $1.avi
mplayer $1.avi

if [ "$2" -e 1]; then
    rm *png *pov
fi
