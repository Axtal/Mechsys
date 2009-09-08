#! /bin/bash
for i in *pov;
do
    povray -V -D +H640 +W640 +FN +I$i -V 2> NULL.txt
done
mencoder "mf://*.png" -mf fps=25 -ovc lavc -lavcopts vcodec=msmpeg4:vbitrate=14400:vhq -o $1
mplayer $1
rm *png *pov
