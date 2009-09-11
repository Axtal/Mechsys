#!/bin/bash

make test_27cubes
./test_27cubes
povray -I test_27cubes.pov
qiv test_27cubes.png
