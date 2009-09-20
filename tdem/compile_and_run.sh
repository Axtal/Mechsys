#!/bin/bash

make test_27cubes
./test_27cubes
povray -I init_test_27cubes.pov
Pov2Avi.sh test_27cubes 1
