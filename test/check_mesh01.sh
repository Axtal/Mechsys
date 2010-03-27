#!/bin/bash

VTUS="mesh01_1tet mesh01_hex_box mesh01_quad mesh01_quad_ring mesh01_tet_box mesh01_tet_hole mesh01_tri"

for vtu in $VTUS; do
    diff $vtu.vtu "$vtu"_py.vtu
done
