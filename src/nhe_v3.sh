#!/bin/sh
ninja -C ../../../forge/clang-mp-11/Ninja/Release install && ../../../bin/dimer NHE xcell.lua ../../data/nhe1_intake_15mM_v3.txt ../../data/nhe1_delta7_15mM_v3.txt 
