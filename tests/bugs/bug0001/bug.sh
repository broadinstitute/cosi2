#!/bin/sh

set -e -o pipefail

./coalescent -p ../../../cosi/tests/bugs/bug0001/test.cosiParams -R ../../../cosi/tests/bugs/bug0001/test.genmap -m -u .001 --output-sim-times --output-end-gens -n 37 --seed 2150193384597205243 -g 1 > /dev/null
