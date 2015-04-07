#!/usr/bin/env bash
# Run one test case.  This script is generated by $srcdir/make_cosi_tests2.py.
set -e -o pipefail -o nounset
$srcdir/runtest.py --test-num 25 --test-dir $srcdir/tests/model011/t_u_0001 --cosi-sim-params '-u .0001' $@
