#!/usr/bin/env bash
# Run one test case.  This script is generated by $srcdir/make_cosi_tests2.py.
set -e -o pipefail -o nounset
$srcdir/tests/runtest.py --test-num 7 --test-dir $srcdir/tests/model003/t  $@
