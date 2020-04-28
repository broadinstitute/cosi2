#!/bin/bash

set -e -o pipefail -x

./configure
make
make check
