#!/usr/bin/env bash

echo "Testing ncmat2endf + njoy + openmc"

#No errors:
set -eux

#Move to empty directory:
test ! -d ./run
mkdir ./run
cd run

#Now do the work, everything should be here:

which njoy
openmc --help
python3 -c 'import openmc'
python3 -c 'import NCrystal'
ncrystal-config -s
ncrystal_ncmat2endf --help
