#!/bin/bash

# get location of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd $DIR

wget https://gitlab.com/unigespc/palabos/-/archive/v2.2.0/palabos-v2.2.0.tar.gz
tar -xf palabos-v2.2.0.tar.gz
rm ./palabos-v2.2.0.tar.gz
