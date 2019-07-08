#!/bin/sh
if [ -d /tmp/CoGAPS ]; then
  echo 'Refusing to use existing directory /tmp/CoGAPS for testing' 1>&2
  exit 2
fi
mkdir /tmp/CoGAPS
../../../script/CoGAPS/CoGAPS_wrapper.R -d ./data.RDS -o /tmp/CoGAPS -y RDS
echo $?
rm -rf /tmp/CoGAPS
