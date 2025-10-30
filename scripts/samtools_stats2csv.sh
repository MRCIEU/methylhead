#!/bin/bash

cat - | grep -e '^SN' \
      | perl -pe 's/^SN[^a-z]+//g' \
      | perl -pe 's/([^:]+):[^0-9]+([0-9a-z+-]+).*/\1,\2/'

