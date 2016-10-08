#!/bin/bash

DICT=$1

cat $DICT | egrep "^@SQ" | cut -f2,3 | sed 's/..://g'

