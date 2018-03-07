#!/bin/bash
egrep "^@SQ" | cut -f-3 | sort  | md5sum - | awk '{print $1}'
