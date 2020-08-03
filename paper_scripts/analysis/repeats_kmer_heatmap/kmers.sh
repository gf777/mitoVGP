#!/bin/bash

for db in meryl/*; do printf "%s\t%s\n" $db $(meryl print $db 2> /dev/null | sort -rnk2 | head -5 | sort -k1 | awk '{printf $1","}' | sed 's/.$//'); done
