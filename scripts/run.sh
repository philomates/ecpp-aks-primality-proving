#!/bin/bash

if [[ $# != 2 ]]; then
  echo "usage: run.sh program numfile"
  exit 1
fi

for n in $(cat $2); do
  echo "Testing $n"
  (time echo $n | $1) 2>&1 
  echo ""
done

