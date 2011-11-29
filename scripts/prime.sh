#!/bin/bash

if [[ $# != 1 ]]; then
  echo "usage: prime.sh digits"
  exit 1
fi

GPRIME=../gprime

if [[ ! -e $GPRIME ]]; then
  echo "make gprime has not been run"
  exit 1
fi

NUM=1
for i in $(seq 2 $1); do
  NUM=${NUM}0
done

$GPRIME $NUM ${NUM}0

