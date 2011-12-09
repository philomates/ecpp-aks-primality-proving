#!/bin/bash

if [[ $# != 1 ]]; then
  echo "usage: format.sh outfile"
  exit 1
fi

DIR=`dirname $0`

cat $1 | grep Testing | awk '{print $2}' > $$.1

cat $1 | grep real | sed 's/[ms]/ /g' | awk '{print ($2 * 60 + $3)}' > $$.2

for n in $(cat $2); do
  $DIR/digits.py $n >> $$.3
done

echo "# Number	Time	Decimal-Digits	Binary-Digits"
paste $$.1 $$.2 $$.3
rm $$.1 $$.2 $$.3

