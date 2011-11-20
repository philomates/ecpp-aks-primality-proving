#/bin/bash

if [[ $# != 2 ]]; then
  echo "usage: format.sh output numfile"
  exit 1
fi

DIR=`dirname $0`

cat $1 | grep real | sed 's/[ms]/ /g' | awk '{print ($2 * 60 + $3)}' > $$.1

for n in $(cat $2); do
  python $DIR/digits.py $n >> $$.2
done

echo "# Number	Time	Decimal-Digits	Binary-Digits"
paste $2 $$.1 $$.2
rm $$.1 $$.2

