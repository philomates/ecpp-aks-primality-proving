#/bin/bash

if [[ $# != 2 ]]; then
  echo "usage: format.sh output numfile"
  exit 1
fi

cat $1 | grep real | awk '{print $2}' > $$.1

for n in $(cat $2); do
  python digits.py $n >> $$.2
done

echo "# Number	Time	Decimal-Digits	Binary-Digits"
paste $2 $$.1 $$.2
rm $$.1 $$.2

