#!/usr/bin/env bash
# Adjusts the read-ends in a read BED by Tn5 offsets
awk_prog='
BEGIN {OFS = FS}
{
  if ($6 == "+") {
      $2 = $2 + 4
  } else if ($6 == "-") {
      $3 = $3 - 5
  }
  print $0
}
'

awk -F $'\t' "$awk_prog" $1
