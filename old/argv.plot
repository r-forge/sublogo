#!/bin/bash
seqfile=$1
R --slave --args "$seqfile" "$2" "$3" "$4" "$5" "$6" "$7" "$8" < argv.plot.R
convert $seqfile.pdf $seqfile.png
