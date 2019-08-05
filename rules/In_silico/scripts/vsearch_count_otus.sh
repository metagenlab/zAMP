#!/bin/bash
_file=$0
_db=$1
_out=$2

if [ -s "$_file" ]
then
    echo "$_file has some data." && \
    vsearch --usearch_global $0 \
        -otutabout $2 \
        -id 1 \
        -strand plus \
        --db $1
        # do something as file has data
else
    echo "$_file is empty." && \
    touch $0
    # do something as file is empty
fi
2> {log}
