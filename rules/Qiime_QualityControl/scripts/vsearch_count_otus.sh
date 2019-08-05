#!/bin/bash
_file={input[samples]}
[ $# -eq 0 ] && {{ echo "Usage: $0 filename"; exit 1; }}
[ ! -f "$_file" ] && {{ echo "Error: $0 file not found."; exit 2; }}

if [ -s "$_file" ]
then
    echo "$_file has some data." && \
    vsearch --usearch_global {input[samples]} \
        -otutabout {output} \
        -id 1 \
        -strand plus \
        --db {input[rep_seq]}
        # do something as file has data
else
    echo "$_file is empty." && \
    touch {input[samples]}
    # do something as file is empty
fi
2> {log}
