#!/bin/bash
echo > /tmp/screen.log
[ -d /content/drive/MyDrive/tmp ] || mkdir /content/drive/MyDrive/tmp
time python recover_initial_state__embedding.py \
2147483647 \
"257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" \
150 \
27 \
--category 1 \
--level 7 \
--verbose 2 \
--block-size 20 \
--threads 1 \
--sieve \
--workout/dim4free-dec 2
