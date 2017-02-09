#!/bin/bash
ls | grep -P [0-9]{6} -o | sort -u

# -P ... perl style regex
# [0-9]{6} ... match patterns that have 6 consecutive strings
# -o ... print only matching strings (instead of whole line)
# sort -u ... sort and print only one of each
