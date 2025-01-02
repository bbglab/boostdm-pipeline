# uses conda environment boostdm-new-pipeline

# usage example: 

# bash scan_errors.sh trace.txt

cat $1 | grep FAILED | cut -f3 | python _scan_errors.py
