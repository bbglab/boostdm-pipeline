# uses conda environment boostdm-new-pipeline

cat $1 | grep FAILED | cut -f3 | python _scan_errors.py
