# uses conda environment boostdm-new-pipeline

cat trace.txt | grep FAILED | cut -f3 | python _scan_errors.py
