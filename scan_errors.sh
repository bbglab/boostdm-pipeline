# uses conda environment boostdm-new-pipeline

cat trace-20240626-38172655.txt | grep FAILED | cut -f3 | python _scan_errors.py
