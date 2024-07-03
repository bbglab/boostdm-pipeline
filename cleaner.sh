# Clean a few cache and logging items

rm -r work/
rm -r .nextflow/
rm timeline-*
rm trace-*
rm .nextflow.log*

find . -maxdepth 4 -type d -name .idea -exec rm -r {} +

find . -maxdepth 4 -type d -name __pycache__ -exec rm -r {} +



