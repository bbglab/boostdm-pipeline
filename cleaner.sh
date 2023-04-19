# Clean a few cache and logging items

rm -r work/
rm -r .nextflow/
rm timeline.html*
rm trace.txt*
rm .nextflow.log*

find . -maxdepth 4 -type d -name .idea -exec rm -r {} +
# rm -r .idea/

find . -maxdepth 4 -type d -name __pycache__ -exec rm -r {} +
#rm -r __pycache__/



