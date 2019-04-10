
#!/bin/bash

paste -d "\n" ./data/*.txt > ./data/resultsFile.csv;

rm ./data/*.txt;

echo "All files merged into resultsFile.csv"

