#!/bin/bash

# Enter the working directory
dir=/data/data2/jab112/2014_mnase_manuscript/scripts/2014_belsky_genes_and_dev/
cd $dir

# Compile the script if needed
javac -d bin/ -sourcepath src/ -cp lib/sam-1.67.jar src/mnasedens/GetDensitySignalAroundFeature.java

# Run the program
java -cp bin:lib/sam-1.67.jar mnasedens.GetDensitySignalAroundFeature
