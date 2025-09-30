#!/bin/bash
# Files from: https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/

# Here is what I did to calculate average local recombination for a specific region (I created a region.bed):

# As a sanity check, you may want to see if the region reported in Hugo’s paper (need to convert to hg38) has a recombination rate close to 0.526
# https://www.nature.com/articles/s41586-020-2818-3#Sec2:~:text=chromosome%C2%A03%20(45%2C859%2C651%E2%80%9345%2C909%2C024%20(hg19))
# Hugo region (hg19): chromosome 3: 45859651 to 45909024 
# Converted to hg38: chromosome 3: 45818159 to 45867532  <------------------

mkdir -p bedGraph

# recomb1000GAvg.bw (based on 1KG)
bigWigToBedGraph data/recomb1000GAvg.bw bedGraph/recomb1000GAvg.bedGraph
bedtools map -a Hugo_region_hg38.bed -b bedGraph/recomb1000GAvg.bedGraph -c 4 -o mean
# chr3    45818159        45867532        0.591216

# recombAvg.bw (based on deCODE)
bigWigToBedGraph data/recombAvg.bw bedGraph/recombAvg.bedGraph
bedtools map -a Hugo_region_hg38.bed -b bedGraph/recombAvg.bedGraph -c 4 -o mean
# chr3    45818159        45867532        0.105390877



