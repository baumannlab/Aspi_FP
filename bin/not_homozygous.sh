#!/bin/bash
for file in ../data/pysam/*.prof;
do
    awk '$8!="Homozygous"' $file > ${file}.not_hom;
done
