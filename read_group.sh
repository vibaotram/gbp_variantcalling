#!/bin/bash

### create sam header @RG from this format @E100033056L1C001R00200000130/1

# f1=$(echo $1 | cut -d' ' -f1)

if [[ $1 == *".gz" ]];
then
    header=$(zcat $1 | head -n 1)
else
    header=$(cat $1 | head -n 1)
fi

if [[ -z $(echo $header | grep : )]]
then
    id=$(echo $header | cut -f1 -d"L" | sed 's/@//')
    sm=$(basename $(dirname $1))
    l=$(echo $header | cut -f2 -d"L" | cut -f1 -d"C")
    echo "@RG\tID:$id.$l\tSM:$sm\tPL:Illumina"
else
    id=$(echo $header | cut -f1 -d":" | sed 's/@//')
    sm=$(basename $(dirname $1))
    l=$(echo $header | cut -f2 -d":")
    echo "@RG\tID:$id.$l\tSM:$sm\tPL:Illumina"
