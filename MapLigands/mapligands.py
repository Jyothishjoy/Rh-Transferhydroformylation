#!/bin/env bash

for lig in `mapLigand.py -ls elements:P,P`; do
    mapLigand.py PP-template.xyz -l 1,2=$lig -o $lig.xyz;
done

