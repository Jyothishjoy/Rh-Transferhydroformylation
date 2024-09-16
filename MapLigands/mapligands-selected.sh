#!/bin/env bash

ligand_names=("AcridinePhos" "AnthaPhos" "Benzoxantphos" "Benzyllnixantphos" "DBFphos" "DPEPhos" "Isopropxantphos" "Nixantphos" "PhosxantPhos" "SixantPhos" "Thixantphos" "TransPhos" "Xanthene-XantPhos" "Xantphos")  

for lig in "${ligand_names[@]}"; do
    mapLigand.py PP-template.xyz -l 1,2="$lig" -o "$lig".xyz
done

