#!/bin/sh

echo ""

# assign inputs
GWAS_DATA=SYounkin_MayoGWAS_09-05-08

for CHR in 22; do
                
    ./prephase.sh $GWAS_DATA $CHR
        
done