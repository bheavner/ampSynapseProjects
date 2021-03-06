#!/bin/bash

echo ""

# assign inputs
GWAS_DATA=SYounkin_MayoGWAS_09-05-08

for CHR in $(seq 1 22); do

    qsub -S /bin/bash -V -cwd -M james.a.eddy@gmail.com -m abe -j y \
        -N chunkify_chr${CHR} \
        shell/chunkify.sh $GWAS_DATA $CHR ;

done
