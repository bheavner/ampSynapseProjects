#!/bin/bash

# assign inputs
# GWAS_DATA="$1"
# CHR="$2"
GWAS_DATA=SYounkin_MayoGWAS_09-05-08
CHR=22

# directories
ROOT_DIR=./
DATA_DIR=${ROOT_DIR}data/
GWAS_DIR=${DATA_DIR}gwas_results/
GENMAP_DIR=${DATA_DIR}haplotypes/
GWAS_HAP_DIR=${GENMAP_DIR}${GWAS_DATA}.phased/
GWAS_IMP_DIR=${GWAS_DIR}${GWAS_DATA}.imputed/

# executables
GTOOL_EXEC=${ROOT_DIR}resources/gtool
QCTOOL_EXEC=${ROOT_DIR}resources/qctool/qctool

# specify additional data files
SAMPLE_FILE=${GWAS_HAP_DIR}${GWAS_DATA}.chr${CHR}.phased.sample

# get list of imputed genotype files for chromosome
CHUNK_LIST=$(ls -d -1 ${GWAS_IMP_DIR}*.* | grep "chr${CHR}.*.imputed$")

# merge all imputed genotype files for chromosome
GEN_FILE="${GWAS_IMP_DIR}${GWAS_DATA}.chr${CHR}.imputed.gen"

echo "Merging all chunk files for chromsome ${CHR}..."
echo
cat $CHUNK_LIST > $GEN_FILE

# create new directory to store results
RESULTS_DIR=${GWAS_DIR}${GWAS_DATA}.imputed.qc/

if [ ! -e "$RESULTS_DIR" ]; then
	mkdir "$RESULTS_DIR"
fi

# perform QC with qctool
QC_FILE=${RESULTS_DIR}${GWAS_DATA}.chr${CHR}.imputed.qc.gen

echo "Performing QC on merged file..."
echo
time $QCTOOL_EXEC -g $GEN_FILE -og $QC_FILE \
	-snp-missing-rate 0.05 -maf 0 1 -info 0.4 1 -hwe 20

# qctool adds an extra columns of NA for some reason; need to remove
TMP_FILE=`mktemp qcgen.XXX`

echo "Removing extraneous first coumn from QC output..."
echo
cut -d " " -f 2- $QC_FILE > $TMP_FILE
mv $TMP_FILE $QC_FILE

# convert to ped/map with gtool
PED_FILE="${RESULTS_DIR}${GWAS_DATA}.chr${CHR}.imputed.ped"
MAP_FILE="${RESULTS_DIR}${GWAS_DATA}.chr${CHR}.imputed.map"

echo "$GTOOL_EXEC -G --g $QC_FILE --s $SAMPLE_FILE --ped $PED_FILE --map $MAP_FILE --phenotype plink_pheno"
echo

echo "Converting gen/sample format to ped/map..."
echo
time $GTOOL_EXEC -G --g $QC_FILE --s $SAMPLE_FILE \
 	--ped $PED_FILE --map $MAP_FILE \
 	--phenotype plink_pheno --chr ${CHR}

# gtool swaps the first 2 columns of the ped file; switch back
TMP_FILE=`mktemp ped.XXX`

echo "Swaping first two columns of gtool generated ped file..."
echo
paste <(cut -f 1-2 $PED_FILE | awk '{OFS = "\t"; t = $1; $1 = $2; $2 = t; print}') <(cut -f 3- $PED_FILE) > $TMP_FILE
mv $TMP_FILE $PED_FILE