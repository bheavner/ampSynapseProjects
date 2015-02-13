# This script downloads the merged count_files produced by running merge_file_counts.R on SNAPR output, transposes the
# dataframe so genes are rows and samples are columns, removes samples to be omitted for QC reasons, and reuploads
# them to synapse with appropriate provenance that this was run.

library(synapseClient)
library(R.utils)
library(edgeR)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

codeFile <- ("https://github.com/bheavner/ampSynapseProjects/blob/master/rnaseqAnalysis/reformat_merged_readcounts.R")

# The files for the first batch are:

# ad_pilot_rnaseq_gene_id_counts.txt.gz ('syn3160436')
# ad_pilot_rnaseq_transcript_id_counts.txt.gz ('syn3160437')

# psp_pilot_rnaseq_gene_id_counts.txt.gz ('syn3160443')
# psp_pilot_rnaseq_transcript_id_counts.txt.gz ('syn3160444')

# mouse_tau_rnaseq_gene_id_counts.txt.gz ('syn3160706')
# mouse_tau_rnaseq_transcript_id_counts.txt.gz ('syn3160709')

adCountFileSynapseIDs <- c('syn3160436', 'syn3160437')
pspCountFileSynapseIDs <- c('syn3160443', 'syn3160444')
mouseCountFileSynapseIDs <- c('syn3160706', 'syn3160709')

# samples to exclude for QC reasons:

adOmit  <- c('5215', 'NA98-300', 'NA02-092')
pspOmit <- c('NA00-002', 'NA00-164', 'NA03-008', 'NA99-042', 'NA04-056', 'NA01-150', 'NA01-058', 'NA00-213', 'NA03-154', 'NA03-291', 'NA04-020')

for (dataset in list(adCountFileSynapseIDs, pspCountFileSynapseIDs, mouseCountFileSynapseIDs)) {
    if (identical(dataset, adCountFileSynapseIDs)) {
        omit <- adOmit
    } else if (identical(dataset, pspCountFileSynapseIDs)) {
        omit <- pspOmit
    } else {
        omit <- NULL
    }

    for (mergedCountFile in dataset) {
        message("Processing SynID: ", mergedCountFile)

        # Download file from Synapse
        originalCountFile <- synGet(mergedCountFile)

        # unzip file and load for processing
        localFilePath <- getFileLocation(originalCountFile)

        if(!file.exists(sub('.gz', '', localFilePath))) {
            gunzip(localFilePath)
        }

        localFilePath <- sub('.gz', '', localFilePath) #trim the .gz suffix

        rawCounts <- read.table(localFilePath, header = TRUE)

        # begin processing - first transpose b/c James did it differently than DGEList expects
        transposedCounts <- t(rawCounts)

        # remove samples that don't pass QC
        keepSamples <- setdiff(colnames(transposedCounts), omit) #probably not the best way..
        transposedCounts <- transposedCounts[, keepSamples]

        # write the data to local dir

        newFileName <- sub('.txt.gz', '', originalCountFile$properties$name)
        newFileName <- paste0(newFileName, "_transposed.txt", sep="")

        write.table(transposedCounts, newFileName, quote = FALSE, sep = "\t", row.names = TRUE)

        # package it up, then create a Synapse object for the output file and upload with provenance

        gzip(newFileName)

        newFileName <- paste0(newFileName, ".gz", sep="")

        parentId <- originalCountFile$properties$parentId

        transposedCountFile <- File(newFileName, parentId = parentId)

        transposedCountFile <- synStore(transposedCountFile,
                                        activityName="Transposed merged readcount file",
                                        used=list(list(name = "reformat_merged_readcounts.R",
                                                       url = codeFile, wasExecuted = T),
                                                  list(entity=originalCountFile,
                                                       wasExecuted=F)))
    }
}