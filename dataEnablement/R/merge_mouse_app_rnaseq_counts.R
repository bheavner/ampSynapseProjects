
# This script uses the merging function saved in merge_file_counts.R to combine
# SNAPR count files for reprocessed pilot APP mouse RNAseq data from Mayo

library(synapseClient)
library(tools)
source("dataEnablement/R/merge_count_files.R")

# point to specific version run
codeFile <- ("https://github.com/bheavner/ampSynapseProjects/blob/eb198e01b658c3ee2d58be414a23504bedbb0c31/dataEnablement/R/merge_mouse_app_rnaseq_counts.R")

# Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# Define paths for required Synapse objects
app_rnaseq_counts <- "syn2875347" # all count files in a zipped directory

# Download files from Synapse
app_count_files <- synGet(app_rnaseq_counts)
app_files_path <- getFileLocation(app_count_files)

# Get name of temporary directory to store unzipped files (same as name of
# original compressed directory)
fileDir <- file_path_sans_ext(basename(app_files_path))
tmpDir <- tempdir()
unzip(app_files_path, exdir = tmpDir)

inputDir <- file.path(tmpDir, fileDir)
prefix <- "mouse_app_rnaseq"

countTypes <- c("gene_id", "transcript_id")

for (countType in countTypes) {
    message(paste("Merging", prefix, "files of count type", countType, "..."))

    # Create the merged file and store the output file path
    merged_file <- create_merged_file(inputDir, countType, prefix)

    # Create a Synapse object for the output file and upload
    merged_file_object <- File(path = merged_file,
                               parentId = app_count_files$properties$parentId)
    merged_file_object <- synStore(merged_file_object,
                                   activityName="Build merged readcount file",
                                   used=list(list(name = "merge_mouse_app_rnaseq_counts.R",
                                                  url = codeFile, wasExecuted = T),
                                             list(entity=originalCountFile,
                                                  wasExecuted=F)))
}

unlink(inputDir, recursive = TRUE)
