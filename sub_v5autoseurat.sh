#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o /dev/null   # Prevent default SGE output
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=10G
#$ -l scratch=10G
#$ -l h_rt=72:00:00


##### LOAD SOFTWARE #####

module load CBI
module load r/4.4.3-gcc13

if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

###### ARGUMENTS ######

## fullpath = full path to the folder containing the cellranger outs folder and where the output version analysis directory will be deposited
## version = what version is this? Analysis folder will be named after this variable
## do_marks = boolean (TRUE or FALSE) indicating if you want FindAllMarkers to be run.
## soupx = boolean (TRUE or FALSE) indicating if you want to run off the soupx corrected matrix. Assumes there is a directory named "soupx_results" in the fullpath folder containing soupx output.

fullpath=$1
basename="${fullpath##*/}" 
version=$2
para="$basename"_"$2"_parameters.csv
do_marks=$3
soupx=$4

dir="${fullpath%/*}"

PROJECT_ROOT="$SGE_O_WORKDIR"
script="${PROJECT_ROOT}/seurat_processing.R"

## LOG directory
LOG_DIR=$fullpath/log
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/autoseuratV5_${JOB_ID}.log"

# Redirect both stdout and stderr to the log file
exec > "$LOG_FILE" 2>&1

echo >&2 "Input/Output Direcotry ${fullpath}"
echo >&2 "autoSeuratV5 Direcotry: ${PROJECT_ROOT}"
echo >&2 "LOG Directory: ${LOG_FILE}"
echo >&2 "Project name: ${basename}"
echo >&2 "Version: ${version}"
echo >&2 "Parameters file: ${para}"
echo >&2 "Run FindAllMarkers? : ${do_marks}"
echo >&2 "Run from Soupx corrected matrix? : ${soupx}"

## move to temp dir and create new folder, copy needed files to TMPDIR

cd $TMPDIR
mkdir $basename


cp -r $dir/$basename/outs/filtered_feature_bc_matrix $basename/
cp $dir/$basename/$para $basename/

if [[ $soupx ]]; then
	gzip $dir/$basename/soupx_results/matrix.mtx
	gzip $dir/$basename/soupx_results/barcodes.tsv
	gzip $dir/$basename/soupx_results/features.tsv
	cp -r $dir/$basename/soupx_results $basename/
fi

Rscript $script $basename $version $para $PROJECT_ROOT $do_marks $soupx

cp -r $basename/$version $dir/$basename/

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
