#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o /wynton/group/reiter/lauren/log_autoS                        #-- output directory (fill in)
#$ -e /wynton/group/reiter/lauren/log_autoS                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=10G
#$ -l scratch=10G
#$ -l h_rt=72:00:00


##### LOAD SOFTWARE #####

module load CBI
module load r/4.3

if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

###### ARGUMENTS ######

## fullpath = full path to the folder containing the FASTQs and where the output version analysis directory will be deposited
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

echo >&2 $SCRIPT_DIR
echo >&2 $PROJECT_ROOT

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
