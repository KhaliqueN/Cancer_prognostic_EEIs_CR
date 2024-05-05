#! /bin/bash
# Downloaded from https://github.com/rnnh/bioinfo-notebook.git
# Modified by Khalique Newaz 

# sh ./sra2count.sh -a ../data/Homo_sapiens.GRCh38.105_exons.gtf -o ../data/featureCounts -f TCGA-AA-3662-11A-01R-1723-07 --removetemp --verbose -p 5 TCGA-AA-3662-11A-01R-1723-07 
# Help/usage text
usage="$(basename "$0") [options] -a|--annotation <annotation_file> \
	-i|--index <index_file> -o|--outdir <output_directory> <SRR ID(s)> \n
\n
This script generates gene count table(s) using featureCounts.\n
\n
Required arguments: \n
\t      -a | --annotation\t     input genome annotation file \n
\t      -f | --bamfile\t     input bam file name \n
\t      -o | --outdir\t     directory to store counts \n
\t      -l | --pairedFlag\t Flag to indicate whether the bam file is paired or not
\t      SRR ID(s)\t\t           Sequence Read Archive Run ID(s) (SRR...) \n
\n
Optional arguments: \n
\t      -h | --help\t\t         show this help text and exit \n
\t      -p | --processors\t	number (n) of processors to use (default: 1) \n
\t      --fastq-dump\t\t        use 'fastq-dump' instead of the 'fasterq-dump'\n
\t      --verbose\t\t           make output of script more verbose\n
\t	--removetemp\t\t	remove read and alignment files once they are\n
\t	\t\t\t  		no longer needed (minimises disk space needed) \n
\t	--log\t\t\t		redirect terminal output to log file
"

# Setting VERBOSE to 0
# This will be changed to "1" if --verbose is given as an argument,
# resulting in more verbose script output
VERBOSE=0

# Setting REMOVETEMP to 0
# This will be changed to "1" if --removetemp is given as an argument,
# resulting in *.fastq, *.fastq.gz, *.sam, *.bam and *.tsv.summary, being
# removed once they are no longer needed to create a featureCounts table
REMOVETEMP=0

# Setting LOG to 0
# This will be changed to "1" if --log is given as an argument,
# resulting in the terminal output from this script being redirected to a log
# file
LOG=0

# Setting default number of PROCESSORS to use
PROCESSORS=1

# Creating an empty variable for SRRs to be downloaded and aligned to genome
SRRs=""

# Print usage instructions if script is called without any arguments
if [ "$1" = "" ] ; then
  echo -e "ERROR: please provide input files. \n"
  echo -e $usage
  exit 1
fi

# Iterating through the input arguments with a while loop
while (( "$#" )); do
	case "$1" in
		-h|--help)
			echo -e $usage
			exit
			;;
		-a|--annotation)
			ANNOTATION=$2
			shift 2
			;;
		-f|--bamfile)
			BAMFILE=$2
			shift 2
			;;
		-o|--outdir)
			OUTDIR=$2
			shift 2
			;;
		-l|--pairedFlag)
			PAIREDFLAG=$2
			shift 2
			;;
		-p|--processors)
			PROCESSORS=$2
			shift 2
			;;
		--verbose)
			VERBOSE=1
			shift
			;;
		--removetemp)
			REMOVETEMP=1
			shift
			;;
		--log)
			LOG=1
			shift
			;;
		--) # end argument parsing
			shift
			break
			;;
		-*|--*) # unsupported flags
			echo -e "ERROR: $1 is an invalid option. \n" >&2
			echo -e $usage
			exit 1
			;;
		*) # preserve SRR ID(s) as positional arguments
			SRRs="$SRRs $1"
			shift
			;;
	esac
done

if [ $LOG -eq "1" ]
then
	# Redirecting terminal output to log file
	exec 3>&1 4>&2
	trap 'exec 2>&4 1>&3' 0 1 2 3
	exec 1>fd_to_fC_$(date +%Y%m%d_%H%M%S).log 2>&1
fi

if [ ! -d $OUTDIR/ ]; then
  mkdir $OUTDIR
fi

# Beginning the main body of the script

echo -e		~~~~~~~~~~~~~  BAM t o F E A T U R E C O U N T S  ~~~~~~~~~~~~~
echo Script started: $(date)

# Loop through the input SRR IDs
for SRR in $SRRs
do
	printf "\n"
	echo ================================================================================
	echo BAM FILE: $SRR
	sleep 1s
	echo ================================================================================
	sleep 1s

	# b=$(basename $SRR | cut -d"." -f1)

	if [ $VERBOSE -eq "1" ]
	then
        	echo Generating count table using featureCounts...
        	sleep 2s
  fi

  ## -p is to indicate whether paired or not

  if [ $PAIREDFLAG -eq "1" ]
  then
		featureCounts -p -O -t exon -g exon_id -F GTF -T $PROCESSORS -a $ANNOTATION -o $OUTDIR/$BAMFILE.tsv $SRR
  fi

  if [ $PAIREDFLAG -eq "0" ]
  then
		featureCounts -O -t exon -g exon_id -F GTF -T $PROCESSORS -a $ANNOTATION -o $OUTDIR/$BAMFILE.tsv $SRR
  fi

	if [ $VERBOSE -eq "1" ]
	then
        	echo Results written to $BAMFILE.tsv
        	sleep 2s

        	echo Head of $BAMFILE.tsv
        	sleep 2s
        	head $OUTDIR/$BAMFILE.tsv
        	sleep 2s

        	echo Tail of $BAMFILE.tsv
        	sleep 2s
        	tail $OUTDIR/$BAMFILE.tsv
        	sleep 2s
  fi


	if [ $REMOVETEMP -eq "1" ]
	then
		echo Removing temporary files...
		rm $OUTDIR/*.tsv.summary
	fi

done

echo Script finished: $(date)


