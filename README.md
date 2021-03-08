# NearVANA
Analyze Illumina sequencing data and identify virus sequence (using BlastN, BlastX or Diamond)

# Installation
- Clone git repository
- Create a dedicated environment conda
conda env create --name NearVANA-env --file NearVANA-env.spec-file.txt
(<NearVANA-env> can be replaced by any name)

# To do only once I - Taxonomy and definition file
NearVANA need some files from NCBI to determine the taxonomy of a hit. This files are named :
- fullnamelineage.dmp
- nucl_gb.accession2taxid
- prot.accession2taxid

They are avalaible on the NCBI ftp. There is a git programm that automatize their download (https://github.com/ayantoine/NCBI-ViralTaxo)

In addition, NearVANA need to acces two file that contains hit definitions.
- NuclAccId2Def.tsv
- ProtAccId2Def.tsv

These file can be empty at begining, they will be completed progressively among analyses.


# To do only once II - Cluster conf file
Cluster conffile contains data to launch job on your cluster. The file must contains the following informations

SDIR : NearVANA script directory
NUCACC: Path to nucl_gb.accession2taxid
NUCDEF: Path to to file with definition for nuc
PROACC: Path to prot.accession2taxid
PRODEF: Path to to file with definition for prot
DBLINEAGE: Path to fullnamelineage.dmp
VIRMINLEN: Path to ViralFamily2MinLen.tsv #cf Prerequisite III

SCALL : keyword to launch a job
SPARAM : list of param to add to a job
MULTICPU : Number of CPU for parallelized task
MULTIMEMORY : Quantity of memory for biggest task
STASKARRAY : keyword to define a array of task
SMAXTASK : keyword to define the maximum amount of parallele task during a task array
SRENAME : keyword to rename a job

SMAXSIMJOB : integer that define the maximum amount of parallele task during a task array
SMAXARRAYSIZE : integer that define the maximum size of array. Use 0 if array haven't limited size on your cluster

STASKID : keyword to ask a task Id during in array
SPSEUDOTASKID : keyword to ask a task Id during the call of an array

VIRNTDB : path to the blast viral database nucleotidic
ALLNTDB : path to the blast database nucelotidic
VIRPTDB : path to the blast viral database proteic
ALLPTDB : path to the blast database proteic
VIRPTDB_DIAMOND : path to the diamond viral database proteic
ALLPTDB_DIAMOND : path to the diamond database proteic

Below is an exemple of Cluster keyword for an SGE cluster. Beware of the space after value of STASKARRAY and SMAXTASK.
SCALL=qsub
SPARAM=-q\ long.q\ -V
STASKARRAY=-t\ 
SMAXTASK=\ -tc\ 
SRENAME=-N

STASKID=$SGE_TASK_ID
SPSEUDOTASKID=\$TASK_ID

Below is an exemple of Cluster keyword for an SLURM cluster.
SCALL=sbatch
SPARAM="--cpus-per-task=1 --mem=25G"
STASKARRAY="--array="
SMAXTASK=%
SRENAME=--job-name

STASKID=$SLURM_ARRAY_TASK_ID
SPSEUDOTASKID=%a


# To do only once III - Reference size for family virus
NearVANA compare contigs size to reference genome to help users to determine interesting fragment. The file must lead the following organization:

Famility MinSize

Data are facultative but, even empty, the file must exist.

# To do for each analysis IV - Analysis arg file
NearVANA conffile contains data that are specific for 1 analysis. The file must contains the following informations

#Warning! Do not change variable name
#Project id
PID=
#Data file
DATA=
#Adaptator file
ADAP=

#Pair-end input (TRUE for pair-end, FALSE for mono fastq) [beware of uppercase!]
PAIREND=
#Metada (TRUE for multiplex, FALSE without multiplex (RNAseq)) [beware of uppercase!]
METADATA=
#Multiplex (TRUE for multiplex, FALSE without multiplex (RNAseq)) [beware of uppercase!]
MULTIPLEX=
#Use unassigned reads into assembly, (TRUE to use unassigned, FALSE to reject unassigned) [beware of uppercase!]
UNASSIGNED=
#Use substraction (TRUE to use unassigned, FALSE to reject unassigned) [beware of uppercase!]
SUBSTRACTION=
#Substractive library (PhiX)
SUBS=

#Cluster Conffile
CONF=

#Use BlastX (TRUE/FALSE, beware uppercase)
BLASTX=
#Use BlastN (TRUE/FALSE, beware uppercase)
BLASTN=
#Use Diamond (TRUE/FALSE, beware uppercase)
DIAMOND=

# To do for each analysis V - Data file
NearVANA conffile contains data that are specific for 1 analysis. The file must contains the following informations

#Liste of plateId, between "()", separed by " "
#Ex: PLATE=(P1 P2 ...)
PLATE=( P1 )
#Details of each plateId, defined as plateId=(R1.fq R2.fq metadata.tsv dodeca.tsv)
#Ex: P1=(R1-001.fq R2-001.fq dodeca-001.tsv metadata-001.tsv)
#Ex: P2=(R1-002.fq R2-002.fq dodeca-002.tsv metadata-002.tsv)
#If fastq aren't pair-end, use only one file and specify FALSE for PAIREND value in Arg file
#If data aren't multiplexed, forget the file and specify FALSE for MULTIPLEX value in Arg file
#If data haven't metada, forget the file and specify FALSE for METADATA value in Arg file
#WARNING: very important: the line must begin with the same value as writed in PLATE content
P1=()

# Prerequisite V - Data files
1) Sequences files
Pair of read, in fastq format, in two separate ".gz" archive

2) Metadata file (facultative)
A tsv file that contains metadata of samples, with columns as following :

PLATE_ID HOST_NAME COORD LOCATION COUNTRY ECO DATE QUANTITY WEIGTH

With:
PLATE_ID: Unique identifier for the sample (ex: S01). No space, no special characters
HOST_NAME: Name of the species sampled
COORD: GPS coordinate (facultative)
LOCATION: City (facultative)
COUNTRY: Country (facultative)
ECO: Ecosystem specification (facultative)
DATE: Date of sampling (facultative)
QUANTITY: Quantity of host (facultative)
WEIGTH: Weight of the sample, in mg (facultative)

Example:
Helicoverpa_1	Helicoverpa armigera	43.705225,3.859843	Les Matelles	France	serre	03/09/2018	17 59

3) Adaptators file
A tsv file that contains the adaptators used during the sequencing

ADAPTATOR1 ADAPTATOR2

Example:
AGATCGGAAGAGCAC	AGATCGGAAGAGCGT

4) Dodeca file
A tsv file that contains data to link SampleId to their DODECA (the nucelotidic tag for the specified SampleId). 
The file is organized as following:

SAMPLE_ID DODECA

Example:
D01	AACAAGACGT	AACAACCTCGGCGG
D02	AACACACTCA	AACCGAGTCGCGAG

5) PhiX.fa (or whatever sequence that must be substracted from the sequencing data)
The fasta file of the element that need to be substracted from the sequencing data. The fasta can contains multiple sequences.

# Output
All files are prefixed with the Prefix of the analysis.
Demultiplexing_Hyper_Distribution.tab : distribution of assigned reads among sample
Demultiplexing_Hyper.tab : tsv file with each read identifier linked to a sample
BlastN_result.xlsx : tsv file with all blastN identification results
BlastX_result.xlsx : tsv file with all blastX identification results
BlastAll_result.xlsx : tsv file with all blastX and BlastN identification results
