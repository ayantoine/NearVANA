# Installation

1. Copy all script in your git folder (assuming ~/Git)
```
git clone https://github.com/ayantoine/NearVANA.git
```

2. Create NearVANA conda environment
```
conda env create --name NearVANA-env --file ~/Git/NearVANA/Env/NearVANA-env.spec-file.txt
```

3. Download and store database in your database folder (assuming ~/Database)
```
bash ~/Git/NearVANA/BuildLocalDB.sh ~/Database
```

Nb: ~/Database/NuclAccId2Def.tsv and ~/Database/ProtAccId2Def.tsv empty at the beginning. They will be completed as the various analyzes of the workflow progress according to the needs of your data. 

# Configuration
## Cluster configuration file
The cluster conffile contains all the cluster information neeeded for the workflow. It must be defined only once time by cluster and will be the same for all analysis on this.

The file template is avalaible at `~/Git/NearVANA/Template/NearVANA.conf.template`

With:\
GITDIR= path to the NearVANA git folder\
LOCALDB= path to the NearVANA Database folder

NUCACC= path to Viral_nucl_gb.accession2taxid.short.tsv (in ~/Database/[date]/ folder)\
NUCDEF= path to NuclAccId2Def.tsv (in ~/Database/ folder)\
PROACC= path to Viral_prot.accession2taxid.short.tsv (in ~/Database/[date]/ folder)\
PRODEF= path to ProtAccId2Def.tsv (in ~/Database/ folder)\
DBLINEAGE= path to fullnamelineage.dmp (in ~/Database/[date]/ folder)\
VIRMINLEN= path to ViralFamily2MinLen.tsv (in ~/Database/ folder)\
VMR= path to VMR_160919_MSL34.tsv (in ~/Database/ folder)

SCALL= Keyword to submit job on the cluster (for slurm: "sbatch" or "sbatch --partition=dgimi-didi")\
SPARAM= Keyword to define ressource for unitary job (for slurm: "--cpus-per-task=1 --mem=25G")\
MULTICPU= Number of CPU for parallelized jobs (ex: 20)\
MULTIMEMORY= Number of memory for multicpu jobs, in G (ex: 100)\
SPARAM_MULTICPU= Keyword to define a parallelized jobs (for slurm: "--cpus-per-task=${MULTICPU} --mem=${MULTIMEMORY}G")\
STASKARRAY= Keyword to define array of jobs (for slurm: "--array=")\
SMAXTASK= Keyword to define the maximum size of an array of jobs (for slurm: %)\
SRENAME= Keyword to rename job Id (for slurm: "--job-name")

SMAXSIMJOB= Number of simultaneous running jobs allowed on the cluster (ex: 30)\
SMAXARRAYSIZE= Number of maximum size for job array in the cluster (ex: 100)

STASKID= Keyword to refering to job taskId (for slurm: $SLURM_ARRAY_TASK_ID)\
SPSEUDOTASKID= Keyword to refering to internal job taskId (for slurm: %a)

VIRNTDB= Path to the Blast genomic viral refseq database on the cluster\
ALLNTDB= Path to the Blast nt database on the cluster\
VIRPTDB= Path to the Blast protein viral refseq database on the cluster\
ALLPTDB= Path to the Blast nr database on the cluster\
VIRPTDB_DIAMOND= Path to the Diamond protein viral refseq database on the cluster\
ALLPTDB_DIAMOND= Path to the Diamond nr database on the cluster

## Analysis configuration file
The analysis argfile contains all the information for the workflorw about the files path and things to do. It must be specifically defined for each analysis.

The file template is avalaible at `~/Git/NearVANA/Template/NearVANA.arg.template`

With:\
PID= Short ID to identify your analysis. All file generated will be prefixed by this PID\
DATA= Path to the data file\
ADAP= Path to the adaptors file

PAIREND= if TRUE, starting data are pair-end\
METADATA= if TRUE, a metadata file will be used\
MULTIPLEX= if TRUE, data are multiplexed\
UNASSIGNED= if TRUE, in case of multiplex, reads not assigned to a sample will be assigned to a theorical sample (unassigned) and be used with all the reads. If FALSE, unassigned reads will be discarded from the analysis.\
SUBSTRACTION= if TRUE, reads matching a substractive library will be discarded (typically PhiX library)\
SUBS= path to the substractive library if needed

CONF= Path to the cluster configuration file

BLASTX= if TRUE, contigs will be identified with BlastX\
BLASTN= if TRUE, contigs will be identified with BlastN\
DIAMOND= if TRUE, contigs will be identified with Diamond

PREFILTER= if TRUE, contigs will be prefiltered against viral refseq, then positive results will be identified with complete refseq. If FALSE, contigs will be directly tested against the compete refseq.

The adaptors file must contains the adaptors used for the sequencing. The file contains only one line, with the the both adaptors separated by a space. Like:\
AGATCGGAAGAGCAC	AGATCGGAAGAGCGT

## Base data file
The file contains the path to the raw data as sequences, metadata and multiplex key.

The file template is avalaible at `~/Git/NearVANA/Template/NearVANA.data.template`

With:\
PLATE= Id for fastq file (or paire of plate for pair-end)\
ID1= for each id, a path to R1.fastq.gz, R2.fastq.gz (if needed), associated metadata file (if needed), associated demultiplex file (if needed)

#### Metadata file
The metada file is a tabulated file, without header, containing following columns:\
Sample: An arbitrary short sample Id\
Species: The sequenced species\
GPS: The GPS coordinates of the sample\
Location: The location of the sample\
Country: The country of the sample\
Biome: A word describing a caracteristic of the sample\
Date: The date of sampling\
Quantity: The quantity of host sequenced\
Weight: The weight of host sequenced (mg)

Missing information can be let empty. Only the Sample shot Id is necessary

#### Demultiplex file
The demultiplex file is a tabulated file, without header, containing following columns:\
Sample: An arbitrary short sample Id, the same as indicated in metadata\
Tag: The nucleotidic sequence associated to the sample.

WARNING: The short sample Id must be exactly the same as the id indicated in the metadata file.




