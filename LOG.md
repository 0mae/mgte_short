# mgte_short


# Set environments

## Set variables and functions

``` {sh}
#| label: set_variables_and_functions
# For PBS
Email=XXX
# Project home directory
Home=$PWD
# print columns
function pri_col () {
    head -n1 $1 | sed 's/\t/\n/g' | awk '{print NR,$1}'
}
# check table header
function hd_check () {
    head $1 | column -t -s "$(printf '\011')"
}
```

## RStudio Server

- Docker image `bio_4.3.2` was used for R programming.

``` {sh}
#| label: build_container
User=XXX
docker image build -t $User/bio_4.3.2 .
```

``` {sh}
#| label: pull_container
# For local
docker pull $User/bio_4.3.2
# For fe1
singularity pull docker://$User/bio_4.3.2
```

- Run RStudio Server
  - http://localhost:8787

``` {sh}
#| label: run_rstudio_server
User=XXX
Pass=XXX
docker container run -p 8787:8787 -v ${PWD}:/home/rstudio -e PASSWORD=$Pass $User/bio_4.3.2
```

## Supercomputer systems

- SuperComputer System (CB202), Institue for Chemical Research, Kyoto
  University
  - https://www.scl.kyoto-u.ac.jp/index_e.html

``` {sh}
#| label: fe1
# For set environment
bash
module load hmmer/3.4
```

- Supercomputer HOKUSAI Bigwaterfall2(HBW2) System, RIKEN
  - https://i.riken.jp/en/supercom/

``` {sh}
#| label: HBW2
# For set environment
project_acc=XXX
```

## Rendering

- `Ctrl + Shift + K` at VSCode

- Or, run following command

``` {sh}
#| label: quarto_render
quarto render LOG.qmd
```

# GTDB r220

# Download databases

## Pfam_A

- RELEASE 37.0

``` {sh}
#| label: download_pfam_a_hmm
mkdir -p db/pfam

# resources_used.vmem=736392kb;resources_used.walltime=00:19:16
qsub -q SMALL -m abe -M $Email -l select=1:ncpus=1:mem=48gb -l walltime=12:00:00 -e qsub_out/download_pfam_a_hmm_e -o qsub_out/download_pfam_a_hmm_o scripts/qsub/download_pfam_a_hmm.sh
{
source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR
cd db/pfam
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/userman.txt
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
}

# Compress PBS output
gzip qsub_out/download_pfam_a_hmm_{e,o}

# release note
cat db/pfam/relnotes.txt | head -n3
<< OUTPUT
PFAM : Multiple alignments and profile HMMs of protein domains
                        RELEASE 37.0
            --------------------------------------
OUTPUT
```

# MgtE

## Protein families

- COGs
  - [COG2239](https://www.ncbi.nlm.nih.gov/research/cog/cog/COG2239/):
    P; MgtE; Mg/Co/Ni transporter MgtE (contains CBS domain)
- Pfam
  - [PF01769](https://www.ebi.ac.uk/interpro/entry/pfam/PF01769/): MgtE;
    Divalent cation transporter
  - [PF03448](https://www.ebi.ac.uk/interpro/entry/pfam/PF03448/):
    MgtE_N; MgtE intracellular N domain
  - [PF00571](https://www.ebi.ac.uk/interpro/entry/pfam/PF00571/): CBS;
    CBS domain

## Collect proteins

``` {sh}
#| label: grep_COG2239; 49,929 hit queries
# resources_used.vmem=1340156kb;resources_used.walltime=00:10:46
qsub -q SMALL -m abe -M $Email -l select=1:ncpus=8:mem=48gb -l walltime=12:00:00 -e qsub_out/grep_COG2239_e -o qsub_out/grep_COG2239_o scripts/qsub/grep_COG2239.sh
{
source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR
bash scripts/bash/grep_cog.sh COG2239 data/tsv/proteins_raw_gtdbtax_cog.tsv.gz | \
  gzip > data/tsv/proteins_raw_gtdbtax_COG2239.tsv.gz
<< cmd
cat <(zcat $2 | head -n1) <(zcat $2 | grep -w $1)
cmd
}
# Compress PBS output
gzip qsub_out/grep_COG2239_{e,o}

# Number of hit queries
zcat data/tsv/proteins_raw_gtdbtax_COG2239.tsv.gz | cut -f 1 | tail -n +2 | sort | uniq | wc -l
<< OUTPUT
49929
OUTPUT
```

``` {sh}
#| label: hmmsearch_PF01769; 57,522 hit queries
cd $Home
mkdir -p data/tsv

# Check accession
zcat db/pfam/Pfam-A.hmm.gz | grep PF01769
<< OUTPUT
ACC   PF01769.21
OUTPUT

# resources_used.vmem=1902784kb;resources_used.walltime=00:02:26
qsub -q SMALL -m abe -M $Email -l select=1:ncpus=8:mem=48gb -l walltime=12:00:00 -e qsub_out/hmmsearch_PF01769_e -o qsub_out/hmmsearch_PF01769_o scripts/qsub/hmmsearch_PF01769.sh
{
source /etc/profile.d/modules.sh
module load hmmer/3.4
cd $PBS_O_WORKDIR
# Make temporal directory
mkdir -p data/tsv/hmmsearch_PF01769
# Extract profile
hmmfetch db/pfam/Pfam-A.hmm.gz PF01769.21 | gzip > db/pfam/PF01769.hmm.gz
# hmmsearch
seq -w 1 501 | xargs -I {} -P 8 -n 1 sh -c "hmmsearch --domtblout data/tsv/hmmsearch_PF01769/tmp_{}.txt db/pfam/PF01769.hmm.gz data/seqs_split/chunk_{}.faa.gz"
}
# Compress PBS output
gzip qsub_out/hmmsearch_PF01769_{e,o}

# Make tsv
seq -w 1 501 | \
  xargs -I {} sh -c "grep -v '^#' data/tsv/hmmsearch_PF01769/tmp_{}.txt" | \
  cat <(head -n3 data/tsv/hmmsearch_PF01769/tmp_001.txt) - | \
  bash scripts/bash/hmmsearch_tbl2tsv.sh | \
  gzip > data/tsv/hmmsearch_PF01769.tsv.gz

# Number of hit queries
zcat data/tsv/hmmsearch_PF01769.tsv.gz | cut -f 1 | tail -n +2 | sort | uniq | wc -l
<< OUTPUT
57522
OUTPUT
```

``` {sh}
#| label: hmmsearch_PF03448; 92,753 hit queries
cd $Home
mkdir -p data/tsv

# Check accession
zcat db/pfam/Pfam-A.hmm.gz | grep PF03448
<< OUTPUT
ACC   PF03448.22
OUTPUT

# resources_used.vmem=1888588kb;resources_used.walltime=00:02:27
qsub -q SMALL -m abe -M $Email -l select=1:ncpus=8:mem=48gb -l walltime=12:00:00 -e qsub_out/hmmsearch_PF03448_e -o qsub_out/hmmsearch_PF03448_o scripts/qsub/hmmsearch_PF03448.sh
{
source /etc/profile.d/modules.sh
module load hmmer/3.4
cd $PBS_O_WORKDIR
# Make temporal directory
mkdir -p data/tsv/hmmsearch_PF03448
# Extract profile
hmmfetch db/pfam/Pfam-A.hmm.gz PF03448.22 | gzip > db/pfam/PF03448.hmm.gz
# hmmsearch
seq -w 1 501 | xargs -I {} -P 8 -n 1 sh -c "hmmsearch --domtblout data/tsv/hmmsearch_PF03448/tmp_{}.txt db/pfam/PF03448.hmm.gz data/seqs_split/chunk_{}.faa.gz"
}
# Compress PBS output
gzip qsub_out/hmmsearch_PF03448_{e,o}

# Make tsv
seq -w 1 501 | \
  xargs -I {} sh -c "grep -v '^#' data/tsv/hmmsearch_PF03448/tmp_{}.txt" | \
  cat <(head -n3 data/tsv/hmmsearch_PF03448/tmp_001.txt) - | \
  bash scripts/bash/hmmsearch_tbl2tsv.sh | \
  gzip > data/tsv/hmmsearch_PF03448.tsv.gz

# Number of hit queries
zcat data/tsv/hmmsearch_PF03448.tsv.gz | cut -f 1 | tail -n +2 | sort | uniq | wc -l
<< OUTPUT
92753
OUTPUT
```

``` {sh}
#| label: hmmsearch_PF00571; 474,937 hit queries
cd $Home
mkdir -p data/tsv

# Check accession
zcat db/pfam/Pfam-A.hmm.gz | grep PF00571
<< OUTPUT
ACC   PF00571.33
OUTPUT

# resources_used.vmem=1884652kb;resources_used.walltime=00:02:32
qsub -q SMALL -m abe -M $Email -l select=1:ncpus=8:mem=48gb -l walltime=12:00:00 -e qsub_out/hmmsearch_PF00571_e -o qsub_out/hmmsearch_PF00571_o scripts/qsub/hmmsearch_PF00571.sh
{
source /etc/profile.d/modules.sh
module load hmmer/3.4
cd $PBS_O_WORKDIR
# Make temporal directory
mkdir -p data/tsv/hmmsearch_PF00571
# Extract profile
hmmfetch db/pfam/Pfam-A.hmm.gz PF00571.33 | gzip > db/pfam/PF00571.hmm.gz
# hmmsearch
seq -w 1 501 | xargs -I {} -P 8 -n 1 sh -c "hmmsearch --domtblout data/tsv/hmmsearch_PF00571/tmp_{}.txt db/pfam/PF00571.hmm.gz data/seqs_split/chunk_{}.faa.gz"
}
# Compress PBS output
gzip qsub_out/hmmsearch_PF00571_{e,o}

# Make tsv
seq -w 1 501 | \
  xargs -I {} sh -c "grep -v '^#' data/tsv/hmmsearch_PF00571/tmp_{}.txt" | \
  cat <(head -n3 data/tsv/hmmsearch_PF00571/tmp_001.txt) - | \
  bash scripts/bash/hmmsearch_tbl2tsv.sh | \
  gzip > data/tsv/hmmsearch_PF00571.tsv.gz

# Number of hit queries
zcat data/tsv/hmmsearch_PF00571.tsv.gz | cut -f 1 | tail -n +2 | sort | uniq | wc -l
<< OUTPUT
474937
OUTPUT
```

# Analysis for low-quality genomes

## MgtE: Protein families

- COGs
  - [COG2239](https://www.ncbi.nlm.nih.gov/research/cog/cog/COG2239/):
    P; MgtE; Mg/Co/Ni transporter MgtE (contains CBS domain)
- Pfam
  - [PF01769](https://www.ebi.ac.uk/interpro/entry/pfam/PF01769/): MgtE;
    Divalent cation transporter
  - [PF03448](https://www.ebi.ac.uk/interpro/entry/pfam/PF03448/):
    MgtE_N; MgtE intracellular N domain
  - [PF00571](https://www.ebi.ac.uk/interpro/entry/pfam/PF00571/): CBS;
    CBS domain

## MgtE: Collect proteins

``` {sh}
#| label: grep_COG2239_lq; 51,471 hit queries
# resources_used.vmem=1340088kb;resources_used.walltime=00:09:41
qsub -q SMALL -m abe -M $Email -l select=1:ncpus=8:mem=48gb -l walltime=12:00:00 -e qsub_out/grep_COG2239_lq_e -o qsub_out/grep_COG2239_lq_o scripts/qsub/grep_COG2239_lq.sh
{
source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR
bash scripts/bash/grep_cog.sh COG2239 data/tsv/lq_proteins_raw_gtdbtax_cog.tsv.gz | \
  gzip > data/tsv/lq_proteins_raw_gtdbtax_COG2239.tsv.gz
<< cmd
cat <(zcat $2 | head -n1) <(zcat $2 | grep -w $1)
cmd
}
# Compress PBS output
gzip qsub_out/grep_COG2239_lq_{e,o}

# Number of hit queries
zcat data/tsv/lq_proteins_raw_gtdbtax_COG2239.tsv.gz | cut -f 1 | tail -n +2 | sort | uniq | wc -l
<< OUTPUT
51471
OUTPUT
```

``` {sh}
#| label: hmmsearch_PF01769_lq; 61,894 hit queries
cd $Home
mkdir -p data/tsv

# Check accession
zcat db/pfam/Pfam-A.hmm.gz | grep PF01769
<< OUTPUT
ACC   PF01769.21
OUTPUT

# resources_used.vmem=1897556kb;resources_used.walltime=00:02:09
qsub -q SMALL -m abe -M $Email -l select=1:ncpus=8:mem=48gb -l walltime=12:00:00 -e qsub_out/hmmsearch_PF01769_lq_e -o qsub_out/hmmsearch_PF01769_lq_o scripts/qsub/hmmsearch_PF01769_lq.sh
{
source /etc/profile.d/modules.sh
module load hmmer/3.4
cd $PBS_O_WORKDIR
# Make temporal directory
mkdir -p data/tsv/hmmsearch_PF01769_lq
# Extract profile
hmmfetch db/pfam/Pfam-A.hmm.gz PF01769.21 | gzip > db/pfam/PF01769.hmm.gz
# hmmsearch
seq -w 1 501 | xargs -I {} -P 8 -n 1 sh -c "hmmsearch --domtblout data/tsv/hmmsearch_PF01769_lq/tmp_{}.txt db/pfam/PF01769.hmm.gz data/seqs_split_lq/chunk_{}.faa.gz"
}
# Compress PBS output
gzip qsub_out/hmmsearch_PF01769_lq_{e,o}

# Make tsv
seq -w 1 501 | \
  xargs -I {} sh -c "grep -v '^#' data/tsv/hmmsearch_PF01769_lq/tmp_{}.txt" | \
  cat <(head -n3 data/tsv/hmmsearch_PF01769_lq/tmp_001.txt) - | \
  bash scripts/bash/hmmsearch_tbl2tsv.sh | \
  gzip > data/tsv/hmmsearch_PF01769_lq.tsv.gz

# Number of hit queries
zcat data/tsv/hmmsearch_PF01769_lq.tsv.gz | cut -f 1 | tail -n +2 | sort | uniq | wc -l
<< OUTPUT
61894
OUTPUT
```

``` {sh}
#| label: hmmsearch_PF03448_lq; 95,461 hit queries
cd $Home
mkdir -p data/tsv

# Check accession
zcat db/pfam/Pfam-A.hmm.gz | grep PF03448
<< OUTPUT
ACC   PF03448.22
OUTPUT

# resources_used.vmem=1889700kb;resources_used.walltime=00:02:03
qsub -q SMALL -m abe -M $Email -l select=1:ncpus=8:mem=48gb -l walltime=12:00:00 -e qsub_out/hmmsearch_PF03448_lq_e -o qsub_out/hmmsearch_PF03448_lq_o scripts/qsub/hmmsearch_PF03448_lq.sh
{
source /etc/profile.d/modules.sh
module load hmmer/3.4
cd $PBS_O_WORKDIR
# Make temporal directory
mkdir -p data/tsv/hmmsearch_PF03448_lq
# Extract profile
hmmfetch db/pfam/Pfam-A.hmm.gz PF03448.22 | gzip > db/pfam/PF03448.hmm.gz
# hmmsearch
seq -w 1 501 | xargs -I {} -P 8 -n 1 sh -c "hmmsearch --domtblout data/tsv/hmmsearch_PF03448_lq/tmp_{}.txt db/pfam/PF03448.hmm.gz data/seqs_split_lq/chunk_{}.faa.gz"
}
# Compress PBS output
gzip qsub_out/hmmsearch_PF03448_lq_{e,o}

# Make tsv
seq -w 1 501 | \
  xargs -I {} sh -c "grep -v '^#' data/tsv/hmmsearch_PF03448_lq/tmp_{}.txt" | \
  cat <(head -n3 data/tsv/hmmsearch_PF03448_lq/tmp_001.txt) - | \
  bash scripts/bash/hmmsearch_tbl2tsv.sh | \
  gzip > data/tsv/hmmsearch_PF03448_lq.tsv.gz

# Number of hit queries
zcat data/tsv/hmmsearch_PF03448_lq.tsv.gz | cut -f 1 | tail -n +2 | sort | uniq | wc -l
<< OUTPUT
95461
OUTPUT
```

``` {sh}
#| label: hmmsearch_PF00571_lq; 443,719 hit queries
cd $Home
mkdir -p data/tsv

# Check accession
zcat db/pfam/Pfam-A.hmm.gz | grep PF00571
<< OUTPUT
ACC   PF00571.33
OUTPUT

# resources_used.vmem=1879092kb;resources_used.walltime=00:01:59
qsub -q SMALL -m abe -M $Email -l select=1:ncpus=8:mem=48gb -l walltime=12:00:00 -e qsub_out/hmmsearch_PF00571_lq_e -o qsub_out/hmmsearch_PF00571_lq_o scripts/qsub/hmmsearch_PF00571_lq.sh
{
source /etc/profile.d/modules.sh
module load hmmer/3.4
cd $PBS_O_WORKDIR
# Make temporal directory
mkdir -p data/tsv/hmmsearch_PF00571_lq
# Extract profile
hmmfetch db/pfam/Pfam-A.hmm.gz PF00571.33 | gzip > db/pfam/PF00571.hmm.gz
# hmmsearch
seq -w 1 501 | xargs -I {} -P 8 -n 1 sh -c "hmmsearch --domtblout data/tsv/hmmsearch_PF00571_lq/tmp_{}.txt db/pfam/PF00571.hmm.gz data/seqs_split_lq/chunk_{}.faa.gz"
}
# Compress PBS output
gzip qsub_out/hmmsearch_PF00571_lq_{e,o}

# Make tsv
seq -w 1 501 | \
  xargs -I {} sh -c "grep -v '^#' data/tsv/hmmsearch_PF00571_lq/tmp_{}.txt" | \
  cat <(head -n3 data/tsv/hmmsearch_PF00571_lq/tmp_001.txt) - | \
  bash scripts/bash/hmmsearch_tbl2tsv.sh | \
  gzip > data/tsv/hmmsearch_PF00571_lq.tsv.gz

# Number of hit queries
zcat data/tsv/hmmsearch_PF00571_lq.tsv.gz | cut -f 1 | tail -n +2 | sort | uniq | wc -l
<< OUTPUT
443719
OUTPUT
```
