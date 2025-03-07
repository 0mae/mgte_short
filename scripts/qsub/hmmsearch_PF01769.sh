source /etc/profile.d/modules.sh
module load hmmer/3.4
cd $PBS_O_WORKDIR
# Make temporal directory
mkdir -p data/tsv/hmmsearch_PF01769
# Extract profile
hmmfetch db/pfam/Pfam-A.hmm.gz PF01769.21 | gzip > db/pfam/PF01769.hmm.gz
# hmmsearch
seq -w 1 501 | xargs -I {} -P 8 -n 1 sh -c "hmmsearch --domtblout data/tsv/hmmsearch_PF01769/tmp_{}.txt db/pfam/PF01769.hmm.gz data/seqs_split/chunk_{}.faa.gz"
