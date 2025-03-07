source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR
bash scripts/bash/grep_cog.sh COG2239 data/tsv/lq_proteins_raw_gtdbtax_cog.tsv.gz | \
  gzip > data/tsv/lq_proteins_raw_gtdbtax_COG2239.tsv.gz
<< cmd
cat <(zcat $2 | head -n1) <(zcat $2 | grep -w $1)
cmd
