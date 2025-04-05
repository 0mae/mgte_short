#!/usr/bin/bash

cat $1 | \
	sed '1d;3d' | \
	#sed '1s/^# //;1s/ name//g;1s/E-value/eval_full/;1s/E-value/eval_best/;1s/score/score_full/;1s/score /score_best /;1s/bias/bias_full/;1s/bias /bias_best /;1s/ of target//' | \
	sed '1s/^# //;1s/ name//g;1s/E-value/eval_full/;1s/c-Evalue/c_eval_dom/;1s/i-Evalue/i_eval_dom/;1s/score/score_full/;1s/score /score_dom /;1s/bias/bias_full/;1s/bias /bias_dom /;1s/ of target//;1s/from/from_hmm/;1s/to/to_hmm/;1s/from /from_ali /;1s/to /to_ali /;1s/from /from_env /;1s/to /to_env /;1s/accession/tacc/;1s/accession/qacc/;1s/#/num_dom/;1s/of/total_dom/'  | \
	sed '/^#/d' | \
	awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}'
