#! /bin/bash

#for sp in Xl Xt Dr Bf Ci Cg
for sp in Gg Ps
do
	echo $sp
	perl get_dN_dS_ratio.pl $sp/yn $sp/$sp.gene $sp/Mm.gene $sp/${sp}_dNdS_res.txt
done	
