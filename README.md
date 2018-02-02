## Analysis details for "Genetic dissection of assortative mating behaviour" ##

Richard M. Merrill, Pasi Rastas, Maria C. Melo, Sarah Barker, John Davey, W. Owen McMillan & Chris D. Jiggins

This repository documents analyses for the manuscript "Genetic dissection of assortative mating behaviour". These bespoke scripts are presented for transparency only and will require editing to be applied to other datasets. 

The repository includes : 

a) Scripts by Pasi Rastas for phasing data etc - contained within `scripts.zip`, detailed below (1)

b) Raw data for QTL analysis, including: 

`complete_preference_data.csv` = which includes courtship data for all 292 males included in the study. This is included for reference only and is not called by any script.

`data_for_Rqtl.csv` = phenotype and genotype data for males with RADseq genotypes in R/qtl format

`IDs_with_pheno_no_RAD.csv` =  phenotype data for individuals with no RADseq genotypes but known B/D locus genotype (determined from forewing colour)

`species_brood_data.csv` = phenotype data for parental species, F1 hybrids and backcross hybrids towards both species, but only includes individuals which performed courtships at least once.

c) Scripts for QTL and associated analysis by Richard Merrill, detailed below (2)

```
1_genome_scan.R
2_permutations_courtship_prop.R
2_permutations_trial_prop.R
3_qtl_significance.R
4_effect_size.R
5_QTL1_sims.R
5_QTL17_sims.R
5_QTL18_sims.R
```

d) Dervived data files for QTL and associated analysis by Richard Merrill, detailed below (2)

```
genome_scan_lod_score.csv
genome_wide_thresholds.csv
permutations_courtship_prop.csv
permutations_trial_prop.csv
QTL_1_sim_results.csv
QTL_17_sim_results.csv
QTL_18_sim_results.csv
qtl_data.csv
qtl_summary.csv
```

### 1) Scripts by Pasi Rastas to phase genotype data for each family and chromosome

scripts contained in `anal_code/phasing`

do phasing and masking of the data
*.filtered.gz is the output of Filtering2 module of LM3
map[1-21].txt ouput of Lep-MAP3 OrderMarker2 with marker number (column 1) mapped to "chr pos" (phased data is outputted by outputPhasedData=1)
snps_mapped.txt lists all snps in "chr pos" format

 create header

```
zcat all.filtered.gz |head -n 4|tail -n 3|./transpose_tab|awk '($3!=0)'|/.transpose_tab|head -n 2|awk -vOFS="\t" -vFS="\t" '{$2=$2 "\tcMPos"; print}' >header.txt
```

 add tabs between segregation patterns

```
for i in `seq 1 21`
do
paste <(cut -f 1-3 map$i.txt) <(cut -f 4- map$i.txt|sed -e 's/\t//g' -e 's/[01]/&\t/g') >map${i}_mapped.txt
done
```


create phased data

```
for i in `cat families.txt` 
do
zcat $i.filtered.gz|sed -e 's/*//g'|awk -f simpleConvert.awk| awk -f phase.awk >phased$i.txt
done

awk '{if (NR==1) print "CHR\tPOS\tcMPos"; else print $0 "\tNA"}' snps_mapped.txt >snps.txt
```

create phased_all.txt

```
awk 'BEGIN{s="paste snps.txt <(cut -f 3- phasedC20.txt)"}{s = s " <(cut -f 3- phased" $1 ".txt)"}END{print s " >phased_all.txt"}' families.txt|bash
```

calculate flips, correspondance of 

```
for i in `seq 1 21`
do
(head -n1 header.txt;(paste <(cut -f 1,2 map${i}_mapped.txt >tmp; awk '(NR==FNR){data[$1,$2]=$0}(NR!=FNR){print data[$1,$2]}' phased_all.txt tmp) map${i}_mapped.txt))|awk '(NR==1){for (i=4;i<=NF;++i) family[i]=$i}(NR>1){for (i=4;i<=NF/2;++i) if ($i==$(i+NF/2)) ++d[family[i], 1]; else if ($i!="-") ++d[family[i], -1]}END{for (i=4;i<=NF/2;++i) print family[i] "\t" d[family[i], 1]+0 "\t" d[family[i], -1]+0}'|awk '($2>$3) {print $1 "\t+"}($2<$3) {print $1 "\t-"}($2==$3) {print $1 "\t?"}' >flips$i.txt
done
for i in `seq 1 21`
do
awk -vFS="\t" -vOFS="\t" '(NR==FNR) {f[NR+3]=$2} (NR!=FNR){for (i=4;i<=NF;++i) if (f[i]=="-") $i=1-$i; print}' flips$i.txt map${i}_mapped.txt >map${i}_phased.txt
done
```

mask non-informative markers

```
zcat ../all.filtered.gz|awk -f inf.awk|awk 'BEGIN{print "#"} 1'|cut -f 3-|paste snps.txt - >inf_all.txt

for i in `seq 1 21`
do 
paste map${i}_phased.txt <(cut -f 1,2 map${i}_phased.txt >tmp; awk '(NR==FNR){data[$1,$2]=$0}(NR!=FNR){print data[$1,$2]}' inf_all.txt tmp)|awk 'BEGIN{size=334}{for (i=1;i<=NF;++i) data[NR,i]=$i}END{for (j=2;j<=NR;++j) for (i=4;i<=size;++i) if (data[j-1, i]!=data[j,i] && data[j-1, i]!="-" && data[j, i]!="-") {k=j-1;while (k >= 1 && data[k,i+size]!=1 && data[k,i+size]!=3) {data[k,i]="-";--k} k=j+1;while (k <= NR && data[k,i+size]!=1 && data[k,i+size]!=3) {data[k,i]="-"; ++k}; if (data[j,i+size]!=1 && data[j,i+size]!=3) data[j,i]="-" } for (j=1;j<=NR;++j) {s=data[j,1];for (i=2;i<=size;++i) s = s "\t" data[j,i]; print s}}' >map${i}_masked.txt; 
done
```

final output, map[1-21]_masked.txt

### 2) Scripts by Richard Merrill for QTL analyses:

### i) Genome-wide QTL scan

```
Rscript 1_genome_scan.R
```

R script to run a binomial GLMM at each genetic (cM) position across the genome (i.e. cydno x melpomene linkage map). GLMM includes id as an individual level random factor to account for overdispersion (see for example, Elston et al 2001, Parasitology). In addition, runs a non-parametric QTL analysis using R/qtl (using function scanone with model = "np"). Requires phenotype and genotype data for backcross hybrids in R/qtl format:  `data_for_Rqtl.csv`. Produces lod scores across the genome: `genome_scan_lod_score.csv`; and summary of qtl peaks:  `qtl_summary`


### ii) Permutations to determine genome-wide-significance 

```
Rscript 2_permutations_courtship_prop.R -20
Rscript 2_permutations_trial_prop.R -20
```

Scripts to run permutations used to determine genome-wide-significance: First (`2_permutations_courtship_prop.R`) for the proportion of courtships directed towards H. melpomene females; and second (`2_permutations_trial_prop.R`) for the proportion of trials in which courtship was initiated towards either a) H. melpomene, or b) H. cydno. Permeates phenotype data over genotypes (and covariates) to provide null distribution of lod scores, i.e. correcting for multiple tests across genome. This is best run over multiple threads (here = 20) as it will take a very long time. Requires phenotype and genotype data for backcross hybrids in R/qtl format:  `data_for_Rqtl.csv`. Produces distribution of lod scores: `permutations_courtship_prop.csv` and `permutations_trial_prop.csv`.




### iii) Test qtl significance

```
Rscript 3_qtl_significance.R 
```

Script to determine significance of individual qtl, using the null distribution of lod scores  generated above through permutation. Also builds a GLMM of the proportion of courtships towards H. melpomene with qtl on chr 1, 17 18 as explanatory factors. Tests whether to retain each qtl in a multiple-qtl model using a) likelihood ratio tests and b) a (very conservative) penalised lod score approach (see Browman and Sen 2009). Requires phenotype and genotype data for backcross hybrids in R/qtl format:  `raw_data/data_for_Rqtl.csv`. Produces table of qtl thresholds for GLMM analyses: `permutations_trial_prop.csv`; and determines P-values from permutations which are added to the qtl summary: `qtl_summary`. Produces `qtl_data.csv`, which is used below.


### iv) Determine effect sizes of qtl

```
Rscript 4_effect_size.R 
```

Script to determines effect sizes from GLMMs (measured as the proportion of the parental difference) and associated 95% confidence intervals (for figures 2 and 3b). Does this for qtl on chromosomes 1, 17 and 18, as well as for individuals without rad-seq data (but with known phenotype). Requires: `qtl_data.csv` generated in `3_qtl_significance.R`; data for hybrids without rad sequence data: `IDs_with_pheno_no_RAD.csv`; data for parental species (includes data for parentals, as well as F1 and backcross hybrids): `species_brood_data.csv`; and, the qtl summary: `qtl_summary.csv`.


### v) Estimate Beavis effect through simulation

```
Rscript 5_QTL1_sims.R
Rscript 5_QTL17_sims.R
Rscript 5_QTL18_sims.R
```

Script to runs simulation (for each qtl) to determine the extent to which the effects of our QTL may be over-estimated due to the Beavis effect. Takes a long time to run. Requires `qtl_data.csv` generated by `4_qtl_significance.R`. Produces: `QTL_1_sim_results.csv` etc. Note that the distribution of 'significant' simulations excludes runs in which the GLMM did not converge.
