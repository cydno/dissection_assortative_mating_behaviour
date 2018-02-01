# R. M. Merrill Oct 2017

# Genome-wide QTL analysis - R scripts to run a binomial glmm at each genetic (cM) position across the genome (i.e. cydno x melpomene linkage map). In addition runs a non-parametric QTL analysis using R/qtl (i.e. scanone with model="np").


# GLMM includes id as an individual level random factor to account for overdispersion (see for example, Elston et al 2001 Parasitology).

# load required packages:
library(qtl)
library(lme4)
library(MASS)

# function for binonial glmm genome scan:

get_lod_score<-function(locus_genotypes, courtships, id, brood, GF) {
  model      <- glmer(courtships ~ locus_genotypes + (1|id),  family = binomial )
  null_model <- glmer(courtships ~ 1 + (1|id),  family = binomial )
  lod_score<-((logLik(model) - logLik(null_model))/log(10))[1]
  return(lod_score)
}

##########################################################################################
# LOAD GENOTYPE/PHENOTYPE DATA AND REMOVE REDUNDANT MARKERS 

# Load data:
masked_data<-read.cross(format = "csv", file = "../raw_data/data_for_Rqtl.csv", genotypes = c("AA","AB", "BB"), alleles = c("A","B"), estimate.map = FALSE, crosstype = "bc")

# Reduce the data to remove redundent ('duplicate') markers i.e. those at same genetic (cM) position. This speeds up analysis considerably, however it should be noted that the 'peak' lod score will fall on mulitiple markers at same genetic, but not physical, position 

dups<-findDupMarkers(masked_data, exact.only = FALSE, adjacent.only = TRUE)
dups_to_drop<-unique(unlist(unique(dups)))

no_dup_data <- masked_data

while (length(dups_to_drop) != 0) {
  dups<-findDupMarkers(no_dup_data, exact.only = FALSE, adjacent.only = TRUE)
  dups_to_drop<-unique(unlist(unique(dups)))
  print(length(dups_to_drop))
  no_dup_data <- drop.markers(no_dup_data, dups_to_drop)
}

reduced_data <- no_dup_data

# Use jittermap because some markers are at the same cM position ...
reduced_data<-jittermap(reduced_data)
# ... and calculate probability of genotypes at each locus.
reduced_data<- calc.genoprob(reduced_data, step = 1)

##########################################################################################
# Pull out genotype and phenotype information

# Calculate genotype probabilies across genome.
pr <- pull.genoprob(reduced_data, omit.first.prob=TRUE)

# Pull out chromsome and position of markers:
chr <- (pull.genoprob(reduced_data,  omit.first.prob=TRUE, include.pos.info=TRUE))$chr
# (Will give warning message)
pos <- (pull.genoprob(reduced_data, omit.first.prob=TRUE, include.pos.info=TRUE))$pos
# (Will give warning message)

# Pull out id, maternal and paternal brood data :
id <- pull.pheno(reduced_data, pheno.col = "id") 

# Pull out phenotype data for GLMM analysis:
court_cp <- pull.pheno(reduced_data, pheno.col = "court_cp")
court_mp <- pull.pheno(reduced_data, pheno.col = "court_mp")
courtships <- cbind(court_mp, court_cp)									

trials_melp <- pull.pheno(reduced_data, pheno.col = "trials_court_melp")
trials_cydno <- pull.pheno(reduced_data, pheno.col = "trials_court_cydno")
total_trials<-pull.pheno(reduced_data, pheno.col = "total_trials")
mp_trials <- cbind(trials_melp, (total_trials-trials_melp))				
cp_trials <-  cbind(trials_cydno, (total_trials-trials_cydno))			


##########################################################################################
##########################################################################################
##########################################################################################
# SET UP DATA-FRAME TO STORE RESULTS 
# Data-frame ("run_data") to store output (chr, pos, lod etc) and ensure lod vector is empty

run_data <- data.frame(chr,
					   pos)

# GLMM analysis results					   
run_data$locus_court<-NA 			# locus id for courtship glmm
run_data$court_lod<-NA  			# lod score for courtship glmm
run_data$mp_trial_lod<-NA  			# lod score for mp trials
run_data$mp_trial_lod<-NA  			# lod score for cp trials

# Non-parametric analysis results
run_data$court_lod_np<-NA  			# lod score non-parametric scan implemented with R/qtl
run_data$mp_trial_lod_np<-NA  		# lod score non-parametric scan implemented with R/qtl
run_data$cp_trial_lod_np<-NA  		# lod score non-parametric scan implemented with R/qtl

##########################################################################################

# RUN GENOME SCANS FOR TOTAL COURTSHIPS 

# GLMM

# Run genome scan and store lod scores: 

for (locus in 1:length(pos)) {
	run_data$locus_court[locus] <- locus
	locus_genotypes <- pr[,locus]
	run_data$court_lod[locus] <- get_lod_score(locus_genotypes, courtships, id )
}

# non-parametric methods with R/qtl

out.court.np <- scanone(reduced_data, pheno.col = "prop_court_mp", model="np")
run_data$court_lod_np <- out.court.np$lod

##########################################################################################
# RUN GENOME SCANS FOR COURTSHIP INITIATION


for (locus in 1:length(pos)) {
	run_data$locus_court[locus] <- locus
	locus_genotypes <- pr[,locus]
	run_data$mp_trial_lod[locus] <- get_lod_score(locus_genotypes, mp_trials, id )
	run_data$cp_trial_lod[locus] <- get_lod_score(locus_genotypes, cp_trials, id )
}
	
# "run_data" will contain lod scores for both courtship and trial genome scans at each cM position .

# set up phenotype (i.e. proportion of trials in which courtship towards i) melpomene or ii) cydno occured)
reduced_data$pheno <- cbind(reduced_data$pheno, prop_mp_trials = (trials_melp / total_trials))
reduced_data$pheno <- cbind(reduced_data$pheno, prop_cp_trials = (trials_cydno / total_trials))

# non-parametric methods with R/qtl
out.mp_trials.np <- scanone(reduced_data, pheno.col = "prop_mp_trials", model="np")
out.cp_trials.np <- scanone(reduced_data, pheno.col = "prop_cp_trials", model="np")

run_data$mp_trial_lod_np <- out.mp_trials.np$lod
run_data$cp_trial_lod_np <- out.cp_trials.np$lod

##########################################################################################

# Write this to a csv file.
write.csv(run_data,"../derived_data/genome_scan_lod_score.csv",  row.names = FALSE)

##########################################################################################
# Produce summary of QTL positions
##########################################################################################

# Gather data for chromosomes with QTL, i.e. 1, 17 and 18

chr1 <- run_data[run_data$chr == 1,]
chr17 <- run_data[run_data$chr == 17,]
chr18 <- run_data[run_data$chr == 18,]


# Get positions of maximum lod score (for GLMM analysis) 

# Build empty data.frame
pheno <- NA
chr <- NA
cM <- NA
lod <- NA
qtl_pos <- data.frame(pheno, chr, cM, lod)
qtl_pos <- qtl_pos[0,]


# Fill data.frame with data for position of max lod scores
for (i in c(1, 17, 18)) {
	
	# for GLMM analses
	max_lod <- max(run_data[run_data$chr == i,]$court_lod)
	pheno <- "courtship"
	chr <- i
	cM <- run_data[run_data$chr == i & run_data$court_lod == max_lod,]$pos
	lod <- run_data[run_data$chr == i & run_data$court_lod == max_lod,]$court_lod
	qtl_pos[nrow(qtl_pos) + 1,] <- c(pheno, chr, cM, lod) 

	max_lod <- max(run_data[run_data$chr == i,]$mp_trial_lod)
	pheno <- "mp_trial"
	chr <- i
	cM <- run_data[run_data$chr == i & run_data$mp_trial_lod == max_lod,]$pos
	lod <- run_data[run_data$chr == i & run_data$mp_trial_lod == max_lod,]$mp_trial_lod
	qtl_pos[nrow(qtl_pos) + 1,] <- c(pheno, chr, cM, lod) 
	
	max_lod <- max(run_data[run_data$chr == i,]$cp_trial_lod)
	pheno <- "cp_trial"
	chr <- i
	cM <- run_data[run_data$chr == i & run_data$cp_trial_lod == max_lod,]$pos
	lod <- run_data[run_data$chr == i & run_data$cp_trial_lod == max_lod,]$cp_trial_lod
	qtl_pos[nrow(qtl_pos) + 1,] <- c(pheno, chr, cM, lod) 

	
	# for non-parametric analyis
	max_lod <- max(run_data[run_data$chr == i,]$court_lod_np)
	pheno <- "courtship_np"
	chr <- i
	cM <- run_data[run_data$chr == i & run_data$court_lod_np == max_lod,]$pos
	lod <- run_data[run_data$chr == i & run_data$court_lod_np == max_lod,]$court_lod_np
	qtl_pos[nrow(qtl_pos) + 1,] <- c(pheno, chr, cM, lod) 


	max_lod <- max(run_data[run_data$chr == i,]$mp_trial_lod_np)
	pheno <- "mp_trial_np"
	chr <- i
	cM <- run_data[run_data$chr == i & run_data$mp_trial_lod_np == max_lod,]$pos
	lod <- run_data[run_data$chr == i & run_data$mp_trial_lod_np == max_lod,]$mp_trial_lod_np
	qtl_pos[nrow(qtl_pos) + 1,] <- c(pheno, chr, cM, lod) 
		
	max_lod <- max(run_data[run_data$chr == i,]$cp_trial_lod_np)
	pheno <- "cp_trial_np"
	chr <- i
	cM <- run_data[run_data$chr == i & run_data$cp_trial_lod_np == max_lod,]$pos
	lod <- run_data[run_data$chr == i & run_data$cp_trial_lod_np == max_lod,]$cp_trial_lod_np
	qtl_pos[nrow(qtl_pos) + 1,] <- c(pheno, chr, cM, lod) 
}

# Write summary as csv
write.csv(qtl_pos,"../derived_data/qtl_summary.csv", row.names = FALSE)
##########################################################################################


