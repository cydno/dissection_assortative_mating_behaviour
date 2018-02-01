# PERM DATA BINOMIAL GLM
# Permeates phenotype data over genotypes (and covariates) to provide null distribution of 
# lod scores ... i.e. correcting for multiple tests across genome.

# Load required packages

library(qtl)
library(lme4)
library(parallel)
library(plyr)

# Read in reduced data, use jittermap to account for individuals with same cM
args<-commandArgs(TRUE)
threads<-ifelse(is.na(args[1]),1,args[1])
message(paste("Running with ",threads,"threads"))
message("Loading data")
pasi_data<-read.cross(format = "csv", file = "../raw_data/data_for_Rqtl.csv", genotypes = c("AA","AB", "BB"), alleles = c("A","B"), estimate.map = FALSE, crosstype = "bc")
pasi_data <-jittermap(pasi_data)


# Remove duplicate markers (i.e. those adjacent markers with that are the same in all individuals)
# Iterations need to carried out over unique combinations of genotypes ... including more 
# simply increases computation time.

message("Removing duplicates")

dups<-findDupMarkers(pasi_data, exact.only = FALSE, adjacent.only = TRUE)
dups_to_drop<-unique(unlist(unique(dups)))

no_dup_data <- pasi_data

while (length(dups_to_drop) != 0) {
  dups<-findDupMarkers(no_dup_data, exact.only = FALSE, adjacent.only = TRUE)
  dups_to_drop<-unique(unlist(unique(dups)))
  print(length(dups_to_drop))
  no_dup_data <- drop.markers(no_dup_data, dups_to_drop)
}

reduced_data <- no_dup_data
reduced_data<- calc.genoprob(reduced_data)
pr <- pull.genoprob(reduced_data, omit.first.prob=TRUE)


# Pull out relevant phenotypes and phenotypes:
message("Extract phenotypes")

trials_cydno <- pull.pheno(reduced_data, pheno.col = "trials_court_cydno") 
trials_melp <- pull.pheno(reduced_data, pheno.col = "trials_court_melp")  
total_trials <- pull.pheno(reduced_data, pheno.col = "total_trials")
id <- pull.pheno(reduced_data, pheno.col = "id") 


# Pull out position of markers

positions <- (pull.genoprob(reduced_data, omit.first.prob=TRUE, include.pos.info=TRUE))$pos

# Set up data frame to store info from simulations i.e. max Lod scores

simulation_number <- NA
max_mp_LN <- NA
max_cp_LN <- NA
perm_data <- data.frame(simulation_number,
					    max_mp_LN,
					    max_cp_LN)
					  
perm_data  <- perm_data [0,]

# Start 1000 interations

get_lod_score<-function(locus_genotypes, all_trials, id) {
  model      <- glmer(all_trials ~ locus_genotypes + (1|id),  family = binomial )
  null_model <- glmer(all_trials ~ 1               + (1|id),  family = binomial )
  lod_score<-((logLik(model) - logLik(null_model))/log(10))[1]
  return(lod_score)
}

get_lod_scores<-function(locus, pr, all_trials_cyd_sample, all_trials_mel_sample, id_sample) {  
  locus_genotypes<-pr[,locus]
  return(c(locus,
           get_lod_score(locus_genotypes, all_trials_cyd_sample, id_sample),
           get_lod_score(locus_genotypes, all_trials_mel_sample, id_sample)
         )
  )
}

# just 20 interactions for now.
 
for (iter in 1:1000) {
  message(paste("Iteration",iter))

    start_time<-Sys.time()
	# Produce random order to be used for phenotype and 'covariates'
	randomise_pheno <- sample(c(1:length(trials_cydno)), size = length(trials_cydno))

	trials_cyd_sample <- trials_cydno[randomise_pheno]
	trials_mel_sample <- trials_melp[randomise_pheno]

	total_trials_sample <- total_trials[randomise_pheno]

	other_trials_cyd_sample <- total_trials_sample - trials_cyd_sample
	other_trials_mel_sample <- total_trials_sample - trials_mel_sample
	
	all_trials_cyd_sample <- cbind(trials_cyd_sample, other_trials_cyd_sample)
	all_trials_mel_sample <- cbind(trials_mel_sample, other_trials_mel_sample)
	
	id_sample<-id[randomise_pheno]
	
	# run glm across each position in the genome i.e. 1 ... to end of genome

	message("Processing positions")
	
	lod_scores<-mclapply(1:length(positions), get_lod_scores, pr, all_trials_cyd_sample, all_trials_mel_sample, id_sample, mc.cores=threads)
	
	run_data<-ldply(lod_scores)
	names(run_data)<-c("marker_number", "cp_LOD_score", "mp_LOD_score")

	message("Find max lod value")
	# find max lod value 
	
	max_mp_LN_value <- max(run_data$mp_LOD_score)
	max_cp_LN_value <- max(run_data$cp_LOD_score)
	
	#add run data tp data frame
	perm_data[iter,]$simulation_number <- iter
	perm_data[iter,]$max_mp_LN <- max_mp_LN_value
	perm_data[iter,]$max_cp_LN <- max_cp_LN_value
	
	# print some progress ...
	print(Sys.time()-start_time)
	print(iter)
	print(quantile(perm_data$max_mp_LN, c(0.90, 0.95)))
	print(quantile(perm_data$max_cp_LN, c(0.90, 0.95)))
}

write.csv(perm_data, "permutations_trial_prop.csv")
##########################################################################################