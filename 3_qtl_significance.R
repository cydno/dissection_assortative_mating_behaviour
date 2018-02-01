# R. Merrill Oct 2017
# Analysis of proportions

# Load_data
run_data <- read.csv("../derived_data/genome_scan_lod_score.csv")
qtl_pos <- read.csv("../derived_data/qtl_summary.csv")
courtship_permutations <- read.csv("../derived_data/permutations_courtship_prop.csv")
trial_permutations <- read.csv("../derived_data/permutations_trial_prop.csv")

# load required packages:
library(qtl)
library(lme4)


##########################################################################################
# Load raw data:

masked_data<-read.cross(format = "csv", file = "../raw_data/data_for_Rqtl.csv", genotypes = c("AA","AB", "BB"), alleles = c("A","B"), estimate.map = FALSE, crosstype = "bc")

# reduce the data to remove redundent ('duplicate') markers (speeds up analysis considerably):
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

trials_melp <- pull.pheno(reduced_data, pheno.col = "trials_court_melp")
trials_cydno <- pull.pheno(reduced_data, pheno.col = "trials_court_cydno")
total_trials<-pull.pheno(reduced_data, pheno.col = "total_trials")
reduced_data$pheno <- cbind(reduced_data$pheno, prop_mp_trials = (trials_melp / total_trials))
reduced_data$pheno <- cbind(reduced_data$pheno, prop_cp_trials = (trials_cydno / total_trials))


##########################################################################################
# Significance of individual QTL
##########################################################################################

# 1) Determine significance thresholds (at alpha = 0.1 and 0.05):

# Make data frame for thresholds.

pheno <- c("courtship", "mp_trial", "cp_trial", "courtship_np", "mp_trial_np", "cp_trial_np") 
alpha_0.1 <- NA
alpha_0.05 <- NA

genome_wide_thresholds <- data.frame(pheno, alpha_0.1, alpha_0.05)

# GLMM thresholds:
genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship",]$alpha_0.05 <- quantile(courtship_permutations$max_LN, 0.95)
genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship",]$alpha_0.1 <- quantile(courtship_permutations$max_LN, 0.9)

genome_wide_thresholds[genome_wide_thresholds$pheno == "mp_trial",]$alpha_0.05 <- quantile(trial_permutations$max_mp_LN, 0.95)
genome_wide_thresholds[genome_wide_thresholds$pheno == "mp_trial",]$alpha_0.1 <- quantile(trial_permutations$max_mp_LN, 0.9)

genome_wide_thresholds[genome_wide_thresholds$pheno == "cp_trial",]$alpha_0.05 <- quantile(trial_permutations$max_cp_LN, 0.95)
genome_wide_thresholds[genome_wide_thresholds$pheno == "cp_trial",]$alpha_0.1 <- quantile(trial_permutations$max_cp_LN, 0.9)

# Non-parametric thresholds:
##########################################################################################
# run permutations
out.perm_courtships<- scanone(reduced_data, pheno.col = "prop_court_mp", n.perm = 10000, model="np")
out.perm_mp_trials.np <- scanone(reduced_data, pheno.col = "prop_mp_trials", n.perm = 10000, model="np")
out.perm_cp_trials.np <- scanone(reduced_data, pheno.col = "prop_cp_trials", n.perm = 10000, model="np")

# and aquire thresholds at alpha = 0.1 and 0.05
genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship_np",]$alpha_0.05 <- summary(out.perm_courtships)[1]
genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship_np",]$alpha_0.1 <- summary(out.perm_courtships) [2]

genome_wide_thresholds[genome_wide_thresholds$pheno == "mp_trial_np",]$alpha_0.05 <- summary(out.perm_mp_trials.np)[1]
genome_wide_thresholds[genome_wide_thresholds$pheno == "mp_trial_np",]$alpha_0.1 <- summary(out.perm_mp_trials.np)[2]

genome_wide_thresholds[genome_wide_thresholds$pheno == "cp_trial_np",]$alpha_0.05 <- summary(out.perm_cp_trials.np)[1]
genome_wide_thresholds[genome_wide_thresholds$pheno == "cp_trial_np",]$alpha_0.1 <- summary(out.perm_cp_trials.np)[2]

# Write to csv:

write.csv(genome_wide_thresholds, "../derived_data/genome_wide_thresholds.csv", row.names = FALSE)

# 2) Determine p-values for individual QTL:

# For GLMM methods:

qtl_pos$P_value <- NA
get_p.value <- function(x,perc) 1 - (ecdf(x)(perc))

qtl_pos[qtl_pos$pheno == "courtship",]$P_value <- get_p.value(courtship_permutations$max_LN, qtl_pos[qtl_pos$pheno == "courtship",]$lod)
qtl_pos[qtl_pos$pheno == "mp_trial",]$P_value <- get_p.value(trial_permutations$max_mp_LN, qtl_pos[qtl_pos$pheno == "mp_trial",]$lod)
qtl_pos[qtl_pos$pheno == "cp_trial",]$P_value <- get_p.value(trial_permutations$max_cp_LN, qtl_pos[qtl_pos$pheno == "cp_trial",]$lod)

# For non-parametric methods:

# prop. courtships:
out.court.np <- scanone(reduced_data, pheno.col = "prop_court_mp", model="np")
summary_court.np <- summary(out.court.np,  perms=out.perm_courtships, pvalues= TRUE)

qtl_pos[qtl_pos$pheno == "courtship_np" & qtl_pos$chr == "1" ,]$P_value <- summary_court.np[summary_court.np$chr == 1,]$pval
qtl_pos[qtl_pos$pheno == "courtship_np" & qtl_pos$chr == "17" ,]$P_value <- summary_court.np[summary_court.np$chr == 17,]$pval
qtl_pos[qtl_pos$pheno == "courtship_np" & qtl_pos$chr == "18" ,]$P_value <- summary_court.np[summary_court.np$chr == 18,]$pval

# prop trials mp
out.mp_trials.np <- scanone(reduced_data, pheno.col = "prop_mp_trials", model="np")
summary_mp_trial.np <-summary(out.mp_trials.np,  perms=out.perm_mp_trials.np, pvalues= TRUE)

qtl_pos[qtl_pos$pheno == "mp_trial_np" & qtl_pos$chr == "1" ,]$P_value <- summary_mp_trial.np[summary_mp_trial.np$chr == 1,]$pval
qtl_pos[qtl_pos$pheno == "mp_trial_np" & qtl_pos$chr == "17" ,]$P_value <- summary_mp_trial.np[summary_mp_trial.np$chr == 17,]$pval
qtl_pos[qtl_pos$pheno == "mp_trial_np" & qtl_pos$chr == "18" ,]$P_value <- summary_mp_trial.np[summary_mp_trial.np$chr == 18,]$pval

# prop trials cp
out.cp_trials.np <- scanone(reduced_data, pheno.col = "prop_cp_trials", model="np")
summary_cp_trial.np <-summary(out.cp_trials.np,  perms=out.perm_cp_trials.np, pvalues= TRUE)

qtl_pos[qtl_pos$pheno == "cp_trial_np" & qtl_pos$chr == "1" ,]$P_value <- summary_cp_trial.np[summary_cp_trial.np$chr == 1,]$pval
qtl_pos[qtl_pos$pheno == "cp_trial_np" & qtl_pos$chr == "17" ,]$P_value <- summary_cp_trial.np[summary_cp_trial.np$chr == 17,]$pval
qtl_pos[qtl_pos$pheno == "cp_trial_np" & qtl_pos$chr == "18" ,]$P_value <- summary_cp_trial.np[summary_cp_trial.np$chr == 18,]$pval

# Write over qtl_summary.csv to include p values.
write.csv(qtl_pos,"../derived_data/qtl_summary.csv", row.names = FALSE)

##########################################################################################
# MODEL COURTSHIPS
##########################################################################################
# GLMM of courtships including multiple QTL

# Make data.frame of i) courtships a phenotypes and ii) genotypes for each male at QTL positions

# Get phenotype data 

id <- pull.pheno(reduced_data, pheno.col = "id")
court_cp <- pull.pheno(reduced_data, pheno.col = "court_cp")
court_mp <- pull.pheno(reduced_data, pheno.col = "court_mp")
trials_court_cydno <- pull.pheno(reduced_data, pheno.col = "trials_court_cydno")
trials_court_melp <- pull.pheno(reduced_data, pheno.col = "trials_court_melp")
total_trials <- pull.pheno(reduced_data, pheno.col = "total_trials")

# Pull out genotypes at QTL
mar1 <- find.marker(reduced_data, 1, as.numeric(qtl_pos[qtl_pos$chr == 1 & qtl_pos$pheno == "courtship",]$cM))
QTL1 <- pull.geno(reduced_data)[,mar1]
if (sum(is.na(QTL1)) != 0) {
	message("Warning! missing chr 1 genotypes ... ")
}
mar17 <- find.marker(reduced_data, 17, as.numeric(qtl_pos[qtl_pos$chr == 17 & qtl_pos$pheno == "courtship",]$cM))
QTL17 <- pull.geno(reduced_data)[,mar17]
if (sum(is.na(QTL17)) != 0) {
	message("Warning! missing chr 1 genotypes ... ")
	QTL17 <- pull.geno(fill.geno(reduced_data))[,"Hmel217014_351457"]
	QTL17
}
mar18 <- find.marker(reduced_data, 18, as.numeric(qtl_pos[qtl_pos$chr == 18 & qtl_pos$pheno == "courtship",]$cM))
QTL18 <- pull.geno(reduced_data)[,mar18]
if (sum(is.na(QTL18)) != 0) {
	message("Warning! missing chr 1 genotypes ... ")
}

mar18b <- find.marker(reduced_data, 18, as.numeric(qtl_pos[qtl_pos$chr == 18 & qtl_pos$pheno == "mp_trial",]$cM))
QTL18b <- pull.geno(reduced_data)[,mar18]
if (sum(is.na(QTL18b)) != 0) {
	message("Warning! missing chr 1 genotypes ... ")
}

# Fill data frame with phentype and genptype data
QTL_data_all <- data.frame(id,
							court_cp,
							court_mp,
							trials_court_cydno,
							trials_court_melp, 
							total_trials,
							QTL1,
							QTL17,
							QTL18,
							QTL18b
							)
														
QTL_data<-QTL_data_all[complete.cases(QTL_data_all),]	
QTL_data$total_court <- QTL_data$court_mp  + QTL_data$court_cp 

QTL_data$geno1 <- NA
QTL_data$geno17 <- NA
QTL_data$geno18 <- NA
QTL_data$geno18b <- NA

QTL_data[QTL_data$QTL1 == 1, ]$geno1 <- "A"
QTL_data[QTL_data$QTL1 == 2, ]$geno1 <- "B"

QTL_data[QTL_data$QTL17 == 1, ]$geno17 <- "A"
QTL_data[QTL_data$QTL17 == 2, ]$geno17 <- "B"

QTL_data[QTL_data$QTL18 == 1, ]$geno18 <- "A"
QTL_data[QTL_data$QTL18 == 2, ]$geno18 <- "B"

QTL_data[QTL_data$QTL18b == 1, ]$geno18b <- "A"
QTL_data[QTL_data$QTL18b == 2, ]$geno18b <- "B"

##########################################################################################
# Build GLMMs to determine significance of removing individual QTL from full model, for
# courtship data

# Response variable
courtships <- cbind(QTL_data$court_mp, QTL_data$court_cp) 

# Specify GLMMs: binomial error structure, proportion of courtships towards melpomene
# as the response ("courtships"), and genotype at loci on chr 1, 17 and/or 18 as fixed 
# effects, and male id as a random effect to account for overdispersion.  

full_model <- glmer(courtships ~ geno1 + geno17 + geno18  +  (1|id),  family = binomial,  data = QTL_data )
no_QTL1 <- glmer(courtships ~  geno17 + geno18  + (1|id),  family = binomial,  data = QTL_data )
no_QTL17 <- glmer(courtships ~ geno1 +  geno18  + (1|id),  family = binomial,  data = QTL_data )
no_QTL18 <- glmer(courtships ~ geno1 +  geno17  + (1|id),  family = binomial,  data = QTL_data )
null_model <- glmer(courtships ~ 1 + (1|id),  family = binomial,  data = QTL_data )

# Likelihood ratio tests. Note that this does not correct for 'multiple testing' across 
# the genome

LRT1 <- anova(full_model, no_QTL1)
LRT17 <- anova(full_model, no_QTL17)
LRT18 <- anova(full_model, no_QTL18)

# Report Chisq values for LRTs

qtl_pos$ChiSq <- NA
qtl_pos[qtl_pos$pheno == "courtship" & qtl_pos$chr == "1" ,] $ChiSq <- LRT1$Chisq[2]
qtl_pos[qtl_pos$pheno == "courtship" & qtl_pos$chr == "17" ,]$ChiSq <- LRT17$Chisq[2]
qtl_pos[qtl_pos$pheno == "courtship" & qtl_pos$chr == "18" ,]$ChiSq <- LRT18$Chisq[2]

# Determine penalised lod score (see Sen and Browman R/qtl book). This utilises the genome- wide significance threshold (determined through permutation) to determine significance, (taking into account 'multiple-testing' across the genome).

# Determine lod score for full GLMM (all three QTL) and models with single QTL removed.

lod_score_full<-((logLik(full_model) - logLik(null_model))/log(10))[1]
lod_score_no1<-((logLik(no_QTL1 ) - logLik(null_model))/log(10))[1]
lod_score_no17<-((logLik(no_QTL17) - logLik(null_model))/log(10))[1]
lod_score_no18<-((logLik(no_QTL18) - logLik(null_model))/log(10))[1]

# Report lod score for full model (i.e. with QTL1, QTL17 and QTL18)
qtl_pos[(nrow(qtl_pos)+1),]$pheno <- "courtship"
qtl_pos[(nrow(qtl_pos)),]$chr <- "all"
qtl_pos[(nrow(qtl_pos)),]$lod <- lod_score_full

# Calculate and report penalised lod after removal of QTL, i.e. reported removal_pLOD
# for chr 1 is the pLOD for a model including QTL17 and QTL18 but NOT QTL1 (but reported 
# removal_pLOD for full model is the pLOD for full model (i.e. QTL1, QTL17 and QTL18). Note
# that the threshold  was determined previously through permutation. 

qtl_pos$removal_pLOD <- NA
qtl_pos[qtl_pos$pheno == "courtship" & qtl_pos$chr == "1" ,]$removal_pLOD <- lod_score_no1 - (genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship",]$alpha_0.05*2)
qtl_pos[qtl_pos$pheno == "courtship" & qtl_pos$chr == "17" ,]$removal_pLOD <- lod_score_no17 - (genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship",]$alpha_0.05*2)
qtl_pos[qtl_pos$pheno == "courtship" & qtl_pos$chr == "18" ,]$removal_pLOD <- lod_score_no18 - (genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship",]$alpha_0.05*2)
qtl_pos[qtl_pos$pheno == "courtship" & qtl_pos$chr == "all" ,]$removal_pLOD <- lod_score_full - (genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship",]$alpha_0.05*3)

# Determine change in penalised lod to for model comaparason.
pLOD_full <- lod_score_full - (genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship",]$alpha_0.05*3)
pLOD_no1  <- lod_score_no1 - (genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship",]$alpha_0.05*2)
pLOD_no17 <- lod_score_no17 - (genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship",]$alpha_0.05*2)
pLOD_no18 <- lod_score_no18 - (genome_wide_thresholds[genome_wide_thresholds$pheno == "courtship",]$alpha_0.05*2)

qtl_pos$change_pLOD <- NA
qtl_pos[qtl_pos$pheno == "courtship" & qtl_pos$chr == "1" ,]$change_pLOD <- pLOD_no1 - pLOD_full
qtl_pos[qtl_pos$pheno == "courtship" & qtl_pos$chr == "17" ,]$change_pLOD <- pLOD_no17 - pLOD_full
qtl_pos[qtl_pos$pheno == "courtship" & qtl_pos$chr == "18" ,]$change_pLOD <- pLOD_no18 - pLOD_full

# Write over qtl_summary.csv to include p values.
write.csv(qtl_pos,"../derived_data/qtl_summary.csv", row.names = FALSE)

# Write QTL_data to .csv (rquired by effect size)
write.csv(QTL_data,"../derived_data/qtl_data.csv", row.names = FALSE)
##########################################################################################
