# R. M. Merrill Oct 2017

# Takes positions from QTL genome scan and determines effect sizes calculated as the proportion
# of difference between parental species (i.e. pure melopmene and cydno males).

# Second part produces effect size plots (proportion of courtships for alternative genotypes.)

###########################################################################################
###########################################################################################

# load required packages:
library(lme4)
library(ggplot2)

###########################################################################################
###########################################################################################

# Data required:

# from qtl_significance.R 
QTL_data <- read.csv("../derived_data/qtl_data.csv")

# Data without RAD genotypes 
without_RAD <- read.csv("../raw_data/IDs_with_pheno_no_RAD.csv")

# & data for parentals (and other individuals)
court_data <- read.csv("../raw_data/species_brood_data.csv")

# & qtl summary
qtl_summary <- read.csv("../derived_data/qtl_summary.csv")
###########################################################################################
# Build GLMM for each qtl indivually and for all 3 QTL
###########################################################################################

# Exactly as for the 'multiple QTL' (penalised LOD score) analysis i.e. GLMM with binomial
# error structure, proportion of courtships towards melpomene as the response, genotype
# at loci on chr 1, 17 and/or 18 as fixed effects, and male id as a random effect to account
# for overdispersion.  

# tidy up QTL_data frame for GLMMs ...

QTL_data$total_court <- QTL_data$court_mp  + QTL_data$court_cp 

QTL_data$geno1 <- NA
QTL_data$geno17 <- NA
QTL_data$geno18 <- NA

QTL_data[QTL_data$QTL1 == 1, ]$geno1 <- "A"
QTL_data[QTL_data$QTL1 == 2, ]$geno1 <- "B"

QTL_data[QTL_data$QTL17 == 1, ]$geno17 <- "A"
QTL_data[QTL_data$QTL17 == 2, ]$geno17 <- "B"

QTL_data[QTL_data$QTL18 == 1, ]$geno18 <- "A"
QTL_data[QTL_data$QTL18 == 2, ]$geno18 <- "B"


# convery QTL loci ('geno1' etc) to factors:

QTL_data$geno1 <- as.factor(QTL_data$geno1)
QTL_data$geno17 <- as.factor(QTL_data$geno17)
QTL_data$geno18 <- as.factor(QTL_data$geno18)

# combine courtships towards melpomene and females to produce response variable
courtships <- cbind(QTL_data$court_mp, QTL_data$court_cp) 

# call GLMMs (with lme4 package)
only1 <- glmer(courtships~ geno1  + (1|id),  family = binomial,  data = QTL_data )
only17 <- glmer(courtships ~ geno17  + (1|id),  family = binomial,  data = QTL_data )
only18 <- glmer(courtships ~ geno18  + (1|id),  family = binomial,  data = QTL_data )
all <- glmer(courtships ~ geno1 + geno17 + geno18  + (1|id),  family = binomial,  data = QTL_data )

# Write over qtl_data.csv to include effects.
write.csv(QTL_data,"../derived_data/qtl_data.csv", row.names = FALSE)

###########################################################################################
# Determine means and 95% confidence intervals for loci individually, and then all three
###########################################################################################

# With thanks to Justin Touchon for the orginal code.
### ... check with Justin that these are 95% c.i.s ... and describe ...

################################### chr1 effects #########################################

newdat<-expand.grid(geno1=levels(QTL_data$geno1), courtships=0)
mm <- model.matrix(terms(only1),newdat)
newdat$courtships <- mm %*% fixef(only1)
pvar1 <- diag(mm %*% tcrossprod(vcov(only1),mm))
newdat <- data.frame(
    newdat
    , plo = newdat$courtships-2*sqrt(pvar1)
    , phi = newdat$courtships+2*sqrt(pvar1)
)
chr1_effects<-newdat
chr1_effects[,2:4] <- mapply(plogis, chr1_effects[,2:4])

################################### chr17 effects ########################################

newdat<-expand.grid(geno17=levels(QTL_data$geno17), courtships=0)
mm <- model.matrix(terms(only17),newdat)
newdat$courtships <- mm %*% fixef(only17)
pvar1 <- diag(mm %*% tcrossprod(vcov(only17),mm))
newdat <- data.frame(
    newdat
    , plo = newdat$courtships-2*sqrt(pvar1)
    , phi = newdat$courtships+2*sqrt(pvar1)
)
chr17_effects<-newdat
chr17_effects[,2:4] <- mapply(plogis, chr17_effects[,2:4])

################################### chr18 effects ########################################

newdat<-expand.grid(geno18=levels(QTL_data$geno18), courtships=0)
mm <- model.matrix(terms(only18),newdat)
newdat$courtships <- mm %*% fixef(only18)
pvar1 <- diag(mm %*% tcrossprod(vcov(only18),mm))
newdat <- data.frame(
    newdat
    , plo = newdat$courtships-2*sqrt(pvar1)
    , phi = newdat$courtships+2*sqrt(pvar1)
)
chr18_effects<-newdat
chr18_effects[,2:4] <- mapply(plogis, chr18_effects[,2:4])

################################### 'all' effects ########################################

newdat<-expand.grid(geno1=levels(QTL_data$geno1), geno17=levels(QTL_data$geno17), geno18=levels(QTL_data$geno18), courtships=0)
mm <- model.matrix(terms(all),newdat)
newdat$courtships <- mm %*% fixef(all)
pvar1 <- diag(mm %*% tcrossprod(vcov(all),mm))
newdat <- data.frame(
    newdat
    , plo = newdat$courtships-2*sqrt(pvar1)
    , phi = newdat$courtships+2*sqrt(pvar1)
)
all_effects<-newdat
all_effects[,4:6] <- mapply(plogis, all_effects[,4:6])


###########################################################################################
# Determine means and 95% confidence intervals for B locus (v. tight linkage to QTL on chr 18)
# for individuals without RADseq (genetic) data 
###########################################################################################

# Using the same methods as above, determine means and confidence intervals for males 
# with red or white forwings (i.e. Bb and bb at the B colour pattern locus, which segragates 
# almost perfectly with the QTL on chr 18). This analysis incoroporates ??? male hybrids
# for which no RADseq (genetic) data is available. As a result these individuals are 
# not included in the whole genome QTL analysis allowing an estimation of effect size 
# independent of determining QTL significance. 

# 'without_RAD' = data.frame with hybrid individuals with no RADseq data

# combine courtships towards melpomene and females to produce response variable
courtships <- cbind(without_RAD$court_mp, without_RAD$court_cp)

# call GLMM (with lme4 package), same as before but with forewing colour as the fixed effect
no_rad <- glmer(courtships ~ forewing  + (1|id),  family = binomial,  data = without_RAD )

# determine means and 95% confidence intervals
newdat<-expand.grid(forewing=levels(without_RAD$forewing), courtships=0)
mm <- model.matrix(terms(no_rad),newdat)
newdat$courtships <- mm %*% fixef(no_rad)
pvar1 <- diag(mm %*% tcrossprod(vcov(no_rad),mm))
newdat <- data.frame(
    newdat
    , plo = newdat$courtships-2*sqrt(pvar1)
    , phi = newdat$courtships+2*sqrt(pvar1)
)
no_rad_effects<-newdat
no_rad_effects[,2:4] <- mapply(plogis, no_rad_effects[,2:4])



###########################################################################################
# Determine means for parental species, i.e. melpomene and cydno males assayed in the 
# same way as the hybrids
###########################################################################################

mp_mean <- plogis(coef(summary(glmer(cbind(court_mp, court_cp) ~  1 + (1|insectary_id),  family = binomial,  data = court_data[court_data$Type == "MP",] )))[ , "Estimate"])
cp_mean <- plogis(coef(summary(glmer(cbind(court_mp, court_cp) ~  1 + (1|insectary_id),  family = binomial,  data = court_data[court_data$Type == "CP",] )))[ , "Estimate"])

###########################################################################################
# Get effect sizes measured as the proportion of the difference between parental species
# explained
###########################################################################################

parental_diff <- mp_mean - cp_mean

QTL1_effect 	<- (chr1_effects$courtships[2] - chr1_effects$courtships[1]) / parental_diff
QTL17_effect 	<- (chr17_effects$courtships[2] - chr17_effects$courtships[1]) / parental_diff
QTL18_effect 	<- (chr18_effects$courtships[2] - chr18_effects$courtships[1]) / parental_diff
all_QTL_effect 	<- (all_effects$courtships[8] - all_effects$courtships[1]) / parental_diff

# 'correct' for fact that measured pref for homozygotes at all three QTL is lower than
# that measured for pure cydno
all_QTL_effect_corrected 	<-(all_effects$courtships[8] - cp_mean) / parental_diff

B_effect 		<- (no_rad_effects$courtships[1] - no_rad_effects$courtships[2]) / parental_diff 	# note that this is other way around because no_rad_effects$courtships[1] = red no_rad_effects$courtships[2] = white

# print 'effect sizes'
QTL1_effect 	
QTL17_effect 	
QTL18_effect 	
all_QTL_effect 
all_QTL_effect_corrected
B_effect 

# Report effect sizes for three QTL and all together.
qtl_summary$prop_parent_difference 	<- NA
qtl_summary[qtl_summary$pheno == "courtship" & qtl_summary$chr == 1,]$prop_parent_difference <- QTL1_effect 
qtl_summary[qtl_summary$pheno == "courtship" & qtl_summary$chr == 17,]$prop_parent_difference <- QTL17_effect 
qtl_summary[qtl_summary$pheno == "courtship" & qtl_summary$chr == 18,]$prop_parent_difference <- QTL18_effect 
qtl_summary[qtl_summary$pheno == "courtship" & qtl_summary$chr == "all",]$prop_parent_difference <- all_QTL_effect 

# Report corrected effect size, for 'all' QTL
qtl_summary$prop_parent_difference_corrected <- NA
qtl_summary[qtl_summary$pheno == "courtship" & qtl_summary$chr == "all",]$prop_parent_difference_corrected <- all_QTL_effect_corrected

# Report effect of preference for colour pattern for non-RAD individuals
qtl_summary[(nrow(qtl_summary) + 1),]$pheno <- "courtship"
levels(qtl_summary$chr) <- c(levels(qtl_summary$chr), "b_locus")
qtl_summary[(nrow(qtl_summary)),]$chr <- "b_locus"
qtl_summary[(nrow(qtl_summary)),]$prop_parent_difference <- B_effect

# Write over qtl_summary.csv to include effects.
write.csv(qtl_summary,"../derived_data/qtl_summary.csv", row.names = FALSE)

###########################################################################################
# Make effect plots
###########################################################################################
#
# Produces plots of proportion of courtships directed towards melpomene females.
# Plot includes individual data points for each male the size of which represents the
# total number of courtships for that individual. This is important because the binomial
# GLMMs are wighted by number of courtship events. Means, and associated confidence intevals
# for each genotype class are generated. Orange and blue dashed lines represent mean values
# for melpomene and cydno males measured in the same way. Black hori
# 


# calculate proportion courtships towards melpomen for each male 
QTL_data$prop_melp <- QTL_data$court_mp / QTL_data$total_court

# make 'position marker' for genotype at each QTL
QTL_data$pos_1 <- 1
QTL_data[QTL_data$geno1 == "B",]$pos_1<- 2

QTL_data$pos_17 <- 1
QTL_data[QTL_data$geno17 == "B",]$pos_17<- 2

QTL_data$pos_18 <- 1
QTL_data[QTL_data$geno18 == "B",]$pos_18<- 2

QTL_data$pos_all <- NA
QTL_data[QTL_data$geno1 == "A" & QTL_data$geno17 == "A" & QTL_data$geno18 == "A",]$pos_all<- 1
QTL_data[QTL_data$geno1 == "B" & QTL_data$geno17 == "B" & QTL_data$geno18 == "B",]$pos_all<- 2

#################################### chr1 plot ###########################################

effect_chr1 <- ggplot(QTL_data, aes(x = pos_1, y = prop_melp)) +
	
	geom_segment(x = -Inf, y = mp_mean, xend = Inf, yend = mp_mean, colour = "orange", size = 1, linetype=2) +
	geom_segment(x = -Inf, y = cp_mean, xend = Inf, yend = cp_mean, colour = "light blue", size = 1, linetype=2) +
	
	geom_jitter(alpha = 0.4, position=position_jitterdodge(jitter.width=0.4, dodge.width=0.85), colour = "dark gray", aes(size=total_court), show.legend = F) +
  	theme_bw() +
  	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  	theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
  	labs(x = "", y = "") +
  	scale_size(name = "") +
  	
  	# add means:
  	annotate("point",size = 3, shape = 15, x = 1, y = chr1_effects$courtships[1]) +
  	annotate("point",size = 3, shape = 15, x = 2, y = chr1_effects$courtships[2]) +
  	
  	geom_segment(x = 1, y = chr1_effects$plo[1], xend = 1, yend = chr1_effects$phi[1], colour = "black", size = 0.25) +
  	geom_segment(x = 1 - 0.05, y = chr1_effects$plo[1], xend = 1 +0.05, yend = chr1_effects$plo[1], colour = "black", size = 0.25) +
  	geom_segment(x = 1 - 0.05, y = chr1_effects$phi[1], xend = 1+0.05, yend = chr1_effects$phi[1], colour = "black", size = 0.25) +
  	
  	geom_segment(x = 2, y = chr1_effects$plo[2], xend = 2, yend = chr1_effects$phi[2], colour = "black", size = 0.25) +
  	geom_segment(x = 2 - 0.05, y = chr1_effects$plo[2], xend = 2 +0.05, yend = chr1_effects$plo[2], colour = "black", size = 0.25) +
  	geom_segment(x = 2 - 0.05, y = chr1_effects$phi[2], xend = 2+0.05, yend = chr1_effects$phi[2], colour = "black", size = 0.25) +
  	
  	geom_segment(x = 1, y = chr1_effects$courtships[1], xend = 2, yend = chr1_effects$courtships[2], colour = "black", linetype=2, size = 0.01) +
  	geom_segment(x = 2.5, y = chr1_effects$courtships[1], xend = 2.5, yend = chr1_effects$courtships[2], colour = "black", size = 1)
  	
  	
#################################### chr17 plot ##########################################
 	  	
effect_chr17 <- ggplot(QTL_data, aes(x = pos_17, y = prop_melp)) +
	
	geom_segment(x = -Inf, y = mp_mean, xend = Inf, yend = mp_mean, colour = "orange", size = 1, linetype=2) +
	geom_segment(x = -Inf, y = cp_mean, xend = Inf, yend = cp_mean, colour = "light blue", size = 1, linetype=2) +
	
	geom_jitter(alpha = 0.4, position=position_jitterdodge(jitter.width=0.4, dodge.width=0.85), colour = "dark gray", aes(size=total_court),show.legend = F) +
  	theme_bw() +
  	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  	theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
  	labs(x = "", y = "") +
  	scale_size(name = "")  +
  	
  	# add means:
  	annotate("point",size = 3, shape = 15, x = 1, y = chr17_effects$courtships[1]) +
  	annotate("point",size = 3, shape = 15, x = 2, y = chr17_effects$courtships[2]) +
  	
  	geom_segment(x = 1, y = chr17_effects$plo[1], xend = 1, yend = chr17_effects$phi[1], colour = "black", size = 0.25) +
  	geom_segment(x = 1 - 0.05, y = chr17_effects$plo[1], xend = 1 +0.05, yend = chr17_effects$plo[1], colour = "black", size = 0.25) +
  	geom_segment(x = 1 - 0.05, y = chr17_effects$phi[1], xend = 1+0.05, yend = chr17_effects$phi[1], colour = "black", size = 0.25) +
  	
  	geom_segment(x = 2, y = chr17_effects$plo[2], xend = 2, yend = chr17_effects$phi[2], colour = "black", size = 0.25) +
  	geom_segment(x = 2 - 0.05, y = chr17_effects$plo[2], xend = 2 +0.05, yend = chr17_effects$plo[2], colour = "black", size = 0.25) +
  	geom_segment(x = 2 - 0.05, y = chr17_effects$phi[2], xend = 2+0.05, yend = chr17_effects$phi[2], colour = "black", size = 0.25)	+
  	
  	geom_segment(x = 1, y = chr17_effects$courtships[1], xend = 2, yend = chr17_effects$courtships[2], colour = "black", linetype=2, size = 0.01) +
  	geom_segment(x = 2.5, y = chr17_effects$courtships[1], xend = 2.5, yend = chr17_effects$courtships[2], colour = "black", size = 1)
  	
#################################### chr18 plot ##########################################

 effect_chr18 <- ggplot(QTL_data, aes(x = pos_18, y = prop_melp)) +
	
	geom_segment(x = -Inf, y = mp_mean, xend = Inf, yend = mp_mean, colour = "orange", size = 1, linetype=2) +
	geom_segment(x = -Inf, y = cp_mean, xend = Inf, yend = cp_mean, colour = "light blue", size = 1, linetype=2) +
	
	geom_jitter(alpha = 0.4, position=position_jitterdodge(jitter.width=0.4, dodge.width=0.85), colour = "dark gray", aes(size=total_court),show.legend = F) +
  	theme_bw() +
  	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  	theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
  	labs(x = "", y = "") +
  	scale_size(name = "") +
  	
  	# add means:
  	annotate("point",size = 3, shape = 15, x = 1, y = chr18_effects$courtships[1]) +
  	annotate("point",size = 3, shape = 15, x = 2, y = chr18_effects$courtships[2]) +
  	
  	geom_segment(x = 1, y = chr18_effects$plo[1], xend = 1, yend = chr18_effects$phi[1], colour = "black", size = 0.25) +
  	geom_segment(x = 1 - 0.05, y = chr18_effects$plo[1], xend = 1 +0.05, yend = chr18_effects$plo[1], colour = "black", size = 0.25) +
  	geom_segment(x = 1 - 0.05, y = chr18_effects$phi[1], xend = 1+0.05, yend = chr18_effects$phi[1], colour = "black", size = 0.25) +
  	
  	geom_segment(x = 2, y = chr18_effects$plo[2], xend = 2, yend = chr18_effects$phi[2], colour = "black", size = 0.25) +
  	geom_segment(x = 2 - 0.05, y = chr18_effects$plo[2], xend = 2 +0.05, yend = chr18_effects$plo[2], colour = "black", size = 0.25) +
  	geom_segment(x = 2 - 0.05, y = chr18_effects$phi[2], xend = 2+0.05, yend = chr18_effects$phi[2], colour = "black", size = 0.25)	+
  	
  	geom_segment(x = 1, y = chr18_effects$courtships[1], xend = 2, yend = chr18_effects$courtships[2], colour = "black", linetype=2, size = 0.01) +
  	geom_segment(x = 2.5, y = chr18_effects$courtships[1], xend = 2.5, yend = chr18_effects$courtships[2], colour = "black", size = 1)
  	 	

#################################### all QTL plot ########################################

effect_all <- ggplot(QTL_data, aes(x = pos_all, y = prop_melp)) +
	geom_segment(x = -Inf, y = mp_mean, xend = Inf, yend = mp_mean, colour = "orange", size = 1, linetype=2) +
	geom_segment(x = -Inf, y = cp_mean, xend = Inf, yend = cp_mean, colour = "light blue", size = 1, linetype=2) +
	
	geom_jitter(alpha = 0.4, position=position_jitterdodge(jitter.width=0.4, dodge.width=0.85), colour = "dark gray", aes(size=total_court)) +
  	theme_bw() +
  	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  	theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
  	labs(x = "", y = "") +
  	scale_size(name = "") +
  	
  	# add means:
  	annotate("point",size = 3, shape = 15, x = 1, y = all_effects$courtships[1]) +
  	annotate("point",size = 3, shape = 15, x = 2, y = all_effects$courtships[8]) +
  	
  	geom_segment(x = 1, y = all_effects$plo[1], xend = 1, yend = all_effects$phi[1], colour = "black", size = 0.25) +
  	geom_segment(x = 1 - 0.05, y = all_effects$plo[1], xend = 1 +0.05, yend = all_effects$plo[1], colour = "black", size = 0.25) +
  	geom_segment(x = 1 - 0.05, y = all_effects$phi[1], xend = 1+0.05, yend = all_effects$phi[1], colour = "black", size = 0.25) +
  	
  	geom_segment(x = 2, y = all_effects$plo[8], xend = 2, yend = all_effects$phi[8], colour = "black", size = 0.25) +
  	geom_segment(x = 2 - 0.05, y = all_effects$plo[8], xend = 2 +0.05, yend = all_effects$plo[8], colour = "black", size = 0.25) +
  	geom_segment(x = 2 - 0.05, y = all_effects$phi[8], xend = 2+0.05, yend = all_effects$phi[8], colour = "black", size = 0.25)	+
  	
  	geom_segment(x = 1, y = all_effects$courtships[1], xend = 2, yend = all_effects$courtships[8], colour = "black", linetype=2, size = 0.01) +
  	geom_segment(x = 2.5, y = 	cp_mean , xend = 2.5, yend = all_effects$courtships[8], colour = "black", size = 1) +
	geom_segment(x = 2.5, y = 	all_effects$courtships[1] , xend = 2.5, yend = cp_mean, colour = "dark gray", size = 1)
	

#################################### no RADseq plot ######################################

# calculate proportion courtships towards melpomene for each male without RADseq data

without_RAD$total_court <- without_RAD$court_mp + without_RAD$court_cp
without_RAD$prop_melp <- without_RAD$court_mp / without_RAD$total_court

## make 'position marker' for genotypes at B locus, i.e. white (bb) and red (Bb) males
without_RAD$pos <- 1
without_RAD[without_RAD$forewing == "RED",]$pos<- 2

effect_no_RAD <- ggplot(without_RAD, aes(x = pos, y = prop_melp)) +
	
	geom_segment(x = -Inf, y = mp_mean, xend = Inf, yend = mp_mean, colour = "orange", size = 1, linetype=2) +
	geom_segment(x = -Inf, y = cp_mean, xend = Inf, yend = cp_mean, colour = "light blue", size = 1, linetype=2) +
	
	geom_jitter(alpha = 0.4, position=position_jitterdodge(jitter.width=0.4, dodge.width=0.85), colour = "dark gray", aes(size=total_court),show.legend = F) +
  	theme_bw() +
  	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  	theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
  	labs(x = "", y = "") +
  	scale_size(name = "") +
  	
  	# add means:
  	annotate("point",size = 3, shape = 15, x = 1, y = no_rad_effects$courtships[2]) +
  	annotate("point",size = 3, shape = 15, x = 2, y = no_rad_effects$courtships[1]) +
  	
  	geom_segment(x = 1, y = no_rad_effects$plo[2], xend = 1, yend = no_rad_effects$phi[2], colour = "black", size = 0.25) +
  	geom_segment(x = 1 - 0.05, y = no_rad_effects$plo[2], xend = 1 +0.05, yend = no_rad_effects$plo[2], colour = "black", size = 0.25) +
  	geom_segment(x = 1 - 0.05, y = no_rad_effects$phi[2], xend = 1+0.05, yend = no_rad_effects$phi[2], colour = "black", size = 0.25) +
  	
  	geom_segment(x = 2, y = no_rad_effects$plo[1], xend = 2, yend = no_rad_effects$phi[1], colour = "black", size = 0.25) +
  	geom_segment(x = 2 - 0.05, y = no_rad_effects$plo[1], xend = 2 +0.05, yend = no_rad_effects$plo[1], colour = "black", size = 0.25) +
  	geom_segment(x = 2 - 0.05, y = no_rad_effects$phi[1], xend = 2+0.05, yend = no_rad_effects$phi[1], colour = "black", size = 0.25)	+
  	
  	geom_segment(x = 1, y = no_rad_effects$courtships[2], xend = 2, yend = no_rad_effects$courtships[1], colour = "black", linetype=2, size = 0.01) +
  	geom_segment(x = 2.5, y = no_rad_effects$courtships[2], xend = 2.5, yend = no_rad_effects$courtships[1], colour = "black", size = 1)
##########################################################################################  	 	
