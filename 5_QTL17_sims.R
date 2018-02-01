# Simulations to compare significant and non-sgnificant effect sizes:
# R. Merrill Nov 2017


# QTL Chr 17
##########################################################################################
# Load equired R packages:

library(qtl)
library(simr)
library(lme4)

##########################################################################################
args<-commandArgs(TRUE)
effect<-ifelse(is.na(args[1]),1,args[1])
effect<-as.numeric(effect)

qtl17_effect <- c(0.208, 0.40, 0.577, 0.749, 0.913, 1.076, 1.238, 1.399, 0.8587)
qtl17_sim_percent <- c(5, 10, 15, 20, 25, 30, 35, 40, "recorded")
effect_data_frame <- data.frame(qtl17_effect, qtl17_sim_percent)

sim_percent <- effect_data_frame[effect_data_frame$qtl17_effect == effect,]$qtl17_sim_percent 
message(paste("Running with effect size:",effect, "relating to", sim_percent , "of parental difference" ))

# Get data

# from qtl_significance.R 
QTL_data <- read.csv("../derived_data/qtl_data.csv")


###########################################################################################
# Determine means for parental species, i.e. melpomene and cydno males assayed in the 
# same way as the hybrids (same as in effect.size.R)
###########################################################################################
court_data <- read.csv("../raw_data/species_brood_data.csv")
mp_mean <- plogis(coef(summary(glmer(cbind(court_mp, court_cp) ~  1 + (1|insectary_id),  family = binomial,  data = court_data[court_data$Type == "MP",] )))[ , "Estimate"])
cp_mean <- plogis(coef(summary(glmer(cbind(court_mp, court_cp) ~  1 + (1|insectary_id),  family = binomial,  data = court_data[court_data$Type == "CP",] )))[ , "Estimate"])
parental_diff <- mp_mean - cp_mean
##########################################################################################


# FUNCTION TO GET LOD SCORE FROM GLMM (SAME AS THAT USED FOR QTL ANALYSIS)
#########
get_lod_score<-function(locus_genotypes, courtships, id) {
  model      <- glmer(courtships ~ locus_genotypes + (1|id),  family = binomial )
  null_model <- glmer(courtships ~ 1 + (1|id),  family = binomial )
  lod_score<-((logLik(model) - logLik(null_model))/log(10))[1]
  return(lod_score)
}


# Function to get effect size from GLMM in form of difference between genotypes in proportion 
# of courtships directed towards melpomene or cydno females 
get_effect<-function(locus_genotypes, courtships, id) {
	model <- glmer(courtships ~ locus_genotypes + (1|id),  family = binomial )
	model_intercept <- coef(summary(model))[1 , "Estimate"]
	#print(model_intercept)
	model_estimate <- coef(summary(model))[2 , "Estimate"]
	#print(model_estimate)
	absolute_effect <- (plogis(model_intercept + model_estimate)) - plogis(model_intercept)
	return(absolute_effect)
}

########

# Function to simulate data and determine effect

simulate_data <- function (glmm_model, num_sim) {
	good_simulation <- "no"
		# Simulate data, test and retrieve LOD score
		simulated_data <- simulate(glmm_model, num_sim)
		options(warn=2)		# n.b options (2) changes warnings to errors or (1) vice-versa
		test.mod1 <- try( lod <- get_lod_score(QTL_data$geno17, cbind(simulated_data$sim_1[,1], simulated_data$sim_1[,2]), QTL_data$id ))
		options(warn=1)
		## ... but discard if model does not converge
		if(class(test.mod1) == "try-error") { good_simulation <- "no" } else { good_simulation <- "yes" }
		lod <- get_lod_score(QTL_data$geno17, cbind(simulated_data$sim_1[,1], simulated_data$sim_1[,2]), QTL_data$id )
		effect <- get_effect(QTL_data$geno17, cbind(simulated_data$sim_1[,1], simulated_data$sim_1[,2]), QTL_data$id )
	out_list<-list(lod, effect, good_simulation) 
	return(out_list)
}

#PRODUCE MODEL FOR COURTSHIPS FROM REAL DATA USING lme4:
rm <- glmer(cbind(court_mp, court_cp) ~ geno17  + (1|id),  family = binomial,  data = QTL_data )
#summary(rm)

model_intercept <- coef(summary(rm))[1 , "Estimate"]
model_estimate <- coef(summary(rm))[2 , "Estimate"]

# Check difference between parents
#plogis(model_intercept + model_estimate) - plogis(model_intercept)



##########################################################################################

# SET-UP DATA.FRAME TO STORE SIMULATION RESULTS

# Set up effect sizes to be tested and number of simulation for each
# Number of simulations per 'treatment': to be changed (test = 10)

# Set up data-frame to store simulation data:
sim <- NA
sim_run <- data.frame(sim)
sim_run$simulated_effect <- NA

sim_run$sim_effect <- NA
sim_run$sim_lod <- NA
#sim_run$simulated_effect <- NA
#sim_run$white_mean <- NA
#sim_run$diff_mean <- NA
sim_run$good_sim <- NA

sim_run <- sim_run[0,]


# RUN SIMULATIONS

start_time<-Sys.time()
#num_sim <- 2

num_sim <- 10000


# Effect size calculated as: (plogis(model_intercept  + *effect*) - plogis(model_intercept ))/parental_diff

#5	0.208
#10	0.40
#15	0.577
#20	0.749
#25	0.913
#30	1.076
#35	1.238
#40	1.399
#recorded	0.8587 



#for (x in 1:length(effect)) {
	fixef(rm)["geno17B"] <- effect
	#print(effect[x])
	for (i in 1:num_sim) {
		if(i %% 500==0) { print(i) }
		sim <- i
		simulated_effect <- as.vector(sim_percent)
		out_list <- simulate_data(rm, 1)
		sim_lod <- out_list[1]
		sim_effect <- out_list[2]
		good_sim <- out_list[3]
		sim_run[nrow(sim_run) + 1,] <- c(sim, simulated_effect, sim_effect, sim_lod,  good_sim) 
	}
#}
#print( Sys.time() - start_time)

sim_run$par_diff <- sim_run$sim_effect / parental_diff

write.csv(sim_run, paste0("../derived_data/beavis_sim_results/QTL_17_sim_results_", as.character(sim_percent), ".csv"),  row.names = FALSE)
