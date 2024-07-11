#SENSITIVITY ANALYSIS with randomly generated binary variable for disease

#binary variables generated using rbinom() with different probabilities for success


#perform 2-sample instrumental variable regression
#Two-Sample Two-Stage Least Squares
#instrument: AB; exposure: MB (indexes, genuses, species); outcome: randomly generated binary variable
# indicating disease incidence
# AB -> MB regression in MB cohort
# AB -> disease regression in EGC cohort (200k; excluding MB cohort)
# AB usage 2004-2014; disease incidence 2014-2022; exclude recent (6 months) AB users
# sensitivity analysis with disease data randomly generated (expect no causal effect of MB on disease)


####DATA + LIBRARIES####
setwd("C:/Users/Administrator/Documents/01_projektid/09_AB_instrument")
load("02_input/andmed.RData")
load("02_input/andmed_200k.RData")

prevotella_bacteroides = readRDS("02_input/metaphlan_PB_ratio.rds")
#set a seed for starting the random number generation
set.seed(38572935)

library(dplyr)
library(Epi)
library(ggplot2)
library(xlsx)


####PREPROCESSING####
exclude = andmed$skood
EGC = andmed_200k
EGC = EGC %>% filter(!skood %in% exclude)
MB = merge(prevotella_bacteroides, andmed, by="skood")
rm(andmed_200k, andmed, exclude)
# two datasets: MB and EGC

#gender-variable
MB$gender = factor(MB$gender, levels = c(0,1), labels = c("Man", "Woman"))

EGC$gender = EGC$gender-1
EGC$gender = factor(EGC$gender, levels = c(0,1), labels = c("Man", "Woman"))

#age
EGC$age = 2014- EGC$birthyear


##generate random binary variables indicating disease incidence
# a vector of different probabilities to be used for generating data
disease_prob = seq(0.025,0.975,0.025)
#names for the randomly generated variables (name contains the probability of success in generated data; e.g. disease proportion)
rand_disease_names = paste("random_", disease_prob, sep="")

#create an empty dataset for the generated data
rand_disease_data = matrix(NA, nrow = nrow(EGC), ncol = length(rand_disease_names))
colnames(rand_disease_data) <- rand_disease_names

for (i in 1:length(rand_disease_names)){
  #for each disease-probablility generate data from binomial distribution with the defined probablility
  rand_disease_data[, rand_disease_names[i]] = rbinom(nrow(EGC),1,disease_prob[i])
}

#merge the randomly generated disease-data with the original EGC dataset
EGC = cbind(EGC, rand_disease_data)

#################
#### FILTERS ####
#################



#which diseases to include
diseases_inc = c(names(EGC[match("random_0.025",names(EGC)):match("random_0.975",names(EGC))]))

#exposures
exposures = "PB_ratio"

#minimum number of cases (for main analysis - does not apply for sub-group analysis)
min_cases = 50

#max number of AB used (max in MB cohort 41)
AB_max = 5

#age exclusion (MB sample age range 23-89)
min_age = 23
max_age = 50


#####################
#### FILTERS END ####
#####################

####APPLY FILTERS TO DATASETS###
# max N of AB used
EGC = EGC %>% filter(AB_n < AB_max+1)
MB = MB %>% filter(AB_n < AB_max+1)

#age
EGC = EGC %>% filter(age >= min_age) %>% filter(age <= max_age)
MB = MB %>% filter(Age_at_MBsample >= min_age) %>% filter(Age_at_MBsample <= max_age)


####ANALYSIS####


for (j in 1:length(exposures)){
  exposure = exposures[j]
  
  #make empty data-frame for storing the results 
  #each exposure will get it's own table with results for each disease as rows
  res_inc = matrix(NA, nrow = length(diseases_inc), ncol = 10)
  colnames(res_inc) <- c("N", "N_ctrl", "beta", "CI_low", "CI_high", "p_value", "exposure", "beta_exp", "p_value_exp", "group")
  rownames(res_inc) <- diseases_inc

  
  #count the number of incident cases (excluding prevalent cases from dataset beforehand)
  #make a for-loop separately for each disease
  for (i in 1:length(diseases_inc)){
    #no filtering based on prevalent-disease
    temporary = EGC
    #count the number of incident cases (when prevalent excluded)
    res_inc[diseases_inc[i], "N"] = table(temporary[diseases_inc[i]])["1"]
    #count the number of controls (when prevalent excluded)
    res_inc[diseases_inc[i], "N_ctrl"] = table(temporary[diseases_inc[i]])["0"]
  }
  
  #make into dataframe
  res_inc = as.data.frame(res_inc)
 
  # step1 - y
  # step2 - x
  
  # step2 of 2-sample-2-stage least squares (association between instrument and exposure; same for each disease)
  step2 = lm(scale(MB[, match(exposure, names(MB))])~AB_n, data=MB)
  
  #for loop for step1 and combining the two steps for each disease
  for (i in 1:length(diseases_inc)){
    
    #no filtering based on prevalent-disease
    temporary = EGC
    
    #pursue step1 (association between instrument and outcome (disease)) 
    step1 = glm(temporary[, match(diseases_inc[i], names(temporary))]~AB_n, family="binomial", data=temporary)
    
    #calculate beta_tsls (two sample two stage least squares) (exposure = microbiome; outcome = disease)
    #beta_exposure_outcome = beta_outcome / beta_exposure
    beta_tsls = step1$coefficients[2]/step2$coefficients[2]
    ##extract information from models to calculate sd for beta_tsls
    #standard deviation of step1 coefficient (instrument on disease beta)
    step1_sd = summary(step1)$coefficients["AB_n","Std. Error"]
    #standard deviation of step2 coefficient (instrument on MB beta)
    step2_sd = summary(step2)$coefficients["AB_n","Std. Error"]
    # step2 coefficient (instrument on MB)
    step2_beta = summary(step2)$coefficients["AB_n","Estimate"]
    
    #sd of tstsls beta based on Pacini and Windmeijer (2016)
    sd_tsls = sqrt((step1_sd**2 + beta_tsls**2 * step2_sd**2)/(step2_beta**2))
    
    #calculate p-value based on z-score
    z_score = beta_tsls/sd_tsls
    p_value = 2*pnorm(abs(z_score), lower.tail = FALSE)
    
    #95% confidence interval for beta_tsls
    CI_low = beta_tsls + 1.96*sd_tsls
    CI_high = beta_tsls - 1.96*sd_tsls
    
    #store values to the results table
    res_inc[diseases_inc[i], "beta"] = beta_tsls
    res_inc[diseases_inc[i], "CI_low"] = CI_low
    res_inc[diseases_inc[i], "CI_high"] = CI_high
    res_inc[diseases_inc[i], "p_value"] = p_value
    
  }
  
  #name of exposure
  res_inc$exposure = exposure
  #beta from step2 (association between instrument and exposure) - strength of the instrument
  res_inc$beta_exp = summary(step2)$coefficients["AB_n","Estimate"]
  #p-value from step2 (association between instrument and exposure) - strength of the instrument
  res_inc$p_value_exp = summary(step2)$coefficients["AB_n","Pr(>|t|)"]
  #group (currently using all data, but later might want to do sub-analyses by gender (or other groupings))
  res_inc$group = "All"
  #names of diseases (have to be same as row-names)
  res_inc$disease = diseases_inc
  
  
  ####export the results table####
  setwd("C:/Users/Administrator/Documents/01_projektid/09_AB_instrument/03_output/04_sensitivity_randbinary/")
  write.xlsx(res_inc, paste(paste("sensitivity_randbinary", exposure, min_cases, AB_max, min_age, max_age, sep="_"), ".xlsx", sep=""))
  
  ####PLOT THE RESULT####
  pdf(file = paste(paste("sensitivity_randbinary", exposure, min_cases, AB_max, min_age, max_age, sep="_"), ".pdf", sep=""), 
      height=15, width=6)
 
  p2 = ggplot(res_inc %>% filter(group=="All"), aes(x=disease, y=beta, ymin=CI_low, ymax=CI_high),asp=1) + 
    geom_pointrange(fatten=2) + 
    geom_hline(yintercept = 0, linetype=2) +
    coord_flip(ylim=c(-1,1)) +
    theme_bw() +
    theme(legend.position="none",aspect.ratio = 2.5)+
    ylab("beta") +
    xlab(" ") + 
    ggtitle(paste("res_inc", "; exposure: ", exposure, "; min cases: ", min_cases, 
                  "; \n max AB: ", AB_max, "; age min: ", min_age, "; age max: ",
                  max_age, "; \n N(EGC): ", nrow(EGC), "; N(MB): ", nrow(MB),
                  "; \n instrument beta: ", round(res_inc$beta_exp, 3), 
                  ", and p-value: ", round(res_inc$p_value_exp, 3), sep=""))
  print(p2)
  
  dev.off()
}
