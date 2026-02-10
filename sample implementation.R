#Sourcing function files 
source("https://raw.githubusercontent.com/Cormy1/Bayesian-Ranking/main/R/bayesrank_utils.R") 

medal.file <- "data/Medalcounts_paris-2024.csv"

#run bayesian ranking - retuns list with mcmc sims along with all reuslts and calculations pertaining to ranking
br.cc <- bayesrank_run(medal_file = medal.file,
                       model = "beta.conditional", #select model to run - "beta.unique", "beta.conditional", beta.conditional.countryspecific
                       control = list(
                         burn_in   = 10000,
                         model_run = 30000,
                         thin      = 30,
                         n_chains  = 4)
) 

#run this function to check chains 
check_convergence(br.cc)
