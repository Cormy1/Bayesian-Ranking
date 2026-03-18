#Sourcing function files 
source("https://raw.githubusercontent.com/Cormy1/Bayesian-Ranking/main/R/bayesrank_utils.R") 

# Data file pre-processed to required structure - country, iso_a3, total_pop_july, competed, Medals.1.team (individual single medallists and team medals), Medals.2,Medals.3,Medals.4 ,... Medals.n, Medals.unique (number of competitors to win medals)
medal.file <- read.csv("data/Medalcounts_paris-2024.csv")
names(medal.file)

#run bayesian ranking - retuns list with mcmc sims along with all results and calculations pertaining to ranking
br.cc <- bayesrank_run(medal_file = medal.file,
                       model = "beta.conditional", #select model to run - "beta.unique", "beta.conditional", beta.conditional.countryspecific
                       control = list(
                         burn_in   = 10000,
                         model_run = 30000, # run longer 
                         thin      = 30,
                         n_chains  = 1 #recommend running 4
                         )
) 

#run this function to check chains 
check_convergence(br.cc)


br.cc$posterior_ranks #posterior distribution of ranks

br.cc$posterior_prob #posterior distribution of medals per capita

br.cc$probs #summary of posterior distribution of medals per capita

br.cc$ranks

br.cc$results #table of bayesian ranks and other common ranking metrics, medal counts population,....

br.cc$sig #matrix of countries with logical TRUE/FALSE indicating wheter significantly different in Bayesian rank , NA if not medal winner

