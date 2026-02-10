
# Silently loading all required packages when sourced
required_packages <- c("tidyverse", "dplyr", "rjags", "coda", "shiny", "ggplot2")

# Function to load packages
load_packages <- function() {
  missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "), 
         "\nInstall with: install.packages(c('", paste(missing, collapse = "', '"), "'))")
  }
  invisible(sapply(required_packages, library, character.only = TRUE, quietly = TRUE))
}

# Auto-load packages when sourced
load_packages()


#' 1) Function to turn the suitably structured dataset into a list object for jags mcmc sampling
make_model_list <- function(medalcounts, unique.model = FALSE) {
  df <- read.csv(medalcounts) %>% filter(competed == TRUE)
  
  m <- df$Medals.unique
  population <- as.numeric(df$total_pop_july)
  n_country <- nrow(df)
  
  medal_data_list <- list(M = m, N = population, n = n_country)
  
  if (!unique.model) {
    # M1 from team medals OR Medals.1
    m1_col <- grep("^(Medals\\.1\\.team|Medals\\.1)$", colnames(df), value = TRUE)[1]
    if (length(m1_col)) medal_data_list$M1 <- as.numeric(df[[m1_col]])
    
    # M2, M3, M4, M5, ... dynamically
    multi_cols <- grep("^Medals\\.[2-9][0-9]*$", colnames(df), value = TRUE)
    for (col in multi_cols) {
      k <- as.numeric(sub("^Medals\\.([0-9]+)$", "\\1", col))
      medal_data_list[[paste0("M", k)]] <- as.numeric(df[[col]])
    }
  }
  
  medal_data_list
}


#' 2) Function to run JAGS mcmc smapling - adjust as needed
jags_run <- function(jags, 
                     model_data.list, 
                     burn.in = 200000, 
                     model.run = 600000, 
                     thin = 20, 
                     n_chains = 4, 
                     model_data_vars = NULL,
                     var_names = c("a", "b", "p", "p1", "p2", "p3", "p4", "q2", "q3", "q4", "loglik"),
                     inits = NULL) {
  
  if (is.null(model_data_vars)) model_data_vars <- names(model_data.list)
  

  existing_vars <- model_data_vars[model_data_vars %in% names(model_data.list)]
  data <- model_data.list[existing_vars]
  
  cat("Using:", paste(existing_vars, collapse = ", "), "\n")
  
  model.jags <- jags.model(textConnection(jags), data = data, n.chains = n_chains, inits = inits)
  update(model.jags, n.iter = burn.in)
  coda.samples(model.jags, variable.names = var_names, n.iter = model.run, thin = thin)
}


# )
# init_mix <- function() {
#   list(
#     u = runif(3, 0.1, 0.9),
#     log_b = runif(3, 5, 12),
#     X1 = runif(1, 0.3, 0.7),
#     X2 = runif(1, 0.3, 0.7),
#     X3 = runif(1, 0.3, 0.7),
#     w = as.numeric(rdirichlet(1, rep(1, 3))),
#     z = sample(1:3, n, replace = TRUE)
#   )
# }

# 3) Function to tidy mcmc chains to go forward with analysis once satisfied convergence has been reached
tidy_mcmc <- function(mcmc, medalcounts) {
  medalcounts <- read.csv(medalcounts) %>% dplyr::filter(competed == TRUE)
  
  # Coerce to mcmc.list
  if (inherits(mcmc, "mcmc")) mcmc <- coda::mcmc.list(mcmc)
  
  # Extract samples - USE INTEGER INDICES NOT LOGICAL
  chain_dfs <- lapply(mcmc, as.data.frame)
  psims <- do.call(rbind, chain_dfs)
  
  cols <- colnames(psims)
  
  # Detect model type
  has_p1 <- any(grepl("^p[1-9][0-9]*\\[", cols))
  has_p <- any(grepl("^p\\[", cols))
  
  if (has_p1) {
    # Conditional model - use integer indices
    p_indices <- grep("^p[1-9][0-9]*\\[", cols)
    post_pc <- matrix(0, nrow(psims), nrow(medalcounts))
    
    for (idx in p_indices) {
      col_name <- cols[idx]
      k <- as.numeric(sub("^p([0-9]+)\\[.*", "\\1", col_name))
      country_idx <- as.numeric(sub(".*\\[([0-9]+)\\]$", "\\1", col_name))
      post_pc[, country_idx] <- post_pc[, country_idx] + k * psims[, idx]
    }
    
  } else if (has_p) {
    # Unique model - INTEGER INDICES ONLY
    p_indices <- grep("^p\\[", cols)
    post_pc <- psims[, p_indices, drop = FALSE]
    
  } else {
    stop("No p parameters found")
  }
  
  colnames(post_pc) <- medalcounts$iso_a3[1:ncol(post_pc)]
  post_pc  # Return as matrix/data.frame
}

# 4) Function to rank each posterior draw 
rank_mcmc <- function(processed_mcmc, # this one
                      medalcounts # for removing non-medal winners for ranking
                      ){
  
  medalcounts <- read.csv(medalcounts)
  
  country_data <- medalcounts  %>%
    filter(competed == TRUE)
  
  non_medals <- as.vector(country_data%>%
                          mutate(medal_winner = Medals.total > 0)%>%
                          filter(medal_winner == FALSE)%>%
                          select(iso_a3))


processed_mcmc[colnames(processed_mcmc) %in% non_medals$iso_a3] <- NA # setting nonmedal winning countries to 0 probability so as not to affect rank
  
  # ranking medal winners
  t(apply(processed_mcmc, 1, function(x) {
    rank(-x, ties.method = "first", na.last = "keep") # bets not average or whatever variant of ths
  }))
}

# 5) Function to summarise posterior draws for estimating mean and median posterior probabilities and credible intervals
post_prob_summary <- function(processed_mcmc) {
  # Calculate summaries
  median_p <- apply(processed_mcmc, 2, median, na.rm = TRUE)
  mean_p <- apply(processed_mcmc, 2, mean, na.rm = TRUE)
  cred_p <- apply(processed_mcmc, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Extract country names from columns
  country_names <- colnames(processed_mcmc)
  
  # Return properly structured dataframe
  data.frame(
    iso_a3 = country_names,
    p_mean = as.numeric(mean_p),
    p_median = as.numeric(median_p),
    p_credlow = as.numeric(cred_p[1, ]),
    p_credhigh = as.numeric(cred_p[2, ]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# 6) Ranking the distribution of ranked posterior draws on mean and median and calculating associated credible intervals
post_rank_rank <- function(post_ranks) {
  post <- as.data.frame(post_ranks) %>% select(where(~ !any(is.na(.x))))
  
  if (ncol(post) == 0) {
    return(data.frame(iso_a3 = character(0)))
  }
  
  average_ranks <- apply(post, 2, mean, na.rm = TRUE)
  median_ranks <- apply(post, 2, median, na.rm = TRUE)
  cred_ranks <- apply(post, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  
  data.frame(
    iso_a3 = names(post),
    rank_mean = as.numeric(average_ranks),
    rank_median = as.numeric(median_ranks),
    rank_credlow = as.numeric(cred_ranks[1, ]),
    rank_credhigh = as.numeric(cred_ranks[2, ]),
    rank_mean_beta = rank(average_ranks),
    rank_median_beta = rank(median_ranks),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# 7) Function to produce a df of all the results and summaries of the Bayesian ranking algorithm
results <- function(medalcounts,
                    probs.bayesian, 
                    ranks.bayesian){
  
  medalcounts <- read.csv(medalcounts)
  country_data <- medalcounts %>%  
    filter(competed == TRUE)%>%
    mutate(medal_type = case_when(
      is.na(Medals.total) ~ "Non-medalist",
      Medals.total == 0 ~ "Non-medalist", 
      Medals.total != Medals.1.team ~ "Multi-medal winners",
      TRUE ~ "Single medal winners"
    ),
    medal_winner = if_else(
      replace_na(Medals.total, 0) > 0,
      "medal_winner",
      "non_medal_winner"
    ))%>%
    select(country,
           games = slug_game,
           iso_a3, 
           population = total_pop_july,
           medal_type,
           medal_winner,
           medal_total = Medals.total,
           medals_1 = Medals.1.team,
           medals_2 = Medals.2, 
           medals_3 = Medals.3, 
           medals_4 = Medals.4, 
           medals_5 = Medals.5, # Wasn't any here but will leave for future proofing 
           medals_1.team = Medals.team) %>%# Recorded and modelled as single medals but were won in a team event
    arrange(-medal_total)%>%
    mutate(rank_total = row_number(), 
           medals.multi.winners = medals_2+medals_3+medals_4+medals_5)
  
  df1 <- left_join(country_data, probs.bayesian, by = "iso_a3")
  
  df2 <- left_join(df1, ranks.bayesian, by = "iso_a3")
  
  df2 %>% mutate(observed_mpm = ifelse(medal_total > 0, (medal_total / population) * 1e6, 0),
                 median_estimate_mpm =  (p_median * 1e6), 
                 mean_estimate_mpm =  (p_mean * 1e6), # just adding in as extra 
                 estimate_mpm_credlow = p_credlow*1e6,
                 estimate_mpm_credhigh = p_credhigh*1e6
  ) %>%
    mutate(rank_pc = rank(-observed_mpm, ties.method = "first", na.last = "keep"))
  
  
  
}


compute_rank_sig_matrix <- function(post_ranks, alpha = 0.05) {
  
  # all.countries <- colnames(post_ranks)
  countries <- colnames(post_ranks)[colSums(!is.na(post_ranks)) > 0]
  post_ranks <- post_ranks[, countries, drop = FALSE]
  
  pairs <- expand_grid(i = seq_along(countries),
                       j = seq_along(countries)) %>%
    filter(i < j) %>%
    mutate(row = countries[i],
           col = countries[j])

  post_rankdiffs <- pairs %>%
    rowwise() %>%
    mutate(prob = mean(post_ranks[, row] < post_ranks[, col])) %>% #proportion of posterior draws where country i is ranked better tahn country j
    ungroup()
  
  M <- matrix(TRUE,
              length(countries), length(countries),
              dimnames = list(countries, countries))
  
  for (k in seq_len(nrow(post_rankdiffs))) {
    i <- post_rankdiffs$row[k]
    j <- post_rankdiffs$col[k]
    p <- post_rankdiffs$prob[k]
    
    sig <- (p > 1 - alpha) | (p < alpha)
    
    M[i, j] <- !sig  # TRUE, not credibly different
    M[j, i] <- !sig
  }
  
  diag(M) <- TRUE
  
  countries <- rownames(M)
  
  # convert to data frame with one row per country
  sig_df <- data.frame(country = countries, M, row.names = NULL)
  
  return(sig_df)
}


#Running bayseina ranking in full
bayesrank_run <- function(medal_file, model, control = list(), alpha = 0.05) {
  ctrl <- modifyList(list(burn_in = 100000, model_run = 300000, thin = 30, n_chains = 4), control)
  
  unique_model <- model == "beta.unique"
  
  # DETECT max medals BEFORE making data list
  df <- read.csv(medal_file) %>% filter(competed == TRUE)
  
  if (!unique_model) {
    medal_cols <- grep("^Medals\\.[2-9][0-9]*$", colnames(df), value = TRUE)
    # ONLY count columns with actual data
    max_medals <- 1  # at least M1
    for (col in medal_cols) {
      if (sum(as.numeric(df[[col]]), na.rm = TRUE) > 0) {
        k <- as.numeric(sub("^Medals\\.([0-9]+)$", "\\1", col))
        max_medals <- max(max_medals, k)
      }
    }
  } else {
    max_medals <- NULL
  }
  
  datalist <- make_model_list(medal_file, unique.model = unique_model)
  jags_txt <- jags_model(model, max_medals = max_medals)
  
  cat("Max medals (non-zero):", max_medals, "\n")
  
  if (unique_model) {
    model_data_vars <- c("M", "N", "n")
    var_names <- c("a", "b", "p", "loglik")
  } else {
    model_data_vars <- c(paste0("M", 1:max_medals), "N", "n")
    var_names <- c("a", "b", "a1", "b1", "p", paste0("p", 1:max_medals), 
                   paste0("q", 2:max_medals), paste0("loglik", 1:max_medals), "loglik")
  }
  
  mcmc <- jags_run(jags_txt, datalist, burn.in = ctrl$burn_in, model.run = ctrl$model_run,
                   thin = ctrl$thin, n_chains = ctrl$n_chains, 
                   model_data_vars = model_data_vars, var_names = var_names)
  
  processed <- tidy_mcmc(mcmc, medal_file)
  post_ranks <- rank_mcmc(processed, medal_file)
  probs <- post_prob_summary(processed)
  ranks <- post_rank_rank(post_ranks)
  res <- results(medal_file, probs, ranks)
  sig <- compute_rank_sig_matrix(post_ranks, alpha)
  
  list(medal_file = medal_file, mcmc = mcmc, processed = processed, 
       post_ranks = post_ranks, probs = probs, ranks = ranks, results = res, sig = sig)
}

check_convergence <- function(bayesian_ranking) {
  
  mcmc<- bayesian_ranking$mcmc
  medal_file <- bayesian_ranking$medal_file
  
  required_libs <- c("shiny", "ggplot2", "dplyr", "coda")
  missing_libs <- required_libs[!sapply(required_libs, requireNamespace, quietly = TRUE)]
  if (length(missing_libs) > 0) {
    stop("Missing: ", paste(missing_libs, collapse = ", "), 
         "\nRun: install.packages(c('", paste(missing_libs, collapse = "', '"), "'))")
  }
  sapply(required_libs, library, character.only = TRUE)
  
  medalcounts <- read.csv(medal_file)
  country_data <- medalcounts %>% filter(competed == TRUE)
  
  if (!inherits(mcmc, "mcmc.list")) {
    mcmc <- as.mcmc.list(lapply(mcmc, as.mcmc))
  }
  n_chains <- length(mcmc)
  
  ui <- fluidPage(
    titlePanel("MCMC Convergence Check"),
    sidebarLayout(
      sidebarPanel(
        selectInput("chain", "Chain:", choices = setNames(1:n_chains, paste("Chain", 1:n_chains))),
        radioButtons("param_type", "Show:", c("Global" = "global", "Country" = "country")),
        conditionalPanel(condition = "input.param_type == 'country'",
                         selectInput("country", "Country:", choices = setNames(1:nrow(country_data), country_data$country))
        ),
        selectInput("param", "Parameter:", choices = NULL),
        uiOutput("rhat_info"),
        uiOutput("country_info"),
        width = 3
      ),
      mainPanel(
        plotOutput("trace", height = "200px"),
        plotOutput("density", height = "300px"),
        plotOutput("autocorr", height = "200px")
      )
    )
  )
  
  server <- function(input, output, session) {
    
    
    axis_limits <- reactive({
      req(input$param)
      chain_vals <- as.numeric(mcmc[[1]][, input$param])
      posterior_range <- range(chain_vals, na.rm = TRUE)
      pad <- (posterior_range[2] - posterior_range[1]) * 0.05
      
      if (grepl("^q2\\[[0-9]+\\]$", input$param)) {
        list(x = c(0, 0.2), y = c(0, 50))
      } else if (grepl("^p\\[[0-9]+\\]$", input$param)) {
        list(x = c(0, 2e-6), y = c(0, 4e6))
      } else {
        list(x = c(posterior_range[1] - pad, posterior_range[2] + pad), y = NULL)
      }
    })
    
    
    observe({
      req(input$param_type)
      all_params <- colnames(mcmc[[1]])
      
      # Get previous param
      prev_param <- input$param
      
      if (input$param_type == "global") {
        params <- all_params[!grepl("\\[[0-9]+\\]$", all_params)]
        # Try to keep same global param
        selected <- if (prev_param %in% params) prev_param else params[1]
        updateSelectInput(session, "param", choices = params, selected = selected)
      } else {
        req(input$country)
        country_idx <- as.numeric(input$country)
        params <- all_params[grepl(paste0("\\[", country_idx, "\\]$"), all_params)]
        
        if (length(params) > 0) {
          # Try to keep same parameter type (p, q2, etc.) for new country
          param_prefix <- sub("\\[.*\\]$", "", prev_param)
          candidate <- paste0(param_prefix, "[", country_idx, "]")
          selected <- if (candidate %in% params) candidate else params[1]
          
          updateSelectInput(session, "param", choices = params, selected = selected)
        }
      }
    })
    
    output$rhat_info <- renderUI({
      req(input$param)
      if (n_chains >= 2) {
        param_mcmc <- lapply(mcmc, function(ch) as.mcmc(ch[, input$param, drop = FALSE]))
        rhat <- round(gelman.diag(as.mcmc.list(param_mcmc))$psrf[1], 3)
        status <- if (rhat < 1.1) "Converged" else "Check convergence"
        bg <- if (rhat < 1.1) "#d4edda" else "#fff3cd"
        border <- if (rhat < 1.1) "#28a745" else "#ffc107"
        
        div(style = paste0("padding: 12px; margin: 10px 0; border-radius: 4px; background: ", bg, 
                           "; border-left: 4px solid ", border),
            h5(strong("Gelman-Rubin R-hat: "), rhat),
            p(status))
      }
    })
    
    output$country_info <- renderUI({
      req(input$param_type == "country", input$country, input$param)
      
      country_num <- as.numeric(input$country)
      if (country_num > nrow(country_data)) return(NULL)
      
      row <- country_data[country_num, ]
      div(style = "padding: 12px; margin: 10px 0; border-radius: 4px; background: #f8f9fa; 
                 border-left: 4px solid #6c757d;",
          h5(strong(row$country)),
          p(paste("Population:", format(row$total_pop_july, big.mark = ","))),
          p(paste("Medals:", row$Medals.total, "(", row$Medals.unique, "reps)")),
          p(paste("Multi-medal:", row$Medals.2 + row$Medals.3 + row$Medals.4 + row$Medals.5)))
    })
    
    chain_vals <- reactive({
      req(input$param, input$chain)
      as.numeric(mcmc[[as.numeric(input$chain)]][, input$param])
    })
    
    output$trace <- renderPlot({
      df <- data.frame(Iter = 1:length(chain_vals()), Val = chain_vals())
      ggplot(df, aes(Iter, Val)) + geom_line(color = "#2E86AB", linewidth = 1) +
        labs(title = paste("Trace:", input$param), y = input$param)+ theme_minimal(14)
    })
    
    output$density <- renderPlot({
      df <- data.frame(Val = chain_vals())
      lims <- axis_limits()
      ggplot(df, aes(Val)) +
        geom_density(fill = "#2E86AB", alpha = 0.7, color = NA) +
        geom_rug() + 
        labs(title = paste("Density:", input$param),
                          subtitle = sprintf("Mean=%.4g (n=%d)", mean(df$Val), nrow(df))) +
        coord_cartesian(xlim = lims$x, ylim = lims$y) + 
        theme_minimal(14)
    })
    
    output$autocorr <- renderPlot({
      acf_vals <- acf(as.mcmc(chain_vals()), plot = FALSE)$acf[-1]
      df <- data.frame(Lag = 1:length(acf_vals), ACF = acf_vals)
      
      ggplot(df, aes(Lag, ACF)) + geom_hline(yintercept = 0.2, color = "red", linetype = 2) +
        geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8, width = 0.8) +
        labs(title = "Autocorrelation", subtitle = "Drop below red line quickly") + theme_minimal(14)
    })
  }
  
  shinyApp(ui, server)
}



# JAGS models

jags_model <- function(model, max_medals = NULL){
  #update this to original RTE brainstorm model
  if (model == "beta.unique"){
    
    return("model {
  for (c in 1:n) {
  
    # Likelihood models for different types of medal winner 
    M[c] ~ dpois(N[c] * p[c])
    
    
    #Prior for probability of being a unique medal winner
    p[c] ~ dbeta(a, b) T(10^-9, 1)
    #log pointwise predictive density
    
    loglik[c] <- logdensity.pois(M[c], N[c]*p[c]) # pointwise log-likelyhood for each country 
  }
  
  # Hyperpriors for beta prior for probability of being a medal winner
  a ~ dunif(0, 1)
  b ~ dunif(10^4, 10^8) 
}")
  }
  else if (model == "beta.conditional") { #This has been done with teh help of AI
    
    # Auto-detect max_medals if not provided
    if (is.null(max_medals)) max_medals <- 4
    
    # Build likelihood section dynamically
    likelihood <- paste0(
      sapply(1:max_medals, function(k) {
        sprintf("    M%d[c] ~ dpois(N[c] * p%d[c])", k, k)
      }), collapse = "\n"
    )
    
    # Build derived probabilities dynamically
    # p1[c] = p[c]*(1-q2)
    # p2[c] = p[c]*q2*(1-q3)
    # p3[c] = p[c]*q2*q3*(1-q4)
    # ...
    derived_probs <- c(
      "    p1[c] = p[c]*(1-q2)",
      sapply(2:(max_medals-1), function(k) {
        q_prod <- paste0("q", 2:k, collapse = "*")
        sprintf("    p%d[c] = p[c]*%s*(1-q%d)", k, q_prod, k+1)
      }),
      sprintf("    p%d[c] = p[c]*%s", max_medals, paste0("q", 2:max_medals, collapse = "*"))
    )
    derived_probs <- paste(derived_probs, collapse = "\n")
    
    # Build loglik sum dynamically
    loglik_lines <- sapply(1:max_medals, function(k) {
      sprintf("loglik%d[c] <- logdensity.pois(M%d[c], N[c]*p%d[c])", k, k, k)
    })
    loglik_lines <- paste("    ", loglik_lines, collapse = "\n")
    loglik_sum <- paste0("loglik[c] <- ", paste0("loglik", 1:max_medals, "[c]", collapse = " + "))
    
    # Build q priors dynamically
    q_priors <- c(
      "  q2 ~ dbeta(a1, b1)",
      sapply(3:max_medals, function(k) sprintf("  q%d ~ dbeta(1, 1)", k))
    )
    q_priors <- paste(q_priors, collapse = "\n")
    
    # Assemble full model
    model_string <- sprintf("model {
  for (c in 1:n) {
  
    # Likelihood models for different types of medal winner 
%s
    
    # Prior for probability of being a medal winner
    p[c] ~ dbeta(a, b) T(10^-9, 1)

    # Derived probabilities of individual medals 
%s

%s
    %s
  }
  
  # Prior probability of transitioning from m medals to m + 1 medals 
%s
  
  # Hyperpriors for beta prior for probability of being a medal winner
  a ~ dunif(0, 1)
  b ~ dunif(10^4, 10^8) 
  
  a1 = 1
  b1 = 1
  f = 1
}", likelihood, derived_probs, loglik_lines, loglik_sum, q_priors)
    
    return(model_string)
  }
  else if (model == "beta.conditional.countryspecific") {
    
    # Auto-detect max_medals if not provided
    if (is.null(max_medals)) max_medals <- 4
    
    # Build likelihood section dynamically (SAME as conditional)
    likelihood <- paste0(
      sapply(1:max_medals, function(k) {
        sprintf("    M%d[c] ~ dpois(N[c] * p%d[c])", k, k)
      }), collapse = "\n"
    )
    
    # Build derived probabilities - ONLY DIFFERENCE: q2 becomes q2[c]
    derived_probs <- c(
      "    p1[c] <- p[c]*(1-q2[c])",  # q2[c] not q2
      sapply(2:(max_medals-1), function(k) {
        q_prod <- paste0("q", 2:k, collapse = "*")
        q_prod <- gsub("q2", "q2[c]", q_prod, fixed = TRUE)  # Replace q2 with q2[c]
        sprintf("    p%d[c] <- p[c]*%s*(1-q%d)", k, q_prod, k+1)
      }),
      {
        q_prod <- paste0("q", 2:max_medals, collapse = "*")
        q_prod <- gsub("q2", "q2[c]", q_prod, fixed = TRUE)
        sprintf("    p%d[c] <- p[c]*%s", max_medals, q_prod)
      }
    )
    derived_probs <- paste(derived_probs, collapse = "\n")
    
    # Build loglik (SAME as conditional)
    loglik_lines <- sapply(1:max_medals, function(k) {
      sprintf("loglik%d[c] <- logdensity.pois(M%d[c], N[c]*p%d[c])", k, k, k)
    })
    loglik_lines <- paste(loglik_lines, collapse = "\n")
    loglik_sum <- paste0("loglik[c] <- ", paste0("loglik", 1:max_medals, "[c]", collapse = " + "))
    
    # Build q priors - q2 is now INSIDE loop, q3+ OUTSIDE
    q_priors <- sapply(3:max_medals, function(k) {
      sprintf("  q%d ~ dbeta(1, 1)", k)
    })
    q_priors <- paste(q_priors, collapse = "\n")
    
    # Assemble model
    model_string <- sprintf("model {
  for (c in 1:n) {
    # Likelihood
%s
    
    # Priors
    p[c] ~ dbeta(a, b) T(10^-10, )
    q2[c] ~ dbeta(a1, b1) T(10^-10, )  # COUNTRY-SPECIFIC

    # Derived probabilities
%s

    # Log-likelihood
%s
    %s
  }
  
  # Global q priors
%s
  
  # Hyperpriors
  a ~ dunif(0, 1)
  b ~ dunif(10^4, 10^8)
  f ~ dunif(0.025, 0.05)
  a1 = f * b1
  b1 ~ dunif(0, 1000)
}", likelihood, derived_probs, loglik_lines, loglik_sum, q_priors)
    
    return(model_string)
  }
  else if(model == "logit-normal"){
    return("model {
  for (i in 1:n) {

    # Likelihood models for different types of medal winners
    M1[i] ~ dpois(N[i] * p1[i])
    M2[i] ~ dpois(N[i] * p2[i])
    M3[i] ~ dpois(N[i] * p3[i])
    M4[i] ~ dpois(N[i] * p4[i])

    # Logit-normal prior on transformed scale for the medal-winning probability:
    theta[i] ~ dnorm(mu, tau)
    p[i] <- 1 / (1 + exp(-theta[i]))

    # Derived probabilities for medal counts
    p1[i] <- p[i] * (1 - X1)
    p2[i] <- p[i] * X1 * (1 - X2)
    p3[i] <- p[i] * X1 * X2 * (1 - X3)
    p4[i] <- p[i] * X1 * X2 * X3
  }

  # Prior for transition probabilities between medal counts
  X1 ~ dbeta(1, 1)
  X2 ~ dbeta(1, 1)
  X3 ~ dbeta(1, 1)

  # Hyperpriors for the logit-normal prior on p[i]
  mu ~ dnorm(-15, 0.3)      
  tau ~ dgamma(0.001, 0.001)   
}")
    
  }
else if(model == "mixture-beta"){
    return(
      "model {
#defining mixture weights for each component
epsilon <- 1.0E-6
  w[1:3] ~ ddirch(alpha_w[]) 
 
# stick breaking for a[i] to try and fix identifiability issues
for (k in 1:3) {
    u[k] ~ dunif(0, 1)
}
  
  a[1] <- u[1]
  a[2] <- u[1] + (1 - u[1]) * u[2]
  a[3] <- u[1] + (1 - u[1]) * u[2] + (1 - (u[1] + (1 - u[1]) * u[2])) * u[3]

# hyperpriors for beta distribution indexed by component of mixture
for(k in 1:3) {
  log_b[k] ~ dnorm(log(1e6), 1) T(0, 20)
  b[k] <- exp(log_b[k]) #trying to make mcmc smpling more efficient and explore entire parameter space
}

  for (i in 1:n) {
  #decide which mixture component to use
  z[i] ~ dcat(w[])
  
    # Likelihood models for different types of medal winner 
    M1[i] ~ dpois(N[i] * p1[i])
    M2[i] ~ dpois(N[i] * p2[i])
    M3[i] ~ dpois(N[i] * p3[i])
    M4[i] ~ dpois(N[i] * p4[i])
    
    
    #Prior for probability of being a medal winner, adding epsilon to a instead of truncation
    p[i] ~ dbeta(epsilon + a[z[i]], b[z[i]]) 

    # Derived probabilities of individual medals 
    p1[i] = p[i]*(1-X1) #probabaility of exactly 1 medal = prob of a medal - prob of 2 or more
    p2[i] = p[i]*X1*(1-X2) #probabaility of exactly 2 medals = prob of 2 or more - prob of 3 or more
    p3[i] = p[i]*X1*X2*(1-X3)
    p4[i] = p[i]*X1*X2*X3

    
  }
  # Prior probability (uninformed) of tranistioning from m medals to m + 1 medals 
    X1 ~ dbeta(1, 1) 
    X2 ~ dbeta(1, 1)
    X3 ~ dbeta(1, 1)
  
  #hyperpriors for weights
   alpha_w[1] <- 1
   alpha_w[2] <- 1
   alpha_w[3] <- 1

}"
    )
  }
else{
  stop("Model type not recognized. Choose `beta`, `logit-normal` or `mixture-beta`.")
}
}
