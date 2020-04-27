plot.covidcast_result <- function(obj, cases)
  viz_melanie(obj, cases)

## COVID viz 
viz_melanie <- function(obj, cases) {
  library(tidyverse)
  result  <- obj$result
  fit     <- obj$extracted
  datList <- obj$config

  summary_fit <- rstan::summary(result)

  outcomeGen <- function(idx) {
    as.data.frame(fit[[idx]]) %>% 
    gather(key = day, value = estim) %>%
    mutate(outcome = idx)
  }

  new_inf <- outcomeGen("new_inf")
  new_sym <- outcomeGen("new_sym")
  new_hos <- outcomeGen("new_hos")
  new_die <- outcomeGen("new_die")

  diag_inf <- outcomeGen("diag_inf")
  diag_sym <- outcomeGen("diag_sym")
  diag_hos <- outcomeGen("diag_hos")
  diag_all <- outcomeGen("diag_all")

  occur_cas <- outcomeGen("occur_cas")
  occur_hos <- outcomeGen("occur_hos")
  occur_die <- outcomeGen("occur_die")

  obs_cas <- select(cases, NEW_COVID_CASE_COUNT) %>%
    mutate(day = seq(1, n(), 1), outcome = "obs_cas") %>% 
    rename(estim = NEW_COVID_CASE_COUNT) %>%
    select(day, estim, outcome) 
  
  obs_hos <- select(cases, HOSPITALIZED_CASE_COUNT) %>%
    mutate(day = seq(1, n(), 1), outcome = "obs_hos") %>% 
    rename(estim = HOSPITALIZED_CASE_COUNT) %>%
    select(day, estim, outcome) 

  obs_die <- select(cases, DEATH_COUNT) %>%
    mutate(day = seq(1, n(), 1), outcome = "obs_die") %>% 
    rename(estim = DEATH_COUNT) %>%
    select(day, estim, outcome) 

  new <- rbind(new_inf, new_sym, new_hos, new_die) %>%
    group_by(day, outcome) %>%
      summarise(median = median(estim), 
      lo = quantile(estim, 0.025), 
      hi = quantile(estim, 0.975)) %>%
    ungroup() %>%
    mutate(day = as.numeric(substr(day, start = 2, stop = 4)) - 10) %>%
    arrange(day)

  diag <- rbind(diag_inf, diag_sym, diag_hos, diag_all) %>%
    group_by(day, outcome) %>%
    summarise(median = median(estim), 
              lo = quantile(estim, 0.025), 
              hi = quantile(estim, 0.975)) %>%
    ungroup() %>%
    mutate(day = as.numeric(substr(day, start = 2, stop = 4)) - 10) %>%
    arrange(day)


  fit_to_data <- rbind(occur_cas, occur_hos, occur_die) %>%
    group_by(day, outcome) %>%
    summarise(median = median(estim), 
              lo = quantile(estim, 0.025), 
              hi = quantile(estim, 0.975)) %>%
    ungroup() %>%
    mutate(day = as.numeric(substr(day, start = 2, stop = 4)) - 10) %>%
    arrange(day)


  input_data <- rbind(obs_cas, obs_hos, obs_die) %>%
                group_by(day, outcome) %>% summarise(median = median(estim), 
                                         lo = quantile(estim, 0.025), 
                                         hi = quantile(estim, 0.975)) %>%
                ungroup() %>% arrange(day)


  ## PLOTS
  ## fitted to data 
  pdfnam <- "Figures_covidcast_4.21.pdf"
  pdf(file=pdfnam,width=8, height=6) 


  ggplot2::ggplot() + 
    geom_point(data = input_data, aes(x = day, y = median, color = outcome)) + 
    geom_line(data = input_data, aes(x = day, y = median, color = outcome), linetype = 2) + 
    geom_line(data = filter(fit_to_data, day >=0), aes(x = day, y = median, color = outcome)) + 
    geom_ribbon(data = filter(fit_to_data, day >= 0), 
                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), alpha=0.3) +
    scale_color_manual(values = c('#a6cee3','#b2df8a','#fb9a99','#1f78b4','#33a02c','#e31a1c'), 
                       labels = c("Reported Cases", "Reported Deaths", "Reported Hospitalizations", 
                                  "Fitted Reported Cases", "Fitted Reported Deaths", "Fitted Reported Hospitalizations"), 
                       guide = guide_legend(override.aes = list(
                       linetype = c(rep("solid", 3), rep("dashed", 3))))) +
    labs(x = "Days Since March 1, 2020",
         y = "Count",
         title = "Observed and Fitted COVID-19 Cases, Hosptialization, Deaths in NYC: March 2 - April 15",
         caption = "Data accessed April 16, 2020",
         color = "")

  ## all cases to data 

  ggplot2::ggplot() + 
    geom_point(data = filter(input_data, outcome == "obs_cas"), 
               aes(x = day, y = median, color = outcome)) + 
    geom_line(data = filter(input_data, outcome == "obs_cas"), 
              aes(x = day, y = median, color = outcome), linetype = 2) + 
    geom_line(data = filter(new, outcome == "new_inf", day >=0), 
              aes(x = day, y = median, color = outcome)) + 
    geom_ribbon(data = filter(new, outcome == "new_inf", day >= 0), 
                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), alpha=0.3) +
    scale_color_manual(values = c('#a6cee3','#fc8d62'), 
                       labels = c("Modeled Cases", "Observed Cases"), 
                       guide = guide_legend(override.aes = list(
                         linetype = c("dashed","solid")))) +
    labs(x = "Days Since Start",
         y = "Count",
         title = "Observed and Fitted COVID-19 Cases, Hosptialization, Deaths: NYC March 2 - April 15",
         subtitle = "with 95% uncertainty intervals",
         caption = "Data accessed April 16, 2020",
         color = "")

  ## comparison: fitted, diagnosed, all cases
  ggplot2::ggplot() + 
    geom_line(data = filter(fit_to_data, outcome == "occur_cas"), 
              aes(x = day, y = median, color = outcome)) + 
        geom_ribbon(data = filter(fit_to_data, outcome == "occur_cas"), 
                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), 
                linetype = 2, alpha=0.2) +
    geom_line(data = filter(diag, outcome == "diag_all"), 
              aes(x = day, y = median, color = outcome)) +
       geom_ribbon(data = filter(diag, outcome == "diag_all"), 
                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), 
                linetype = 2, alpha=0.2) +
    geom_line(data = filter(new, outcome == "new_inf"), 
              aes(x = day, y = median, color = outcome)) + 
        geom_ribbon(data = filter(new, outcome == "new_inf"), 
                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), 
                linetype = 2, alpha=0.2) +
    scale_color_manual(values = c('#fc8d62','#66c2a5','#8da0cb'), 
                       breaks = c("new_inf", "diag_all", "occur_cas"),
                       labels = c("new_inf" = "Modeled All Cases",
                                  "diag_all" = "Modeled All Diagnosed", 
                                  "occur_cas" = "Fitted Reported Cases"))+
    labs(x = "Days Since March 1, 2020",
         y = "Count",
         title = "Modele dTotal, Diagnosed, and Reported COVID-19 Cases: NYC March 2 - April 15",
         subtitle = "with 95% uncertainty intervals",
         caption = "Data accessed April 16, 2020",
         color = "")

  # plotting: all, diagnosed, reported outcomes
  ggplot2::ggplot() +
    geom_line(data = fit_to_data, 
              aes(x = day, y = median, color = outcome)) + 
    geom_line(data = filter(diag, outcome != "diag_inf"), 
              aes(x = day, y = median, color = outcome)) +
    geom_line(data = new, 
              aes(x = day, y = median, color = outcome)) + 
    scale_color_manual(values = c('#08519c','#3182bd', '#6baed6','#bdd7e7',
                                  '#238b45', '#74c476','#bae4b3',
                                  '#6a51a3','#9e9ac8','#cbc9e2'), 
                       breaks = c("new_inf", "new_sym", "new_hos", "new_die", 
                                  "diag_all", "diag_sym", "diag_hos", 
                                  "occur_cas", "occur_hos", "occur_die"),
                       labels = c("Modeled All Cases", "Modeled Symptomatic", 
                                  "Modeled Hospitalizable", "Modeled Deaths", 
                                  "Modeled All Diagnosed", "Modeled Diagnosed at Symptomatic", 
                                  "Modeled Diagnosed at Hospitalizable", 
                                  "Fitted Cases", "Fitted Hospitalizations", 
                                  "Fitted Deaths")) +
    labs(x = "Days Since March 1, 2020",
         y = "Count",
         title = "Incidence Outcomes Compared to Data: NYC March 2 - April 15",
         subtitle = "with 95% uncertainty intervals",
         caption = "Data accessed April 16, 2020",
         color = "")


  new.par <- par(mfrow=c(2,3))
  mod <- fit[["p_sym_if_inf"]]
  mn <- as.numeric(datList[["pri_p_sym_if_inf_a"]])
  ss <- as.numeric(datList[["pri_p_sym_if_inf_b"]])
  hist(mod, main = "P(Progression to Symptomatic)", xlab = "", ylab = "")
  zz <- range(c(mod,0,1))
  x_pts <- seq(zz[1],zz[2],length.out=101)
  y_pts <- dbeta(x_pts,mn,ss)
  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
  
  mod <- fit[["p_hos_if_sym"]]
  mn <- as.numeric(datList[["pri_p_hos_if_sym_a"]])
  ss <- as.numeric(datList[["pri_p_hos_if_sym_b"]])
  hist(mod, main = "P(Progression to Hospitalizable)", xlab = "", ylab = "")
  zz <- range(c(mod,0,1))
  x_pts <- seq(zz[1],zz[2],length.out=101)
  y_pts <- dbeta(x_pts,mn,ss)
  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
  
  mod <- fit[["p_die_if_hos"]]
  mn <- as.numeric(datList[["pri_p_die_if_hos_a"]])
  ss <- as.numeric(datList[["pri_p_die_if_hos_b"]])
  hist(mod, main = "P(Progression to Death)", xlab = "", ylab = "")
  zz <- range(c(mod,0,1))
  x_pts <- seq(zz[1],zz[2],length.out=101)
  y_pts <- dbeta(x_pts,mn,ss)
  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
  
  mod <- fit[["p_diag_if_inf"]]
  mn <- as.numeric(datList[["pri_p_diag_if_inf_a"]])
  ss <- as.numeric(datList[["pri_p_diag_if_inf_b"]])
  hist(mod, main = "P(Diagnosis at Infectious)", xlab = "", ylab = "")
  zz <- range(c(mod,0,1))
  x_pts <- seq(zz[1],zz[2],length.out=101)
  y_pts <- dbeta(x_pts,mn,ss)
  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
  
  mod <- fit[["p_diag_if_sym"]]
  mn <- as.numeric(datList[["pri_p_diag_if_sym_a"]])
  ss <- as.numeric(datList[["pri_p_diag_if_sym_b"]])
  hist(mod, main = "P(Diagnosis at Symptomatic)", xlab = "", ylab = "")
  zz <- range(c(mod,0,1))
  x_pts <- seq(zz[1],zz[2],length.out=101)
  y_pts <- dbeta(x_pts,mn,ss)
  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
  
  mod <- fit[["p_diag_if_hos"]]
  mn <- as.numeric(datList[["pri_p_diag_if_hos_a"]])
  ss <- as.numeric(datList[["pri_p_diag_if_hos_b"]])
  hist(mod, main = "P(Diagnosis at Hospitalizable)", xlab = "", ylab = "")
  zz <- range(c(mod,0,1))
  x_pts <- seq(zz[1],zz[2],length.out=101)
  y_pts <- dbeta(x_pts,mn,ss)
  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)


  library("EpiEstim")
  library("lubridate")

  data1 <- group_by(new_inf, day) %>% summarise(cases = median(estim)) %>% 
    mutate(day = substr(day, start = 2, stop =4)) %>%
    arrange(as.numeric(day)) %>%
    mutate(I = round(cases), 
           dates = seq(as.Date("2020/02/21"), as.Date("2020/04/15"), by =1))

  # assuming a mean SI of 4.1 days with a sd of 2.1 days
  lognorm <- c(rlnorm(n = 100000, meanlog = 1.41, sdlog = 0.75))
  lognorm2 <- table(round(lognorm))
  lognorm_dist <- c(0, prop.table(lognorm2))

  config <- make_config(list(t_start = seq(8, nrow(data1)-4), 
                             t_end = seq(11, nrow(data1)-1),
                             si_distr = lognorm_dist))

  R_est  <- estimate_R(data1,
                       method = "non_parametric_si",
                       config = config)

  print(estimate_R_plots(R_est, what = "R") + 
    labs(title = "Preliminary R Effective Estimate, NYC", 
         subtitle = "Assuming mean SI of 4.1 days (sd 2.1 days)") )

  dev.off()
  ############


  ## comparison: diagnosed to fitted 
  ggplot2::ggplot() + 
    geom_line(data = diag, aes(x = day, y = median, color = outcome)) + 
    geom_line(data = filter(fit_to_data, outcome != "occur_die"), 
              aes(x = day, y = median, color = outcome), linetype = 2) + 
    scale_color_manual(values = c('#a6cee3','#fb9a99', '#cab2d6', '#fdbf6f', '#1f78b4','#e31a1c'), 
                       labels = c("Modeled All Diagnosed", "Modeled Diagnosed at Hospitalization", 
                                  "Modeled Diagnoed at Infectious", "Modeled Diagnosed at Symptomatic",
                                  "Fitted Cases", "Fitted Hospitalizations"), 
                       guide = guide_legend(override.aes = list(
                         linetype = c(rep("solid", 4), rep("dashed", 2))))) +
    labs(x = "Days Since March 2, 2020",
         y = "Count",
         title = "Modeled (Diagnosed) and Fitted COVID-19 Cases in NYC, March 2 - April 15",
         color = "")

  ## comparison: modeled (all) to fitted input 
  ggplot2::ggplot() + 
    geom_line(data = filter(new, outcome != "new_sym"), aes(x = day, y = median, color = outcome)) + 
    geom_line(data = fit_to_data, aes(x = day, y = median, color = outcome), linetype = 2) + 
    scale_color_manual(values = c('#a6cee3','#b2df8a','#fb9a99','#1f78b4','#33a02c','#e31a1c'), 
                       labels = c("Modeled All Deaths", "Modeled All Hospitalizations", "Modeled All Infections", 
                                  "Fitted Cases", "Fitted Deaths", "Fitted Hospitalizations"), 
                       guide = guide_legend(override.aes = list(
                         linetype = c(rep("solid", 3), rep("dashed", 3))))) +
    labs(x = "Days Since March 2, 2020",
         y = "Count",
         title = "Modeled (All) and Fitted COVID-19 Cases in NYC, March 2 - April 15",
         color = "")

  ###
  ## SAME but with uncertainty 

  ## comparison: diagnosed to fitted 
  ggplot2::ggplot() + 
    geom_line(data = diag, aes(x = day, y = median, color = outcome)) + 
      geom_ribbon(data = diag,
                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), 
                linetype = 2, alpha=0.3) +
    geom_line(data = filter(fit_to_data, outcome != "occur_die"), 
              aes(x = day, y = median, color = outcome), linetype = 2) + 
      geom_ribbon(data = filter(fit_to_data, outcome != "occur_die"), 
                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), linetype = 2, alpha=0.3) +
    scale_color_manual(values = c('#a6cee3','#fb9a99', '#cab2d6', '#fdbf6f', '#1f78b4','#e31a1c'), 
                       labels = c("Modeled All Diagnosed", "Modeled Diagnosed at Hospitalization", 
                                  "Modeled Diagnoed at Infectious", "Modeled Diagnosed at Symptomatic",
                                  "Fitted Cases", "Fitted Hospitalizations"), 
                       guide = guide_legend(override.aes = list(
                         linetype = c(rep("solid", 4), rep("dashed", 2))))) +
    labs(x = "Days Since March 2, 2020",
         y = "Count",
         title = "Modeled (Diagnosed) and Fitted COVID-19 Cases in NYC, March 2 - April 15",
         color = "")

  ## comparison: modeled (all) to fitted input 
  ggplot2::ggplot() + 
    geom_line(data = filter(new, outcome != "new_sym"), aes(x = day, y = median, color = outcome)) + 
        geom_ribbon(data = filter(new, outcome != "new_sym"),
                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), linetype = 2, alpha=0.3) +
    geom_line(data = fit_to_data, aes(x = day, y = median, color = outcome)) + 
        geom_ribbon(data = fit_to_data, 
                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), linetype = 2, alpha=0.3) +
    scale_color_manual(values = c('#a6cee3','#b2df8a','#fb9a99','#1f78b4','#33a02c','#e31a1c'), 
                       labels = c("Modeled All Deaths", "Modeled All Hospitalizations", "Modeled All Infections", 
                                  "Fitted Cases", "Fitted Deaths", "Fitted Hospitalizations"), 
                       guide = guide_legend(override.aes = list(
                         linetype = c(rep("solid", 6))))) +
    labs(x = "Days Since March 2, 2020",
         y = "Count",
         title = "Modeled (All) and Fitted COVID-19 Cases in NYC, March 2 - April 15",
         color = "")
}




