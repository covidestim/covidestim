#' @export
#' @rdname Rt.estimate
RtEst <- function(...) UseMethod('RtEst')

#' Calculate the Effective Reprouction Number from CovidCast output
#'
#' Creates a figure with Rt estimates over time and 95% CI and undelrying df. 
#'
#' @param cc The result of calling \code{\link{run}}. An object of class
#'   \code{covidcast_result}.
#'   
#' @param sample_fraction The fraction of iterations to sample for the Rt
#' calualation; defaults to 2/3 samples after warm-up.  
#' 
#' @param window The size of the moving window over which Rt should be 
#' estimated; defaults to 7 days. 
#' 
#' @param mean.si The mean serial interval to use in the Rt estimate; defaults
#' to 4.7 days. 
#' 
#' @param std.si The standard deviation of the serial interval; defaults to 
#' 2.9 days.  
#'
#' @param renderPDF A logical scalar. If true, will render a pdf to the working
#'   directory with name \code{covidcast_Rt.pdf}.
#'
#' @return an df with estimates
#'
#' @export
Rt.estimate <- function(cc, 
                        sample_fraction = (2/3), 
                        window = 7, 
                        mean.si = 4.7,
                        std.si = 2.9,
                        pdf.name = covidcast_Rt.pdf)
  
fit         = cc$extracted 
tot_iter    = cc$iter
warm        = cc$warmup
iter        = (tot_iter - warm) * cc$chains
n_sample    = sample_fraction * iter
day_start   = 2
day_end     = window + day_start - 1
ndb         = cc$config$N_days_before

sample_iter <- sample(iter, n_sample)

inf <- as.data.frame(fit[["new_inf"]]) %>% 
             gather(key = day, value = I) %>%
             mutate(day = as.numeric(substr(day, start = 2, stop = 4))) %>%
             filter(day > ndb) %>%
             group_by(day) %>% 
             mutate(iter = seq(1,iter)) %>% 
             ungroup() 

# number of moving windows in which Rt will be estimated     
# total day - window works if 1st day is day 2 of data
windows <- length(unique(inf$day)) - window 

# create matrices to store output of epi estim call 
# where each row is value for a window, each column is estimates from an iter
mn <- matrix(nrow = windows, ncol = n_sample)
lo <- matrix(nrow = windows, ncol = n_sample)
hi <- matrix(nrow = windows, ncol = n_sample)


#' @import EpiEstim
# now we generate estimates of the mean, upper, and lower bounds of Rt from 
# each sampled iteration of case data. 
for(i in seq_along(sample_iter)) {
  inf2 <- filter(inf, iter == sample_iter[i]) 
  
  config <- EpiEstim::make_config(list(
                             t_start = seq(day_start, 
                                           nrow(inf2)-(day_end-day_start)), 
                             t_end = seq(day_end, 
                                         nrow(inf2)),
                             mean_si = mean.si, 
                             std_si = std.si))
  out <- EpiEstim::estimate_R(inf2$I,
                    method = "parametric_si",
                    config = config)
  R_est <- out[["R"]]
  
  mn[,i] <- R_est$`Mean(R)`
  lo[,i] <- R_est$`Quantile.0.025(R)`
  hi[,i] <- R_est$`Quantile.0.975(R)`
}

# convert matriaces into df, calculate bounds on each estimate
make_est_df <- function(x){
    as.data.frame(x) %>% 
    rowid_to_column() %>% 
    gather(key = "iter", 
           value = "est", 
           2:points) %>%
    group_by(rowid) %>%
    summarise(mn = mean(est), 
              lo = quantile(est, 0.025), 
              hi = quantile(est, 0.975))
  
}

Rt    <- make_est_df(mn)
lower <- make_est_df(lo) 
upper <- make_est_df(hi) 

# combine: 
Rt_df <- as.data.frame(cbind(Rt$mn, 
                             lower$lo, 
                             upper$hi)) %>% 
          mutate(day = seq(1,nrow(Rt),1)) 

#' @import ggplot2
plota <- ggplot2::ggplot() +
              geom_hline(aes(yintercept = 1, color = "red"), 
                         size = 0.5, show.legend = FALSE) +
              geom_line(data = Rt_df, aes(x = day, y = V1)) + 
              geom_ribbon(data = Rt_df, 
                          aes(x = day, y = V1, ymin=V2, ymax=V3), alpha=0.3) +
              scale_x_continuous(breaks = seq(1,n_days,10)) +
              labs(y = "Rt", 
                   x = "Days Since Start", 
                   title = "Effective Reproduction Number Estimate") 

pdf(file = pdfname, width = 8, height = 6) 
print(plota)
dev.off()

return(Rt_df)

