### Correlating population estimates with AFM case data - Northern America

pacman::p_load(
  ape
  , reshape2
  , plyr
  , dplyr
  , lubridate
  , janitor
  , flextable
  , forecast
  , tibble
  , cowplot
  , dispRity
  , ggnewscale
  , tidyr
  , tibble
  , ggplot2
  , coda
  , bdskytools
  , beastio
  , pammtools
  , forcats
  , ggpubr
  , readr
  , magrittr
)

# AFM case data
USA_cases <- read.csv("./data/USA_AFM.csv") # AFM case data from cdc.gov



## Trajectories 
traj_B <- read_tsv("./samples/samp_3mo1/D68_B.traj")

traj_B %<>% 
  filter(Sample>=40000000) %>% # Burnin (states)
  filter(variable=="N") # Only keep popsize estimates

ages <- seq(0,40,length.out=1001)

# pull HPDs
gridded_traj_B <- traj_B %>%
  group_by(Sample, type) %>%
  reframe(N=approx(age, value, ages, method="constant", f=0, yright=0)$y,
          age=ages) %>%
  group_by(type,age) %>%
  summarize(low=quantile(N,0.025),
            high=quantile(N,0.975),
            med=quantile(N,0.5))

# plot the median and HPD
popestimate_B <- ggplot(gridded_traj_B %>%
                          mutate(Date=ymd("2022-12-10")-age*365) %>%
                          mutate(Location=recode(factor(type),
                                                 "0"="North America",
                                                 "1"="Europe")),
                        aes(Date, N, col=Location, fill=Location, y=med, ymin=low, ymax=high)) +
  geom_ribbon(alpha=0.5) +
  geom_line() + ylab("Population size") +
  scale_x_date(date_breaks="4 years", limits = c(as.Date("1996-01-01"), NA), date_minor_breaks= "1 year",  date_labels="%Y") +
  scale_fill_manual(values = c("mediumseagreen", "gold2"), na.value = NA, na.translate = F, name = "Region") +
  scale_color_manual(values = c("mediumseagreen", "gold2"), na.value = NA, na.translate = F, name = "Region") +
  ggtitle("Clade B + C") +
  theme_bw()


traj_A <- read_tsv("./samples/samp_3mo1/D68_A.traj")

traj_A %<>% 
  filter(Sample>=40000000) %>% # Burnin (states)
  filter(variable=="N") # Only keep popsize estimates

ages <- seq(0,40,length.out=1001)

gridded_traj_A <- traj_A %>%
  group_by(Sample, type) %>%
  reframe(N=approx(age, value, ages, method="constant", f=0, yright=0)$y,
          age=ages) %>%
  group_by(type,age) %>%
  summarize(low=quantile(N,0.025),
            high=quantile(N,0.975),
            med=quantile(N,0.5))


popestimate_A <- ggplot(gridded_traj_A %>%
                          mutate(Date=ymd("2022-12-10")-age*365) %>%
                          mutate(Location=recode(factor(type),
                                                 "0"="Northern America",
                                                 "1"="Western Europe")),
                        aes(Date, N, col=Location, fill=Location, y=med, ymin=low, ymax=high)) +
  geom_ribbon(alpha=0.5) +
  geom_line() + ylab("Population size") +
  scale_x_date(date_breaks="4 years", limits = c(as.Date("1996-01-01"), NA), date_minor_breaks= "1 year",  date_labels="%Y") +
  scale_fill_manual(values = c("mediumseagreen", "gold2"), na.value = NA, na.translate = F, name = "Region") +
  scale_color_manual(values = c("mediumseagreen", "gold2"), na.value = NA, na.translate = F, name = "Region") +
  ggtitle("Clade A + D") +
  theme_bw()

# plot by superclade
estimated_traj <- ggarrange(popestimate_A, popestimate_B, ncol = 1, common.legend = T, legend = "bottom")
#ggsave("./plots/estimated_traj.png", estimated_traj, height = 4, width = 6.5, units = "in", dpi = 700)


# group by region instead ----

gridded_traj_NA_B <- traj_B %>%
  filter(type == "0") %>%
  group_by(Sample) %>%
  reframe(N=approx(age, value, ages, method="constant", f=0, yright=0)$y,
          age=ages) %>%
  group_by(age) %>%
  summarize(low=quantile(N,0.025),
            high=quantile(N,0.975),
            med=quantile(N,0.5))%>%
  mutate(clade = "B")


gridded_traj_NA_A <- traj_A %>%
  filter(type == "0") %>%
  group_by(Sample) %>%
  reframe(N=approx(age, value, ages, method="constant", f=0, yright=0)$y,
          age=ages) %>%
  group_by(age) %>%
  summarize(low=quantile(N,0.025),
            high=quantile(N,0.975),
            med=quantile(N,0.5)) %>%
  mutate(clade = "A")

gridded_traj_NA <- bind_rows(gridded_traj_NA_A, gridded_traj_NA_B)

gridded_traj_EU_B <- traj_B %>%
  filter(type == "1") %>%
  group_by(Sample) %>%
  reframe(N=approx(age, value, ages, method="constant", f=0, yright=0)$y,
          age=ages) %>%
  group_by(age) %>%
  summarize(low=quantile(N,0.025),
            high=quantile(N,0.975),
            med=quantile(N,0.5))%>%
  mutate(clade = "B")


gridded_traj_EU_A <- traj_A %>%
  filter(type == "1") %>%
  group_by(Sample) %>%
  reframe(N=approx(age, value, ages, method="constant", f=0, yright=0)$y,
          age=ages) %>%
  group_by(age) %>%
  summarize(low=quantile(N,0.025),
            high=quantile(N,0.975),
            med=quantile(N,0.5)) %>%
  mutate(clade = "A")

gridded_traj_EU <- bind_rows(gridded_traj_EU_A, gridded_traj_EU_B)


# Time series cross-correlation
gridded_traj_EU %<>%
  mutate(Date=ymd("2022-12-10")-age*365) 

gridded_traj_NA %<>%
  mutate(Date=ymd("2022-12-10")-age*365)

USA_cases_grp <- USA_cases %>%
  mutate(Date = as.Date(Month, format = "%m/%d/%Y"),
         group_date = floor_date(Date, "month")) %>%
  group_by(group_date, .add = T)%>%
  summarise(group_case= sum(Cases))

# convert to monthly to match AFM data
grid_NA_grp <- gridded_traj_NA %>%
  mutate(group_date = floor_date(Date, "month"))%>%
  group_by(group_date, clade)%>%
  summarise(med= sum(med))

corr_data_NA <- left_join(grid_NA_grp, USA_cases_grp, by = "group_date")

corLag_A_NA <-ccf(corr_data_NA[corr_data_NA$clade == "A",]$med,corr_data_NA[corr_data_NA$clade == "A",]$group_case,lag.max=6, na.action = na.pass)

corLag_B_NA <-ccf(corr_data_NA[corr_data_NA$clade == "B",]$med,corr_data_NA[corr_data_NA$clade == "B",]$group_case,lag.max=6, na.action = na.pass)


# convert to monthly to match AFM data
grid_EU_grp <- gridded_traj_EU %>%
  mutate(group_date = floor_date(Date, "month"))%>%
  group_by(group_date, clade)%>%
  summarise(med= sum(med))

corr_data_EU <- left_join(grid_EU_grp, USA_cases_grp, by = "group_date")

corLag_A_EU <-ccf(corr_data_EU[corr_data_EU$clade == "A",]$med,corr_data_EU[corr_data_EU$clade == "A",]$group_case,lag.max=6, na.action = na.pass)

corLag_B_EU <-ccf(corr_data_EU[corr_data_EU$clade == "B",]$med,corr_data_EU[corr_data_EU$clade == "B",]$group_case,lag.max=6, na.action = na.pass)


# only 122 months have case data so this is the true sample

A_EU_corr <- do.call(cbind.data.frame, corLag_A_EU) %>%
  select(acf, lag, n.used) %>%
  ggplot(mapping = aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2 / sqrt(122), color = "#f54b42", linetype = 2) +
  geom_hline(yintercept = -(2 / sqrt(122)), color = "#f54b42",  linetype = 2) +
  geom_col(width = .1) +
  scale_y_continuous(limits = c(-1,1))+
  labs(title = "A + D clade Europe"
       , y = "Cross-correlation"
       , x = "Lag (months)") +
  theme_bw()

B_EU_corr<- do.call(cbind.data.frame, corLag_B_EU) %>%
  select(acf, lag, n.used) %>%
  ggplot(mapping = aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2 / sqrt(122), color = "#4287f5", linetype = 2) +
  geom_hline(yintercept = -(2 / sqrt(122)), color = "#4287f5",  linetype = 2) +
  geom_col(width = .1) +
  scale_y_continuous(limits = c(-1,1))+
  labs(title = "B + C clade Europe"
       , y = "Cross-correlation"
       , x = "Lag (months)") +
  theme_bw()

A_NA_corr <- do.call(cbind.data.frame, corLag_A_NA) %>%
  select(acf, lag, n.used) %>%
  ggplot(mapping = aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2 / sqrt(122), color = "#f54b42", linetype = 2) +
  geom_hline(yintercept = -(2 / sqrt(122)), color = "#f54b42",  linetype = 2) +
  geom_col(width = .1) +
  scale_y_continuous(limits = c(-1,1))+
  labs(title = "A + D clade North America"
       , y = "Cross-correlation"
       , x = "Lag (months)") +
  theme_bw()

B_NA_corr<- do.call(cbind.data.frame, corLag_B_NA) %>%
  select(acf, lag, n.used) %>%
  ggplot(mapping = aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2 / sqrt(122), color = "#4287f5", linetype = 2) +
  geom_hline(yintercept = -(2 / sqrt(122)), color = "#4287f5",  linetype = 2) +
  geom_col(width = .1) +
  scale_y_continuous(limits = c(-1,1))+
  labs(title = "B + C clade North America"
       , y = "Cross-correlation"
       , x = "Lag (months)") +
  theme_bw()


traj_AFM_corr <- ggarrange(A_EU_corr, B_EU_corr, A_NA_corr, B_NA_corr, ncol = 2, nrow = 2)

title_traj_AFM_corr <- annotate_figure(traj_AFM_corr, top = text_grob("Correlation of reconstructed trajectories with reported AFM"))


#figure 5
us_correlation <- ggarrange(A_NA_corr, B_NA_corr, ncol = 2)
ggsave("./plots/fig5.png", us_correlation, height = 2, width = 6, units = "in", dpi = 700)


# Sliding window cross-correlation ----
# sliding window - B
xy <- cbind(corr_data_NA[corr_data_NA$clade=="B",]$med, corr_data_NA[corr_data_NA$clade=="B",]$group_case)
xy <- xy[complete.cases(xy),]

start_date <- as.Date("2014-08-01")
window_size <- 24
max_lag <- 6
step_size <- 6
lags <- -max_lag:max_lag
times <- seq(1, nrow(xy) - window_size + 1, by = step_size)
dates <- start_date + months(times - 1)  # convert index to dates

# Initialize empty list to collect results
results <- list()

# Compute sliding window cross-correlation
for (i in seq_along(times)) {
  t <- times[i]
  xy_win <- xy[t:(t + window_size - 1),]
  ccf_vals <- ccf(xy_win[,1], xy_win[,2], lag.max = max_lag, plot = FALSE)$acf[,1,1]
  results[[length(results) + 1]] <- data.frame(
    date = dates[i],
    time = t,
    lag = lags,
    correlation = ccf_vals
  )
}

# Combine and tidy
ccf_df <- bind_rows(results)

# Plot heatmap with ggplot2
b_slide <- ggplot(ccf_df, aes(x = date, y = lag, fill = correlation, width = 185)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "goldenrod2", midpoint = 0, limits = c(-1,1)) +
  scale_x_date(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))+
  labs(
    title = "Sliding Window Cross-Correlation",
    x = "Time (2-year window start)",
    y = "Lag (months)",
    fill = "cross-correlation"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# sliding window - A
xy <- cbind(corr_data_NA[corr_data_NA$clade=="A",]$med, corr_data_NA[corr_data_NA$clade=="A",]$group_case)
xy <- xy[complete.cases(xy),]

start_date <- as.Date("2014-08-01")
window_size <- 24
max_lag <- 6
step_size <- 6
lags <- -max_lag:max_lag
times <- seq(1, nrow(xy) - window_size + 1, by = step_size)
dates <- start_date + months(times - 1)  # convert index to dates

# Initialize empty list to collect results
results <- list()

# Compute sliding window cross-correlation
for (i in seq_along(times)) {
  t <- times[i]
  xy_win <- xy[t:(t + window_size - 1),]
  ccf_vals <- ccf(xy_win[,1], xy_win[,2], lag.max = max_lag, plot = FALSE)$acf[,1,1]
  results[[length(results) + 1]] <- data.frame(
    date = dates[i],
    time = t,
    lag = lags,
    correlation = ccf_vals
  )
}

# Combine and tidy
ccf_df <- bind_rows(results)

# Plot heatmap with ggplot2
a_slide <- ggplot(ccf_df, aes(x = date, y = lag, fill = correlation, width = 185)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "goldenrod2", midpoint = 0, , limits = c(-1,1)) +
  scale_x_date(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(
    title = "Sliding Window Cross-Correlation",
    x = "Time (2-year window start)",
    y = "Lag (months)",
    fill = "Cross-correlation"
  ) +
  theme_bw()+
  theme(legend.position = "bottom") 

ggsave("./plots/traj_AFM_correlation.png", traj_AFM_corr, height = 4, width = 6.5, units = "in", dpi = 700)
us_slide_correlation <- ggarrange(A_NA_corr, B_NA_corr, a_slide, b_slide, align = "v", common.legend = T, legend = "bottom", labels = "AUTO")
ggsave("./plots/fig5.png", us_slide_correlation, height = 5.3, width = 8, units = "in")


# plot trajectories against AFM

popestimate_NA <- ggplot(gridded_traj_NA %>%
                           mutate(clade=recode(factor(clade),
                                               "A"="A+D",
                                               "B"="B+C")),
                         aes(Date, N, col=clade, fill=clade, y=med, ymin=low, ymax=high)) +
  geom_rect(aes(xmin = as.Date("2014-01-01"), xmax = as.Date("2015-01-01"), ymin = -Inf, ymax = Inf), fill = "grey85", color = NA)+
  geom_rect(aes(xmin = as.Date("2016-01-01"), xmax = as.Date("2017-01-01"), ymin = -Inf, ymax = Inf), fill = "grey85", color = NA)+
  geom_rect(aes(xmin = as.Date("2018-01-01"), xmax = as.Date("2019-01-01"), ymin = -Inf, ymax = Inf), fill = "grey85", color = NA)+
  geom_ribbon(alpha=0.5, linewidth = 0.2) +
  geom_line() + ylab("Population size") +
  scale_x_date(expand = expansion(), date_breaks="4 years", limits = c(as.Date("1994-01-01"), as.Date("2024-01-01")), date_minor_breaks= "1 year",  date_labels="%Y") +
  scale_fill_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  scale_color_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  coord_cartesian(xlim = c(as.Date("1996-01-01"), as.Date("2025-01-01"))) +
  ggtitle("North America") +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), plot.margin = unit(c(0.01,0.1,0.1,0.1), "cm"), panel.grid.minor.y = element_blank(), axis.text = element_text(size = 12))

popestimate_EU <- ggplot(gridded_traj_EU %>%
                           mutate(clade=recode(factor(clade),
                                               "A"="A+D",
                                               "B"="B+C")),
                         aes(Date, N, col=clade, fill=clade, y=med, ymin=low, ymax=high)) +
  geom_rect(aes(xmin = as.Date("2014-01-01"), xmax = as.Date("2015-01-01"), ymin = -Inf, ymax = Inf), fill = "grey85", color = NA)+
  geom_rect(aes(xmin = as.Date("2016-01-01"), xmax = as.Date("2017-01-01"), ymin = -Inf, ymax = Inf), fill = "grey85", color = NA)+
  geom_rect(aes(xmin = as.Date("2018-01-01"), xmax = as.Date("2019-01-01"), ymin = -Inf, ymax = Inf), fill = "grey85", color = NA)+
  geom_ribbon(alpha=0.5, linewidth = 0.2) +
  geom_line() + ylab("Population size") +
  scale_x_date(expand = expansion(), date_breaks="4 years", limits = c(as.Date("1994-01-01"), as.Date("2024-01-01")), date_minor_breaks= "1 year",  date_labels="%Y") +
  scale_fill_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  scale_color_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  coord_cartesian(xlim = c(as.Date("1996-01-01"), as.Date("2025-01-01"))) +
  ggtitle("Europe") +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), plot.margin = unit(c(0.01,0.1,0.1,0.1), "cm"), panel.grid.minor.y = element_blank(), axis.text = element_text(size = 12))

estimated_traj_loc <- ggarrange(popestimate_NA, popestimate_EU, ncol = 1, common.legend = T, legend = "bottom")
#ggsave("./plots/estimated_traj_loc_3mo.png", estimated_traj_loc, height = 4, width = 6.5, units = "in", dpi = 700)
estimated_traj_loc

AFM_USA <- ggplot()+
  geom_rect(aes(xmin = 2014, xmax = 2015, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_rect(aes(xmin = 2016, xmax = 2017, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_rect(aes(xmin = 2018, xmax = 2019, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_col(data = USA_cases,aes(x = decimal_date(as.Date(Month, format = "%m/%d/%Y")), y = Cases), fill = "navy")+
  scale_x_continuous(expand = expansion(), breaks = seq(1996, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
  scale_y_continuous(expand = expansion(add = c(0,5)),  breaks = seq(0,90, by = 30))+
  coord_cartesian(xlim = c(1996, 2025)) +
  theme_bw() +
  labs(title = "AFM cases - USA (cdc.gov) - Older data on inset graphs"
       , y = "Cases"
       , x = "Date") +
  theme(axis.title.x = element_blank(), plot.margin = unit(c(0.01,0.1,0.1,0.1), "cm"), axis.title.y = element_text(margin = margin(0,22,0,0)), axis.text = element_text(size = 12))

#Inset some earlier data
## Ayscue et al. 2014 - Calif.
CA_2014 <- data.frame(Cases = c(1,0,3,3,1,2,0,0,0,1,0,1,0,1,0,0,0,0,2,0,2,2,2,2,0), Month = seq(as.Date("2012-06-01"), as.Date("2014-06-01"), by = "month"))

## Cortese et al. 2020 Midwest/Mountain West USA
MW_2020 <- data.frame(Cases = c(1,1,1,0,0,1,2,0,2,1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,2,0,0,2,0,0,1,0,1,0,0,0,0,0,0,0,2,0,0,1,1,1,1,0,0,0,0,1,0,2,1,0,0,13,2),
                      Month = seq(as.Date("2005-01-01"), as.Date("2014-11-01"), by = "2 months"))

MW_plot <- ggplot() +
  geom_rect(aes(xmin = 2014, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_col(data = MW_2020 ,aes(x = decimal_date(as.Date(Month, format = "%m/%d/%Y")), y = Cases), fill = "navy")+
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
  scale_y_continuous(expand = expansion(add = c(0,1)),  breaks = seq(0,15, by = 5))+
  theme_bw() +
  theme(axis.title = element_blank(), plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"), axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.minor.y = element_blank(), axis.text = element_text(size = 11))

CA_plot <- ggplot() +
  geom_rect(aes(xmin = 2014, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_col(data = CA_2014 ,aes(x = decimal_date(as.Date(Month, format = "%m/%d/%Y")), y = Cases), fill = "navy")+
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
  scale_y_continuous(expand = expansion(add = c(0,2)),  breaks = seq(0,15, by = 5))+
  theme_bw() +
  theme(axis.title = element_blank(), plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"), axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.minor.y = element_blank(), axis.text = element_text(size = 11))


tst_inset <- ggdraw() +
  draw_plot(AFM_USA) +
  draw_plot(MW_plot, x = .3405, y = .55, height = .27, width = 0.348) +
  draw_plot(CA_plot, x = .581, y = .415, height = .135, width = 0.0935) +
  draw_label("Cortese et al. 2020", x = 0.24, y = 0.72, size = 10 ) +
  draw_label("Midwest/Mountain west", x = 0.24, y = 0.63, size = 7 ) +
  
  draw_label("Ayscue et al. 2014", x = 0.48 , y = 0.525, size = 10 ) +
  draw_label("California", x = 0.48, y = 0.445, size = 7 ) 

ggsave("./plots/AFM_USA.png", tst_inset, height = 1.6, width = 6.5, units = "in", dpi = 700)

case_traj <- plot_grid(tst_inset,estimated_traj_loc, ncol = 1, rel_heights = c(1,2.1))
case_traj

ggsave("./plots/fig4.svg", case_traj, height = 5.5, width = 6.5, units = "in", dpi = 700)


#Europe cases AFM - helfferich et al. 2025 (eurosurv) (annual)
EUR_AFM <- data.frame(Cases = c(37,5,39,13,3,6,15,12), group_date = seq(as.Date("2016-01-01"), as.Date("2023-01-01"), by = "year"))
ggplot(EUR_AFM)+
  geom_col(aes(y=Cases, x=group_date), fill = "gold2", width = 355, position = position_nudge(x = 184))+
  scale_x_date(expand = c(0,0), date_breaks = "1 year", date_labels = "%Y")+
  scale_y_continuous(expand = c(0,0), limits = c(0,40))+
  labs(title = "AFM cases - Europe (Helfferich et al. 2025 Eurosurv.)"
       , y = "AFM Cases"
       , x = "Date") +
  theme_bw()

ggsave("./plots/AFM_europe_suppl.png", dpi = 700, width = 6, height = 3, units = "in")

