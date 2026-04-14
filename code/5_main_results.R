# Main manuscript results plots

pacman::p_load(
  ape
  , reshape2
  , plyr
  , dplyr
  , lubridate
  , janitor
  , flextable
  , cowplot
  , dispRity
  , ggnewscale
  , gginnards
  , tidyr
  , readr
  , tibble
  , ggplot2 
  , coda
  , bdskytools
  , beastio
  , pammtools
  , forcats
  , ggpubr
  , magrittr
  , readr
  ,fmcmc

  )


## universal options ##

file_root <- "/Volumes/dnj13/home/resub_bdmm/samp_3mo2/" # where to find your beast outputs
burnin = 0.5 # burnin for time plots (proportion)
save_plots = F # should outputs be saved?
# folder to save plots

## Figure 3 : Rt trajectories --------------------------------------------------
# metadata
load(paste0(file_root, "A_sub_3mo_2.rdata"))
A_3mo1 <- A_meta_strat_date_subs
load(paste0(file_root, "B_sub_3mo_2.rdata"))
B_3mo1 <- B_meta_strat_date_subs

# logfile
tr_3m01 <- beastio::readLog(paste0(file_root, "log_run1.log"), burnin=burnin)
tr_3m02 <- beastio::readLog(paste0(file_root, "log_run2.log"), burnin=burnin)
tr_3m03 <- beastio::readLog(paste0(file_root, "log_run3.log"), burnin=burnin)
tr_3m04 <- beastio::readLog(paste0(file_root, "log_run4.log"), burnin=burnin)
tr_3m05 <- beastio::readLog(paste0(file_root, "log_run5.log"), burnin=burnin)
tr_3m06 <- beastio::readLog(paste0(file_root, "log_run6.log"), burnin=burnin)
tr_3m07 <- beastio::readLog(paste0(file_root, "log_run7.log"), burnin=burnin)
tr_3m08 <- beastio::readLog(paste0(file_root, "log_run8.log"), burnin=burnin)

tr_3mo1 <- append_chains(tr_3m01, tr_3m02, tr_3m03, tr_3m04, tr_3m05, tr_3m06, tr_3m07, tr_3m08)


A_3mo1$clade <- factor(A_3mo1$clade, levels = c("D","A","Other"))
B_3mo1$clade <- factor(B_3mo1$clade,  levels = c("B","C","Other"))


## split A/B and NA/EU
Rt_A_NA_3m <- beastio::getLogFileSubset(tr_3mo1, "ReSPEpi.Ai\\d+_North.America")
Rt_A_EU_3m <- beastio::getLogFileSubset(tr_3mo1, "ReSPEpi.Ai\\d+_Europe")

times_A_3m <- beastio::getLogFileSubset(tr_3mo1, "ReSPEpi.Ai\\d+_endtime")
origin_A_3m <- beastio::getLogFileSubset(tr_3mo1, "originBDMMPrime.A")

# convert from relative to origin to relative to MRSD 
# seems origin is given backwards in time from the MRSD and changepoints forward in time from the origin?
reltimes_A_3m <- -origin_A_3m + times_A_3m

# relative to MRSD - match the inputs!
times_A_3m <- 2021.981 + reltimes_A_3m

Rt_B_NA_3m <- beastio::getLogFileSubset(tr_3mo1, "ReSPEpi.Bi\\d+_North.America")
Rt_B_EU_3m <- beastio::getLogFileSubset(tr_3mo1, "ReSPEpi.Bi\\d+_Europe")

times_B_3m <- beastio::getLogFileSubset(tr_3mo1, "ReSPEpi.Bi\\d+_endtime")
origin_B_3m <- beastio::getLogFileSubset(tr_3mo1, "originBDMMPrime.B")

reltimes_B_3m <- -origin_B_3m + times_B_3m
times_B_3m <- 2022.94 + reltimes_B_3m

# HPDs
Rt_A_NA_hpd_3m    <- t(beastio::getHPDMedian(Rt_A_NA_3m))
Rt_A_EU_hpd_3m    <- t(beastio::getHPDMedian(Rt_A_EU_3m))

Rt_B_NA_hpd_3m    <- t(beastio::getHPDMedian(Rt_B_NA_3m))
Rt_B_EU_hpd_3m    <- t(beastio::getHPDMedian(Rt_B_EU_3m))


# plot A North America ----
A_NA_std <- Rt_A_NA_hpd_3m %>%
  t() %>%
  as.data.frame() %>%
  add_column(grid_end = c(times_A_3m[1,], 2024)) %>%
  filter(grid_end > 2021.981-mean(origin_A_3m), grid_end < 2023.2)

A_NA_plot <- ggplot() +
  #  geom_rect(aes(xmin = 2014, xmax = 2015, ymin = -Inf, ymax = Inf), fill = "grey85")+
  #  geom_rect(aes(xmin = 2016, xmax = 2017, ymin = -Inf, ymax = Inf), fill = "grey85")+
  #  geom_rect(aes(xmin = 2018, xmax = 2019, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_stepribbon(data = A_NA_std, aes(ymin = lower, ymax = upper, x = grid_end), alpha = 0.5,  direction = "vh") +
  geom_hline(yintercept =1, color = "#f54b42") +
  geom_step(data = A_NA_std, aes(x = grid_end, y = med), direction = "vh") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994.2, 2024.2, by = 4), minor_breaks = seq(1994.2, 2024.2, by = 1), labels = c("March\n1994", "March\n1998", "March\n2002", "March\n2006", "March\n2010", "March\n2014", "March\n2018", "March\n2022") )+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(0,2.1)) +
  theme_bw()+
  labs(title = "A + D clade North America"
       , y = "Rt"
       , x = "Date")+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0.1,0.1,0,0.1), "cm"), axis.text = element_text(size = 12))

A_NA_samps <- ggplot()+
  # geom_rect(aes(xmin = 2014, xmax = 2015, ymin = -Inf, ymax = Inf), fill = "grey85")+
  # geom_rect(aes(xmin = 2016, xmax = 2017, ymin = -Inf, ymax = Inf), fill = "grey85")+
  # geom_rect(aes(xmin = 2018, xmax = 2019, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_point(data = A_3mo1[A_3mo1$region == "North America",], aes(x = date, y = fct_rev(clade)), alpha = 0.3, color = "#f54b42") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994.2, 2024.2, by = 4), minor_breaks = seq(1994.2, 2024.2, by = 1), labels = c("March\n1994", "March\n1998", "March\n2002", "March\n2006", "March\n2010", "March\n2014", "March\n2018", "March\n2022") )+
  scale_y_discrete(drop = F)+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(1,3)) +
  theme_bw() +
  theme(axis.title = element_blank(), plot.margin = unit(c(0,0.1,0.1,0.1), "cm"), axis.text = element_text(size = 12))

A_NA_plot <- plot_grid(A_NA_plot, A_NA_samps, ncol = 1, rel_heights = c(6.5,3.5), align = "v")

# plot A Europe ----
A_EU_std <- Rt_A_EU_hpd_3m %>%
  t() %>%
  as.data.frame() %>%
  add_column(grid_end = c(times_A_3m[1,], 2024)) %>%
  filter(grid_end > 2021.981-mean(origin_A_3m), grid_end < 2023.2)

A_EU_plot <- ggplot() +
  geom_stepribbon(data = A_EU_std, aes(ymin = lower, ymax = upper, x = grid_end), alpha = 0.5,  direction = "vh") +
  geom_hline(yintercept =1, color = "#f54b42") +
  geom_step(data = A_EU_std, aes(x = grid_end, y = med), direction = "vh") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994.2, 2024.2, by = 4), minor_breaks = seq(1994.2, 2024.2, by = 1), labels = c("March\n1994", "March\n1998", "March\n2002", "March\n2006", "March\n2010", "March\n2014", "March\n2018", "March\n2022") )+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(0,2.1)) +
  theme_bw()+
  labs(title = "A + D clade Europe"
       , y = "Rt"
       , x = "Date")+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0.1,0.1,0,0.1), "cm"), axis.text = element_text(size = 12))

A_EU_samps <- ggplot()+
  geom_point(data = A_3mo1[A_3mo1$region == "Europe",], aes(x = date, y = fct_rev(clade)), alpha = 0.3, color = "#f54b42") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994.2, 2024.2, by = 4), minor_breaks = seq(1994.2, 2024.2, by = 1), labels = c("March\n1994", "March\n1998", "March\n2002", "March\n2006", "March\n2010", "March\n2014", "March\n2018", "March\n2022") )+
  scale_y_discrete(drop = F)+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(1,3)) +
  theme_bw() +
  theme(axis.title = element_blank(), plot.margin = unit(c(0,0.1,0.1,0.1), "cm"), axis.text = element_text(size = 12))

A_EU_plot <- plot_grid(A_EU_plot, A_EU_samps, ncol = 1, rel_heights = c(6.5,3.5), align = "v")

# plot B North America ----

B_NA_std <- Rt_B_NA_hpd_3m %>%
  t() %>%
  as.data.frame() %>%
  add_column(grid_end = c(times_B_3m[1,], 2024)) %>%
  filter(grid_end > 2022.94-mean(origin_B_3m), grid_end < 2023.2)

B_NA_plot <- ggplot() +
  geom_stepribbon(data = B_NA_std, aes(ymin = lower, ymax = upper, x = grid_end), alpha = 0.5,  direction = "vh") +
  geom_hline(yintercept =1, color = "#4287f5") +
  geom_step(data = B_NA_std, aes(x = grid_end, y = med), direction = "vh") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994.2, 2024.2, by = 4), minor_breaks = seq(1994.2, 2024.2, by = 1), labels = c("March\n1994", "March\n1998", "March\n2002", "March\n2006", "March\n2010", "March\n2014", "March\n2018", "March\n2022") )+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(0,2.1)) +
  theme_bw()+
  labs(title = "B + C clade North America"
       , y = "Rt"
       , x = "Date")+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0.1,0.1,0,0.1), "cm"), axis.text = element_text(size = 12))

B_NA_samps <- ggplot()+
  geom_point(data = B_3mo1[B_3mo1$region == "North America",], aes(x = date, y = fct_rev(clade)), alpha = 0.3, color = "#4287f5") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994.2, 2024.2, by = 4), minor_breaks = seq(1994.2, 2024.2, by = 1), labels = c("March\n1994", "March\n1998", "March\n2002", "March\n2006", "March\n2010", "March\n2014", "March\n2018", "March\n2022") )+
  scale_y_discrete(drop = F)+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(1,3)) +
  theme_bw() +
  theme(axis.title = element_blank(), plot.margin = unit(c(0,0.1,0.1,0.1), "cm"), axis.text = element_text(size = 12))

B_NA_plot <- plot_grid(B_NA_plot, B_NA_samps, ncol = 1, rel_heights = c(6.5,3.5), align = "v")

# plot B Europe ----

B_EU_std <- Rt_B_EU_hpd_3m %>%
  t() %>%
  as.data.frame() %>%
  add_column(grid_end = c(times_B_3m[1,], 2024)) %>%
  filter(grid_end > 2022.94-mean(origin_B_3m), grid_end < 2023.2)

B_EU_plot <- ggplot() +
  geom_stepribbon(data = B_EU_std, aes(ymin = lower, ymax = upper, x = grid_end), alpha = 0.5,  direction = "vh") +
  geom_hline(yintercept =1, color = "#4287f5") +
  geom_step(data = B_EU_std, aes(x = grid_end, y = med), direction = "vh") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994.2, 2024.2, by = 4), minor_breaks = seq(1994.2, 2024.2, by = 1), labels = c("March\n1994", "March\n1998", "March\n2002", "March\n2006", "March\n2010", "March\n2014", "March\n2018", "March\n2022") )+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(0,2.1)) +
  theme_bw()+
  labs(title = "B + C clade Europe"
       , y = "Rt"
       , x = "Date") +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0.1,0.1,0,0.1), "cm"), axis.text = element_text(size = 12))

B_EU_samps <- ggplot()+
  geom_point(data = B_3mo1[B_3mo1$region == "Europe",], aes(x = date, y = fct_rev(clade)), alpha = 0.3, color = "#4287f5") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994.2, 2024.2, by = 4), minor_breaks = seq(1994.2, 2024.2, by = 1), labels = c("March\n1994", "March\n1998", "March\n2002", "March\n2006", "March\n2010", "March\n2014", "March\n2018", "March\n2022") )+
  scale_y_discrete(drop = F)+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(1,3)) +
  theme_bw() +
  theme(axis.title = element_blank(), plot.margin = unit(c(0,0.1,0.1,0.1), "cm"), axis.text = element_text(size = 12))

B_EU_plot <- plot_grid(B_EU_plot, B_EU_samps, ncol = 1, rel_heights = c(6.5,3.5), align = "v")

# plot
d68_sky <- plot_grid(A_NA_plot, B_NA_plot, A_EU_plot, B_EU_plot, ncol = 1)
d68_sky

if(save_plots == T){
  ggsave("../BDMM_prime_revised/plots/Rt_3mo1.svg", d68_sky, height = 9.5, width = 7, units = "in", dpi = 700)
  ggsave("../BDMM_prime_revised/plots/Rt_3mo1.png", d68_sky, height = 10, width = 8, units = "in", dpi = 700)
}


# export as a csv
A_NA_std$clade = "A+D"
A_NA_std$region = "North America"
B_NA_std$clade = "B+C"
B_NA_std$region = "North America"
A_EU_std$clade = "A+D"
A_EU_std$region = "Europe"
B_EU_std$clade = "B+C"
B_EU_std$region = "Europe"

write.csv(bind_rows(A_NA_std, B_NA_std, A_EU_std, B_EU_std), file = "table_S3.csv", row.names = F)


## Figure 4 : Population trajectories ------------------------------------------
# read in trajectory files sampled in beast
traj_B <- read_tsv(paste0(file_root, "D68_B.1.traj"))%>%
  filter(variable == "N") 

  burnt = round(burnin*length(traj_B$Sample)) # can calculate the trajmin as a burnin 
  traj_B = traj_B[(burnt+1):length(traj_B$Sample),]  
  
traj_A <- read_tsv(paste0(file_root, "D68_A.1.traj")) %>%
  filter(variable == "N") 

burnt = round(burnin*length(traj_A$Sample)) # can calculate the trajmin as a burnin
traj_A = traj_A[(burnt+1):length(traj_A$Sample),]  
# extract across a grid of times
ages <- seq(0,40,length.out=1001)

gridded_traj_B <- traj_B %>%
  group_by(Sample, type) %>%
  reframe(N=approx(age, value, ages, method="constant", f=0, yright=0)$y,
          age=ages) %>%
  group_by(type,age) %>%
  dplyr::summarize(low=quantile(N,0.025),
                   high=quantile(N,0.975),
                   med=quantile(N,0.5))

gridded_traj_A <- traj_A %>%
  group_by(Sample, type) %>%
  reframe(N=approx(age, value, ages, method="constant", f=0, yright=0)$y,
          age=ages) %>%
  group_by(type,age) %>%
  dplyr::summarize(low=quantile(N,0.025),
                   high=quantile(N,0.975),
                   med=quantile(N,0.5))


# plot the trajectories split by location and colored by clade
gridded_traj_A$clade <- "A+D"
gridded_traj_B$clade <- "B+C"

gridded_all <- bind_rows(gridded_traj_A, gridded_traj_B)

## grey vertical stripes on final plot
stripes <- matrix(c("1994-01-01", "1995-01-01",
                    "1996-01-01", "1997-01-01",
                    "1998-01-01", "1999-01-01",
                    "2000-01-01", "2001-01-01",
                    "2002-01-01", "2003-01-01",
                    "2004-01-01", "2005-01-01",
                    "2006-01-01", "2007-01-01",
                    "2008-01-01", "2009-01-01",
                    "2010-01-01", "2011-01-01",
                    "2012-01-01", "2013-01-01",
                    "2014-01-01", "2015-01-01",
                    "2016-01-01", "2017-01-01",
                    "2018-01-01", "2019-01-01",
                    "2020-01-01", "2021-01-01",
                    "2022-01-01", "2023-01-01",
                    "2024-01-01", "2025-01-01"),ncol=2,byrow=TRUE) %>%
  data.frame()%>%
  set_colnames(c("Start", "End")) %>%
  mutate(across(everything(), as.Date))


gridded_all %<>%
  mutate(Date=ymd("2022-12-10")-age*365) %>%
  mutate(Location=recode(factor(type),
                         "1"="North America",
                         "0"="Europe"))


plot1 <- ggplot(gridded_all,
       aes(Date, N, col=clade, fill=clade, y=med)) +
    #  aes(Date, N, col=clade, fill=clade, y=med, ymin=low, ymax=high)) +
 #geom_ribbon(alpha=0.5) +
  geom_line() + ylab("Population size") +
  scale_x_date(date_breaks="4 years", limits = c(as.Date("1995-12-01"), NA), date_minor_breaks= "1 year",  date_labels="%Y") +
  scale_fill_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  scale_color_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  scale_y_continuous(expand = expansion(0,0))+
  facet_grid(rows = vars(Location))+
  theme_bw()

plot1 <- append_layers(plot1, geom_rect(data=stripes,inherit.aes=FALSE, aes(xmin=Start, xmax=End, ymin=-Inf, ymax=+Inf), fill='grey95'), position = "bottom")

plot2 <- ggplot(gridded_all,
              #  aes(Date, N, col=clade, fill=clade, y=med)) +
    aes(Date, N, col=clade, fill=clade, y=med, ymin=low, ymax=high)) +
  geom_ribbon(alpha=0.5) +
  geom_line() + ylab("Population size") +
  scale_x_date(date_breaks="4 years", limits = c(as.Date("1997-12-01"), NA), date_minor_breaks= "1 year",  date_labels="%Y") +
  scale_fill_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  scale_color_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  scale_y_continuous(expand = expansion(0,0))+
  coord_cartesian(xlim = c(as.Date("1998", format = "%Y"), as.Date("2025", format = "%Y"))) +
  facet_grid(rows = vars(Location))+
  theme_bw()+
  theme(legend.position = "bottom")


plot2 <- append_layers(plot2, geom_rect(data=stripes,inherit.aes=FALSE, aes(xmin=Start, xmax=End, ymin=-Inf, ymax=+Inf), fill='grey95'), position = "bottom")
plot2

if(save_plots == T){
ggsave("../BDMM_prime_revised/plots/pop_sky_bdmm_6mo.png", plot2, height = 5, width = 7, units = "in", dpi = 700)
}

# plot the AFM cases
# AFM cases - CDC
USA_cases <- read.csv("../BDMM_prime_revised/data/AFM_data_cdc.csv") # AFM case data from cdc.gov

USA_cases$Location <- "AFM Cases - USA"

hack_facet <- 
gridded_all %>%
  bind_rows(., USA_cases)

  hack_facet$Location <-  factor(hack_facet$Location, levels = c("AFM Cases - USA", "North America", "Europe"))

  case_plot <- ggplot(hack_facet,
         aes(Date, N, col=clade, fill=clade, y=med, ymin=low, ymax=high)) +
  geom_ribbon(alpha=0.5, linewidth = 0.5) +
  geom_line() + ylab("Population size") +
  scale_x_date(date_breaks="4 years", limits = c(as.Date("1997-12-01"), NA), date_minor_breaks= "1 year",  date_labels="%Y") +
  scale_fill_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  scale_color_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  scale_y_continuous(expand = expansion(0,0))+
  coord_cartesian(xlim = c(as.Date("2000", format = "%Y"), as.Date("2024", format = "%Y"))) +
  geom_col(aes(x = as.Date(Month, format = "%m/%d/%Y"), y = Cases), fill = "navyblue", show.legend = F, width = 31)+
  facet_grid(rows = vars(Location), scales = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom")
    


traj_afm <- append_layers(case_plot, geom_rect(data=stripes,inherit.aes=FALSE, aes(xmin=Start, xmax=End, ymin=-Inf, ymax=+Inf), fill='grey95'), position = "bottom")
traj_afm

if(save_plots == T){
ggsave("../BDMM_prime_revised/plots/AFM_traj_3mo1.svg", traj_afm, height = 5, width = 7, units = "in", dpi = 700)
ggsave("../BDMM_prime_revised/plots/AFM_traj_3mo1.png", traj_afm, height = 5, width = 7, units = "in", dpi = 700)
}


## Figure 5 : Correlations -----------------------------------------------------

# Time series cross-correlation
USA_cases_grp <- USA_cases %>%
  mutate(Date = as.Date(Month, format = "%m/%d/%Y"),
         group_date = floor_date(Date, "month")) %>%
  group_by(group_date, .add = T)%>%
  summarise(group_case= sum(Cases))

# convert to monthly to match AFM data
grid_grp <- gridded_all %>%
  mutate(group_date = floor_date(Date, "month"))%>%
  group_by(group_date, clade, Location)%>%
  summarise(med= sum(med))


corr_data_NA <- left_join(grid_grp %>% filter(Location == "North America"), USA_cases_grp, by = "group_date")

corLag_A_NA <-ccf(corr_data_NA[corr_data_NA$clade == "A+D",]$med,corr_data_NA[corr_data_NA$clade == "A+D",]$group_case,lag.max=6, na.action = na.pass, pl = F)

corLag_B_NA <-ccf(corr_data_NA[corr_data_NA$clade == "B+C",]$med,corr_data_NA[corr_data_NA$clade == "B+C",]$group_case,lag.max=6, na.action = na.pass, pl = F)


corr_data_EU <- left_join(grid_grp %>% filter(Location == "Europe"), USA_cases_grp, by = "group_date")

corLag_A_EU <-ccf(corr_data_EU[corr_data_EU$clade == "A+D",]$med,corr_data_EU[corr_data_EU$clade == "A+D",]$group_case,lag.max=6, na.action = na.pass, pl = F)

corLag_B_EU <-ccf(corr_data_EU[corr_data_EU$clade == "B+C",]$med,corr_data_EU[corr_data_EU$clade == "B+C",]$group_case,lag.max=6, na.action = na.pass, pl = F)


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
us_correlation <- ggarrange(A_NA_corr, B_NA_corr, ncol = 2)

# sliding window

# sliding window - B
xy <- cbind(corr_data_NA[corr_data_NA$clade=="B+C",]$med, corr_data_NA[corr_data_NA$clade=="B+C",]$group_case)
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
xy <- cbind(corr_data_NA[corr_data_NA$clade=="A+D",]$med, corr_data_NA[corr_data_NA$clade=="A+D",]$group_case)
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


us_slide_correlation <- ggarrange(A_NA_corr, B_NA_corr, a_slide, b_slide, align = "v", common.legend = T, legend = "bottom", labels = "AUTO")
if(save_plots == T){
ggsave("../BDMM_prime_revised/plots/slide_corr.png", us_slide_correlation, height = 5.3, width = 8, units = "in")
ggsave("../BDMM_prime_revised/plots/slide_corr.svg", us_slide_correlation, height = 5.3, width = 8, units = "in")
}
us_slide_correlation

# migration posterior difference
mig1 <- data.frame(beastio::getLogFileSubset(tr_3mo1, "migrationRateSPEpi"))
mig1$diff <- mig1$migrationRateSPEpi.AEurope_to_North.America -  mig1$migrationRateSPEpi.ANorth.America_to_Europe
(1/nrow(mig1))*sum(mig1$diff)
quantile(mig1$diff, c(0.025, 0.975))

table(mig1$migrationRateSPEpi.AEurope_to_North.America > mig1$migrationRateSPEpi.ANorth.America_to_Europe)


## Europe correlations - only have incomplete annual data ----
#Europe cases D68 - helfferich et al. 2025 (eurosurv)
EUR_AFM <- data.frame(Cases = c(37,5,39,13,3,6,15,12), Month = seq(as.Date("2016-01-01"), as.Date("2023-01-01"), by = "year"))
EUR_68 <- data.frame(Cases = c(15,0,12,7,2,3,8,1), Month = seq(as.Date("2016-01-01"), as.Date("2023-01-01"), by = "year"))

ggplot()+
  geom_col(data = EUR_AFM,aes(y=Cases, x=Month), fill = "gold2", alpha = 0.3, width = 365, position = position_nudge(x = 184))+
  geom_col(data = EUR_68,aes(y=Cases, x=Month), fill = "gold2", alpha = 1, width = 365, position = position_nudge(x = 184))+
  scale_y_continuous(expand = c(0,0), limits = c(0,40))+
  scale_x_date(expand = c(0,0), date_breaks = "1 year", date_labels = "%Y")+
  labs(title = "AFM cases - Europe (Helfferich et al. 2025 Eurosurv.)"
       , y = "AFM Cases"
       , x = "Date") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank())

ggsave("./plots/AFM_europe_suppl_D68.png", dpi = 700, width = 6, height = 3, units = "in")
ggsave("./plots/AFM_europe_suppl_D68.svg", dpi = 700, width = 6, height = 3, units = "in")

# correlation of AFM with predicted cases - Europe
# convert to monthly to match AFM data
grid_grp_yr <- gridded_all %>%
  mutate(group_date_yr = floor_date(Date, "year"))%>%
  group_by(group_date_yr, clade, Location)%>%
  summarise(med= sum(med))

EUR_AFM$group_date_yr <- EUR_AFM$Month

corr_data_EU <- left_join(grid_grp_yr %>% filter(Location == "Europe"), EUR_AFM, by = "group_date_yr")

corLag_A_EU <-ccf(corr_data_EU[corr_data_EU$clade == "A+D",]$med,corr_data_EU[corr_data_EU$clade == "A+D",]$Cases,lag.max=2, na.action = na.pass, pl = F)

corLag_B_EU <-ccf(corr_data_EU[corr_data_EU$clade == "B+C",]$med,corr_data_EU[corr_data_EU$clade == "B+C",]$Cases,lag.max=2, na.action = na.pass, pl = F)


A_EU_corr <- do.call(cbind.data.frame, corLag_A_EU) %>%
  select(acf, lag, n.used) %>%
  ggplot(mapping = aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2 / sqrt(7), color = "#f54b42", linetype = 2) +
  geom_hline(yintercept = -(2 / sqrt(7)), color = "#f54b42",  linetype = 2) +
  geom_col(width = .1) +
  scale_y_continuous(limits = c(-1,1))+
  labs(title = "A + D clade Europe"
       , y = "Cross-correlation"
       , x = "Lag (years)") +
  theme_bw()

B_EU_corr<- do.call(cbind.data.frame, corLag_B_EU) %>%
  select(acf, lag, n.used) %>%
  ggplot(mapping = aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2 / sqrt(7), color = "#4287f5", linetype = 2) +
  geom_hline(yintercept = -(2 / sqrt(7)), color = "#4287f5",  linetype = 2) +
  geom_col(width = .1) +
  scale_y_continuous(limits = c(-1,1))+
  labs(title = "B + C clade Europe"
       , y = "Cross-correlation"
       , x = "Lag (years)") +
  theme_bw()


# 2020 cutoff
corr_data_EU <- left_join(grid_grp_yr %>% filter(Location == "Europe"), EUR_AFM, by = "group_date_yr")
corr_data_EU %<>% filter(group_date_yr<as.Date("2020-01-01"))

corLag_A_EU <-ccf(corr_data_EU[corr_data_EU$clade == "A+D",]$med,corr_data_EU[corr_data_EU$clade == "A+D",]$Cases,lag.max=2, na.action = na.pass, pl = F)

corLag_B_EU <-ccf(corr_data_EU[corr_data_EU$clade == "B+C",]$med,corr_data_EU[corr_data_EU$clade == "B+C",]$Cases,lag.max=2, na.action = na.pass, pl = F)


A_EU_corr_pre <- do.call(cbind.data.frame, corLag_A_EU) %>%
  select(acf, lag, n.used) %>%
  ggplot(mapping = aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2 / sqrt(4), color = "#f54b42", linetype = 2) +
  geom_hline(yintercept = -(2 / sqrt(4)), color = "#f54b42",  linetype = 2) +
  geom_col(width = .1) +
  scale_y_continuous(limits = c(-1,1))+
  labs(title = "A + D clade Europe pre-2020"
       , y = "Cross-correlation"
       , x = "Lag (years)") +
  theme_bw()

B_EU_corr_pre<- do.call(cbind.data.frame, corLag_B_EU) %>%
  select(acf, lag, n.used) %>%
  ggplot(mapping = aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2 / sqrt(4), color = "#4287f5", linetype = 2) +
  geom_hline(yintercept = -(2 / sqrt(4)), color = "#4287f5",  linetype = 2) +
  geom_col(width = .1) +
  scale_y_continuous(limits = c(-1,1))+
  labs(title = "B + C clade Europe pre-2020"
       , y = "Cross-correlation"
       , x = "Lag (years)") +
  theme_bw()

corr_data_EU <- left_join(grid_grp_yr %>% filter(Location == "Europe"), EUR_AFM, by = "group_date_yr")
corr_data_EU %<>% filter(group_date_yr>=as.Date("2020-01-01"))

corLag_A_EU <-ccf(corr_data_EU[corr_data_EU$clade == "A+D",]$med,corr_data_EU[corr_data_EU$clade == "A+D",]$Cases,lag.max=2, na.action = na.pass, pl = F)

corLag_B_EU <-ccf(corr_data_EU[corr_data_EU$clade == "B+C",]$med,corr_data_EU[corr_data_EU$clade == "B+C",]$Cases,lag.max=2, na.action = na.pass, pl = F)


A_EU_corr_post <- do.call(cbind.data.frame, corLag_A_EU) %>%
  select(acf, lag, n.used) %>%
  ggplot(mapping = aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2 / sqrt(3), color = "#f54b42", linetype = 2) +
  geom_hline(yintercept = -(2 / sqrt(3)), color = "#f54b42",  linetype = 2) +
  geom_col(width = .1) +
  scale_y_continuous(limits = c(-1,1))+
  labs(title = "A + D clade Europe post-2020"
       , y = "Cross-correlation"
       , x = "Lag (years)") +
  theme_bw()

B_EU_corr_post<- do.call(cbind.data.frame, corLag_B_EU) %>%
  select(acf, lag, n.used) %>%
  ggplot(mapping = aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2 / sqrt(3), color = "#4287f5", linetype = 2) +
  geom_hline(yintercept = -(2 / sqrt(3)), color = "#4287f5",  linetype = 2) +
  geom_col(width = .1) +
  scale_y_continuous(limits = c(-1,1))+
  labs(title = "B + C clade Europe post-2020"
       , y = "Cross-correlation"
       , x = "Lag (years)") +
  theme_bw()