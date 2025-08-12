# Plot Rt trajectories for each superclade and location

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
  , tidyr
  , tibble
  , ggplot2 
  , coda
  , bdskytools
  , beastio
  , pammtools
  , forcats
  , ggpubr
)


# metadata
load("./samples/samp_3mo1/A_sub_3mo_1.rdata")
A_3mo1 <- A_meta_strat_date_subs
load("./samples/samp_3mo1/B_sub_3mo_1.rdata")
B_3mo1 <- B_meta_strat_date_subs

# logfile
tr_3m01 <- beastio::readLog("./samples/samp_3mo1/dated_rt.log", burnin=0.4)

load( "./allseq_labels/allseq.rdata")
all_meta$date <- all_meta$V3

all_meta$plot <- revalue(all_meta$assignment, c("B3-lin1" = "B",
                                                "C" = "C",
                                                "B3" = "B",
                                                "B1" = "B",
                                                "B3-lin2" = "B",
                                                "B+C" = "Other",
                                                "B" = "B",
                                                "B2" = "B",
                                                "B3-US18" = "B",
                                                "B3-EU19" = "B",
                                                "D-19" = "D",
                                                "B3-EU18" = "B",
                                                "A+D" = "Other",
                                                "A1" = "A",
                                                "D" = "D",
                                                "A" = "A"))
all_meta$names <- gsub(" ", "-", all_meta$names)

A_3mo1 <- left_join(A_3mo1[, c("newlabs", "dategroup")], all_meta, by = join_by("newlabs"=="names"))
B_3mo1 <- left_join(B_3mo1[, c("newlabs", "dategroup")], all_meta, by = join_by("newlabs"=="names"))
A_3mo1$plot <- factor(A_3mo1$plot)
B_3mo1$plot <- factor(B_3mo1$plot)

## split A/B and NA/EU
Rt_A_NA_3m <- beastio::getLogFileSubset(tr_3m01, "ReSPEpi.Ai\\d+_Northern.America")
Rt_A_EU_3m <- beastio::getLogFileSubset(tr_3m01, "ReSPEpi.Ai\\d+_Western.Europe")

times_A_3m <- beastio::getLogFileSubset(tr_3m01, "ReSPEpi.Ai\\d+_endtime")
origin_A_3m <- beastio::getLogFileSubset(tr_3m01, "originBDMMPrime.A")

# convert from relative to origin to relative to MRSD 
# seems origin is given backwards in time from the MRSD and changepoints forward in time from the origin?
reltimes_A_3m <- -origin_A_3m + times_A_3m

# relative to MRSD - match the inputs!
times_A_3m <- 2021.981 + reltimes_A_3m

Rt_B_NA_3m <- beastio::getLogFileSubset(tr_3m01, "ReSPEpi.Bi\\d+_Northern.America")
Rt_B_EU_3m <- beastio::getLogFileSubset(tr_3m01, "ReSPEpi.Bi\\d+_Western.Europe")

times_B_3m <- beastio::getLogFileSubset(tr_3m01, "ReSPEpi.Bi\\d+_endtime")
origin_B_3m <- beastio::getLogFileSubset(tr_3m01, "originBDMMPrime.B")

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
  filter(grid_end > 2021.98-mean(origin_A_3m), grid_end < 2022.2)

A_NA_plot <- ggplot() +
  geom_rect(aes(xmin = 2014, xmax = 2015, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_rect(aes(xmin = 2016, xmax = 2017, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_rect(aes(xmin = 2018, xmax = 2019, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_stepribbon(data = A_NA_std, aes(ymin = lower, ymax = upper, x = grid_end), alpha = 0.5,  direction = "vh") +
  geom_hline(yintercept =1, color = "#f54b42") +
  geom_step(data = A_NA_std, aes(x = grid_end, y = med), direction = "vh") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(0,2)) +
  theme_bw()+
  labs(title = "A + D clade North America"
       , y = "Rt"
       , x = "Date")+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0.1,0.1,0,0.1), "cm"), axis.text = element_text(size = 12))

A_NA_samps <- ggplot()+
  geom_rect(aes(xmin = 2014, xmax = 2015, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_rect(aes(xmin = 2016, xmax = 2017, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_rect(aes(xmin = 2018, xmax = 2019, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_point(data = A_3mo1[A_3mo1$region == "North America",], aes(x = date, y = fct_rev(plot)), alpha = 0.3, color = "#f54b42") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
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
  filter(grid_end > 2021.98-mean(origin_A_3m), grid_end < 2022.2)

A_EU_plot <- ggplot() +
  geom_stepribbon(data = A_EU_std, aes(ymin = lower, ymax = upper, x = grid_end), alpha = 0.5,  direction = "vh") +
  geom_hline(yintercept =1, color = "#f54b42") +
  geom_step(data = A_EU_std, aes(x = grid_end, y = med), direction = "vh") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(0,2)) +
  theme_bw()+
  labs(title = "A + D clade Western Europe"
       , y = "Rt"
       , x = "Date")+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0.1,0.1,0,0.1), "cm"), axis.text = element_text(size = 12))

A_EU_samps <- ggplot()+
  geom_point(data = A_3mo1[A_3mo1$region == "Europe",], aes(x = date, y = fct_rev(plot)), alpha = 0.3, color = "#f54b42") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
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
  filter(grid_end > 2021.98-mean(origin_B_3m), grid_end < 2022.2)

B_NA_plot <- ggplot() +
  geom_rect(aes(xmin = 2014, xmax = 2015, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_rect(aes(xmin = 2016, xmax = 2017, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_rect(aes(xmin = 2018, xmax = 2019, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_stepribbon(data = B_NA_std, aes(ymin = lower, ymax = upper, x = grid_end), alpha = 0.5,  direction = "vh") +
  geom_hline(yintercept =1, color = "#4287f5") +
  geom_step(data = B_NA_std, aes(x = grid_end, y = med), direction = "vh") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(0,2)) +
  theme_bw()+
  labs(title = "B + C clade North America"
       , y = "Rt"
       , x = "Date")+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0.1,0.1,0,0.1), "cm"), axis.text = element_text(size = 12))

B_NA_samps <- ggplot()+
  geom_rect(aes(xmin = 2014, xmax = 2015, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_rect(aes(xmin = 2016, xmax = 2017, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_rect(aes(xmin = 2018, xmax = 2019, ymin = -Inf, ymax = Inf), fill = "grey85")+
  geom_point(data = B_3mo1[B_3mo1$region == "North America",], aes(x = date, y = fct_rev(plot)), alpha = 0.3, color = "#4287f5") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
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
  filter(grid_end > 2021.98-mean(origin_B_3m), grid_end < 2022.2)

B_EU_plot <- ggplot() +
  geom_stepribbon(data = B_EU_std, aes(ymin = lower, ymax = upper, x = grid_end), alpha = 0.5,  direction = "vh") +
  geom_hline(yintercept =1, color = "#4287f5") +
  geom_step(data = B_EU_std, aes(x = grid_end, y = med), direction = "vh") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(0,2)) +
  theme_bw()+
  labs(title = "B + C clade Europe"
       , y = "Rt"
       , x = "Date") +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0.1,0.1,0,0.1), "cm"), axis.text = element_text(size = 12))

B_EU_samps <- ggplot()+
  geom_point(data = B_3mo1[B_3mo1$region == "Europe",], aes(x = date, y = fct_rev(plot)), alpha = 0.3, color = "#4287f5") +
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
  scale_y_discrete(drop = F)+
  coord_cartesian(xlim = c(1996, 2025), ylim = c(1,3)) +
  theme_bw() +
  theme(axis.title = element_blank(), plot.margin = unit(c(0,0.1,0.1,0.1), "cm"), axis.text = element_text(size = 12))

B_EU_plot <- plot_grid(B_EU_plot, B_EU_samps, ncol = 1, rel_heights = c(6.5,3.5), align = "v")

# plot ----
d68_sky <- plot_grid(A_NA_plot, A_EU_plot, B_NA_plot, B_EU_plot, ncol = 1)
d68_sky


