# pull in all the metafiles and plot samples
pacman::p_load(tidyverse,
               ggpubr,
               grid,
               ggplot2)

# samples used in each run
load("./samples/samp_3mo1/A_sub_3mo_1.rdata")
A_3mo1 <- A_meta_strat_date_subs
load("./samples/samp_3mo1/B_sub_3mo_1.rdata")
B_3mo1 <- B_meta_strat_date_subs

all_3mo1 <- rbind(A_3mo1, B_3mo1) %>%
  mutate(dategroup = as.Date(dategroup),
         region = as.factor(region))

plt_3mo1 <- ggplot(all_3mo1)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F)+
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(limits = c(as.Date("1992-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(region ~ color)+
  theme_bw()+
  labs(x = NULL,
       y = NULL,
       title = "Proportional sample 1") +
  theme(panel.grid.minor = element_blank())

load("./samples/samp_3mo2/A_sub_3mo_2.rdata")
A_3mo2 <- A_meta_strat_date_subs
load("./samples/samp_3mo2/B_sub_3mo_2.rdata")
B_3mo2 <- B_meta_strat_date_subs

all_3mo2 <- rbind(A_3mo2, B_3mo2) %>%
  mutate(dategroup = as.Date(dategroup),
         region = as.factor(region))

plt_3mo2 <- ggplot(all_3mo2)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F, width = 93)+
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(limits = c(as.Date("1992-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(region ~ color)+
  theme_bw()+
  labs(x = NULL,
       y = NULL,
       title = "Proportional sample 2") +
  theme(panel.grid.minor = element_blank())


load("./samples/samp_6mo1/A_sub_6mo_1.rdata")
A_6mo1 <- A_meta_strat_date_subs
load("./samples/samp_6mo1/B_sub_6mo_1.rdata")
B_6mo1 <- B_meta_strat_date_subs

all_6mo1 <- rbind(A_6mo1, B_6mo1) %>%
  mutate(dategroup = as.Date(dategroup),
         region = as.factor(region))

plt_6mo1 <-ggplot(all_6mo1)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F, width = 186)+
  scale_y_continuous(expand = c(0,0), breaks = c(0, 20, 40)) +
  scale_x_date(limits = c(as.Date("1992-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(region ~ color)+
  theme_bw()+
  labs(x = NULL,
       y = NULL,
       title = "Uniform sample (January) 1") +
  theme(panel.grid.minor = element_blank())



load("./samples/samp_6mo2/A_sub_6mo_2.rdata")
A_6mo2 <- A_meta_strat_date_subs
load("./samples/samp_6mo2/B_sub_6mo_2.rdata")
B_6mo2 <- B_meta_strat_date_subs

plt_6mo2 <- ggplot(all_6mo2)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F, width = 186)+
  scale_y_continuous(expand = c(0,0), breaks = c(0, 20, 40)) +
  scale_x_date(limits = c(as.Date("1992-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(region ~ color)+
  theme_bw()+
  labs(x = NULL,
       y = NULL,
       title = "Uniform sample (January) 2") +
  theme(panel.grid.minor = element_blank())

load("./samples/samp_6mo1_off/A_sub_6mo_offset_1.rdata")
A_6mooff1 <- A_meta_strat_date_subs
load("./samples/samp_6mo1_off/B_sub_6mo_offset_1.rdata")
B_6mooff1 <- B_meta_strat_date_subs

all_6mooff1 <- rbind(A_6mooff1, B_6mooff1) %>%
  mutate(dategroup = as.Date(dategroup),
         region = as.factor(region))

plt_6mooff1 <- ggplot(all_6mooff1)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F, width = 186)+
  scale_y_continuous(expand = c(0,0), breaks = c(0, 20, 40)) +
  scale_x_date(limits = c(as.Date("1992-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(region ~ color)+
  theme_bw()+
  labs(x = NULL,
       y = NULL,
       title = "Uniform sample (March) 1") +
  theme(panel.grid.minor = element_blank())

load("./samples/samp_6mo2_off/A_sub_6mo_offset_2.rdata")
A_6mooff2 <- A_meta_strat_date_subs
load("./samples/samp_6mo2_off/B_sub_6mo_offset_2.rdata")
B_6mooff2 <- B_meta_strat_date_subs

all_6mooff2 <- rbind(A_6mooff2, B_6mooff2) %>%
  mutate(dategroup = as.Date(dategroup),
         region = as.factor(region))

plt_6mooff2 <- ggplot(all_6mooff2)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F, width = 186)+
  scale_y_continuous(expand = c(0,0), breaks = c(0, 20, 40)) +
  scale_x_date(limits = c(as.Date("1992-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#f54b42", "#4287f5"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(region ~ color)+
  theme_bw()+
  labs(x = NULL,
       y = NULL,
       title = "Uniform sample (March) 2") +
  theme(panel.grid.minor = element_blank())

allsamps <- ggarrange(plt_3mo1, plt_3mo2, plt_6mo1, plt_6mo2, plt_6mooff1, plt_6mooff2, ncol = 2, nrow = 3, labels = "AUTO")

allsamps <- annotate_figure(allsamps, left = textGrob("Number of sequences", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                            bottom = textGrob("Date", gp = gpar(cex = 1)))


ggsave("./plots/samplingschemes_suppl.png", dpi = 700, height = 9, width = 8, units = "in")
