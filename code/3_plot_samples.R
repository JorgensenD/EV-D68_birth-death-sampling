# pull in all the metafiles and plot samples
pacman::p_load(tidyverse,
               ggpubr,
               grid,
               ggplot2,
               tabyl, 
               flextable)
# samples used in each run
load("./samples/samp_3mo1/A_sub_3mo_1.rdata")
A_3mo1 <- A_meta_strat_date_subs
load("./samples/samp_3mo1/B_sub_3mo_1.rdata")
B_3mo1 <- B_meta_strat_date_subs

all_3mo1 <- rbind(A_3mo1, B_3mo1) %>%
  mutate(dategroup = as.Date(dategroup),
         newregion = as.factor(newregion))

plt_3mo1 <- ggplot(all_3mo1)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F)+
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(limits = c(as.Date("1996-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#4287f5", "#f54b42"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(newregion ~ color)+
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
         newregion = as.factor(newregion))

plt_3mo2 <- ggplot(all_3mo2)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F, width = 93)+
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(limits = c(as.Date("1996-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#4287f5", "#f54b42"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(newregion ~ color)+
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
         newregion = as.factor(newregion))

plt_6mo1 <-ggplot(all_6mo1)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F, width = 186)+
  scale_y_continuous(expand = c(0,0), breaks = c(0, 20, 40)) +
  scale_x_date(limits = c(as.Date("1996-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#4287f5", "#f54b42"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(newregion ~ color)+
  theme_bw()+
  labs(x = NULL,
       y = NULL,
       title = "Uniform sample (January) 1") +
  theme(panel.grid.minor = element_blank())



load("./samples/samp_6mo2/A_sub_6mo_2.rdata")
A_6mo2 <- A_meta_strat_date_subs
load("./samples/samp_6mo2/B_sub_6mo_2.rdata")
B_6mo2 <- B_meta_strat_date_subs

all_6mo2 <- rbind(A_6mo2, B_6mo2) %>%
  mutate(dategroup = as.Date(dategroup),
         newregion = as.factor(newregion))

plt_6mo2 <- ggplot(all_6mo2)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F, width = 186)+
  scale_y_continuous(expand = c(0,0), breaks = c(0, 20, 40)) +
  scale_x_date(limits = c(as.Date("1996-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#4287f5", "#f54b42"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(newregion ~ color)+
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
         newregion = as.factor(newregion))

plt_6mooff1 <- ggplot(all_6mooff1)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F, width = 186)+
  scale_y_continuous(expand = c(0,0), breaks = c(0, 20, 40)) +
  scale_x_date(limits = c(as.Date("1996-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#4287f5", "#f54b42"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(newregion ~ color)+
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
         newregion = as.factor(newregion))

plt_6mooff2 <- ggplot(all_6mooff2)+
  geom_bar(aes(x = dategroup , fill = color), show.legend = F, width = 186)+
  scale_y_continuous(expand = c(0,0), breaks = c(0, 20, 40)) +
  scale_x_date(limits = c(as.Date("1996-01-01"), NA), date_breaks = "8 years", date_labels = "%Y")+
  scale_fill_manual(values = c("#4287f5", "#f54b42"), na.value = NA, na.translate = F, name = "Clade") +
  facet_grid(newregion ~ color)+
  theme_bw()+
  labs(x = NULL,
       y = NULL,
       title = "Uniform sample (March) 2") +
  theme(panel.grid.minor = element_blank())

allsamps <- ggarrange(plt_3mo1, plt_3mo2, plt_6mo1, plt_6mo2, plt_6mooff1, plt_6mooff2, ncol = 2, nrow = 3, labels = "AUTO")

allsamps <- annotate_figure(allsamps, left = textGrob("Number of sequences", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                            bottom = textGrob("Date", gp = gpar(cex = 1)))


ggsave("./plots/samplingschemes_suppl.png", dpi = 700, height = 9, width = 8, units = "in")



# violin plot of the distribution of sequence lengths
algn_retained <- read.FASTA("./allseq_labels/algn_retained.fasta")
lengthdata <- data.frame(names = names(algn_retained), ungapped_lengths = sapply(algn_retained, function(x) { sum(x != as.raw(04))}))
# add all metadata using the sequence names
allmeta_used$label_clean <- gsub("/", "", allmeta_used$newlabs)
lengthdata <- left_join(lengthdata, allmeta_used, by = c("names" = "label_clean"))


lengthdata %>%
  filter(region %in% c("Europe", "North America")) %>%
  ggplot(aes(x = ungapped_lengths, y =region)) +
  geom_violin(aes(color = region))+
  geom_point(aes(color = region), position = position_jitter(seed = 1, width = 0.2)) +
  scale_x_continuous(limits = c(300,950), breaks = seq(300,1000, by = 100)) +
  labs(x = "Available VP1 sequence lengths (nucleotides)") +
  scale_color_manual(labels = c("Europe", "North America"), values = c("gold2", "mediumseagreen"), na.value = NA, na.translate = F, name = "Region")+
  theme_bw()+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) 

ggsave("./plots/seqlengths_suppl.png", width = 6, height = 3, units = "in")

# table of included sequences by country and clade

meta_table_S1 <- allmeta_used %>%
  filter(newregion %in% c("Europe", "North America")) %>%
  dplyr::select(country, color) %>%
  tabyl(country, color) %>%
  adorn_totals("row")%>%
  flextable() %>%
  add_header_row(top = TRUE, values = c("Country", "Superclade", "")) %>%
  set_header_labels(country = "", `B/C` = "B/C", `A/D` = "A/D") %>%
  merge_at(i = 1, j = 2:3, part = "header") %>%
  border_remove() %>%  
  theme_booktabs() %>%
  flextable::align(align = "center", j = c(2:3), part = "header") %>%
  bold(i= 1, bold = TRUE, part = "header") %>%
  bold(i = 19, bold = TRUE, part = "body") %>%
  merge_at(i = 1:2, j = 1, part = "header") %>%
  vline(part = "all", j = 1) %>%
  width(j=1, width = 1.5) %>%
  color(color = "#f54b42",i=2, j=2, part = "header")%>%
  color(color = "#4287f5",i=2, j=3, part = "header")
  
save_as_docx("Number of sequences by superclade" = meta_table_S1, path = "./plots/tableS1.docx")

  





