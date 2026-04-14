pacman::p_load(
  ape,
  treeio,
  ggtree,
  ggplot2,
  lubridate,
  phytools,
  plyr,
  dplyr,
  tidyr
)

# load sequences - not included on GitHub but can be downloaded via the included accession numbers
algn_VP1_300 <- algn <- read.FASTA("./data/algn_VP1_300_dec25.fasta", type = "DNA")

# load metadata
metadata <- read.csv("./data/metafile_dec25.csv")

# load tree - manually colored by superclade in figtree
tr_ann <- read.beast("./data/color_tree.tree")

# assign a group to each sequence (in metadata file)
treemetafile <- tr_ann@data[is.na(tr_ann@data$label),]
treemetafile$nodelabel <- tr_ann@phylo$tip.label
treemetafile$`!color` <- as.factor(treemetafile$`!color`) 
# these colors depend on what you assign in figtree
treemetafile$color <- revalue(treemetafile$`!color`, c("#f0231f" = "A/D", "#bf0063" = "A/D", "#ba7918" = "A/D", "#0c91f0" = "B/C", "#27b8b3" = "B/C", "#0d1bb8" = "B/C"))
treemetafile$clade <- revalue(treemetafile$`!color`, c("#f0231f" = "Other", "#bf0063" = "A", "#ba7918" = "D", "#0c91f0" = "Other", "#27b8b3" = "C", "#0d1bb8" = "B"))
treemetafile$node <- as.numeric(treemetafile$node)

# order to match the alignment and output two alignments for each super-clade
allmeta <- left_join(treemetafile, metadata, by = c("nodelabel"="accession"))

# keep only US + Canada - check for other countries
allmeta %>% 
  filter(region =="North America") %>%
  distinct(country)


# make new regions
allmeta %<>%
  mutate(newregion = case_when(
    region == "North America" ~ "North America",
    region == "Europe" ~ "Europe",
    TRUE ~ "External"
  ))

allmeta$newlabs <- paste(allmeta$nodelabel, gsub(" ","-", allmeta$newregion), allmeta$color, gsub(" ","-", allmeta$country), round(as.numeric(allmeta$V3),3), sep = "_")
allmeta$date <- as.numeric(allmeta$date)

names(algn_VP1_300) <- allmeta[match(names(algn_VP1_300), allmeta$nodelabel),]$newlabs
# more descriptive sequence names to run in BEAST later
write.FASTA(algn_VP1_300, "./allseq_labels/allseq.fasta")
save(allmeta, file = "./allseq_labels/allseq.rdata")
write.csv(allmeta%>%select(nodelabel, color, clade, dec.date, country, V3, region, newregion, newlabs), file = "./allseq_labels/metafile.csv")

# check for majorly ambiguous dates - we won't use these for analysis
ambig_retain <- allmeta %>%
  # search for the ambiguous date format
  filter(grepl("\\[", dec.date)) %>%
  # split this into two columns
  mutate(year_range_clean = gsub("\\[|\\]", "", dec.date)) %>%
  separate(year_range_clean, into = c("year_start", "year_end"), sep = ":", convert = TRUE) %>%
  filter(year_end - year_start <= 0.17) # remove ones which are highly variable - with the fast clock rate large ranges will not be easily resolved in the model.
 ## NOW DROP ANYTHING WITH 2mo+ VARIABILITY.


# drop the over-ambiguous ones
  allmeta %<>%
  anti_join(allmeta %>%
  # search for the ambiguous date format
  filter(grepl("\\[", dec.date)) %>%
  # split this into two columns
  mutate(year_range_clean = gsub("\\[|\\]", "", dec.date)) %>%
  separate(year_range_clean, into = c("year_start", "year_end"), sep = ":", convert = TRUE) %>%
  filter(year_end - year_start > 0.17))
  
# Drop post-2022 sequences
allmeta %<>% filter(date<2023, country != "Bermuda")

allmeta_used <- allmeta
save(allmeta_used, file = "./allseq_labels/allseq_retained.rdata")
algn_retained <- algn_VP1_300[allmeta_used$newlabs]
# rename to remove slash for iqtree to treetime
names(algn_retained) <- gsub("/", "", names(algn_retained))

# write tipdates for treetime
tipdates <- allmeta_used %>%
  select(newlabs, date) %>% # select correct columns
  dplyr::rename(accession = newlabs, 
                dec.date = date
  )
tipdates$accession <- gsub("/","", tipdates$accession)


write.csv(tipdates, "./allseq_labels/tipdates_retained.csv", row.names = F, quote = F)

write.FASTA(algn_retained, "./allseq_labels/algn_retained.fasta")

# subset the two regions and two clades - make fasta and metafile for each
B_data <- allmeta %>%
  filter((newregion == "Europe" | newregion == "North America") & color == "B/C")

A_data <- allmeta %>%
  filter((newregion == "Europe" | newregion == "North America") & color == "A/D")

# proportional subsampling approach (run twice) ---- 

B_meta_strat_date_subs <- B_data[-which.max(B_data$date),] %>%
  group_by(newregion) %>%
  dplyr::arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date), "3 months")) %>% 
  group_by(dategroup, .add = T) %>%
  group_modify(~ .x %>% slice_sample(n = max(1, ceiling(0.4 * nrow(.x))))) # 40% of available sequences maximum

# keep most recent sequence always - helps to align the dates
B_meta_strat_date_subs <- rbind(B_meta_strat_date_subs, B_data[which.max(B_data$date),]) 

B_subseqs_new <- algn_VP1_300[B_meta_strat_date_subs$newlabs]
# write sequence
write.FASTA(B_subseqs_new, "./samples/samp_3mo2/B_sub_3mo_2.fasta")
# write metadata
save(B_meta_strat_date_subs, file = "./samples/samp_3mo2/B_sub_3mo_2.rdata")

# check the included dates
range(B_meta_strat_date_subs$date) 

# B only - for supplement
B_only_meta <- B_meta_strat_date_subs %>%
  filter(clade == "B")

B_only_seqs <- algn_VP1_300[B_only_meta$newlabs]
# write sequence
write.FASTA(B_only_seqs, "./samples/B_only/B_only.fasta")
# write metadata
save(B_only_meta, file = "./samples/B_only/B_only.rdata")

A_meta_strat_date_subs <- A_data[-which.max(A_data$date),] %>%
  group_by(newregion) %>%
  dplyr::arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date), "3 months")) %>%
  group_by(dategroup, .add = T)%>%
  group_modify(~ .x %>% slice_sample(n = max(1, ceiling(0.4 * nrow(.x)))))

A_meta_strat_date_subs <- rbind(A_meta_strat_date_subs, A_data[which.max(A_data$date),])

A_subseqs_new <- algn_VP1_300[A_meta_strat_date_subs$newlabs]
write.FASTA(A_subseqs_new, "./samples/samp_3mo2/A_sub_3mo_2.fasta")
save(A_meta_strat_date_subs, file = "./samples/samp_3mo2/A_sub_3mo_2.rdata")

range(A_meta_strat_date_subs$date)

# Uniform subsampling approach ----

B_meta_strat_date_subs <- B_data[-which.max(B_data$date),] %>%
  group_by(newregion) %>%
  dplyr::arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date), "6 months")) %>%
  group_by(dategroup, .add = T)%>%
  slice_sample(n = 50)
B_meta_strat_date_subs <- rbind(B_meta_strat_date_subs, B_data[which.max(B_data$date),])


B_subseqs_new <- algn_VP1_300[B_meta_strat_date_subs$newlabs]
write.FASTA(B_subseqs_new, "./samples/samp_6mo1/B_sub_6mo_1.fasta")
save(B_meta_strat_date_subs, file = "./samples/samp_6mo1/B_sub_6mo_1.rdata")

range(B_meta_strat_date_subs$date)

A_meta_strat_date_subs <- A_data[-which.max(A_data$date),] %>%
  group_by(newregion) %>%
  dplyr::arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date), "6 months")) %>%
  group_by(dategroup, .add = T)%>%
  slice_sample(n = 50)
A_meta_strat_date_subs <- rbind(A_meta_strat_date_subs, A_data[which.max(A_data$date),])

range(A_meta_strat_date_subs$date)

A_subseqs_new <- algn_VP1_300[A_meta_strat_date_subs$newlabs]
write.FASTA(A_subseqs_new, "./samples/samp_6mo1/A_sub_6mo_1.fasta")
save(A_meta_strat_date_subs, file = "./samples/samp_6mo1/A_sub_6mo_1.rdata")


# offset the uniform sample by 3 months ----

B_meta_strat_date_subs <- B_data[-which.max(B_data$date),] %>%
  group_by(newregion) %>%
  dplyr::arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date) %m+% months(3), "6 months")) %>%
  group_by(dategroup, .add = T)%>%
  slice_sample(n = 50)
B_meta_strat_date_subs <- rbind(B_meta_strat_date_subs, B_data[which.max(B_data$date),])


B_subseqs_new <- algn_VP1_300[B_meta_strat_date_subs$newlabs]
write.FASTA(B_subseqs_new, "./samples/samp_6mo1_off/B_sub_6mo_offset_1.fasta")
save(B_meta_strat_date_subs, file = "./samples/samp_6mo1_off/B_sub_6mo_offset_1.rdata")

range(B_meta_strat_date_subs$date)


A_meta_strat_date_subs <- A_data[-which.max(A_data$date),] %>%
  group_by(newregion) %>%
  dplyr::arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date) %m+% months(3), "6 months")) %>%
  group_by(dategroup, .add = T)%>%
  slice_sample(n = 50)
A_meta_strat_date_subs <- rbind(A_meta_strat_date_subs, A_data[which.max(A_data$date),])

range(A_meta_strat_date_subs$date)

A_subseqs_new <- algn_VP1_300[A_meta_strat_date_subs$newlabs]
write.FASTA(A_subseqs_new, "./samples/samp_6mo1_off/A_sub_6mo_offset_1.fasta")
save(A_meta_strat_date_subs, file = "./samples/samp_6mo1_off/A_sub_6mo_offset_1.rdata")



# Plot the full data after removing ambiguities (Fig. 1 A+B)
region_sub <- allmeta_used[allmeta_used$newregion == "North America" | allmeta_used$newregion == "Europe" ,]

seqs_country <- tabyl(region_sub, country , color)


coords <- data.frame(
  country = c(
    "Belgium", "Canada", "Denmark", "Finland",
    "France", "Germany", "Iceland", "Ireland", "Italy",
    "Netherlands", "Norway", "Russia", "Slovenia", "Spain",
    "Sweden", "Switzerland", "United Kingdom", "USA"
  ),
  lon = c(
    4.47,   -106.35,   9.50,   25.75,
    2.21,    10.45,  -19.02,  -8.00,   12.57,
    5.29,     8.47,  33,   14.99,   -3.75,
    18.64,     8.23,   -2.55,  -98.58
  ),
  lat = c(
    50.50,    56.13,   56.00,   61.92,
    46.23,    51.16,   64.96,   53.14,   41.87,
    52.13,    60.47,   58.5,   46.15,   40.46,
    60.13,    46.82,   54.70,   39.83
  )
)
# Merge coordinates with data
plot_data <- left_join(seqs_country, coords, by = "country")
plot_data$tot = plot_data$`A/D`+plot_data$`B/C`

plot_NA <- ggplot() +
  borders("world", colour = "gray80", fill = "white") +
  geom_scatterpie(data = plot_data,
                  aes(x = lon, y = lat, r = sqrt(tot)*0.18),
                  cols = c("A/D", "B/C")) +
  scale_fill_manual(values = c("A/D" = "#f54b42", "B/C" = "#4287f5")) +
  coord_fixed(xlim = c(-140,-57), ylim = c(25, 70)) +
  theme_void()+
  geom_scatterpie_legend(sqrt(plot_data$tot)*0.2, x=-135, y=30, labeller = function(x) round(x^2/0.18^2, 0), n=4)+
  theme(legend.position = "bottom", panel.background = element_rect(fill = "gray95", color = NA), plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))

shared_legend <- get_legend(plot_NA)
plot_NA <- plot_NA + theme(legend.position = "none")

plot_EU <- ggplot() +
  borders("world", colour = "gray80", fill = "white") +
  geom_scatterpie(data = plot_data,
                  aes(x = lon, y = lat, r = sqrt(tot)*0.18),
                  cols = c("A/D", "B/C"), show.legend =F) +
  scale_fill_manual(values = c("A/D" = "#f54b42", "B/C" = "#4287f5")) +
  coord_fixed(xlim = c(-10,35), ylim = c(25, 70)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "gray95", color = NA), plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))


#plot_grid(plot_NA, plot_EU, rel_widths = c(1.83,1)), shared_legend, rel_heights = c(1,0.1), ncol = 1)

map_included <- plot_grid(plot_NA, plot_EU, rel_widths = c(1.83,1))
ggsave("./plots/sequences_available_map.svg", map_included, dpi = 700, width = 9, height = 4, units = "in")

allmeta_used_region <- allmeta_used %>%
  filter(newregion %in% c("North America", "Europe"))


# sequences through time
samp_all_plot <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = 2014, ymin = -Inf, ymax = Inf, alpha = 0.2) +
  geom_histogram(data = allmeta_used_region, aes(x=date, y = ..count.., fill = color), color = "black", linewidth = 0.2, breaks = brks_off) +
  theme_bw() +
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
  geom_vline(xintercept = 2014)+
  scale_y_continuous(expand = expansion(add = c(0,5)))+
  scale_fill_manual(values = c("A/D" = "#f54b42", "B/C" = "#4287f5")) +
  ylab("Genetic Sequences") +
  xlab("Date") +
  guides(fill = guide_legend(position = "bottom", title = "Clade"))+
  facet_grid(color~region)+
  theme(axis.text=element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14), strip.text = element_text(size=14), axis.title = element_text(size=14), strip.text.y = element_blank(), strip.background.y = element_blank())


zoom_plot <- allmeta_used %>%
  filter(newregion %in% c("North America", "Europe"), color == "A/D")

zoom_a <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = 2014, ymin = -Inf, ymax = Inf, alpha = 0.2) +
  geom_histogram(data = zoom_plot, aes(x=date, y = ..count.., fill = color), color = "black", linewidth = 0.2, breaks = brks_off, show.legend = F) +
  theme_bw() +
  scale_x_continuous(expand = expansion(), breaks = seq(1994, 2024, by = 4), minor_breaks = seq(1994, 2024, by = 1))+
  geom_vline(xintercept = 2014)+
  scale_y_continuous(expand = expansion(add = c(0,5)))+
  scale_fill_manual(values = c("A/D" = "#f54b42")) +
  ylab("Genetic Sequences") +
  xlab("Date") +
  guides(fill = guide_legend(position = "bottom", title = "Clade"))+
  facet_grid(~region)+
  theme(axis.text=element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14), strip.text = element_text(size=14), axis.title = element_text(size=14), strip.text.x = element_blank(), strip.text.y = element_blank(), strip.background.y = element_blank())


ggsave("./plots/sequences_available_bar.svg", samp_all_plot, dpi = 700, width = 9, height = 4, units = "in")
ggsave("./plots/sequences_A_zoom.svg", zoom_a, dpi = 700, width = 9, height = 2, units = "in")

