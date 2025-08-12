pacman::p_load(
  ape,
  treeio,
  ggtree,
  ggplot2,
  lubridate,
  phytools,
  plyr,
  dplyr
)

# load sequences - not included on GitHub but can be downloaded via the included accession numbers
algn_VP1_300 <- algn <- read.FASTA("./data/algn_VP1_300.fasta", type = "DNA")

# load metadata
metadata <- read.csv("./data/metafile.csv")

# load tree - manually colored by superclade in figtree
tr_ann <- read.beast("./data/treetime/treesplit_analysis.tree")

# assign a group to each sequence (in metadata file)
treemetafile <- tr_ann@data[is.na(tr_ann@data$label),]
treemetafile$nodelabel <- tr_ann@phylo$tip.label
treemetafile$`!color` <- as.factor(treemetafile$`!color`) 
# these colors depend on what you assign in figtree
treemetafile$color <- revalue(treemetafile$`!color`, c("#e62320" = "B/C", "#0582e8" = "A/D"))
treemetafile$node <- as.numeric(treemetafile$node)

# order to match the alignment and output two alignments for each super-clade
allmeta <- left_join(treemetafile, metadata, by = c("nodelabel"="accession"))
allmeta$newlabs <- paste(allmeta$nodelabel, gsub(" ","-", allmeta$region), allmeta$color, gsub(" ","-", allmeta$country), allmeta$V3, sep = "_")
allmeta$date <- as.numeric(allmeta$date)

names(algn_VP1_300) <- allmeta[match(names(algn_VP1_300), allmeta$nodelabel),]$newlabs
# more descriptive sequence names to run in BEAST later
write.FASTA(algn_VP1_300, "./allseq_labels/allseq.fasta")
save(allmeta, file = "./allseq_labels/allseq.rdata")

# subset the two regions and two clades - make fasta and metafile for each
B_data <- allmeta %>%
  filter((region == "Europe" | region == "North America") & color == "B/C")

A_data <- allmeta %>%
  filter((region == "Europe" | region == "North America") & color == "A/D")

# proportional subsampling approach (run twice) ---- 

B_meta_strat_date_subs <- B_data[-which.max(B_data$date),] %>%
  group_by(region) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date), "3 months")) %>% 
  group_by(dategroup, .add = T) %>%
  group_modify(~ .x %>% slice_sample(n = max(1, ceiling(0.45 * nrow(.x))))) # 45% of available sequences maximum

# keep most recent sequence always - helps to align the dates
B_meta_strat_date_subs <- rbind(B_meta_strat_date_subs, B_data[which.max(B_data$date),]) 

B_subseqs_new <- algn_VP1_300[B_meta_strat_date_subs$newlabs]
# write sequences
write.FASTA(B_subseqs_new, "./samples/samp_3mo1/B_sub_3mo_1.fasta")
# write metadata
save(B_meta_strat_date_subs, file = "./samples/samp_3mo1/B_sub_3mo_1.rdata")

# check the included dates
range(B_meta_strat_date_subs$date) 

A_meta_strat_date_subs <- A_data[-which.max(A_data$date),] %>%
  group_by(region) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date), "3 months")) %>%
  group_by(dategroup, .add = T)%>%
  group_modify(~ .x %>% slice_sample(n = max(1, ceiling(0.45 * nrow(.x)))))

A_meta_strat_date_subs <- rbind(A_meta_strat_date_subs, A_data[which.max(A_data$date),])

A_subseqs_new <- algn_VP1_300[A_meta_strat_date_subs$newlabs]
write.FASTA(A_subseqs_new, "./samples/samp_3mo1/A_sub_3mo_1.fasta")
save(A_meta_strat_date_subs, file = "./samples/samp_3mo1/A_sub_3mo_1.rdata")

range(A_meta_strat_date_subs$date)

# Uniform subsampling approach ----

B_meta_strat_date_subs <- B_data[-which.max(B_data$date),] %>%
  group_by(region) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date), "6 months")) %>%
  group_by(dategroup, .add = T)%>%
  slice_sample(n = 50)
B_meta_strat_date_subs <- rbind(B_meta_strat_date_subs, B_data[which.max(B_data$date),])


B_subseqs_new <- algn_VP1_300[B_meta_strat_date_subs$newlabs]
write.FASTA(B_subseqs_new, "./samples/samp_6mo1/B_sub_6mo_1.fasta")
save(B_meta_strat_date_subs, file = "./samples/samp_6mo1/B_sub_6mo_1.rdata")


A_meta_strat_date_subs <- A_data[-which.max(A_data$date),] %>%
  group_by(region) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date), "6 months")) %>%
  group_by(dategroup, .add = T)%>%
  slice_sample(n = 50)
A_meta_strat_date_subs <- rbind(A_meta_strat_date_subs, A_data[which.max(A_data$date),])


A_subseqs_new <- algn_VP1_300[A_meta_strat_date_subs$newlabs]
write.FASTA(A_subseqs_new, "./samples/samp_6mo1/A_sub_6mo_1.fasta")
save(A_meta_strat_date_subs, file = "./samples/samp_6mo1/A_sub_6mo_1.rdata")

# first and last sample of each subset
range(A_meta_strat_date_subs$date)
# 1999.15 2021.98
range(B_meta_strat_date_subs$date)
# 1996.69 2022.94


# offset the uniform sample by 3 months ----

B_meta_strat_date_subs <- B_data[-which.max(B_data$date),] %>%
  group_by(region) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date) %m+% months(3), "6 months")) %>%
  group_by(dategroup, .add = T)%>%
  slice_sample(n = 50)
B_meta_strat_date_subs <- rbind(B_meta_strat_date_subs, B_data[which.max(B_data$date),])


B_subseqs_new <- algn_VP1_300[B_meta_strat_date_subs$newlabs]
write.FASTA(B_subseqs_new, "./samples/samp_6mo1_off/B_sub_6mo_offset_1.fasta")
save(B_meta_strat_date_subs, file = "./samples/samp_6mo1_off/B_sub_6mo_offset_1.rdata")


A_meta_strat_date_subs <- A_data[-which.max(A_data$date),] %>%
  group_by(region) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(dategroup = lubridate::floor_date(date_decimal(date) %m+% months(3), "6 months")) %>%
  group_by(dategroup, .add = T)%>%
  slice_sample(n = 50)
A_meta_strat_date_subs <- rbind(A_meta_strat_date_subs, A_data[which.max(A_data$date),])


A_subseqs_new <- algn_VP1_300[A_meta_strat_date_subs$newlabs]
write.FASTA(A_subseqs_new, "./samples/samp_6mo1_off/A_sub_6mo_offset_1.fasta")
save(A_meta_strat_date_subs, file = "./samples/samp_6mo1_off/A_sub_6mo_offset_1.rdata")


