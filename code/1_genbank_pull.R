## genbank D68 - updated nextstrain data
## VP1 300 - any sequences with 300 bases of VP1 - don't need to be complete

## API pull
# BiocManager::install("msa")
devtools::install_github("https://github.com/lpmor22/read.gb") # adapted version with a fix

pacman::p_load(rentrez,
               ape,
               xml2,
               read.gb,
               numbers,
               tidyverse,
               zoo,
               janitor,
               anytime,
               countrycode)


## API pull
reference_sequence <- read.table("./data/reference.txt", sep = "\n", header = F)
# search for EV-D68 accessions
newseqs <- entrez_search(db = "nucleotide", term = "txid42789[organism:exp] AND 1970:2024[PDAT]", retmax = 10000)

## need to break requests into chunks
split_into_chunks <- function(vec, chunk_size) {
  split(vec, ceiling(seq_along(vec) / chunk_size))
}
uid_chunks <- split_into_chunks(newseqs$ids, 200)


# pull accession numbers
accession_numbers <- data.frame(V1 = unname(unlist(
  lapply(uid_chunks, function(chunk) {
    summaries <- entrez_summary(db="nucleotide", id=chunk)
    sapply(summaries, function(x) sub("\\..*$", "", x$accessionversion))  # remove version
  })
)))
print(accession_numbers)

# add my reference sequence
new_genbank <- rbind(accession_numbers, reference_sequence)

# extract accessions from both
accns_to_pull <- new_genbank$V1
# split into chunks to pull
accn_chunks <- split_into_chunks(accns_to_pull, 200)

# Function to fetch and write one chunk
fetch_and_write_chunk <- function(chunk) {
  whist <- entrez_post(db = "nuccore", id = chunk)
  recs <- entrez_fetch(db = "nuccore", web_history = whist, rettype = "fasta", retmax = length(chunk))
  cat(recs, file = "./data/D68_new_VP1.fasta", append = TRUE)
  invisible(NULL) 
  # Sys.sleep(0.5)  # Optional pause
}

# pull the sequences
lapply(accn_chunks, fetch_and_write_chunk)

# pull metadata
genbank_chunks <- split_into_chunks(new_genbank$V1, 200)

lapply(seq_along(genbank_chunks), function(i) {
  chunk <- genbank_chunks[[i]]
  
  whist <- entrez_post(db = "nuccore", id = chunk)
  tryCatch({
    recs <- entrez_fetch(db = "nuccore", web_history = whist, rettype = "gb", retmode = "text", retmax = length(chunk))
    if (nchar(recs) > 0) {
      cat(recs, file = "./data/D68_new_VP1.gb", append = TRUE)
    }
  }, error = function(e) {
    cat("Failed on chunk", i, ":", conditionMessage(e), "\n")
  })
})


meta_newnew <- read.gb("./data/D68_new_VP1.gb", Type = "full")

# Some manual fixes were needed:
# - UTR annotations - changed 5' to prime5 and 3' to prime3
# - ampersand in some address fields
# - need to change geo_loc_name to country

# Pull this metadata into R
feature_list <- extract.gb(meta_newnew, "source")
transposeList <- lapply(feature_list,function(x){
  t(x$source) %>%
    data.frame %>%
    row_to_names(row_number=1)
})
meta_genbank <- do.call(plyr::rbind.fill, transposeList) %>%
  add_column(accn = names(transposeList), .before = 1)



# align them all with mafft
## "/opt/homebrew/bin/mafft" --thread -1 --retree 2 --reorder "./data/D68_new_vp1.fasta" > "./data/D68_new_VP1_algn.fasta"

# trim them manually to a known VP1 with aliview

# treetime needs tip dates - rename sequences to just their accession and then pull the dates from the metafile
algn_VP1_300 <- read.FASTA("./data/D68_new_VP1_algn.fasta")
names(algn_VP1_300) <- sapply(strsplit(names(algn_VP1_300), "\\."), `[`, 1)

# remove short sequences
algn_VP1_300 <- algn_VP1_300[lengths(del.gaps(algn_VP1_300))>300]
# remove reference sequence - makes a nicer tree
algn_VP1_300 <- algn_VP1_300[names(algn_VP1_300)!="AY426531"]

# file with decimal dates - use ranges where not complete
meta_genbank$date <- round(decimal_date(anytime::anydate(meta_genbank$collection_date)),3)
# leaves NAs where only month is given, covert these to a range of values for treetime to use
#### THIS IS NOT A COMPLETE CATCH-ALL AS COULD HAVE ONLY YEARS REPORTED IN THE FUTURE ####
meta_genbank[is.na(meta_genbank$date),]$date <- paste0("[",
                                                       round(decimal_date(as.Date(paste("01-",meta_genbank[is.na(meta_genbank$date),]$collection_date, sep=""), format = "%d-%b-%Y")),3),
                                                       ":",
                                                       round(decimal_date(ceiling_date(as.Date(paste("01-",meta_genbank[is.na(meta_genbank$date),]$collection_date, sep=""), format = "%d-%b-%Y"), "month") -1),3),
                                                       "]"
                                                       )

meta_genbank$country <- sapply(strsplit(meta_genbank$country, "\\:"), `[`, 1)

# tipdate file for treetime
tipdates <- meta_genbank %>%
  select(accn, date) %>% # select correct columns
  dplyr::rename(accession = accn, 
         dec.date = date
         )

# drop the ones without a VP1??

write.csv(tipdates, "./data/matched_data/all_tipdates.csv", row.names = F, quote = F)

# check for sequences not in the metafile
table(names(algn_VP1_300) %in% tipdates$accession)
names(algn_VP1_300[-which(names(algn_VP1_300) %in% tipdates$accession)])

write.FASTA(algn_VP1_300, "./data/matched_data/algn_VP1_300.fasta")

# run iqtree from matched data file
# iqtree -s algn_VP1_300.fasta -m HKY -nt AUTO
# run treetime with this tree (wont auto detect iqtree for some reason) (treetime is in a conda env - source activate treetime)
# treetime --aln algn_VP1_300.fasta --dates all_tipdates.csv --clock-filter 0 --stochastic-resolve --covariation  --time-marginal 'false' --confidence --tree algn_VP1_300.fasta.treefile --outdir treetime
# pull in the dates from treetime

tt_dates <- read.table("./data/matched_data/treetime/dates.tsv")

# want V3 and trim to the names included in the FASTA file for BEAST tree dating
avg_dates <- left_join(data.frame(names(algn_VP1_300)), tt_dates[,c("V1", "V3")], by = c("names.algn_VP1_300." = "V1"))
write.table(avg_dates, file = "./data/matched_data/beast_dates_jun25.txt", row.names = F, col.names = F, quote = F, sep = "\t")

# resolve treetime tree with ape to use as starting tree
dated_tree <- read.tree("./data/matched_data/treetime/timetree.nexus")
# pull the tree object from this
starting_tree <- multi2di(dated_tree[[6]])
write.nexus(starting_tree, file = "./data/matched_data/treetime/random_resolve.nexus")



# downsample --------------------------------------------------------------
# need sample dates and locations to take a time-and-location-stratified sample

metadata <- meta_genbank %>%
  select(accn, date, country) %>% # select correct columns
  rename(accession = accn, 
         dec.date = date,
         country = country
  ) %>% # name to match the other meta file
  bind_rows(select(with_acsn, accession, dec.date, country),.) %>% # join basic metadata from genbank and nextstrain
  left_join(tt_dates[,c("V1", "V3")], by = c("accession" = "V1")) %>% # add the dates inferred in treetime
  drop_na() %>% # remove rows which are not in the treetime tree
  mutate(region = countrycode(country,"country.name", "region23"),
         region = recode(region, `Northern Europe` = "Europe", `Southern Europe` = "Europe", `Western Europe` = "Europe", `Northern America` = "North America")) # Combine Europe

write.csv(metadata, "./data/matched_data/metafile.csv")

# meta table
tabyl(metadata$region) %>%  adorn_totals("row") %>% adorn_pct_formatting()

# Clades are manually added in figtree aligning to the clades provided on nextstrain






