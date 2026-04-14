# Inject tip date distributions and operators into BEAST XML.
# METADATA should contain the input ranges and the tip names.
pacman::p_load(
  dplyr,
  tidyr,
  glue
)

# pull names of the tips which are variable
load("./samples/samp_6mo2_off/A_sub_6mo_offset_2.rdata")
load("./samples/samp_6mo2_off/B_sub_6mo_offset_2.rdata")

# dec.date has the range
A_ambig <- A_meta_strat_date_subs %>%
  # search for the ambiguous date format
  filter(grepl("\\[", dec.date)) %>%
  # split this into two columns
  mutate(year_range_clean = gsub("\\[|\\]", "", dec.date)) %>%
  separate(year_range_clean, into = c("year_start", "year_end"), sep = ":", convert = TRUE)

# write the prior XML tags
clade_A_priors <- glue('
<distribution id="{A_ambig$nodelabel}.prior" spec="beast.base.evolution.tree.MRCAPrior" monophyletic="true" tipsonly="true" tree="@Tree.t:D68_A">
  <taxonset id="{A_ambig$nodelabel}" spec="TaxonSet">
    <taxon id="{A_ambig$newlabs}" spec="Taxon"/>
  </taxonset>
  <AlmostUniform id="Uniform.{1:nrow(A_ambig)}" name="distr" lower="{A_ambig$year_start}" upper="{A_ambig$year_end}"/>
</distribution>
')

B_ambig <- B_meta_strat_date_subs %>%
  # search for the ambiguous date format
  filter(grepl("\\[", dec.date)) %>%
  # split this into two columns
  mutate(year_range_clean = gsub("\\[|\\]", "", dec.date)) %>%
  separate(year_range_clean, into = c("year_start", "year_end"), sep = ":", convert = TRUE)

# write the prior XML tags   
clade_B_priors <- glue('
<distribution id="{B_ambig$nodelabel}.prior" spec="beast.base.evolution.tree.MRCAPrior" monophyletic="true" tipsonly="true" tree="@Tree.t:D68_B">
  <taxonset id="{B_ambig$nodelabel}" spec="TaxonSet">
    <taxon id="{B_ambig$newlabs}" spec="Taxon"/>
  </taxonset>
  <AlmostUniform id="Uniform.{(1:+nrow(B_ambig))+nrow(A_ambig)}" name="distr" lower="{B_ambig$year_start}" upper="{B_ambig$year_end}"/>
</distribution>
')
priors <- c(clade_A_priors, clade_B_priors)

# write the operator XML tags
A_ops <- glue('<operator id="TipDatesRandomWalker.{A_ambig$nodelabel}.prior" spec="TipDatesRandomWalker" taxonset="@{A_ambig$nodelabel}" tree="@Tree.t:D68_A" weight="0.1" windowSize="1.0"/>')
B_ops <- glue('<operator id="TipDatesRandomWalker.{B_ambig$nodelabel}.prior" spec="TipDatesRandomWalker" taxonset="@{B_ambig$nodelabel}" tree="@Tree.t:D68_B" weight="0.1" windowSize="1.0"/>')
ops <- c(A_ops, B_ops)

# write the log XML tags
A_log <- glue('<log idref="{A_ambig$nodelabel}.prior"/>')
B_log <- glue('<log idref="{B_ambig$nodelabel}.prior"/>')
logs <- c(A_log, B_log)


# Inject these into our XML template
xml0 = readLines("./samples/samp_6mo2_off/template_update_trajlog.xml")
xml1 = gsub(xml0, pattern = 'INSERT_TIP_DISTN', replacement = paste0(priors, collapse = "\n"))
xml2 = gsub(xml1, pattern = 'INSERT_TIP_OPERATORS', replacement = paste0(ops, collapse = "\n"))
xml3 = gsub(xml2, pattern = 'INSERT_LOG', replacement = paste0(logs, collapse = "\n"))

writeLines(xml3,"./samples/samp_6mo2_off/inject_dates_trajlog.xml")

# Trial it in beast


