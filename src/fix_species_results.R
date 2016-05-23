
library(data.table)

species.results <- readRDS(
  "data/fiveacc/deseq2/wald_species/contrast_results.Rds")

species.results[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]

# save
saveRDS(species.results, "data/species_results.Rds")

quit(save = "no", status = 0)
