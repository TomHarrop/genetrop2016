library(ggplot2)
library(data.table)

source("src/column_plots.R")

theme_slide <- theme_grey(base_size = 16, base_family = "Lato") +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

######################
# PANICLE PHENOTYPES #
######################

# 20 accessions branch # vs. spikelet number
rd2 <- as.data.table(read.table("data/OsOgObOrPTRAPdata.txt",
                                sep = "\t",
                                header = TRUE,
                                stringsAsFactors = FALSE,
                                dec = ","))

# set legend labels
rd2[, Species := plyr::revalue(
  Origin, c(Ob = "O. barthii", Os = "O. sativa",
            Og = "O. glaberrima", Or = "O. rufipogon"))]

# add number of accessions per species
rd2[, n_acc := length(unique(Bar_Code)), by = "Species"]
setkey(rd2, "Species", "Bar_Code", "Sowing_nb", "Repet_nb", "Plant_nb",
       "Panicle_nb")

# select data to plot
plot.data <- unique(rd2[, .(
  Species, Bar_Code, Sowing_nb, Repet_nb, Plant_nb, Panicle_nb, Sp_nb, TA_nb,
  n_acc)])

# add accession # to label
plot.data[, guide_label := paste0(Species, " (", n_acc, ")")]

# draw the plot
ptypes2 <- ggplot(plot.data, aes(x = TA_nb, y = Sp_nb)) +
  theme_slide +
  theme(legend.text = element_text(face = "italic")) +
  geom_point(aes(fill = guide_label), colour = NA, size = 3, alpha = 0.5,
             shape = 21) +
  scale_fill_brewer(palette = "Set1", guide = guide_legend(title = NULL)) +
  geom_smooth(method = lm, se = FALSE, colour = "black", size = 0.5) +
  ylab("Number of spikelets") + xlab("Number of secondary branches")

# 5 accessions diversity
rawData <- data.table(read.csv('data/Phenotype_Panicle_corrected.csv',
                               stringsAsFactors = FALSE))

# remove "Echelle" row
rawData <- rawData[!grep("Echelle", file_name)]

rawData[, Accession := unlist(strsplit(
  file_name, split = "_", fixed = TRUE))[1], by = file_name]
rawData[, Species := plyr::revalue(
  toupper(Accession),
  c(B88 = "O. barthii", IR64 = "O. sativa",
    NIP = "O. sativa", TOG5681 = "O. glaberrima",
    W1654 = "O. rufipogon"))]
rawData[Accession == "tog5681", Accession := "Tog5681"]
rawData[Accession == "Nip", Accession := "Nipponbare"]
rawData[, Indiv := unlist(
  strsplit(file_name, split = "_", fixed = TRUE))[2], by = file_name]
rawData[, Panicle := unlist(
  strsplit(file_name, split = "_", fixed = TRUE))[3], by = file_name]

# select data
setkey(rawData, "Accession", "Indiv", "Panicle")
pd.wide <- unique(rawData[, .(
  Species, Accession, Indiv, Panicle,
  "Spikelets" = Sp_nb,
  "Primary branches" = SA_nb,
  "PB length (mm)" = Total.SA.length,
  "Rachis length (mm)" = PA_length * 10,
  "Secondary branches" = TA_nb,
  "SB length (mm)" = Total.TA.length)])
pd <- melt(pd.wide, id.vars = c("Species", "Accession", "Indiv", "Panicle"))

# order accessions and species
pd[, Accession := factor(
  Accession,
  levels = c("IR64", "Nipponbare", "W1654", "Tog5681", "B88"))]
pd[, Species := factor(
  Species,
  levels = c("O. sativa", "O. rufipogon", "O. glaberrima", "O. barthii"))]

# plot phenotypes
acc.pheno <- ggplot(pd, aes(x = Accession, y = value, colour = Species)) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  theme_slide +
  scale_colour_brewer(palette = "Set1", guide = guide_legend(title=NULL)) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.text = element_text(face = "italic"),
        legend.position = "top",
        axis.text.x	= element_text(angle = 45, hjust = 1)) +
  geom_point(size = 2, alpha = 0.5, position = position_jitter(width = 0.4))

###############
# Basic genes #
###############

tpm <- readRDS("data/fiveacc/tpm/tpm_with_calls.Rds")
stage.table <- readRDS(
  "data/fiveacc/deseq2/wald_species/stage_results_table.Rds")

# lots of "sig" genes
n.sig.stage <- stage.table[padj < 0.05 & abs(log2FoldChange) > log(1.5, 2),
                           length(unique(gene))]

# just highlight a few genes
setorder(stage.table, padj, na.last = TRUE)
stage.genes <- stage.table[1:10, gene]

# plot
plot.data <- tpm[stage.genes]
plot.data[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
plot.data[is.na(symbol), symbol := gene]
plot.data[, Species := factor(plyr::revalue(
  substr(sample, 1, 1),
  c(I = "O. sativa indica",
    J = "O. sativa japonica",
    R = "O. rufipogon",
    G = "O. glaberrima",
    B = "O. barthii")),
  levels = c("O. rufipogon", "O. sativa indica", "O. sativa japonica",
             "O. barthii", "O. glaberrima"))]

stage.genes <- ggplot(
  plot.data,
  aes(x = stage, y = tpm, colour = Species, group = Species, shape = call)) +
  facet_wrap(~symbol, nrow = 2, scales= "free_y") +
  scale_colour_brewer(palette = "Set1", guide = guide_legend(title=NULL)) +
  scale_shape(guide = FALSE) +
  theme_slide + xlab(NULL) +
  ylab("Transcripts per million") +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10, face = "italic"),
        legend.text = element_text(face = "italic")) +
  geom_smooth(method = lm, se = FALSE) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

########################
# Structural variation #
########################

species.results <- readRDS(
  "data/fiveacc/deseq2/wald_species/contrast_results.Rds")
ReturnFourGenes <- function(species.a, species.b,
                            species.results = species.results) {
  result.table <- species.results[
    grepl(species.a, contrast) & grepl(species.b, contrast)]
  setorder(result.table, padj, na.last = TRUE)
  result.table[1:4]
}

# get the number of sig per contrast
spc <- species.results[padj < 0.05, .(
  n_sig = length(unique(gene))), by = contrast]
setkey(spc, contrast)

# fix the tpm data
tpm.clean <- copy(tpm)
tpm.clean[, Species := factor(plyr::revalue(
  species,
  c(I = "O. sativa indica",
    J = "O. sativa japonica",
    R = "O. rufipogon",
    G = "O. glaberrima",
    B = "O. barthii")),
  levels = c("O. rufipogon", "O. sativa indica", "O. sativa japonica",
             "O. barthii", "O. glaberrima"))]

species.results[, unique(contrast)]
pal <- RColorBrewer::brewer.pal(5, "Set1")

# ruf vs. ind
genes <- species.results[contrast == "indica.rufipogon"]
setorder(genes, padj, na.last = TRUE)
gene.names <- genes[1:4, gene]
pd <- tpm.clean[gene %in% gene.names & species %in% c("I", "R")]
pd[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
pd[is.na(symbol), symbol := gene]
struct1 <- ggplot(pd, aes(
  x = stage, y = tpm, colour = Species, shape = call, group = Species)) +
  theme_slide + theme(
    strip.text.x = element_text(size = 10, face = "italic"),
    legend.text = element_text(face = "italic")) +
  scale_colour_manual(values = c(pal[1], pal[2]),
                      guide = guide_legend(title = NULL)) +
  ylab("TPM") + xlab(NULL) +
  scale_shape(guide = FALSE) +
  facet_grid(. ~ symbol, drop = TRUE, scales = "free_y") +
  geom_smooth(method = lm, se = FALSE) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

# ruf vs. sat
genes <- species.results[contrast == "japonica.rufipogon"]
setorder(genes, padj, na.last = TRUE)
gene.names <- genes[1:4, gene]
pd <- tpm.clean[gene %in% gene.names & species %in% c("J", "R")]
pd[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
pd[is.na(symbol), symbol := gene]
struct2 <- ggplot(pd, aes(
  x = stage, y = tpm, colour = Species, shape = call, group = Species)) +
  theme_slide + theme(
    strip.text.x = element_text(size = 10, face = "italic"),
    legend.text = element_text(face = "italic")) +
  scale_colour_manual(values = c(pal[1], pal[3]),
                      guide = guide_legend(title = NULL)) +
  ylab("TPM") + xlab(NULL) +
  scale_shape(guide = FALSE) +
  facet_grid(. ~ symbol, drop = TRUE, scales = "free_y") +
  geom_smooth(method = lm, se = FALSE) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

# ind vs. sat
genes <- species.results[contrast == "japonica.indica"]
setorder(genes, padj, na.last = TRUE)
gene.names <- genes[1:4, gene]
pd <- tpm.clean[gene %in% gene.names & species %in% c("J", "I")]
pd[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
pd[is.na(symbol), symbol := gene]
struct3 <- ggplot(pd, aes(
  x = stage, y = tpm, colour = Species, shape = call, group = Species)) +
  theme_slide + theme(
    strip.text.x = element_text(size = 10, face = "italic"),
    legend.text = element_text(face = "italic")) +
  scale_colour_manual(values = c(pal[2], pal[3]),
                      guide = guide_legend(title = NULL)) +
  ylab("TPM") + xlab(NULL) +
  scale_shape(guide = FALSE) +
  facet_grid(. ~ symbol, drop = TRUE, scales = "free_y") +
  geom_smooth(method = lm, se = FALSE) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

# ind vs. sat
genes <- species.results[contrast == "barthii.glaberrima"]
setorder(genes, padj, na.last = TRUE)
gene.names <- genes[1:4, gene]
pd <- tpm.clean[gene %in% gene.names & species %in% c("B", "G")]
pd[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
pd[is.na(symbol), symbol := gene]
struct4 <- ggplot(pd, aes(
  x = stage, y = tpm, colour = Species, shape = call, group = Species)) +
  theme_slide + theme(
    strip.text.x = element_text(size = 10, face = "italic"),
    legend.text = element_text(face = "italic")) +
  scale_colour_manual(values = c(pal[4], pal[5]),
                      guide = guide_legend(title = NULL)) +
  ylab("TPM") + xlab(NULL) +
  scale_shape(guide = FALSE) +
  facet_grid(. ~ symbol, drop = TRUE, scales = "free_y") +
  geom_smooth(method = lm, se = FALSE) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

#######################
# Domestication genes #
#######################

int.results.table <- readRDS(
  "data/fiveacc/deseq2/wald_domestication/results_table.Rds")
genes <- int.results.table[padj < 0.05, unique(gene)]
dom.genes <- ColumnPlot(genes) + theme_slide +
  theme(axis.text.y = element_text(face = "italic"),
        strip.text.x = element_text(face = "italic"))

dom.continent <- readRDS(
  "data/fiveacc/deseq2/wald_domestication_by_continent/results_table.Rds")

# for africa, take the top 20 in each direction
af.table <- dom.continent[padj < 0.05 & domestication == "africa"]
af.table <- af.table[order(-abs(log2FoldChange))]
af1 <- af.table[log2FoldChange > 0][1:20, unique(gene)]
af2 <- af.table[log2FoldChange < 0][1:20, unique(gene)]

# for asia, we have room to show all
as <- dom.continent[padj < 0.05 & domestication == "asia", unique(gene)]

# draw plots
as.genes <- ColumnPlot(as, species = "asian") + theme_slide +
  theme(axis.text.y = element_text(face = "italic"),
        strip.text.x = element_text(face = "italic"))
af1.genes <- ColumnPlot(af1, species = "african") + theme_slide +
  theme(axis.text.y = element_text(face = "italic"),
        strip.text.x = element_text(face = "italic"))
af2.genes <- ColumnPlot(af2, species = "african") + theme_slide +
  theme(axis.text.y = element_text(face = "italic"),
        strip.text.x = element_text(face = "italic"))