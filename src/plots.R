library(ggplot2)
library(data.table)

source("src/column_plots.R")

theme_slide <- theme_grey(base_size = 16, base_family = "Lato") +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

#############
# LOAD DATA #
#############

# LMD lib stats
lmd.quant <- readRDS("data/lmd/quantStats/libStats.Rds")

# LMD expression values
tpm.lmd.wide <- data.table(readRDS("data/lmd/tpm/tpm.Rds"),
                           keep.rownames = TRUE)
tpm.lmd.long <- melt(tpm.lmd.wide, id.vars = "rn",variable.name = "library",
                     value.name = "Transcripts per million")
setkey(tpm.lmd.long, rn, library)

# LMD expression calls
lmd.expression.wide <- data.table(readRDS(
  'data/lmd/expressedGenes/expGenTT.Rds'))
lmd.expression <- melt(lmd.expression.wide, id.vars = "id",
                       variable.name = "library", value.name = "call")
setkey(lmd.expression, id, library)
lmd.tpm <- lmd.expression[tpm.lmd.long]
setnames(lmd.tpm, "id", "gene")

# add stage
lmd.tpm[, stage := substr(library, start = 1, stop = 2)]
old <- c("n1", "n2", "n3", "n4")
new <- c("RM", "PBM", "ePBM/\nAM", "SM")
lmd.tpm[, stage := factor(
  plyr::mapvalues(stage, from = old, to = new), levels = new)]

# gsea data
gsea <- readRDS('data/lmd/gsea/gseaTable.Rds')
famCat <- readRDS("data/tfdb_fam_cat.Rds")

# 5acc expression values
tpm <- readRDS("data/fiveacc/tpm/tpm_with_calls.Rds")
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

# homeobox data
hb <- readRDS("data/lmd/homeobox/plotData.long.Rds")

# 20 accessions branch # vs. spikelet number
rd2 <- as.data.table(read.table("data/OsOgObOrPTRAPdata.txt",
                                sep = "\t",
                                header = TRUE,
                                stringsAsFactors = FALSE,
                                dec = ","))
# 5 accessions diversity
rawData <- data.table(read.csv('data/Phenotype_Panicle_corrected.csv',
                               stringsAsFactors = FALSE))

# 5acc seq stats
stats.5ac <- readRDS("data/fiveacc/mappingStats/starLogs.Rds")

# basic genes
stage.table <- readRDS(
  "data/fiveacc/deseq2/wald_species/stage_results_table.Rds")

# structural variation
species.results <- readRDS(
  "data/species_results.Rds")

# domestication genes
int.results.table <- readRDS(
  "data/fiveacc/deseq2/wald_domestication/results_table.Rds")
dom.continent <- readRDS(
  "data/fiveacc/deseq2/wald_domestication_by_continent/results_table.Rds")

######################
# PANICLE PHENOTYPES #
######################

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
cols.acc <- RColorBrewer::brewer.pal(5, "Set1")[c(2, 1, 5, 4)]
acc.pheno <- ggplot(pd, aes(x = Accession, y = value, fill = Species)) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  theme_slide +
  scale_fill_manual(values = cols.acc, guide = guide_legend(title=NULL)) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.text = element_text(face = "italic"),
        legend.position = "top",
        axis.text.x	= element_text(angle = 45, hjust = 1)) +
  geom_point(size = 2, alpha = 0.5, shape = 21, colour = NA,
             position = position_jitter(width = 0.4))

#############
# lmd stats #
#############

lmd.quant.table <- lmd.quant[, .(
  Library, Sample, `Number of panicles`, `Dissected area (mm²)`, `Yield (ng)`,
  RIN, `Reads (M)`, `Reads in genes (M)`, `Detected genes`, `rRNA reads (M)`
)]
lmd.quant.table[, Sample := factor(
  plyr::mapvalues(Sample, "ePBM/SBM", "ePBM/AM"),
  levels = c("RM", "PBM", "ePBM/AM", "SM"))]

# lmd.pd <- melt(
#   lmd.quant.table,
#   id.vars = c("Library", "Sample", "rRNA reads (M)"))
# lmd.pd <- lmd.pd[variable %in% c("Yield (ng)", "RIN")]
# lmd.ribosomal <- ggplot(lmd.pd,
#        aes(x = value, y = `rRNA reads (M)`)) +
#   facet_grid(. ~ variable, scales = "free_x") +
#   theme_slide + 
#   xlab(NULL) +
#   scale_fill_brewer(palette = "Set1", guide = guide_legend(title = NULL)) +
#   geom_smooth(colour = "black", size = 0.5, se = TRUE, span = 1) +
#   geom_point(aes(fill = Sample), colour = NA, shape = 21, size = 3,
#              alpha = 0.8)

lmd.ribosomal <- ggplot(lmd.quant.table,
                        aes(x = `Yield (ng)`, y = `rRNA reads (M)`/`Reads (M)`)) +
  theme_slide +
  xlab("RNA input (ng)") + ylab("Proportion of rRNA reads") +
  scale_fill_brewer(palette = "Set1", guide = guide_legend(title = NULL)) +
  geom_smooth(colour = "black", size = 0.5, se = TRUE, span = 1) +
  geom_point(aes(fill = Sample), colour = NA, shape = 21, size = 3,
             alpha = 0.8)

############
# clusters #
############

expressionMatrix <- readRDS('data/lmd/mfuzz/expressionMatrix.Rds')
c1 <- readRDS('data/lmd/mfuzz/c1.Rds')
memCutoff <- 0.5

# get clusters and membership
cluster <- data.table(id = names(c1$cluster), Cluster = c1$cluster,
                      Membership = apply(c1$membership, 1, max), key = "id")

# get standardised VST counts
exprs <- data.table(Biobase::exprs(expressionMatrix), keep.rownames = TRUE)
exprs[,id := rn][, rn := NULL]
setkey(exprs, "id")

plotData.wide <- exprs[cluster[Membership > memCutoff]]
plotData.wide[, number := length(id), by = Cluster]
plotData.wide[, label := factor(paste0("Cluster ", Cluster,
                                       "\n(", number, " genes)"))]

# relevel clusters
centres.wide <- data.table(c1$centers)

# re-order the cluster for the plot
centres.wide[, Cluster := paste("Cluster", 1:nrow(centres.wide))]

# find the changes between RM and PBM and PBM and SM for each cluster
centres.wide[, c("n1n2", "n2n4") :=
               list(PBM - RM,
                    SM - PBM)]
# divide these changes into categories
centres.wide[, c("dn1n2", "dn2n4") :=
               list(c('dec', 'small', 'inc')[cut(
                 n1n2, breaks = c(-Inf, -0.5, 0.5, Inf))],
                 c('dec', 'small', 'inc')[cut(
                   n2n4, breaks = c(-Inf, -1, 1, Inf))])]               

# first, show gradual increase / decrease
centres.wide[dn1n2 == dn2n4, cOrder := c(1,2)[order(RM)]]

# next, big changes in n1n2 but small in n2n4
centres.wide[dn2n4 == 'small', cOrder := c(3,4)[order(RM)]]

# small changes in n1n2, then large change
centres.wide[dn1n2 == 'small', cOrder := c(5,6)[order(RM)]]

# complex patterns 
centres.wide[!dn1n2 == dn2n4 & !dn1n2 == "small" & !dn2n4 == "small",
             cOrder := c(7,8)[order(SM)]]

# order any leftovers on RM
if (any(is.na(centres.wide[, cOrder]))) {
  orderMax <- max(centres.wide[,cOrder], na.rm = TRUE)
  centres.wide[is.na(cOrder), cOrder := c((orderMax + 1):nrow(centres.wide))]
}

# relevel the clusters by cOrder
plotData.wide[, label :=
                factor(label, levels = levels(label)[order(centres.wide[,cOrder])])]

# add label to centres.wide
setkey(centres.wide, 'cOrder')
centres.wide[, label := plotData.wide[, levels(label)]]
centres.wide[, label := factor(label, levels = label)]

# make long
plotData <- reshape2::melt(plotData.wide,
                           id.vars = c("id", "Cluster", "Membership", "label", "number"),
                           variable.name = "Stage",
                           value.name = "Scaled read count")

# fix stage label
plotData[, Stage := plyr::mapvalues(Stage, "ePBM.SBM", "ePBM/AM")]

# add centres to plot
centres <- reshape2::melt(centres.wide, id.vars = 'label',
                          measure.vars = c("RM", "PBM", "ePBM.SBM", "SM"),
                          variable.name = 'Stage',
                          value.name = "Scaled read count")
centres[, Stage := plyr::mapvalues(Stage, "ePBM.SBM", "ePBM/AM")]

# set up heatscale
heatscale <- RColorBrewer::brewer.pal(n = 6, name = "YlOrRd")

# main cluster figure
mfuzzClusters <- ggplot(
  plotData, aes(x = Stage, y = `Scaled read count`,
                colour = Membership, group = id)) +
  theme_slide +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.key.size = grid::unit(8, "point")) +
  xlab(NULL) +
  scale_colour_gradientn(colours = heatscale, limits = c(0, 1),
                         breaks = seq(0, 1, 0.2), guide = FALSE) +
  geom_line(alpha = 0.8) +
  geom_line(data = centres, mapping = aes(group = 1),
            colour = "black", alpha = 0.5) +
  facet_wrap("label", ncol = 4)

# number of genes per family
family.clusters.wide <- readRDS("data/lmd/mfuzz/familyClusters.Rds")
family.clusters <- melt(
  family.clusters.wide, id.vars = c("Fam", "n", "nexpr", "total"),
  variable.name = "Cluster", value.name = "Number of genes")
pd <- family.clusters[total > 0]
pd[`Number of genes` > 0, tile_lab := `Number of genes`]
pd[, mnpc := max(`Number of genes`), by = Fam]

# order the y-axis
pd[, Cluster := factor(Cluster, levels = rev(levels(Cluster)))]

# order the x-axis
x.ord.frame <- data.frame(family.clusters.wide[, .(
  Fam, C1, C2, C3, C4, C5, C6, C7, C8)],
  row.names = "Fam")
x.ord.mat <- as.matrix(x.ord.frame)
hc <- hclust(dist(x.ord.mat, method = "minkowski"), method = "ward.D2")
tfs.ord <- rownames(x.ord.frame)[hc$order]
pd[, Fam := factor(Fam, levels = tfs.ord)]

cluster.tfs <- ggplot(pd[`Number of genes` > 0 & mnpc > 1],
                      aes(y = Cluster, fill = `Number of genes`, x = Fam)) +
  theme_slide +
  xlab(NULL) + ylab(NULL) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colours = heatscale, guide = FALSE) +
  geom_raster() +
  geom_text(aes(label = tile_lab))

##############
# ALOG genes #
##############

mads.svp <- data.table(gene = c(
  "LOC_Os02g52340", "LOC_Os03g08754", "LOC_Os06g11330", "LOC_Os02g07030",
  "LOC_Os06g46030", "LOC_Os10g33780"), key = "gene")
mads.svp[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
mads.svp[is.na(symbol), symbol := gene]

# add expression values
pd.alog <- lmd.tpm[mads.svp]
cols <- RColorBrewer::brewer.pal(3, "Set1")[c(2,1)]
alog <- ggplot(pd.alog, aes(x = stage, y = `Transcripts per million`,
                            group = symbol, colour = call)) +
  theme_slide + xlab(NULL) +
  theme(strip.text = element_text(face = "italic")) +
  scale_colour_manual(values = cols, guide = FALSE) +
  facet_wrap(~symbol, nrow = 2, scales = "free_y") +
  geom_smooth(se = FALSE, colour = "grey40", size = 0.5) +
  geom_point(size = 3, position = position_jitter(width = 0.3), alpha = 0.8)

#######
# spl #
#######

spl <- data.table(gene = c(
  "LOC_Os02g04680", "LOC_Os06g49010", "LOC_Os09g31438", "LOC_Os07g32170",
  "LOC_Os08g39890", "LOC_Os01g46984"), key = "gene")

spl[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
spl[is.na(symbol), symbol := gene]

# add expression values
pd.spl <- lmd.tpm[spl]

# order panels by msuId
pd.spl[, symbol := factor(symbol, levels = unique(symbol))]

cols <- RColorBrewer::brewer.pal(3, "Set1")[c(2,1)]
spl.plot <- ggplot(pd.spl, aes(x = stage, y = `Transcripts per million`,
                            group = symbol, colour = call)) +
  theme_slide + xlab(NULL) +
  theme(strip.text = element_text(face = "italic")) +
  scale_colour_manual(values = cols, guide = FALSE) +
  facet_wrap(~symbol, nrow = 2, scales = "free_y") +
  geom_smooth(se = FALSE, colour = "grey40", size = 0.5) +
  geom_point(size = 3, position = position_jitter(width = 0.3), alpha = 0.8)

#################
# TF Enrichment #
#################

# format some labels
setnames(gsea, old = "Test statistic", new = "Test\nstatistic")
gsea[, Stage := plyr::mapvalues(Stage, "ePBM/SBM", "ePBM/\nAM")]

# separate by TF / other proteins
setkey(famCat, 'Family')
setkey(gsea, 'rn')
gsea[, Category := famCat[gsea][,Category]]
gsea[, Category := plyr::mapvalues(Category, from = c("TF", "Other"),
                                   to = c("Transcription factors",
                                          "Other regulators"))]
gsea[, Category := factor(Category, levels = c("Transcription factors",
                                               "Other regulators"))]

# reverse order of y-axis
gsea[, Stage := factor(Stage, levels = rev(levels(Stage)))]

heatscale <- rev(RColorBrewer::brewer.pal(6, "RdBu"))
gsea.plot <- ggplot(gsea, aes(x = rn, y = Stage, fill = `Test\nstatistic`)) +
  theme_slide +
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradientn(colours = heatscale, guide = FALSE) +
  facet_grid(. ~ Category, scales = "free_x", space = "free_x") +
  geom_raster()

############
# homeobox #
############

# fix labels
hb[, Stage := plyr::mapvalues(Stage, "ePBM/SBM", "ePBM\n/AM")]

# reverse y-axis
pd.hb <- hb[grepl("HD-ZIP III", class) | grepl("HD-ZIP IV", class)]
pd.hb[, Stage := factor(Stage, levels = rev(levels(Stage)))]

# re-scale x-axis
pd.hb[, symbol := NULL]
pd.hb[, symbol := oryzr::LocToGeneName(msuId)$symbols, by = msuId]
pd.hb[is.na(symbol), symbol := msuId]
mat.frame <- data.frame(
  dcast(pd.hb, symbol ~ Stage, value.var = "Scaled reads"), row.names = "symbol")
mat <- as.matrix(mat.frame)
hc <- hclust(dist(mat, method = "minkowski"), method = "ward.D2")
ord <- rownames(mat)[hc$order]
pd.hb[, symbol := factor(symbol, levels = rev(ord))]

# plot
heatscale <- RColorBrewer::brewer.pal(6, "YlOrRd")
hdzip <- ggplot(pd.hb, aes(x = symbol, y = Stage, fill = `Scaled reads`)) +
  theme_slide + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_fill_gradientn(colours = heatscale) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  facet_grid(. ~ class, scales = "free_x", space = "free_x") +
  geom_raster()

##############
# 5acc stats #
##############

stats.5ac[, Species := factor(plyr::revalue(
  substr(Library, 1, 1),
  c(I = "O. sativa indica",
    J = "O. sativa japonica",
    R = "O. rufipogon",
    G = "O. glaberrima",
    B = "O. barthii")),
  levels = c("O. rufipogon", "O. sativa indica", "O. sativa japonica",
             "O. barthii", "O. glaberrima"))]

statplot.5ac <- ggplot(
  stats.5ac,
  aes(x = `Number of input reads`/1e6, y = `Number of reads in genes`/1e6,
      group = Species, colour = Species)) +
  theme_slide +
  theme(strip.text.x = element_text(face = "italic")) +
  coord_fixed() +
  xlab("Input reads (M)") + ylab("Reads in genes (M)") +
  facet_wrap(~ Species, nrow = 2) +
  scale_colour_brewer(palette = "Set1", guide = FALSE) +
  geom_hline(aes(yintercept = 20), size = 0.5, linetype = 2, colour = "grey20") +
  geom_smooth(method = "lm", se = FALSE, size = 0.5, alpha = 0.5) +
  geom_point()

###############
# Basic genes #
###############

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

# get the number of sig per contrast
spc <- species.results[padj < 0.05, .(
  n_sig = length(unique(gene))), by = contrast]
setkey(spc, contrast)

# # test vis
# spc[, s1 := unlist(strsplit(contrast, "\\."))[1], by = contrast]
# spc[, s2 := unlist(strsplit(contrast, "\\."))[2], by = contrast]
# 
# hspt <- RColorBrewer::brewer.pal(6, "YlOrRd")
# ggplot(spc, aes(x = s2, y = s1, label = n_sig, size = n_sig, fill = n_sig)) +
#   scale_size(range = c(10,20)) +
#   scale_fill_gradientn(colours = hspt) +
#   geom_point(mapping = aes(y = s1 - 0.5, x = s2), shape = 21) +
#   geom_text(size = 16)
  
pal.struc <- RColorBrewer::brewer.pal(5, "Set1")

# ruf vs. ind
genes <- species.results[contrast == "indica.rufipogon" & !is.na(symbol)]
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
  scale_colour_manual(values = c(pal.struc[1], pal.struc[2]),
                      guide = guide_legend(title = NULL)) +
  ylab("TPM") + xlab(NULL) +
  scale_shape(guide = FALSE) +
  facet_wrap(~ symbol, nrow = 1, scales = "free_y") +
  geom_smooth(method = lm, se = FALSE) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

# ruf vs. sat
genes <- species.results[contrast == "japonica.rufipogon" & !is.na(symbol)]
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
  scale_colour_manual(values = c(pal.struc[1], pal.struc[3]),
                      guide = guide_legend(title = NULL)) +
  ylab("TPM") + xlab(NULL) +
  scale_shape(guide = FALSE) +
  facet_wrap(~ symbol, nrow = 1, scales = "free_y") +
  geom_smooth(method = lm, se = FALSE) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

# ind vs. sat
genes <- species.results[contrast == "japonica.indica" & !is.na(symbol)]
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
  scale_colour_manual(values = c(pal.struc[2], pal.struc[3]),
                      guide = guide_legend(title = NULL)) +
  ylab("TPM") + xlab(NULL) +
  scale_shape(guide = FALSE) +
  facet_wrap(~ symbol, nrow = 1, scales = "free_y") +
  geom_smooth(method = lm, se = FALSE) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

# ind vs. sat
genes <- species.results[contrast == "barthii.glaberrima" & !is.na(symbol)]
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
  scale_colour_manual(values = c(pal.struc[4], pal.struc[5]),
                      guide = guide_legend(title = NULL)) +
  ylab("TPM") + xlab(NULL) +
  scale_shape(guide = FALSE) +
  facet_wrap(~ symbol, nrow = 1, scales = "free_y") +
  geom_smooth(method = lm, se = FALSE) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

#######################
# Domestication genes #
#######################

genes <- int.results.table[padj < 0.05, unique(gene)]
dom.genes <- ColumnPlot(genes) + theme_slide +
  theme(axis.text.y = element_text(face = "italic"),
        strip.text.x = element_text(face = "italic"))

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