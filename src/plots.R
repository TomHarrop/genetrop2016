library(ggplot2)
library(data.table)

theme_slide <- theme_grey(base_size = 16, base_family = "Lato") +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

rawData <- data.table(read.csv('data/Phenotype_Panicle_corrected.csv',
                               stringsAsFactors = FALSE))

# remove "Echelle" row
rawData <- rawData[!grep("Echelle", file_name)]

rawData[, Accession := unlist(strsplit(file_name, split = "_", fixed = TRUE))[1], by = file_name]
rawData[, Species := plyr::revalue(toupper(Accession),
                                   c(B88 = "O. barthii", IR64 = "O. sativa indica",
                                     NIP = "O. sativa japonica", TOG5681 = "O. glaberrima",
                                     W1654 = "O. rufipogon")),
        by = Accession]
rawData[Accession == "tog5681", Accession := "Tog5681"]
rawData[, Indiv := unlist(strsplit(file_name, split = "_", fixed = TRUE))[2], by = file_name]
rawData[, Panicle := unlist(strsplit(file_name, split = "_", fixed = TRUE))[3], by = file_name]

setkey(rawData, "Species", "Indiv", "Panicle")
pd <- unique(rawData[, .(Species, Indiv, Panicle, TA_nb, Sp_nb)])
ptypes <- ggplot(pd, aes(x = TA_nb, y = Sp_nb)) +
  theme_slide +
  theme(legend.text = element_text(face = "italic")) +
  geom_smooth(method = lm, se = FALSE, colour = "black", size = 1, alpha = 0.5) +
  geom_point(aes(colour = Species), size = 3, alpha = 0.8) +
  scale_color_brewer(palette = "Set1", guide = guide_legend(title = NULL)) +
  ylab("Number of spikelets") + xlab("Number of secondary branches")


rd2 <- as.data.table(read.table("data/OsOgObOrPTRAPdata.txt",
                                sep = "\t",
                                header = TRUE,
                                stringsAsFactors = FALSE,
                                dec = ","))

rd2[, Species := plyr::revalue(Origin, c(Ob = "O. barthii", Os = "O. sativa",
                                         Og = "O. glaberrima", Or = "O. rufipogon"))]

setkey(rd2, "Species", "Bar_Code", "Sowing_nb", "Repet_nb", "Plant_nb", "Panicle_nb")

pd2 <- unique(rd2[, .(Species, Bar_Code, Sowing_nb, Repet_nb, Plant_nb, Panicle_nb, Sp_nb, TA_nb)])
ptypes2 <- ggplot(pd2, aes(x = TA_nb, y = Sp_nb)) +
  theme_slide +
  theme(legend.text = element_text(face = "italic")) +
  geom_point(aes(fill = Species), colour = NA, size = 3, alpha = 0.6, shape = 21) +
  scale_fill_brewer(palette = "Set1", guide = guide_legend(title = NULL)) +
  geom_smooth(method = lm, se = FALSE, colour = "black", size = 0.5) +
  ylab("Number of spikelets") + xlab("Number of secondary branches")





