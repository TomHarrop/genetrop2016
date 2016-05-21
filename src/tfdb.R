#!/usr/bin/env Rscript

# download to tempfile
tmp <- tempfile(fileext = ".tab")
call <- paste(c('curl --data-urlencode "_submit_check=1"',
       '--data-urlencode "DOWNLOAD=list of transcription factors"',
       '--data-urlencode "SUBMIT=Download"',
       'http://plntfdb.bio.uni-potsdam.de/v3.0/export.php > ', tmp),
       collapse = " ")
system(call)

# read into R
tfdbRaw <- data.table(read.table(tmp, header = T, sep = '\t',
                                 stringsAsFactors = F))

# sort by Protein.ID and Family
setkey(tfdbRaw, 'Species', 'Family', 'Protein.ID')

# remove duplicates
tfdbRaw <- unique(tfdbRaw)

# extract family and category
famCat <- unique(tfdbRaw["Oryza sativa subsp. japonica",.(Family, Category)])

# save
saveRDS(famCat, "data/tfdb_fam_cat.Rds")

quit(save = "no", status = 0)
