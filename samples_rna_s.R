rna.seq.all <- read.table("~/Documents/Pilot 3/gene_read_counts_BIOS_and_LLD_passQC.tsv")
pilot.3.ids <- read.table("~/Documents/Pilot 3/sc_pilot_3.txt")
sample.locations <- read.table("~/Documents/Pilot 3/sample_locations.txt", header = T)
gene.counts.names <- read.table("~/Documents/Pilot 3/gene.counts.tsv.gz.head.txt")[1,]
tpm.combined <- read.table("~/Documents/Pilot 3/tpm_combined_BIOS_all_files_CONTAINS_DUPLICATES.tsv")[1,]

sample.locations[1:5,]


pilot.3.ids


pilot.3.ids

colnames(rna.seq.all) <- as.character(unlist(rna.seq.all[1,]))
rna.seq.all <- rna.seq.all[-1,]

in.rna.seq <- colnames(rna.seq.all) %in% pilot.3.ids[,1]
in.rna.seq[1:10]

rna.seq.pilot3 <- rna.seq.all[,in.rna.seq]
rna.seq.pilot3[,]

pilot.3.ids[pilot.3.ids[,1] %in% colnames(rna.seq.all),1]

pilot.3.ids[,1] %in% sample.locations$sample_name

# Sample name and file id of all samples
sample.locations[sample.locations$file_id %in% pilot.3.ids[,1],]

# These sample names differ from file_id &
sample.names <- sample.locations[sample.locations$file_id %in% pilot.3.ids[,1] & !sample.locations$sample_name %in% pilot.3.ids[,1],]$sample_name
# None in expression matrix:(
pilot.3.ids[sample.names %in% colnames(rna.seq.all),1]

pilot.3.ids[,1] %in% gene.counts.names
pilot.3.ids[,1] %in% tpm.combined

sample.locations[sample.locations$sample_name %in% pilot.3.ids[,1],]

colnames(rna.seq.all)
dim(rna.seq.all)

grep("lldeep", colnames(rna.seq.all), value = T)

colnames(rna.seq.all)
