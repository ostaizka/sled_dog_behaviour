
### R
library(dada2)
packageVersion("dada2") # 1.22 when this was put together

setwd("/raw_data")

samples <- scan("samples", what="character")
samples2 <- scan("samples2", what="character")

samples %in% samples2

#count.table <- read.table(paste(projectdir,"DADA2/ASVs_counts_",code,".tsv",sep=""),header=TRUE,row.names=1)


########Pool1##########
setwd("/Pool1")
samples <- scan("samples", what="character")
## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier(perhaps it needs to be moved)
list.files()

raw_path <- file.path("Users", "ostaizkaaizpurua", "Projects","GreenDogs","raw_data", fsep="/")
trimmed_path1 <- file.path("Users", "ostaizkaaizpurua", "Projects","GreenDogs","Pool1", fsep="/")
trimmed_path2 <- file.path("Users", "ostaizkaaizpurua", "Projects","GreenDogs","Pool2", fsep="/")

codeDADA1="P1"
projectdir="/Users/ostaizkaaizpurua/Projects/GreenDogs/DADA2/Seq1/"

# one holding the file names of all the forward reads and one with the reverse
plotQualityProfile(reverse_reads[3])

forward_reads <- paste0(samples, "_1_trimmed.fq")
reverse_reads <- paste0(samples2, "_2_trimmed.fq")
filtered_forward_reads <- paste0(samples, "_1_filtered.fq")
filtered_reverse_reads <- paste0(samples2, "_2_filtered.fq")

forward_reads <- file.path(trimmed_path1, paste0(samples, "_1_trimmed.fq"))
reverse_reads <- file.path(trimmed_path1, paste0(samples, "_2_trimmed.fq"))
filtered_forward_reads <- file.path(trimmed_path1, paste0(samples, "_1_filtered.fq"))
filtered_reverse_reads <- file.path(trimmed_path1, paste0(samples, "_2_filtered.fq"))

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=c(240,220), compress=TRUE, multithread=TRUE)

#when some input samples had no read that pass the filter
#terminal: ls *_1_filtered.fq | sed "s/_1_filtered\.fq//g" > samples_fil
samples_fil <- scan("samples_fil", what="character")
filtered_forward_reads <- paste0(samples_fil, "_1_filtered.fq")
filtered_reverse_reads <- paste0(samples_fil, "_2_filtered.fq")

#Generating an error model of our data by learning the specific error-signature of our dataset
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)
#Dereplication
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples_fil
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples_fil
saveRDS(derep_forward, paste(projectdir,"derep_forward_",codeDADA1,".RData",sep=""))
saveRDS(derep_reverse, paste(projectdir,"derep_reverse_",codeDADA1,".RData",sep=""))

#Inferring ASVs
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, multithread=TRUE)
#### Merging forward and reverse reads ####
merged_amplicons_5 <- mergePairs(dada_forward, derep_forward, dada_reverse,
                                 derep_reverse, minOverlap=5, verbose=TRUE)
merged_amplicons_default <- mergePairs(dada_forward, derep_forward, dada_reverse,
                                       derep_reverse, verbose=TRUE)
saveRDS(merged_amplicons_5,paste(projectdir,"merged_amplicons_5_",codeDADA1,".RData",sep=""))
saveRDS(merged_amplicons_default,paste(projectdir,"merged_amplicons_default",codeDADA1,".RData",sep=""))

#maxMismatch is by default 0 but you can make it say 1 or 2 if you’re finding that a lot of your forward and reverse reads are not merging
#minOverlap option has been change because they were sequenced at 250 bp
# Construct sequence table and remove chimeras
seqtab_5 <- makeSequenceTable(merged_amplicons_5)
saveRDS(seqtab_5, paste(projectdir,"seqtab_5_",codeDADA1,".rds",sep="")) 
seqtab_default <- makeSequenceTable(merged_amplicons_default)
saveRDS(seqtab_default,paste(projectdir,"seqtab_default_",codeDADA1,".rds",sep=""))

########Pool1 reverse##########
codeDADA1Rev="P1_rev"

forward_reads_rev1 <- paste0(samples, "_1rev_trimmed.fq")
reverse_reads_rev1 <- paste0(samples, "_2rev_trimmed.fq")
filtered_forward_reads_rev1 <- paste0(samples, "_1rev_filtered.fq")
filtered_reverse_reads_rev1 <- paste0(samples, "_2rev_filtered.fq")

#plotQualityProfile(filtered_forward_reads_rev1[1:4])
#plotQualityProfile(filtered_reverse_reads_rev1[1:4])
filtered_out_rev1 <- filterAndTrim(forward_reads_rev1, filtered_forward_reads_rev1, reverse_reads_rev1, filtered_reverse_reads_rev1, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=c(220,240), compress=TRUE, multithread=TRUE)

#when some input samples had no read that pass the filter
ls *_1rev_filtered.fq | sed "s/_1rev_filtered\.fq//g" > samples1_fil
samples1_fil <- scan("samples1_fil", what="character")
filtered_forward_reads_rev1 <- paste0(samples1_fil, "_1rev_filtered.fq")
filtered_reverse_reads_rev1 <- paste0(samples1_fil, "_2rev_filtered.fq")

#Generating an error model of our data by learning the specific error-signature of our dataset
err_forward_reads_rev1 <- learnErrors(filtered_forward_reads_rev1, multithread=TRUE)
err_reverse_reads_rev1 <- learnErrors(filtered_reverse_reads_rev1, multithread=TRUE)
#Dereplication
derep_forward_rev1 <- derepFastq(filtered_forward_reads_rev1, verbose=TRUE)
names(derep_forward_rev1) <- samples1_fil
derep_reverse_rev1 <- derepFastq(filtered_reverse_reads_rev1, verbose=TRUE)
names(derep_reverse_rev1) <- samples1_fil
saveRDS(derep_forward_rev1, paste(projectdir,"derep_forward_",codeDADA1Rev,".RData",sep=""))
saveRDS(derep_reverse_rev1, paste(projectdir,"derep_reverse_",codeDADA1Rev,".RData",sep=""))
#Inferring ASVs
dada_forward_rev1 <- dada(derep_forward_rev1, err=err_forward_reads_rev1, multithread=TRUE)
dada_reverse_rev1 <- dada(derep_reverse_rev1, err=err_reverse_reads_rev1, multithread=TRUE)
#### Merging forward and reverse reads ####
merged_amplicons_5_rev1 <- mergePairs(dada_forward_rev1, derep_forward_rev1, dada_reverse_rev1,
                                      derep_reverse_rev1, minOverlap=5, verbose=TRUE)
merged_amplicons_default_rev1 <- mergePairs(dada_forward_rev1, derep_forward_rev1, dada_reverse_rev1,
                                            derep_reverse_rev1, verbose=TRUE)
saveRDS(merged_amplicons_5_rev1, paste(projectdir,"merged_amplicons_5_",codeDADA1Rev,".RData",sep=""))
saveRDS(merged_amplicons_default_rev1, paste(projectdir,"merged_amplicons_default_",codeDADA1Rev,".RData",sep=""))

#maxMismatch is by default 0 but you can make it say 1 or 2 if you’re finding that a lot of your forward and reverse reads are not merging
#minOverlap option has been change because they were sequenced at 250 bp
# Construct sequence table and remove chimeras
seqtab_5_rev1 <- makeSequenceTable(merged_amplicons_5_rev1)
saveRDS(seqtab_5_rev1, paste(projectdir,"seqtab_5_",codeDADA1Rev,".rds",sep="")) 
seqtab_default_rev1 <- makeSequenceTable(merged_amplicons_default_rev1)
saveRDS(seqtab_default_rev1, paste(projectdir,"seqtab_default_",codeDADA1Rev,".rds",sep="")) 
#Generating a count table
table(nchar(getSequences(seqtab_5)))

#save the working environment
save.image(paste(projectdir,"all_Dada_final.RData",sep=""))

##########Pool2#################
setwd("/Users/ostaizkaaizpurua/Projects/GreenDogs/Pool2")
codeDADA2="P2"
samples2 <- scan("samples2", what="character")

# one holding the file names of all the forward reads (only one direction)
forward_reads2 <- paste0(samples2, "_1_trimmed.fq")
# and one with the reverse
reverse_reads2 <- paste0(samples2, "_2_trimmed.fq")

# and variables holding file names for the forward and reverse
filtered_forward_reads2 <- paste0(samples2, "_1_filtered.fq")
filtered_reverse_reads2 <- paste0(samples2, "_2_filtered.fq")


#plotQualityProfile(filtered_reverse_reads[1:4])
#plotQualityProfile(filtered_forward_reads[1:4])
#filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=c(223,223), compress=TRUE, multithread=TRUE)
#try different truncLen values
filtered_out2 <- filterAndTrim(forward_reads2, filtered_forward_reads2, reverse_reads2, filtered_reverse_reads2, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=c(240,220), compress=TRUE, multithread=TRUE)

#when some input samples had no read that pass the filter
ls *_1_filtered.fq | sed "s/_1_filtered\.fq//g" > samples2_fil
samples2_fil <- scan("samples2_fil", what="character")
filtered_forward_reads2 <- paste0(samples2_fil, "_1_filtered.fq")
filtered_reverse_reads2 <- paste0(samples2_fil, "_2_filtered.fq")

#Generating an error model of our data by learning the specific error-signature of our dataset
err_forward_reads2 <- learnErrors(filtered_forward_reads2, multithread=TRUE)
err_reverse_reads2 <- learnErrors(filtered_reverse_reads2, multithread=TRUE)
#Dereplication
derep_forward2 <- derepFastq(filtered_forward_reads2, verbose=TRUE)
names(derep_forward2) <- samples2
derep_reverse2 <- derepFastq(filtered_reverse_reads2, verbose=TRUE)
names(derep_reverse2) <- samples2
saveRDS(derep_forward2, paste(projectdir,"derep_forward_",codeDADA2,".RData",sep=""))
saveRDS(derep_reverse2, paste(projectdir,"derep_reverse_",codeDADA2,".RData",sep=""))
#derep_reverse2 <- readRDS(paste(projectdir,"DADA2/derep_reverse_",codeDADA2,".RData",sep=""))
#Inferring ASVs
dada_reverse2 <- dada(derep_reverse2, err=err_reverse_reads2, multithread=TRUE)
dada_forward2 <- dada(derep_forward2, err=err_forward_reads2, multithread=TRUE)
#### Merging forward and reverse reads ####
merged_amplicons2_5 <- mergePairs(dada_forward2, derep_forward2, dada_reverse2,
                                  derep_reverse2, minOverlap=5, verbose=TRUE)
merged_amplicons2_default <- mergePairs(dada_forward2, derep_forward2, dada_reverse2,
                                        derep_reverse2, verbose=TRUE)
saveRDS(merged_amplicons2_5,paste(projectdir,"merged_amplicons_5_",codeDADA2,".RData",sep=""))
saveRDS(merged_amplicons2_default, paste(projectdir,"merged_amplicons_default",codeDADA2,".RData",sep=""))

#maxMismatch is by default 0 but you can make it say 1 or 2 if you’re finding that a lot of your forward and reverse reads are not merging
#minOverlap option has been change because they were sequenced at 250 bp
# Construct sequence table and remove chimeras
seqtab2_5 <- makeSequenceTable(merged_amplicons2_5)
saveRDS(seqtab2_5,paste(projectdir,"seqtab_5_",codeDADA2,".rds",sep="")) 
seqtab2_default <- makeSequenceTable(merged_amplicons2_default)
saveRDS(seqtab2_default, paste(projectdir,"seqtab_default_",codeDADA2,".rds",sep="")) 

#################Pool2 reverse##########################

#samples2 <- scan("samples2", what="character")
codeDADA2Rev="P2_rev"

forward_reads2_rev1 <- paste0(samples2, "_1rev_trimmed.fq")
reverse_reads2_rev1 <- paste0(samples2, "_2rev_trimmed.fq")
filtered_forward_reads2_rev1 <- paste0(samples2, "_1rev_filtered.fq")
filtered_reverse_reads2_rev1 <- paste0(samples2, "_2rev_filtered.fq")

#plotQualityProfile(filtered_reverse_reads[1:4])
#plotQualityProfile(filtered_forward_reads[1:4])

filtered_out2_rev1 <- filterAndTrim(forward_reads2_rev1, filtered_forward_reads2_rev1, reverse_reads2_rev1, filtered_reverse_reads2_rev1, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=c(220,240), compress=TRUE, multithread=TRUE)

#when some input samples had no read that pass the filter
ls *_1rev_filtered.fq | sed "s/_1rev_filtered\.fq//g" > samples2_filt_rev
samples2_filt_rev <- scan("samples2_filt_rev", what="character")
filtered_forward_reads2_rev1 <- paste0(samples2_filt_rev, "_1rev_filtered.fq")
filtered_reverse_reads2_rev1 <- paste0(samples2_filt_rev, "_2rev_filtered.fq")

#Generating an error model of our data by learning the specific error-signature of our dataset
err_forward_reads2_rev1 <- learnErrors(filtered_forward_reads2_rev1, multithread=TRUE)
err_reverse_reads2_rev1 <- learnErrors(filtered_reverse_reads2_rev1, multithread=TRUE)
#Dereplication
derep_forward2_rev1 <- derepFastq(filtered_forward_reads2_rev1, verbose=TRUE)
names(derep_forward2_rev1) <- samples2
derep_reverse2_rev1 <- derepFastq(filtered_reverse_reads2_rev1, verbose=TRUE)
names(derep_reverse2_rev1) <- samples2
saveRDS(derep_forward2_rev1, paste(projectdir,"derep_forward_",codeDADA2Rev,".RData",sep=""))
saveRDS(derep_reverse2_rev1, paste(projectdir,"derep_reverse_",codeDADA2Rev,".RData",sep=""))
#Inferring ASVs
dada_forward2_rev1 <- dada(derep_forward2_rev1, err=err_forward_reads2_rev1, multithread=TRUE)
dada_reverse2_rev1 <- dada(derep_reverse2_rev1, err=err_reverse_reads2_rev1, multithread=TRUE)
#### Merging forward and reverse reads ####
merged_amplicons2_5_rev1 <- mergePairs(dada_forward2_rev1, derep_forward2_rev1, dada_reverse2_rev1,
                                       derep_reverse2_rev1, minOverlap=5, verbose=TRUE)
merged_amplicons2_default_rev1 <- mergePairs(dada_forward2_rev1, derep_forward2_rev1, dada_reverse2_rev1,
                                             derep_reverse2_rev1, verbose=TRUE)
saveRDS(merged_amplicons2_5_rev1, paste(projectdir,"merged_amplicons_5_",codeDADA2Rev,".RData",sep=""))
saveRDS(merged_amplicons2_default_rev1, paste(projectdir,"merged_amplicons_default_",codeDADA2Rev,".RData",sep=""))

#maxMismatch is by default 0 but you can make it say 1 or 2 if you’re finding that a lot of your forward and reverse reads are not merging
#minOverlap option has been change because they were sequenced at 250 bp
# Construct sequence table and remove chimeras
seqtab2_5_rev1 <- makeSequenceTable(merged_amplicons2_5_rev1)
saveRDS(seqtab2_5_rev1, paste(projectdir,"seqtab_5_",codeDADA2Rev,".rds",sep=""))
seqtab_default2_rev1 <- makeSequenceTable(merged_amplicons2_default_rev1)
saveRDS(seqtab_default2_rev1, paste(projectdir,"seqtab_default_",codeDADA2Rev,".rds",sep=""))

#save the working environment
save.image(paste(projectdir,"all_Dada_final.RData",sep=""))
############### Merge multiple runs (when we have different runs) ####################
codeDADA1
codeDADA1Rev
codeDADA2
codeDADA2Rev
codeGeneral="all"

st1 <- readRDS(paste(projectdir,"seqtab_default_",codeDADA1,".rds",sep=""))
st2 <- readRDS(paste(projectdir,"seqtab_default_",codeDADA1Rev,".rds",sep=""))
st3 <- readRDS(paste(projectdir,"seqtab_default_",codeDADA2,".rds",sep=""))
st4 <- readRDS(paste(projectdir,"seqtab_default_",codeDADA2Rev,".rds",sep=""))
rownames(st2) <- paste0(rownames(st2),"_rev")
rownames(st4) <- paste0(rownames(st4),"_rev")
st.all <- mergeSequenceTables(st1, st2, st3, st4)


#Chimera identification
seqtab.nochim <- removeBimeraDenovo(st.all, multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,paste(projectdir,"seqtab.nochim_",codeGeneral,".RData",sep=""))
#seqtab.nochim <- readRDS("/Users/jtk712/Documents/Projects/GreenDogs/DADA2/Seq1/seqtab.nochim_all.RData")

#Assigning taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/Users/ostaizkaaizpurua/Documents/Silva_DADA2/silva_nr_v138_train_set.fa.gz", tryRC=T)
saveRDS(taxa, paste(projectdir,"taxa_",codeGeneral,".RData",sep=""))
#taxa <- readRDS("/Users/jtk712/Documents/Projects/GreenDogs/DADA2/Seq1/taxa_all.RData")


####Summary####
# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons_default, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
summary_tab




#Phylogenetic tree
library(DECIPHER)
library(phangorn)
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

-----
#Extracting the standard goods from DADA2
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

#making and writing out a fasta of our final ASV seqs:
projectdir="/Users/ostaizkaaizpurua/Projects/GreenDogs/DADA2/"
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste(projectdir,"Results/Tables/ASVs_",codeGeneral,".fa",sep=""))

#count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste(projectdir,"Results/Tables/ASVs_counts_",codeGeneral,".tsv",sep=""), sep="\t", quote=F, col.names=NA)

asv_tab1 <- t(seqtab.nochim)
write.table(asv_tab1, paste(projectdir,"Results/Tables/ASVs_counts_names_",codeGeneral,".tsv",sep=""), sep="\t", quote=F, col.names=NA)


#tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, paste(projectdir,"Results/Tables/ASVs_taxonomy_",codeGeneral,".tsv",sep=""), sep="\t", quote=F, col.names=NA)

-----
#save the working environment
save.image(paste(projectdir,"DADA2/all_Dada_final",codeGeneral,".RData",sep=""))


#Phyloseq object
library(phyloseq)
asv_tab1 <- t(seqtab.nochim)
asv_tab1_table <- as.data.frame(asv_tab1)
#change the name and merge samples with same name
names(asv_tab1_table) = gsub(pattern = "_rev", replacement = "", x = names(asv_tab1_table)) 
merge_asv_tab1 <- t(rowsum(t(asv_tab1_table), group = colnames(asv_tab1_table), na.rm = T))
count_phy <- otu_table(merge_asv_tab1, taxa_are_rows=T)
count_phy1 <- otu_table(merge_asv_tab1, taxa_are_rows=F)
#make the object
ps1 <- phyloseq(count_phy,tax_table(taxa),phy_tree(fitGTR$tree))
ps1 <- phyloseq(count_phy,tax_table(taxa))

#change taxa names to ASVs
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))
#Check the change
rank_names(ps1)
table(tax_table(ps1)[, "Phylum"], exclude = NULL)

#add metadata to the phyloseq
projectdirGD="/Users/ostaizkaaizpurua/Projects/GreenDogs/"
#asv_fasta <- c(rbind(asv_headers, asv_seqs))
#write(asv_fasta, paste(projectdirGD,"Results/Tables/ASVs_",codeGeneral,".fa",sep=""))


metadata1 <- read.csv(paste(workingdir,"physeq_filtered1.csv",sep=""),header=TRUE,sep=";", row.names=1)
metadata <- read.table(paste(workingdir,"meta_dog.csv",sep=""),header=TRUE,sep=";",row.names=1)

metadata <- read.table(paste(projectdirGD,"meta_dog.csv",sep=""),header=TRUE,sep=";",row.names=1)
hierarchy_all <- tibble::rownames_to_column(metadata1, "Samples")
sample_tab <- sample_data(metadata1)
ps1_all=merge_phyloseq(ps1, sample_tab)
data.frame(ps1_all@sam_data) %>%
  filter(Working_status=="Lead") %>%
  filter(Extraction=="Saliva")

physeq_filtered_allL=merge_phyloseq(physeq_filtered_new, sample_tab)

############# DECONTAM ################
library(decontam)
sample_data(ps1_all)$is.neg <- sample_data(ps1_all)$Sample_type == "Control"
contamdf.prev <- isContaminant(ps1_all, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))
# getting IDs of those identified as likely contaminants
contam_asvs <- row.names(contamdf.prev[contamdf.prev$contaminant == TRUE, ])
# We can see this by looking at their taxonomic designations in our tax table:
contaminants <- ps1_all@tax_table[row.names(ps1_all@tax_table) %in% contam_asvs, ]
asv.cont <- rownames(contaminants)
#How many ASV contaminants in each sample?
samples <- colnames(ps1_all@otu_table)
contaminants <- c()
for (sample in samples){
  print(sample)
  one.column <- ps1_all@otu_table[,sample]
  if (colSums(one.column) > 0){
    s.depth <- colSums(one.column)
    asv.incolumn <- one.column[one.column > 0]
    total.asvs <- nrow(asv.incolumn)
    if (length(intersect(rownames(asv.incolumn), contam_asvs)) > 0){
      asv.contam <- asv.incolumn[rownames(asv.incolumn) %in%  contam_asvs,]
      contaminant.reads <- colSums(asv.contam)
      contaminant.asvs <- nrow(asv.contam)
      contaminants <- rbind(contaminants,c(s.depth,contaminant.reads,total.asvs,contaminant.asvs))
    }else{
      contaminants <- rbind(contaminants,c(s.depth,0,total.asvs,0))  
    }
  }else{
    contaminants <- rbind(contaminants,c(0,0))
  }
}

row.names(contaminants) <- samples
colnames(contaminants) <- c("total.reads", "contaminant.reads", "total.ASVs","contaminant.ASVs")
as.data.frame(contaminants)
detach("package:phangorn", unload=TRUE)
#detach("package:microbiome", unload=TRUE)
contaminants <- transform(contaminants, cont.read.percetage = contaminant.reads/total.reads*100)
write.table(contaminants,paste(projectdirR,"Tables/contaminants_Phylum.tsv",sep=""), quote = FALSE, sep = ",")

#summary_contaminants <- as.data.frame(contaminants)
#summary_contaminants_filt <- subset(summary_contaminants, total.reads >= 5000)
#summary_contaminants_filt <- subset(summary_contaminants_filt, cont.read.percetage <= 3)

saveRDS(ps1_all,"/Users/ostaizkaaizpurua/Projects/GreenDogs/DADA2/Seq1/physeq_Phylum_contaminants.RData")


# making new phyloseq object without ASV contaminats
physeq <- subset_samples(ps1_all, Sample_type == "Sample")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
goodTaxa <- setdiff(taxa_names(physeq), contam_asvs)
physeq_no_contam <- prune_taxa(goodTaxa, physeq)


#### FILTERING ####
#check the sequencing depth of lead_dogs
all_lead_check <- subset_samples(physeq_no_contam, Working_status == "Lead")
all_lead_check@sam_data
colSums(all_lead_check@otu_table)
#Filter the table
library(hilldiv)
count.all <- data.frame(physeq_no_contam@otu_table) 
count.all_filtdepth <- depth_filt(count.all,8000)
count.all_copyfilt <- copy_filt(count.all_filtdepth,0.0001)
count_filtered <- otu_table(count.all_copyfilt, taxa_are_rows=T)
otu_table(physeq_no_contam) <- count_filtered
physeq_filtered <- prune_taxa(taxa_sums(physeq_no_contam)>0, physeq_no_contam)
physeq_filtered_keep <- subset_samples(physeq_filtered, Keep == "Keep")
physeq_filtered_keep <- prune_taxa(taxa_sums(physeq_filtered_keep)>0, physeq_filtered_keep)
data.frame(physeq_filtered_keep@sam_data) %>%
  filter(Working_status=="Lead")#check that all leaders are there


#TAXONOMY filtering: Remove ASVs identified as Eukaryota or Archaea, and Bacteria lacking Phylum info in the TAXONOMY file
physeq_bacteria <- subset_taxa(physeq_filtered_keep, Kingdom == "Bacteria")
physeq_others <- subset_taxa(physeq_filtered_keep, Kingdom != "Bacteria")#35 Eukaryota
physeq_phylum <- subset_taxa(physeq_bacteria, !is.na(Phylum))
physeq_phylumNA <- subset_taxa(physeq_bacteria, is.na(Phylum))#6
physeq_class <- subset_taxa(physeq_phylum, !is.na(Class))
physeq_classNA <- subset_taxa(physeq_phylum, is.na(Class))#4

# Removing chloroplast and mitocondria
taxonomy <- physeq_phylum@tax_table
is.mitochondria <- physeq_phylum@tax_table[,"Family"] %in% "Mitochondria"
taxonomy <- taxonomy[!is.mitochondria,]
is.chloroplast <- taxonomy[,"Order"] %in% "Chloroplast"
taxonomy <- taxonomy[!is.chloroplast,]
taxonomy.m <- as.matrix(taxonomy)
tax_table(physeq_phylum) <- taxonomy.m
physeq_phylum_filtered <- prune_taxa(taxa_sums(physeq_phylum)>0, physeq_phylum)
tree <- force.ultrametric(tree, method = "extend")

phy_tree(physeq_phylum_filtered) <- tree

physeq_phylum_clean=physeq_phylum_filtered
#Save
asv.all.final <- data.frame(physeq_phylum@otu_table)
write.table(asv.all.final,paste(projectdir,"Results/Tables/asv.all.phylum.tsv",sep=""))

taxonomy.all.final <- data.frame(physeq_phylum@tax_table)
write.table(taxonomy.all.final,paste(projectdir,"Results/Tables/taxonomy.all.phylum.tsv",sep=""))

sample.all.final <- data.frame(physeq_phylum@sam_data)
write.table(sample.all.final,paste(projectdir,"Results/Tables/sample.all.phylum.tsv",sep=""))

saveRDS(physeq_phylum_filtered,paste(projectdirSeq, "physeq_phylum_filtered_tree.RData",sep=""))
physeq_phylum_filtered <- readRDS(paste(projectdirSeq, "physeq_phylum_filtered.RData",sep=""))

saveRDS(physeq_class,"/Users/ostaizkaaizpurua/Projects/GreenDogs/DADA2/Seq1/physeq_Class.RData")

####Rarefaction curves####
library(ranacapa)
curves <- ggrare(physeq_phylum, step = 1000, color = "Extraction", se = FALSE)
pdf(paste(projectdir,"Results/Figures/rarecurve_together.pdf",sep=""),width=8, height=5)
curves + scale_color_manual(values=c("lightcoral", "steelblue3"))
dev.off()

###Number of samples units
#Get samples
metadata <- as.data.frame(physeq_phylum@sam_data)
asv.table <- as.data.frame(physeq_phylum@otu_table)
F_samples = row.names(metadata[metadata$Extraction == "Faecal",])
S_samples = row.names(metadata[metadata$Extraction == "Saliva",])
#Subset countable
countable_F <- asv.table[,colnames(asv.table) %in% F_samples]
countable_F[countable_F > 0] <- 1
countable_S <- asv.table[,colnames(asv.table) %in% S_samples]
countable_S[countable_S > 0] <- 1
#Add them to a list
countlist <- list(Faecal=countable_F,Saliva=countable_S)
library(iNEXT)
out <- iNEXT(countlist, q=0, datatype="incidence_raw")
pdf(paste(projectdir,"Results/Figures/sample_units.pdf",sep=""),width=8, height=5)
ggiNEXT(out)
dev.off()

physeq_faecal <- subset_samples(physeq_phylum, Extraction == "Faecal")
physeq_faecal <- prune_taxa(taxa_sums(physeq_faecal)>0, physeq_faecal)
physeq_faecal <- subset_samples(physeq_faecal, Working_status %in% c("Follower", "Lead"))
physeq_faecal <- prune_taxa(taxa_sums(physeq_faecal)>0, physeq_faecal)

physeq_saliva <- subset_samples(physeq_phylum, Extraction == "Saliva")
physeq_saliva <- prune_taxa(taxa_sums(physeq_saliva)>0, physeq_saliva)

#Faecal outliers
count_fecal <- data.frame(physeq_faecal@otu_table)
samples_fecal <- data.frame(physeq_faecal@sam_data)
samples_fecal$Pack=factor(samples_fecal$Pack)

hierarchy_faecal <- hierarchy_all[which(hierarchy_all[,1] %in% colnames(count_fecal)),]

#Check outliers
pslog <- transform_sample_counts(physeq_phylum, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "Extraction") +
  labs(col = "Pack") +
  coord_fixed(sqrt(evals[2] / evals[1]))

##Removing outlier from the phyloseq
Samples_toRemove <- c("Fecal1.FTA_02_06","Fecal1.FTA_02_02","Fecal8.FTA_05_02","Saliva2.S_04_01","Saliva2.S_04_10")
physeq_phylum_clean <- prune_samples(!(sample_names(physeq_phylum) %in% Samples_toRemove), physeq_phylum)

#Save tables and objects
saveRDS(physeq_phylum_clean,"/Users/ostaizkaaizpurua/Projects/GreenDogs/DADA2/Seq1/physeq_phylum_clean.RData")
plot_tree (physeq_phylum_clean)

count_dogs <- data.frame(physeq_phylum_clean@otu_table)
write.table(count_dogs, paste(projectdirDada,"Results/Tables/ASVs_counts_clean_final.tsv",sep=""), sep="\t", quote=F, col.names=NA)

samples_dogs <- data.frame(physeq_phylum_clean@sam_data)
write.table(samples_dogs, paste(projectdirDada,"Results/Tables/Samples_clean_final.tsv",sep=""), sep="\t", quote=F, col.names=NA)

tax_dogs <- data.frame(physeq_phylum_clean@tax_table)
write.table(tax_dogs, paste(projectdirDada,"Results/Tables/Taxonomy_clean_final.tsv",sep=""), sep="\t", quote=F, col.names=NA)

tree_clean = phy_tree(physeq_phylum_clean)
ape::write.tree(tree_clean, "/Users/ostaizkaaizpurua/Projects/GreenDogs/DADA2/Seq1/tree1_clean.tree")
tree <- ape::read.tree(paste(projectdirSeq, "tree1_clean.tree",sep=""))

save.image(paste(projectdirSeq,"all_Dada_final.RData",sep=""))

