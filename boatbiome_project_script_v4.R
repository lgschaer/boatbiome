#packages used
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.8")
library(dada2)
library(phyloseq)
library(csv)
library(tidyverse)
library(vegan)
#install.packages("FSA")
library(FSA)
#install.packages("remotes")
#remotes::install_github("opisthokonta/tsnemicrobiota")
library(remotes)
library(tsnemicrobiota)

#-----MAKING ASV AND TAXA TABLES-----#

#path to fastq files
path <- "/home/lgschaer/old/boatbiome/April-Boat-Microbiome/Boat-Microbiome/04082019-FASTQ" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
head(sample.names)

#visualize quality profiles
plotQualityProfile(fnFs[1:2])           #forward reads
plotQualityProfile(fnRs[1:2])           #reverse reads

#place filtered files in "filtered", a subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

#standard filtering parameters: 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(295, 175),
                     maxN=0, maxEE=Inf, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

read.stats <- out %>%
  as_tibble() %>%
  summarize(
    max.in = max(reads.in),
    min.in = min(reads.in),
    mean.in = mean(reads.in),
    median.in = median(reads.in),
    sum.in = sum(reads.in),
    count.in = sum(ifelse(reads.in > 0, 1, 0)),
    under1000.in = sum(ifelse(reads.in < 1000, 1, 0)),
    max.out = max(reads.out),
    min.out = min(reads.out),
    mean.out = mean(reads.out),
    median.out = median(reads.out),
    sum.out = sum(reads.out),
    count.out = sum(ifelse(reads.out > 0, 1, 0)),
    under1000.out = sum(ifelse(reads.out < 1000, 1, 0))
  )
View(t(read.stats))

#save csvs
write.csv(out, "/home/lgschaer/old/boatbiome/April-Boat-Microbiome/Boat-Microbiome/out.csv")
write.csv(read.stats, "/home/lgschaer/old/boatbiome/April-Boat-Microbiome/Boat-Microbiome/read.stats.csv")

#learn about error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plot errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

#name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inferance
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#inspecting dada-class object
dadaFs[[1]]
dadaRs[[1]]

#merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/home/lgschaer/old/Oil_Genomes/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save sequence table and taxa table
write_rds(seqtab.nochim, "/home/lgschaer/old/boatbiome/April-Boat-Microbiome/Boat-Microbiome/seqtab.rds")     #sequence table
write_rds(taxa.print, "/home/lgschaer/old/boatbiome/April-Boat-Microbiome/Boat-Microbiome/taxa.rds")          #taxa table


#-----MAKING A PHYLOSEQ OBJECT-----#

#load data into R for phyloseq analysis

#load sample data
sdata <- as.csv("/home/lgschaer/old/boatbiome/phyloseq_input/boatbiome_sampledata.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)

sdata2 <- sdata %>% 
  filter(Sample_Project == "Boat_Microbiome") %>%                                                               #filter to only include samples from this project
  filter(datatype == "water"|datatype == "swab"|datatype == "bilge"|datatype == "dock"|datatype == "air")%>%    #remove quality control blanks
  as_tibble() %>%                                                   #change to tibble format
  #group_by(sampletype) %>%                                          #group by sampletype
  mutate(                                                           #add columns with sample number and significance data
    a = ifelse(sampletype == "air", "AB", "X"),
    b = ifelse(sampletype == "bilge", "A", a),
    c = ifelse(sampletype == "boatback", "C", b),
    d = ifelse(sampletype == "boatside", "NS", c),
    e = ifelse(sampletype == "dock", "D", d),
    sig = ifelse(sampletype == "water", "BCD", e),
    v = ifelse(sampletype == "air", "n=14", "X"),
    w = ifelse(sampletype == "bilge", "n=11", v),
    x = ifelse(sampletype == "boatback", "n=6", w),
    y = ifelse(sampletype == "boatside", "n=7", x),
    z = ifelse(sampletype == "dock", "n=9", y),
    n = ifelse(sampletype == "water", "n=131", z)
  ) %>%
  select(-c("a", "b", "c", "d", "e", "v", "w", "x", "y", "z"))     #remove extra columns added during mutate
head(sdata2)                                                       #view data to make sure everything is OK

sdata3 <- sdata2 %>%                                               #make data to use in phyloseq object with Sample_ID as rownames
  column_to_rownames(var = "Sample_ID")
head(sdata3)                                                       #view data to make sure everything is OK

#load sequence table
sequence_table <- readRDS("/home/lgschaer/old/boatbiome/April-Boat-Microbiome/Boat-Microbiome/seqtab.rds")
colnames(sequence_table) <- NULL                                   #remove column names
sequence_table <- as.data.frame(sequence_table)                    #convert to data frame format
sequence_table[1:4,1:4]                                            #view a portion to make sure it looks how we expect
sequence_table <- rownames_to_column(sequence_table, var = "Sample_ID")    #change the rownames to a column with header "Sample_ID"
sequence_table[1:4,1:4]                                            #view again
class(sequence_table)                                              #check class, should be data frame

sequence_table <- sequence_table %>%
  as_tibble()%>%
  mutate(Sample_ID = str_replace_all(Sample_ID, "-", "_"))%>%      #change Sample IDs to match metadata
  as.data.frame()                                                  #change back to data frame format
sequence_table[1:4,1:4]                                            #view to make sure all is as expected

#join metadata and sequence table to filter samples
seq_table <- sdata2 %>%                                            #start with filtered metadata
  full_join(sequence_table, by = "Sample_ID")                      #left join sequence table by matching up rownames. This will drop all items from sequence table that are not included in the metadata file (samples from other projects, quality control blanks)
seq_table[1:10,1:10]                                               #view to make sure it is what we expect


#subset to only include sampletypes we want
subset_all <- seq_table %>%
  subset(Sample_Project=="Boat_Microbiome")%>%                     #make a subset to only include the samples we want (this should have already been done by using left_join command, but it doesn't hurt to double check!)
  filter(datatype == "water"|datatype == "swab"|datatype == "bilge"|datatype == "dock"|datatype == "air")%>%            #filter to only include desired sample types
  select(-c("Sample_Project", "datatype", "sampletype", "date", "filtertype", "n", "sig")) %>%                          #remove metadata columns from sequence table after desired filtering is done.
  column_to_rownames(var = "Sample_ID")                                                                                 #change the "Sample_ID column to rownames to get sequence table in the right format for phyloseq
subset_all[1:10,1:10]                                              #view a portion of sequence table to make sure all is well
class(subset_all)                                                  #class should be data frame

#make nonzero subset to remove all columns with zero taxa counts
subset_all <- as.matrix(subset_all)                                #change to matrix format
m <- (colSums(subset_all, na.rm=TRUE) != 0)                        #T if colSum is not 0, F otherwise
nonzero <- subset_all[, m]                                         #all the non-zero columns


#load taxa table
taxa_table <- readRDS("/home/lgschaer/old/boatbiome/April-Boat-Microbiome/Boat-Microbiome/taxa.rds")
taxa_table <- as.matrix(taxa_table)                                #change to matrix format
taxa_table[1:5,1:5]                                                #view to make sure everything looks good

#make phyloseq object
samdata = sample_data(sdata3)                                      #define sample data
colnames(nonzero) <- NULL                                          #remove column names from "nonzero"
seqtab = otu_table(nonzero, taxa_are_rows = FALSE)                 #define sequence table
taxtab = tax_table(taxa_table)                                     #define taxa table
rownames(taxtab) <- NULL                                           #remove rownames from taxa table

phyloseq_object_all = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
phyloseq_object_all

#sequencing depth before rarefying
bilge <- subset_samples(phyloseq_object_all, sampletype == "bilge")
median(sample_sums(bilge))
air <- subset_samples(phyloseq_object_all, sampletype == "air")
median(sample_sums(air))
hull <- subset_samples(phyloseq_object_all, sampletype == "boatside")
median(sample_sums(hull))
transom<- subset_samples(phyloseq_object_all, sampletype == "boatback")
median(sample_sums(transom))
dock<- subset_samples(phyloseq_object_all, sampletype == "dock")
median(sample_sums(dock))
water<- subset_samples(phyloseq_object_all, sampletype == "water")
median(sample_sums(water))


#sample counts before rarefying
sample_counts <- sample_data(phyloseq_object_all) %>%
  group_by(sampletype) %>%
  mutate(Count = 1) %>%
  summarise(SumCount = sum(Count))
sample_counts

#normalize data
#Delete samples with a mean of less than 1000
samplesover1000_all <- subset_samples(phyloseq_object_all, sample_sums(phyloseq_object_all) > 1000)

#Check if there are OTUs with no counts, if so how many?
any(taxa_sums(samplesover1000_all) == 0)
sum(taxa_sums(samplesover1000_all) == 0)

#Prune OTUs with no counts 
prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)
any(taxa_sums(prune_samplesover1000_all) == 0)

#make sure seed is set the same each time, set to 81 here
rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))


#sequencing depth after rarefying
bilge <- subset_samples(rarefy_samplesover1000_all, sampletype == "bilge")
median(sample_sums(bilge))
air <- subset_samples(rarefy_samplesover1000_all, sampletype == "air")
median(sample_sums(air))
hull <- subset_samples(rarefy_samplesover1000_all, sampletype == "boatside")
median(sample_sums(hull))
transom<- subset_samples(rarefy_samplesover1000_all, sampletype == "boatback")
median(sample_sums(transom))
dock<- subset_samples(rarefy_samplesover1000_all, sampletype == "dock")
median(sample_sums(dock))
water<- subset_samples(rarefy_samplesover1000_all, sampletype == "water")
median(sample_sums(water))


#assign colors to samples
sample_colors <- c("air" = "lightgoldenrod", "bilge" = "blue", "boatback" = "darkred", "boatside" = "darkolivegreen4", "dock" = "darkorange", "water" = "lightblue")

sample_types <- c("bilge", "boatback", "boatside", "water", "air", "dock")

sample_labels <- c("bilge" ="Bilge", "water" = "Water", "boatback" = "Transom", "boatside" = "Hull", "dock" = "Dock", "air" = "Air")

#filter out eukaryotes and mitochondria
justbacteria <- rarefy_samplesover1000_all %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
justbacteria

#saveRDS file
saveRDS(justbacteria, file = "/home/lgschaer/old/boatbiome/bb-phyloseq.rds")

#-----ALPHA DIVERISTY-----#

#Violin plot of alpha diversity Observed OTUs and Shannon Diversity (with color)
justbacteria %>%                                                     #phyloseq object
  plot_richness(
    x = "sampletype",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = sampletype), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  geom_text(aes(label = n, y = 1), vjust = 1, size = 5)+
  geom_text(aes(label = sig, y = 1), vjust = -25, size = 6)+
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = sample_colors)+                         #set fill colors
  scale_x_discrete(                                                  #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

#Violin plot of alpha diversity Observed OTUs and Shannon Diversity (no color)
justbacteria %>%                                                     #phyloseq object
  plot_richness(
    x = "sampletype",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  geom_text(aes(label = n, y = 1), vjust = 1, size = 5)+
  geom_text(aes(label = sig, y = 1), vjust = -25, size = 6)+
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = sample_colors)+                         #set fill colors
  scale_x_discrete(                                                  #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

#Get alpha diversity data into data frame format
richness <- justbacteria %>%
  estimate_richness(measures = c("Observed", "Shannon")) %>%           #specify which measures
  rownames_to_column(var = "Sample_ID") %>%                                #add column name to SampleID column
  separate(Sample_ID, into = c("X", "Sample_ID"), sep = "X")%>%
  as_tibble() %>%
  select(-c("X"))
head(richness)

alphadiv <- richness %>%
  left_join(sdata2, by = "Sample_ID")
head(alphadiv)

#finding means and medians
rich_summary <- richness %>%
  left_join(sdata2, by = "Sample_ID") %>%                                       #join to meta data
  group_by(sampletype) %>%                                                      #group by sample type
  summarize(                                                                    #add columns for stat summaries
    maxObs = max(Observed),
    maxShan = max(Shannon),
    minObs = min(Observed),
    minShan = min(Shannon),
    rangeObs = maxObs - minObs,
    meanObserved = mean(Observed),
    medianObserved = median(Observed),
    rangeShan = maxShan - minShan,
    meanShannon = mean(Shannon),
    medianShannon = median(Shannon)
  )
head(rich_summary)                                                            #check that we have the columns we want

#-----ALPHA DIVERSITY STATISTICS-----#

#check assumptions for ANOVA --> non-normal distribution, will do Kuruskal-Wallis test instead
histogram(richness$Observed)
histogram(richness$Shannon)

#power analysis to determime whether we have sufficent power use an ANOVA with non-normally distributed data.
#count samples
sample_counts <- sample_data(justbacteria) %>%
  group_by(sampletype) %>%
  mutate(Count = 1)%>%
  summarise(SumCount = sum(Count)) %>%
  mutate(countpercategory = mean(SumCount))
sample_counts

pwr.anova.test(k=6, f=0.23, sig.level=.05, power=.8)

#Kruskal-Wallis Test
set.seed(81)

#Observed
kruskal.test(Observed ~ sampletype, data = alphadiv) 

#Shannon
kruskal.test(Shannon ~ sampletype, data = alphadiv) 

#Dunn test (post hoc)

##Observed
dunnO <- dunnTest(Observed ~ sampletype,
         data=alphadiv,
         method="bh")
dunnO

dunnO <- dunnO$res
View(dunnO)

##Shannon
dunnS <- dunnTest(Shannon ~ sampletype,
                 data=alphadiv,
                 method="bh")
dunnS

dunnS <- dunnS$res
View(dunnS)

#-----BETA DIVERSITY-----#

#ordination
distance <- ordinate(
  physeq = justbacteria, 
  method = "PCoA", 
  distance = "bray"
)
summary(distance)
distance

#t-SNE plot
head(sdata3)

tsne <- tsne_phyloseq(justbacteria, distance = "bray", perplexity = 25, dimensions = 2,
                      precomputed_distance = NULL, pseudocounts = 1, verbose = 1,
                      rng_seed = 81, philr_options = list(), control = list())
summary(tsne)

justbacteria

#tSNE Plot
plot_tsne_phyloseq(justbacteria, tsne, color = "sampletype", shape = "filtertype") +
  geom_point(aes(color = sampletype, fill = sampletype), color = "black", size = 5, show.legend = TRUE) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23))+
  scale_fill_manual(values = sample_colors,
                    #breaks = sample_types,
                    labels = sample_labels) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(color = FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command



#-----BETA DIVERSTIY STATISTICS-----#

#PERMANOVA
set.seed(81)

#subset phyloseq object, all samples by datatype
ab <- subset_samples(justbacteria, sampletype %in% c("air", "bilge"))
abb <- subset_samples(justbacteria, sampletype %in% c("boatback", "air"))
abs <- subset_samples(justbacteria, sampletype %in% c("boatside", "air"))
ad <- subset_samples(justbacteria, sampletype %in% c("air", "dock"))
aw <- subset_samples(justbacteria, sampletype %in% c("air", "water"))
bbb <- subset_samples(justbacteria, sampletype %in% c("boatback", "bilge"))
bbbs <- subset_samples(justbacteria, sampletype %in% c("boatback", "boatside"))
bbw <- subset_samples(justbacteria, sampletype %in% c("boatback", "water"))
bsb <- subset_samples(justbacteria, sampletype %in% c("boatside", "bilge"))
bsw <- subset_samples(justbacteria, sampletype %in% c("boatside", "water"))
db <- subset_samples(justbacteria, sampletype %in% c("dock", "bilge"))
dbb <- subset_samples(justbacteria, sampletype %in% c("boatback", "dock"))
dbs <- subset_samples(justbacteria, sampletype %in% c("boatside", "dock"))
dw <- subset_samples(justbacteria, sampletype %in% c("dock", "water"))
wb <- subset_samples(justbacteria, sampletype %in% c("bilge", "water"))

# Calculate bray curtis distance matrix, all samples
ab_bray <- phyloseq::distance(ab, method = "bray")
abb_bray <- phyloseq::distance(abb, method = "bray")
abs_bray <- phyloseq::distance(abs, method = "bray")
ad_bray <- phyloseq::distance(ad, method = "bray")
aw_bray <- phyloseq::distance(aw, method = "bray")
bbb_bray <- phyloseq::distance(bbb, method = "bray")
bbbs_bray <- phyloseq::distance(bbbs, method = "bray")
bbw_bray <- phyloseq::distance(bbw, method = "bray")
bsb_bray <- phyloseq::distance(bsb, method = "bray")
bsw_bray <- phyloseq::distance(bsw, method = "bray")
db_bray <- phyloseq::distance(db, method = "bray")
dbb_bray <- phyloseq::distance(dbb, method = "bray")
dbs_bray <- phyloseq::distance(dbs, method = "bray")
dw_bray <- phyloseq::distance(dw, method = "bray")
wb_bray <- phyloseq::distance(wb, method = "bray")

# make a data frame from the sample_data, all samples
ab.sam <- data.frame(sample_data(ab))
abb.sam <- data.frame(sample_data(abb))
abs.sam <- data.frame(sample_data(abs))
ad.sam <- data.frame(sample_data(ad))
aw.sam <- data.frame(sample_data(aw))
bbb.sam <- data.frame(sample_data(bbb))
bbbs.sam <- data.frame(sample_data(bbbs))
bbw.sam <- data.frame(sample_data(bbw))
bsb.sam <- data.frame(sample_data(bsb))
bsw.sam <- data.frame(sample_data(bsw))
db.sam <- data.frame(sample_data(db))
dbb.sam <- data.frame(sample_data(dbb))
dbs.sam <- data.frame(sample_data(dbs))
dw.sam <- data.frame(sample_data(dw))
wb.sam <- data.frame(sample_data(wb))

# Adonis test, all samples
adonis(ab_bray ~ sampletype, data = ab.sam)
adonis(abb_bray ~ sampletype, data = abb.sam)
adonis(abs_bray ~ sampletype, data = abs.sam)
adonis(ad_bray ~ sampletype, data = ad.sam)
adonis(aw_bray ~ sampletype, data = aw.sam)
adonis(bbb_bray ~ sampletype, data = bbb.sam)
adonis(bbbs_bray ~ sampletype, data = bbbs.sam)
adonis(bbw_bray ~ sampletype, data = bbw.sam)
adonis(bsb_bray ~ sampletype, data = bsb.sam)
adonis(bsw_bray ~ sampletype, data = bsw.sam)
adonis(db_bray ~ sampletype, data = db.sam)
adonis(dbb_bray ~ sampletype, data = dbb.sam)
adonis(dbs_bray ~ sampletype, data = dbs.sam)
adonis(dw_bray ~ sampletype, data = dw.sam)
adonis(wb_bray ~ sampletype, data = wb.sam)

#-----SOURCE TRACKER----#

#Pepare OTU table for analysis in SourceTracker

#remove otu table from phyloseq object
otutab <- (otu_table(phyloseq_object_all))
otutab[,5:1][1:5,]

#transpose
otu_t <- t(otutab)
otu_t[,5:1][1:5,]

#write table
write.table(otu_t, "/home/lgschaer/old/boatbiome/st_input/otu_t.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

#Analyze output from SourceTracker

#import data
sample_info <- as.csv("/home/lgschaer/old/boatbiome/phyloseq_input/boatbiome_sampledata.csv", row.names=1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sample_info)

t_sample_info <- sample_info %>%
  mutate(SampleID = Sample_ID)%>%
  as_tibble()
head(t_sample_info)

#Ran sourcetracker2 five times, will average the results to show variation
#Run 1
run1 <- read.delim("/home/lgschaer/old/boatbiome/st_input/mixingproportions_run1/mixing_proportions.txt")
run1 <- run1 %>%
  mutate(Air = air, Dock = dock, Water = water)%>%
  select(-c("air", "dock", "water"))%>%
  gather("Air", "Dock", "Water", "Unknown", key = source, value = proportion1)
head(run1)

#Run 2
run2 <- read.delim("/home/lgschaer/old/boatbiome/st_input/mixingproportions_run2/mixing_proportions.txt")
run2 <- run2 %>%
  mutate(Air = air, Dock = dock, Water = water)%>%
  select(-c("air", "dock", "water"))%>%
  gather("Air", "Dock", "Water", "Unknown", key = source, value = proportion2)
head(run2)

#Run 3
run3 <- read.delim("/home/lgschaer/old/boatbiome/st_input/mixingproportions_run3/mixing_proportions.txt")
run3 <- run3 %>%
  mutate(Air = air, Dock = dock, Water = water)%>%
  select(-c("air", "dock", "water"))%>%
  gather("Air", "Dock", "Water", "Unknown", key = source, value = proportion3) 
head(run3)

#Run 4
run4 <- read.delim("/home/lgschaer/old/boatbiome/st_input/mixingproportions_run4/mixing_proportions.txt")
run4 <- run4 %>%
  mutate(Air = air, Dock = dock, Water = water)%>%
  select(-c("air", "dock", "water"))%>%
  gather("Air", "Dock", "Water", "Unknown", key = source, value = proportion4)
head(run4)

#Run 5
run5 <- read.delim("/home/lgschaer/old/boatbiome/st_input/mixingproportions_run5/mixing_proportions.txt")
run5 <- run5 %>%
  mutate(Air = air, Dock = dock, Water = water)%>%
  select(-c("air", "dock", "water"))%>%
  gather("Air", "Dock", "Water", "Unknown", key = source, value = proportion5)
head(run5)

#Join data from the five runs together
wda_prop <- run1 %>%
  left_join(run2, by = c("SampleID", "source")) %>%
  left_join(run3, by = c("SampleID", "source")) %>%
  left_join(run4, by = c("SampleID", "source")) %>%
  left_join(run5, by = c("SampleID", "source")) %>%
  left_join(t_sample_info, by = "SampleID") %>%
  filter(Sample_Project == "Boat_Microbiome") %>%
  as.tibble() %>%
  gather(proportion1, proportion2, proportion3, proportion4, proportion5, key = "run", value = "proportion") %>%
  group_by(source, sampletype) %>%
  summarise(
    sum_prop = sum(proportion),
    mean_prop = mean(proportion),
    sd = sd(proportion)) %>%
  group_by(sampletype) %>%
  mutate(perc_prop = (sum_prop/sum(sum_prop))*100) 
View(wda_prop)

#COLOR FIGURE
#make a custom color vector

sample_colors2 <- c("Air" = "lightgoldenrod", "Dock" = "darkorange", "Water" = "lightblue", "Unknown" = "grey")

View(wda_prop)

#water, air, dock as sources
ggplot(data = wda_prop, aes(x = sampletype, y = sum_prop, fill = source))+
  geom_col(position = "fill", color = "black", size = 0.2)+
  theme_classic()+
  ylab("Proportion")+
  xlab(NULL)+
  scale_fill_manual(values = sample_colors2)+                         #set fill colors
  scale_x_discrete(                                                  #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+ 
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

#BW FIGURE

#make black and white color vector

blkwht <- c("black", "grey30", "grey", "white")

#water, air, dock as sources
ggplot(data = wda_prop, aes(x = sampletype, y = sum_prop, fill = source))+
  geom_col(position = "fill", color = "black", size = 0.2)+
  theme_classic()+
  ylab("Proportion")+
  xlab(NULL)+
  scale_fill_manual(values = blkwht)+                          #set fill colors
  scale_x_discrete(                                                  #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+ 
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))


#-----MAKING TAXA PLOTS-----#

#phylum abundance
phylumabundance <- justbacteria %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Class) 
head(phylumabundance)


#mostcommon taxa
taxa <- phylumabundance %>%
  select(Phylum, Class, sampletype, Abundance, date) %>%          #select columns we want
  group_by(Phylum, Class, sampletype, date) %>%                   #group by Phylum, Class, Sample type and Date
  summarize(
    avg_abundance = mean(Abundance)                               #calculate the average abundance of each taxa
  ) %>%
  ungroup()%>%                                                    #ungroup data
  mutate(
    Phylum = as.character(Phylum),                                #change Phylum and Class columns to character vectors
    Class = as.character(Class),
    Taxa = ifelse(Phylum == "Proteobacteria", Class, Phylum),     #Create "Taxa" column which shows Class for Proteobacteria, Phylum for all other phyla
    Legend = ifelse(avg_abundance < 0.01, "<1%", Taxa),           #Label taxa present at low abundance "< 1%" for both Taxa and Class
  )
head(taxa)

#vector of colors for taxa plot
phylum_colors <- c(
  "black",     "orchid1",   "darkgreen",     "yellow",     "darkcyan",        "cyan",  
  "magenta",   "blue",      "red",           "white",      "tan4",        
  "orange",    "green",     "mediumpurple1", "grey77",     "palegoldenrod",   "firebrick",  
  "purple4",   "coral1",    "lightblue",     "grey47"           
)  

#Taxa plot by date faceted for each sample type

#prepare data
taxa2 <- taxa %>%
  mutate(Sampletype = as.factor(sampletype),
         Date = str_replace_all(date, "^6", "6-"),
         Date = str_replace_all(Date, "^7", "7-"))%>%
  group_by(Legend, Date, Sampletype)%>%
  summarise(avg_abundance2 = sum(avg_abundance))
head(taxa2)

#new facet label names for sampletype variable
type.labs <- c("Bilge", "Transom", "Hull", "Water", "Air", "Dock")
names(type.labs) <- c("bilge", "boatback", "boatside", "water", "air", "dock")

#plot
ggplot(taxa2)+
  geom_col(mapping = aes(x = Date, y = avg_abundance2, fill = Legend), color = "black", position = "fill", show.legend = TRUE)+
  facet_wrap(~Sampletype, labeller = labeller(Sampletype = type.labs))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  theme_minimal()+
  theme(strip.text.x = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 16),
        axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        title = element_text(size = 25))


#taxa plot grouped by sampletype

#prepare data
taxa3 <- taxa %>%
  group_by(Legend, sampletype)%>%
  summarise(avg_abundance2 = sum(avg_abundance))

#plot
ggplot(taxa3)+
  geom_col(mapping = aes(x = sampletype, y = avg_abundance2, fill = Legend), color = "black", position = "fill", show.legend = TRUE)+
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL)+
  theme_minimal()+
  scale_x_discrete(                                                  #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+  
  theme(legend.text = element_text(size = 16),
        axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))

#make data frame to summarize abundances of each phyla.

phylaSum <- justbacteria %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  select(sampletype, date, Phylum, Abundance) %>%
  filter(Abundance >= 0.01) %>%
  group_by(sampletype, Phylum) %>%
  summarise(
    meanPhy = round(mean(Abundance),2),
    sd = round(sd(Abundance),2)
  ) %>%
  mutate(sd = ifelse(is.na(sd), 0, sd))%>%
  unite(phySd, meanPhy, sd, sep = "+/-", remove = TRUE)%>%
  spread(key = "sampletype", value = "phySd", fill = "NA")
View(phylaSum)

write_csv(phylaSum, "/home/lgschaer/old/boatbiome/phylaSum.csv")

#-----COEFICIENT OF VARIATION----#
head(sdata2)


cvar <- justbacteria %>%                                      #phyloseq object
  estimate_richness(                                           #compare diversity of datatype
    measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var = "Sample_ID") %>%
  mutate(Sample_ID = str_replace_all(Sample_ID, "X", ""))%>%
  left_join(sdata2, by = "Sample_ID")%>%
  as.tibble()%>%
  mutate(std = sd(Shannon))%>%
  group_by(sampletype, date, filtertype) %>%
  mutate(
    avg = mean(Shannon),
    CV = std/avg
  ) %>%
  ungroup()%>%
  transmute(Sample_ID, sampletype, CV)
head(cvar)

cvar.sum <- cvar %>%
  group_by(sampletype) %>%
  summarize(
    maxCV = max(CV),
    meanCV = mean(CV),
    minCV = min(CV)
  )
head(cvar.sum)

#make coefficient of variation plot (color version)
ggplot(cvar, aes(x = sampletype, y = CV, fill = sampletype))+
  scale_fill_manual(values = sample_colors)+
  scale_x_discrete(limits = sample_types, labels = sample_labels)+ 
  geom_point(show.legend = FALSE)+
  geom_boxplot(show.legend = FALSE)+
  theme_classic()+
  xlab(NULL)+
  ylab("Coefficient of Variation")+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))

#make coefficient of variation plot (black and white)
ggplot(cvar, aes(x = sampletype, y = CV))+
  scale_x_discrete(limits = sample_types, labels = sample_labels)+ 
  geom_point(show.legend = FALSE)+
  geom_boxplot(show.legend = FALSE)+
  theme_classic()+
  xlab(NULL)+
  ylab("Coefficient of Variation")+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))
