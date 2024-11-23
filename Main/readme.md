------------------------------------------------------------------------

title: "Sample assembly"

author: "Thomas Stocker"

date: "11-11-2024"

bibliography: references.bib

------------------------------------------------------------------------

<<<<<<< Updated upstream:readme.md
# **Welcome to the sample assembly branch:**
=======
# **Experimental design**
>>>>>>> Stashed changes:Main/readme.md

**Welcome to the Sample assembly branch**

This branch deals with the assembly of the samples but also the automated renaming within R for both the species and or samples - this greatly saves time - instead of manual processing. But before we get into any coding, here is the brief description of the samples.

## **Samples:**

As discussed in the "Main" branch, the samples are labeled EF\*\*\*L-1 to EF\*\*\*L-14, with EF\*\*\* representing the cow's identifier, and the L-1 to L-14 denoting the primer used for amplification and the pooling method employed. This is visually explained in Supplementary figure 2 (below and in the "Main" branch").

### **Supplmentary figure 2 \|** **Sample identifiers and corresponding sample names:**

<<<<<<< Updated upstream:readme.md
![](images/clipboard-285113729.png)
=======
# **Laboratory methods:**
>>>>>>> Stashed changes:Main/readme.md

## **Loading the necessary libraries:**

For this branch the following libraries are required:\

```         
library(readxl)
library(dada2)
library(DECIPHER)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(stringr) 
library(readr)
library(tidyverse)
library(openxlsx)
```

Now that the libraries have been loaded we can proceed into sample processing.

## **Sample processing and assigning taxonomy:**

<<<<<<< Updated upstream:readme.md
This section of the branch involves the assembly of samples - assembling the read library and aligning it to a reference library through assignTaxonomy(). The first step involves loading the sequences into R, by calling path your folder with your sequences. Change "Nemabiome_sequences" to your file location for the sequences. The second piece of code involves identifying the forward and reverse reads, which are labelled here by the end of the file name.
=======
# **Bioinformatic preparations:**
>>>>>>> Stashed changes:Main/readme.md

```         
path <-"<path to sequence files>"

<<<<<<< Updated upstream:readme.md
fwd_files <- sort(list.files(path, pattern = "R1", full.names = TRUE)) 
rev_files <- sort(list.files(path, pattern = "R2", full.names = TRUE))
```
=======
#### **Supplementary figure 3 \| Bioinformatic pooling sample identifiers:**
>>>>>>> Stashed changes:Main/readme.md

Now we place all our names into vector so that the code knows what to look for: Then we apply said vector to the original files. Now our files have their proper names.

```         
# to cutsomise this code jist replace _ with whatever file name differentiating between forward and reverse. This code tells R to ignore the _R1/R2
samples = str_extract(basename(fwd_files), "^[^_]+")


names(fwd_files) <- samples
names(rev_files) <- samples
```

Now we move onto trimming of reads based on the primer binding sites, for this we used the forward NEMA1 primer and the reverse NEMA2 primer as these two primers allowed us to amplify ASVs that were identical for both primer sets (Supplementary figure 4). We then identified these primers in R, to adjust this to your dataset, just replace the fwd_primer/ rev_primer with the primers of interest.

```         
#NEMA1 Primer
fwd_primer <- "ACGTCTGGTTCAGGGTTGTT"
fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer)))

<<<<<<< Updated upstream:readme.md
# NEMA2 Primer
rev_primer <- "ATGCTTAAGTTCAGCGGGTA"
rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))
```
=======
# **Analyses of the data:**
>>>>>>> Stashed changes:Main/readme.md

### **Supplementary figure 4 \| Primer binding sites for the rDNA of *Haemonchus contortus* (ON677958):**

<<<<<<< Updated upstream:readme.md
The red arrows represents the NEMA1 primer and the blue arrows represent the NEMA2 primers. Due to the nested location of the primers we can amplify identical ASVs by trimming for the forward NEMA1 primer and the reverse NEMA2 primer.
=======
# spare graphs - included in body of text:
>>>>>>> Stashed changes:Main/readme.md

![](images/Screenshot 2024-11-12 112120.png)

We then count and trim the primers with the code below.

```         
# this is just a function to identify the number of primers present in the sequences - it does not need changing.

count_primers <- function(primer, filename) {
  num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = FALSE)
  return(sum(num_hits > 0))
}

# this will tell us how many primers are there:
count_primers(fwd_primer, fwd_files[[1]])
count_primers(rev_primer, rev_files[[1]])

cutadapt <- path.expand("<path to your cutadapt location>/cutadapt.exe")

# Make sure it works
system2(cutadapt, args = "--version")
# this will tell you the version and serves as an assurance to know it succesfully loaded. 

# this requires no editing:
cut_dir <- file.path(path, "cutadapt")
if (!dir.exists(cut_dir)) dir.create(cut_dir)

fwd_cut <- file.path(cut_dir, basename(fwd_files))
rev_cut <- file.path(cut_dir, basename(rev_files))

names(fwd_cut) <- samples
names(rev_cut) <- samples

# this just creates log files which is good practice.
cut_logs <- path.expand(file.path(cut_dir, paste0(samples, ".log")))

# this trims the forward and reverse primers and discards any untrimmed reads since they are most likely sequencing errors without primer binding sites!
cutadapt_args <- c("-g", fwd_primer, "-a", rev_primer_rev, 
                   "-G", rev_primer, "-A", fwd_primer_rev,
                   "-n", 2, "--discard-untrimmed")

# this applies cutadapt to our samples - no need for editing!
for (i in seq_along(fwd_files)) {
  system2(cutadapt, 
          args = c(cutadapt_args1,
                   "-o", fwd_cut[i], "-p", rev_cut[i], 
                   fwd_files[i], rev_files[i]),
          stdout = cut_logs[i])  }

# this tells us if we got a result!
head(list.files(cut_dir))
tail(list.files(cut_dir))

dir(path="<path to cutadapt>/cutadapt") 
```

The next step is to do quality control on the reads. Since removing reads with less than 200bp can reduce the error rates @pfeifferSystematicEvaluationError2018, we decided to remove any reads below 200bp, through the function truncLen=c(200, 200). This can be changed to any desired number as long as there is enough base pairs left to produce an overlap.

```         
plotQualityProfile(fwd_cut[3:4]) + ggtitle("Forward")

plotQualityProfile(rev_cut[3:4]) + ggtitle("Reverse")

 # Same as for the cutting we create an output directory to store the filtered files
filt_dir <- file.path(path, "filtered")
if (!dir.exists(filt_dir)) dir.create(filt_dir)

fwd_filt <- file.path(filt_dir, basename(fwd_files))
rev_filt <- file.path(filt_dir, basename(rev_files))

names(fwd_filt) <- samples
names(rev_filt) <- samples

# everything is default apart from the truncLen=c(200, 200) which says that the forward and reverse reads get truncated at 200bp and anything below that is discarded:
filtered_out <- filterAndTrim(
  fwd = fwd_cut, 
  filt = fwd_filt,
  rev = rev_cut,
  filt.rev = rev_filt,
  maxEE = c(2, 5), 
  truncQ = 2, 
  rm.phix = TRUE, 
  compress = TRUE, 
  multithread = FALSE, truncLen=c(200,200)
  )  
```

Then plot the error profiles to double check the data all looks appropriate, then incorporate it into the data:

```         
err_fwd <- learnErrors(fwd_filt, multithread = FALSE)
err_rev <- learnErrors(rev_filt, multithread = FALSE)

plotErrors(err_fwd, nominalQ = TRUE)

dada_fwd <- dada(fwd_filt, err = err_fwd, multithread = FALSE)
dada_rev <- dada(rev_filt, err = err_rev, multithread = FALSE)
```

Now just to assembly the data, there is no need for editing:

```         
mergers <- mergePairs(
  dadaF = dada_fwd,
  dadaR = dada_rev,
  derepF = fwd_filt,
  derepR = rev_filt,
  maxMismatch = 1, 
  verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab) 

seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
dim(seqtab_nochim)

table(nchar(getSequences(seqtab_nochim)))
```

This function returns the number of ASVs:

```         
getN <- function(x) sum(getUniques(x))
```

This provides us with a detailed account of how many reads were lost and at what step!

```         
track <- cbind(
  filtered_out, 
  sapply(dada_fwd, getN), 
  sapply(dada_rev, getN), 
  sapply(mergers, getN), 
  rowSums(seqtab_nochim))
  
colnames(track) <- c("raw", "filtered", "denoised_fwd", "denoised_rev", "merged", "no_chim")
rownames(track) <- samples  
head(track)

write.csv(track, file = "nemabiome_track.csv")
```

Now for the preparation of the dataset:

```         
TSeqtab_nochim <- as.data.frame(t(seqtab_nochim))
TSeqtab_nochim$variant<-1:nrow(TSeqtab_nochim)
Tseqtab_nochim<-as.matrix(Tseqtab_nochim)
print(Tseqtab_nochim)

rownames(Tseqtab_nochim) <- paste0("ASV", 1:nrow(Tseqtab_nochim))
colnames(Tseqtab_nochim) <- paste0("Sample", 1:ncol(Tseqtab_nochim))
Tseqtab_nochim
```

This piece of code assigns the taxonomy/ species/ locus

```         
Taxa <- assignTaxonomy(seqtab_nochim,"C:/Users/tsto3616/Downloads/NEMA_USYD_ITS2_Database_assignTax_v3.fa",taxLevels = c("Locus", "Genus", "Species", "Subsp/Strain"), verbose = TRUE, minBoot = 50, multithread=FALSE, tryRC=TRUE) 
print(Taxa)

tax.augmented <- data.frame(Taxa, (Tseqtab_nochim), stringsAsFactors=FALSE)

write.xlsx(tax.augmented, file="ASVs_raw_masterfiles.xlsx")
```

## **Automated dataframe manipulation for a cleaner dataset!**

THis involves the appropriate variables being assigned to the dataset in seperate columns (eg. pooling, primer, cow) and the cleaning of the ASV labels so that the names are consistent with scientific species names (eg *Haemonchus contortus*).
