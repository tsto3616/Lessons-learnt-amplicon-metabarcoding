------------------------------------------------------------------------

title: "Sample assembly"

author: "Thomas Stocker"

date: "11-11-2024"

#if needed: bibliography: "references.bib"

------------------------------------------------------------------------

# **Welcome to the Sample assembly branch**

This branch deals with the assembly of the samples but also the automated renaming within R for both the species and or samples - this greatly saves time - instead of manual processing. But before we get into any coding, here is the brief description of the samples.

## **Samples:**

As discussed in the "Main" branch, the samples are labeled EF\*\*\*L-1 to EF\*\*\*L-14, with EF\*\*\* representing the cow's identifier, and the L-1 to L-14 denoting the primer used for amplification and the pooling method employed. This is visually explained in Supplementary figure 2 (below and in the "Main" branch").

### **Supplmentary figure 2 \|** **Sample identifiers and corresponding sample names:**

![](images/clipboard-2518261799.png)

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
```

Now that the libraries have been loaded we can proceed into sample processing.

## **Sample processing and assigning taxonomy:**

This section of the branch involves the assembly of samples - assembling the read library and aligning it to a reference library through assignTaxonomy(). The first step involves loading the sequences into R, by calling path your folder with your sequences. Change "Nemabiome_sequences" to your file location for the sequences. The second piece of code involves identifying the forward and reverse reads, which are labelled here by the end of the file name. For each sample there is a

```         
path <-"Nemabiome_sequences"

fwd_files <- sort(list.files(path, pattern = "R1", full.names = TRUE)) 
rev_files <- sort(list.files(path, pattern = "R2", full.names = TRUE))
```
