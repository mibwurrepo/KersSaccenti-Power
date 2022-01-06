
### HMP DATA
### https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/introduction.html 

library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling
library(ape)
library(vegan)

#setwd("")

## Upload data
ps.ng.tax <- read_phyloseq(otu.file = "humanmicrobiome.biom", 
                           taxonomy.file = NULL, 
                           metadata.file = "metadata_table.csv", 
                           type = "biom") 
treefile_p1 <- read.tree("humanmicrobiome.tree")
ps.ng.tax <- merge_phyloseq(ps.ng.tax, treefile_p1)
ps.ng.tax

## Remove Domain NA
ps.ng.tax <- subset_taxa(ps.ng.tax, Domain !="NA")  
ps.ng.tax   

## Subset data

## ran every environment separate
ps.ng.tax_GUTSKIN <- subset_samples(ps.ng.tax, scientific_name == "human gut metagenome" |scientific_name == "human skin metagenome"  ) 

#sum(meta$scientific_name == 'human vaginal metagenome')  #86
#sum(meta$scientific_name == 'human gut metagenome') #168
#sum(meta$scientific_name == 'human skin metagenome') #69 
#sum(meta$scientific_name == 'human oral metagenome') #150

## prune taxa
ps.ng.tax_GUTSKIN <- prune_taxa(taxa_sums(ps.ng.tax_GUTSKIN) > 0, ps.ng.tax_GUTSKIN)
ps.ng.tax_GUTSKIN

####################################
### calculate alpha 
#################################### 

ps1.adiv <- estimate_richness(ps.ng.tax_GUTSKIN, measures = c("Chao1", "Shannon", "Simpson"))
ps1.metadata <- as(sample_data(ps.ng.tax_GUTSKIN), "data.frame")

## calculate power 
set.seed(12345)
ps1.metadata$Shannon <- ps1.adiv$Shannon

y = ps1.metadata$Shannon 
group = ps1.metadata$scientific_name

group1 = 'human gut metagenome'
group2 = 'human skin metagenome'

idx1 = which(group==group1)
idx2 = which(group==group2)

COUNT = NULL;
K = 5:69

MaxRep =  2000
for(j in 1 : length(K)){
  
  k = K[j]
  
  PVAL = NULL
  X2 = NULL
  for(i in 1 : MaxRep){
    samp1 = sample(length(idx1),k)
    samp2 = sample(length(idx2),k)
    
    idx1_samp = idx1[samp1]
    idx2_samp = idx2[samp2]
    
    idx_samp = cbind(idx1_samp, idx2_samp)
    
    kruskal.observed_type <- kruskal.test(y[idx_samp] ~ group[idx_samp])
    
    PVAL[i] = kruskal.observed_type$p.value
    X2[i] = kruskal.observed_type$statistic
    
  }
  
  COUNT[j] = length(which(PVAL < 0.01))/MaxRep*100
  
  
}

COUNT 

####################################
### Calculate beta 
####################################

## ran distrance metrics seperate 
dist.bc <- phyloseq::distance(ps.ng.tax_GUTSKIN, method = "bray")  
#dist.j <- phyloseq::distance(ps.ng.tax_GUTSKIN, method = "jaccard", binomal=TRUE) 
#dist.uf <- phyloseq::distance(ps.ng.tax_GUTSKIN, method = "unifrac") 
#dist.wuf <- phyloseq::distance(ps.ng.tax_GUTSKIN, method = "wunifrac")

metadata <- as (sample_data(ps.ng.tax_GUTSKIN), "data.frame")
#metadata <- metadata[order(metadata$sample_type),]

x = as.matrix(dist.bc)   ##Change Distance!! 
y = as.dist(x)

#adonis2(dist.bc~ scientific_name, data = metadata, perm=9999)
adonis2(y ~ scientific_name, data = metadata, perm=9999) 

DistMat0 = as.matrix(y)

idx1 = 1:69  #which(group==group1)    # 1:69, adapt to dataset
idx2 = 70:138 #which(group==group2)   # 70:138

COUNT = NULL;

K = c(5,15,20,25,30,35)  # adapt to dataset
MaxRep = 10  #100 
for(j in 1 : length(K)){
  
  k = K[j]
  
  PVAL = NULL
  X2 = NULL
  for(i in 1 : MaxRep){
    
    samp1 = sample(length(idx1),k)
    samp2 = sample(length(idx2),k)
    
    idx1_samp = idx1[samp1]
    idx2_samp = idx2[samp2]
    
    idx_samp = c(idx1_samp, idx2_samp)
    
    IN = DistMat0
    DistMat = IN[-idx_samp, -idx_samp]
    
    QQ <- as.dist(DistMat)
    print(dim(DistMat))
    MM = metadata
    MMM = MM[-idx_samp,]
    RR <- adonis2(QQ ~ scientific_name, data = MMM , perm=9999) 
    
    PVAL[i] = unlist(RR)[16]  #p-value
    X2[i] = unlist(RR)[13]   #R2
    
  }
  
  COUNT[j] = length(which(PVAL < 0.01))/MaxRep*100
  
  
}

COUNT

