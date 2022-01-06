

library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling
library(ape)
library(vegan)


## LOAD DATA FROM EDO PC
setwd("...")
data.directory = ("...")
treefile_p1 <- read.tree("../humanmicrobiome.tree")

K = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)

MaxRep = 100 #100

alpha = 0.01

#Simulated data differences
simul.perc = c(1,2,5,10,25,50,75)

#Loop on different data sets
for(d in 1:length(simul.perc)){
  
  # Data file in
  data_file_name = sprintf('phylo_simulated_data_deletion_9Feb2021_%i.Rdata',simul.perc[d])
  
  file_in =  file.path(data.directory,data_file_name)
  
  print(file_in)
  
  # results file out
  file_name_out = sprintf('Power_Beta_simulation1_%i.Rdata',simul.perc[d])
  
  #load data
  load(file_in)
  
  ##### LOAD Original 
  ps.ng.tax <- read_phyloseq(otu.file = "humanmicrobiome.biom", 
                             taxonomy.file = NULL, 
                             metadata.file = "metadata_table.csv", 
                             type = "biom") 
  
  treefile_p1 <- read.tree("humanmicrobiome.tree")
  ps.ng.tax <- merge_phyloseq(ps.ng.tax, treefile_p1)
  ps.ng.tax 
  
  ## remove NA
  ps3.subset <- subset_taxa(ps.ng.tax, Domain !="NA")  #remove no bacterial DNA  #6249-6090 taxa 
  ps3.subset   # 4704 taxa 474 samples       
  ps3 <- ps3.subset 
  
  ## extract data 
  meta <- meta(ps3)
  ps.otu <- as.data.frame(ps3@otu_table)
  ps.tree <- ps3@phy_tree
  
  ## subset data e.g. gut
  ps.ng.tax_GUT <- subset_samples(ps3, scientific_name == "human gut metagenome" )
  ps.ng.tax_GUT <- prune_taxa(taxa_sums(ps.ng.tax_GUT) > 0, ps.ng.tax_GUT)
  ps.ng.tax_GUT
  
  ### original
  ps.otu <- as.data.frame(ps.ng.tax_GUT@otu_table)
  ps.tree <- ps.ng.tax_GUT@phy_tree
  tax <- as.data.frame(ps.ng.tax_GUT@tax_table)
  
  tax <- ps.ng.tax_GUT@tax_table
  ps.otu <- ps.ng.tax_GUT@otu_table
  
  colnames(ps.otu) <- paste0("Sample", 1:ncol(ps.otu))
  #head(ps.otu)
  
  OTU = otu_table(ps.otu, taxa_are_rows = TRUE)
  TAX = tax_table(tax)
  
  physeq = phyloseq(OTU, TAX, ps.tree)
  #physeq = phyloseq(OTU, TAX, sampledata)
  
  ### SIM data
  taxS <- ps.ng.tax_GUT.simulated@tax_table
  ps.otuS <- ps.ng.tax_GUT.simulated@otu_table
  
  OTUS = otu_table(ps.otuS, taxa_are_rows = TRUE)
  TAXS = tax_table(taxS)
  
  physeq_Sim = phyloseq(OTUS, TAXS, ps.tree)
  
  physeqM <- merge_phyloseq(physeq, physeq_Sim)
  
  
  sampledataSim = sample_data(data.frame(
    Type = rep(c("Original","Simulated"), each=169),
    Depth = 1:338,
    row.names=sample_names(physeqM),
    stringsAsFactors=FALSE
  ))
  sampledataSim
  
  
  physeqTotal2 <- merge_phyloseq(physeqM, sampledataSim)
  
  
  #Can also store all POWER calculation in one matrix
  POWER = matrix(NA, nrow = 4, ncol = length(K))
  
  rownames(POWER) = c('Bray-Curtis', 'Jaccard','UnWeigthed Unifrac', 'Weigthed Unifrac')
  colnames(POWER) = K
  
  
  for(j in 1 : length(K)){
    
    k = K[j]
    
    PVAL = matrix(NA, nrow = 4, ncol = MaxRep )
    rownames(PVAL) = c('Bray-Curtis', 'Jaccard','UnWeigthed Unifrac', 'Weigthed Unifrac')
    colnames(PVAL) = c(1:MaxRep)
    
    
    for(i in 1 : MaxRep){
      
      print(sprintf('Sample size k = %i ; Sampling %i of %i',k,i,MaxRep))
      
      samp1 = sample(1:169,k) #select k sample from the controls
      samp2 = sample(170:338,k) #select k sample from the simulated
      
      
      ## Make temporary metadata and add a new metadata item
      sampling.labels = rep('OUT', 338); #by default all samples are out
      sampling.labels[c(samp1,samp2)] =  'IN' #Asssign IN labels to the k + k samples selected to be in
      
      sampledataSim.k = sample_data(data.frame(
        Type = rep(c("Original","Simulated"), each=169),
        Depth = c(1:338),
        Sampling = sampling.labels,
        row.names=sample_names(physeqM),
        stringsAsFactors=FALSE
      ))
      
      
      #### Merge dataset with metadata containing info about the samples to be selected
      physeqTotal2.k <- merge_phyloseq(physeqM, sampledataSim.k)
      
      #Subsampling
      phylo.k.all = subset_samples(physeqTotal2.k,Sampling == 'IN') #select the k + k samples 
      
      ### remove zeros
      phylo.k <- prune_taxa(taxa_sums(phylo.k.all) > 0, phylo.k.all)
      
      ## Calculate distances
      dist.b <- distance(phylo.k, method = "bray")   
      dist.j <- distance(phylo.k, method = "jaccard", binary=TRUE)   
      dist.uf <- distance(phylo.k, method = "unifrac")   
      dist.wuf <- distance(phylo.k, method = "wunifrac")   
      
      x.b = as.matrix(dist.b)  
      x.j = as.matrix(dist.j)  
      x.uf = as.matrix(dist.uf)  
      x.wuf = as.matrix(dist.wuf)  
      
      y.b = as.dist(as.matrix(dist.b) )
      y.j = as.dist(as.matrix(dist.j) )
      y.uf = as.dist(as.matrix(dist.uf) )
      y.wuf = as.dist(as.matrix(dist.wuf) )
      
      # Make lables of the sample group
      Groups = rep(c('Contols','Simulated'),each = k)
      
      # Make temporary metadata to pass to adonis
      MMM <- data.frame(Groups)
      
      RR.b <- adonis2(y.b ~ Groups, data = MMM , perm=9999)
      RR.j <- adonis2(y.j ~ Groups, data = MMM , perm=9999) 
      RR.uf <- adonis2(y.uf ~ Groups, data = MMM , perm=9999) 
      RR.wuf <- adonis2(y.wuf ~Groups, data = MMM , perm=9999) 
      
      
      PVAL[1,i] = RR.b$`Pr(>F)`[1]
      PVAL[2,i] = RR.j$`Pr(>F)`[1]
      PVAL[3,i] = RR.uf$`Pr(>F)`[1]
      PVAL[4,i] = RR.wuf$`Pr(>F)`[1]
      
    }
    
    #Calculate POWER
    
    POWER[1,j] = length(which(PVAL[1,] <  alpha))/MaxRep*100
    POWER[2,j] = length(which(PVAL[2,] <  alpha))/MaxRep*100
    POWER[3,j] = length(which(PVAL[3,] <  alpha))/MaxRep*100
    POWER[4,j] = length(which(PVAL[4,] <  alpha))/MaxRep*100
    
  }
  
  save(POWER, file = file_name_out)
}
