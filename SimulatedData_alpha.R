

library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling
library(ape)
library(vegan)
library(picante)

## LOAD DATA
data.directory = ('..')
treefile_p1 <- read.tree('..')
load('..')

ps.ng.tax_GUT.simulated.ctrl <- merge_phyloseq(ps.ng.tax_GUT.simulated.ctrl, treefile_p1)

### Merge the ps objects and add labels
ps.otu <- as.data.frame(ps.ng.tax_GUT.simulated.ctrl@otu_table)
ps.tree <- ps.ng.tax_GUT.simulated.ctrl@phy_tree
tax <- as.data.frame(ps.ng.tax_GUT.simulated.ctrl@tax_table)

tax <- ps.ng.tax_GUT.simulated.ctrl@tax_table
ps.otu <- ps.ng.tax_GUT.simulated.ctrl@otu_table

colnames(ps.otu) <- paste0("Sample", 1:ncol(ps.otu))

OTU = otu_table(ps.otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)

physeq = phyloseq(OTU, TAX, ps.tree)



#Simulated data differences
#simul.perc = c(1,2,5,10,15,20,25)
#simul.perc = c(15,20,25)
#simul.perc = c(2,5,10, 15, 20)
simul.perc = c(2,5)


#K = c(5:10,15, 20, 25, 30) #I guess this is enough, need to based on the other calucaltions we have done

#depend on the data set
K = c(5,10,15,20,25,30,35,40,45,50)

MaxRep = 100 #2000

alpha = 0.01

#Loop on different data sets
for(d in 1:length(simul.perc)){
  
  # Data file in
  data_file_name = sprintf('phylo_simulated_simulated_data_reduction_16March2021_%i.Rdata',simul.perc[d])
  
  file_in =  file.path(data.directory,data_file_name)
  
  print(file_in)
  
  # results file out
  file_name_out = sprintf('Power_Alpha_simulation2_%i.Rdata',simul.perc[d])
  #load data
  load(file_in)

  ### SIM data
  taxS <- ps.ng.tax_GUT.simulated@tax_table
  ps.otuS <- ps.ng.tax_GUT.simulated@otu_table
  
  OTUS = otu_table(ps.otuS, taxa_are_rows = TRUE)
  TAXS = tax_table(taxS)
  
  physeq_Sim = phyloseq(OTUS, TAXS, ps.tree)
  
  physeqM <- merge_phyloseq(physeq, physeq_Sim)
  
  ## label 
  sampledataSim = sample_data(data.frame(
    Type = rep(c("Original","Simulated"), each=169),
    Depth = c(1:338),
    Sampling = rep('IN', 338),
    row.names=sample_names(physeqM),
    stringsAsFactors=FALSE
  ))
  
  #### Merge datasets
  physeqTotal2 <- merge_phyloseq(physeqM, sampledataSim)
  
  ############ CAL BETA
  
  metadata <- as(sample_data(physeqTotal2), "data.frame")
  
  set.seed(123456)
  
  #Can also store all POWER calculation in one matrix
  POWER = matrix(NA, nrow = 4, ncol = length(K))
  
  rownames(POWER) = c('Chao1', 'Shannon','Simpson', 'Phylogetic D')
  colnames(POWER) = K
  
  
  for(j in 1 : length(K)){
    
    k = K[j]
    
    PVAL = matrix(NA, nrow = 4, ncol = MaxRep )
    rownames(PVAL) = c('Chao1', 'Shannon','Simpson', 'Phylogetic D')
    colnames(PVAL) = c(1:MaxRep)
    
    
    for(i in 1 : MaxRep){
      
      print(sprintf('Data set %i of %i (perc = %i): Sample size k = %i ; Sampling %i of %i',d,length(simul.perc),simul.perc[d],k,i,MaxRep))
      
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
      
      ps.otu.k <- as.data.frame(phylo.k@otu_table)
      ps.tree.k <- phylo.k@phy_tree
      meta.k <- meta(phylo.k)
      
      ## Calculate ALPHA
      
      PD <- pd(t(ps.otu.k), ps.tree.k,include.root=T) 
      y.pd.k = PD$PD
      
      ps1.adiv <- estimate_richness(phylo.k, measures = c("Chao1", "Shannon", "Simpson"))
      ps1.metadata <- as(sample_data(phylo.k), "data.frame")
      
      y.Shannon.k = ps1.adiv$Shannon  
      y.Chao1.k  = ps1.adiv$Chao1  
      y.Simpson.k  = ps1.adiv$Simpson 
      
      # Make lables of the sample group
      Groups = rep(c('Contols','Simulated'),each = k)
      
      # test
      RR.Shannon <- kruskal.test(y.Shannon.k ~ Groups)
      RR.Chao1 <- kruskal.test(y.Chao1.k ~ Groups)
      RR.Simpson <- kruskal.test(y.Simpson.k ~ Groups)
      RR.PD <- kruskal.test(y.pd.k ~ Groups)
      
      PVAL[1,i] = RR.Chao1$p.value
      PVAL[2,i] = RR.Shannon$p.value
      PVAL[3,i] = RR.Simpson$p.value
      PVAL[4,i] = RR.PD$p.value
      
    }
    
    #Calculate POWER
    
    POWER[1,j] = length(which(PVAL[1,] <  alpha))/MaxRep*100
    POWER[2,j] = length(which(PVAL[2,] <  alpha))/MaxRep*100
    POWER[3,j] = length(which(PVAL[3,] <  alpha))/MaxRep*100
    POWER[4,j] = length(which(PVAL[4,] <  alpha))/MaxRep*100
    
  }
  
  save(POWER, file = file_name_out)
  
}
