#################### UNIVERSIDADE FEDERAL DA PARAÍBA ##########################
################ DEPARTAMENTO DE SISTEMÁTICA E ECOLOGIA #######################
######################     FILTERING SNPs TUTORIAL     ########################

#Jeronymo Dalapicolla, 2023: Guia de genômica de populações aplicada a mamíferos Neotropicais



##I. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT ----
rm(list=ls())

##II. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! ----
# IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H

##III. INSTALL AND LOAD THE PACKAGES ----
#Basic Packages for installation:
if (!require('remotes'))      install.packages('remotes');           library('remotes')
if (!require('BiocManager'))  install.packages('BiocManager');       library('BiocManager')
if (!require('pacman'))       install.packages('pacman');            library('pacman')
if (!require('devtools'))     install.packages('devtools');          library('devtools')

#From Github or BiocManager:
if (!require('LEA'))          BiocManager::install("LEA");           library('LEA')
if (!require('qvalue'))       BiocManager::install("qvalue");        library('qvalue')

#From CRAN R:
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('reshape2'))     install.packages("reshape2");          library('reshape2')



##IV. AUXILIARY FUNCTIONS ----
create_dir = function(dir_names){
  for (i in 1:length(dir_names)){
    if (dir.exists(file.path(getwd(),dir_names[i])) == FALSE) {dir.create(dir_names[i], recursive = T, showWarnings = T) }
    else (print (paste0(dir_names[i], " has already been created. Be careful you can overwrite files")))}
  
}

PlotK <- function(snmfproject){
  
  Mk1 <- mean(cross.entropy(snmfproject, K=1))
  SDk1 <- sd(cross.entropy(snmfproject, K=1))
  Mk2 <- mean(cross.entropy(snmfproject, K=2))
  SDk2 <- sd(cross.entropy(snmfproject, K=2))
  Mk3 <- mean(cross.entropy(snmfproject, K=3))
  SDk3 <- sd(cross.entropy(snmfproject, K=3))
  Mk4 <- mean(cross.entropy(snmfproject, K=4))
  SDk4 <- sd(cross.entropy(snmfproject, K=4))
  Mk5 <- mean(cross.entropy(snmfproject, K=5))
  SDk5 <- sd(cross.entropy(snmfproject, K=5))
  Mk6 <- mean(cross.entropy(snmfproject, K=6))
  SDk6 <- sd(cross.entropy(snmfproject, K=6))
  Mk7 <- mean(cross.entropy(snmfproject, K=7))
  SDk7 <- sd(cross.entropy(snmfproject, K=7))
  Mk8 <- mean(cross.entropy(snmfproject, K=8))
  SDk8 <- sd(cross.entropy(snmfproject, K=8))
  Mk9 <- mean(cross.entropy(snmfproject, K=9))
  SDk9 <- sd(cross.entropy(snmfproject, K=9))
  Mk10 <- mean(cross.entropy(snmfproject, K=10))
  SDk10 <- sd(cross.entropy(snmfproject, K=10))
  
  CE <- data.frame(K=c(1:10), Mean = c(Mk1,Mk2,Mk3,Mk4,Mk5,Mk6,Mk7,Mk8,Mk9,Mk10),
                   SD = c(SDk1,SDk2,SDk3,SDk4,SDk5,SDk6,SDk7,SDk8,SDk9,SDk10))
  
  library(ggplot2)
  
  ggplot(CE, aes(x=K, y=Mean)) + 
    geom_segment(aes(x=K[1],y= Mean[1], xend =K[10] , yend =Mean[10]), size = 0.5, linetype = 'dashed')+
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.2)+
    geom_line() + 
    geom_point(size=4, shape=21, fill="red", color="darkred") + xlab("Number of ancestral populations") + ylab("Cross-entropy")+
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, vjust=0.5), axis.text = element_text(size=15, face ="bold" , color = "black"), axis.title.x = element_text(size=15, face="bold", color="black"),axis.title.y = element_text(size=15, face="bold", color="black")) +
    scale_x_continuous(breaks = seq(0,10, by=2))
}

Best.run <- function(nrep, optimalK, p1, p2, p3, p4, p5, p6){
  ce1 = LEA::cross.entropy(p1, K = optimalK) # get the cross-entropy of each run for optimal K
  ce2 = LEA::cross.entropy(p2, K = optimalK) # get the cross-entropy of each run for optimal K
  ce3 = LEA::cross.entropy(p3, K = optimalK) # get the cross-entropy of each run for optimal K
  ce4 = LEA::cross.entropy(p4, K = optimalK) # get the cross-entropy of each run for optimal K
  ce5 = LEA::cross.entropy(p5, K = optimalK) # get the cross-entropy of each run for optimal K
  ce6 = LEA::cross.entropy(p6, K = optimalK) # get the cross-entropy of each run for optimal K
  AllProjects <- rbind(ce1, ce2, ce3, ce4, ce5, ce6)
  rownames(AllProjects) <- NULL
  AllProjects <- as.data.frame(AllProjects)
  AllProjects$Project <- c(rep(1, nrep), rep(2,nrep), rep(3, nrep), rep(4, nrep), rep(5, nrep), rep(6, nrep))
  Best_project <- AllProjects[AllProjects[, 1]==min(AllProjects[, 1]), 2]
  Best_runs <- AllProjects[AllProjects$Project==Best_project, ]
  Best_runs$Nrun <- 1:nrow(Best_runs)
  Best_run <- Best_runs[Best_runs[, 1]==min(Best_runs[, 1]), 3]
  print(paste0("Best run is: ", "project = ", Best_project, " run = ", Best_run)) 
}



##V. CONVERT FILES USING PLINK AND TERMINAL ----
#Install PLINK: https://www.cog-genomics.org/plink/
#Open the terminal and move to PLINK folder:
#cd PLINK/
#./plink --vcf [INPUT_NAME.vcf] --recode --out [OUTPUT_NAME] --allow-extra-chr

#Convert PED to LFMM and GENO to use them inside R:
#lfmm_input = ped2lfmm("Inputs/ped/INPUT_NAME.ped", "Inputs/lfmm/INPUT_NAME.lfmm", force = T)
#geno_input = ped2geno("Inputs/ped/INPUT_NAME.ped", output.file = "Inputs/geno/INPUT_NAME.geno", force = T)



###1. CONVERT FILES IN R ----
#A. Load VCF in R
vcf = read.vcfR("Inputs/vcf/steerei_Neutral.vcf")
#check file
vcf

#B. Convert VCF to Genotypes
genotypes= t(extract.gt(vcf, element = "GT", mask = FALSE, as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = TRUE))
#check:
genotypes[1:5,1:5]
row.names(genotypes)

#C. Convert Genotypes to LFMM:
dim(genotypes)
((sum(is.na(genotypes)))/(dim(genotypes)[1]*dim(genotypes)[2]))*100
#5.17%  #Amount of missing data

genotypes[is.na(genotypes)] <- 9 #The missing genotypes have to be encoded with the value 9
genotypes[1:10,1:10]
write.lfmm(genotypes,"Inputs/lfmm/steerei_Neutral.lfmm")


lfmm_input = read.lfmm("Inputs/lfmm/steerei_Neutral.lfmm")
geno_input = lfmm2geno("Inputs/lfmm/steerei_Neutral.lfmm", output.file = "Inputs/geno/steerei_Neutral.geno", force = TRUE)



###### 2. RUN ANALYSIS:
#B. Create folders for alpha values and copy .geno object in each folder:
alpha_values = c(10, 100, 500, 1000, 2000, 4000)
geno_input
for (i in alpha_values){
  path = paste0("./Results/sNMF/Alpha", i, "n")
  if (dir.exists(file.path(getwd(), path)) == FALSE)
  {dir.create(path, recursive = T, showWarnings = F)} else (print (paste0(path, " has already been created. Be careful with overwritting")))
  file.copy(geno_input, path)
}


#C. Set parameters to run SNMF (LEA) using different alpha.
K = c(1:10) # K to be tested
replications = 10 # number of replication by K
ploidy = 2 # species ploidy
CPU = 4 #Number of cores


#D. Run sNMF (LEA) using different alpha.
set.seed(123)
geno_name = "steerei_Neutral"
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("Results/sNMF/Alpha", i,"n/", geno_name, ".geno")
  pro_snmf = snmf(path, K = K, rep = replications, alpha = i, entropy = T, ploidy = ploidy , project = "new", CPU= CPU)
  assign(paste0("project_snmf", loop, "n"), pro_snmf)
}


#E. To load the SNMF projects in a new R session (after quitting R).
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("Results/sNMF/Alpha", i,"n/", geno_name, ".snmfProject")
  pro_snmf = load.snmfProject(path)
  assign(paste0("project", loop, "n"), pro_snmf)
}

#F. summary of the project
summary(project1n)
summary(project2n)
summary(project3n)
summary(project4n)
summary(project5n)
summary(project6n)

#G. View Cross-Entropy plots
PlotK(project1n) 
PlotK(project2n) 
PlotK(project3n) 
PlotK(project4n) 
PlotK(project5n) 
PlotK(project6n) 


#H. Save graphs of cross-Entropy plot with standard deviation error bars
for (i in alpha_values){
  pdf(paste0("./Results/sNMF/Cross_Entropy_sNMF_Alpha_",  i, "n.pdf"), onefile = F)
  path = paste0("Results/sNMF/Alpha", i,"n/",  geno_name, ".snmfProject")
  print(PlotK(load.snmfProject(path)))
  dev.off()
}


#I. Select optimal K value. If K = 1, so you have only one genetic cluster
optimal_K = 3 #We will use 3 to increase number of groups
optimal_K = 2


#J. Select best run (lowest cross-entropy)
best_run = Best.run(nrep=10, optimalK=optimal_K, p1=project1n, p2=project2n, p3=project3n, p4=project4n, p5=project5n, p6=project6n)
#load the best project
best_run_split = scan(text = best_run, what = "")
path_best_run = paste0("Results/sNMF/Alpha", alpha_values[as.numeric(best_run_split[6])],"n/",  geno_name, ".snmfProject")
#set the values
project = load.snmfProject(path_best_run)
run=as.numeric(best_run_split[9])



#L. Add admixture coefficient and replace the population ID to vcf file
#create cols names
coeff = c()
for (i in 1:optimal_K){coeff[i] = paste0("Adx_Coeff_", i)}
coeff

Qmat =
  project %>%
  Q(., run=run, K=optimal_K) %>%
  as.data.frame() %>%
  setNames(coeff) %>%
  mutate(PopID_snmf = apply(., 1, which.max))
head(Qmat)


#C. Load information on samples
#order of samples in vcf
sample_name_vcf = colnames(vcf@gt)[-1]

geo_data = read.csv("Inputs/metadata/coord_steerei.csv") 
head(geo_data)
tail(geo_data)

#select rows in the same order than vcf
genomic_order = as.data.frame(sample_name_vcf)
names(genomic_order) = "Sample_ID"

geo_data_reorder = semi_join(geo_data, genomic_order)
geo_data_reorder = inner_join(genomic_order, geo_data_reorder)
head(geo_data_reorder)

#verify the order:
identical(as.character(sample_name_vcf), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(sample_name_vcf), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(geo_data_reorder$Sample_ID), as.character(sample_name_vcf))

#data frame final:
geo_data_final = cbind(geo_data_reorder, Qmat)
head(geo_data_final)

#H. Structure-like barplot for the Q-matrix
#define colors and name for populations, in this case according to Monteiro et al. 2021:
geo_data_final[geo_data_final$PopID_snmf == 1,] 
geo_data_final[geo_data_final$PopID_snmf == 2,]
geo_data_final[geo_data_final$PopID_snmf == 3,] 

colors_snmf = c('red', 'blue', "yellow") #color for the genetic cluster
labels_snmf = c("POP1 - Madeira River", "POP2 - Bolivia", "POP3 - Juruá River") #name for the genetic cluster



#create a dataframe (df) for ggplot2, all localities
coeff
#col 'pop' is my Locality information
df = geo_data_final %>%
  dplyr::select(Sample_ID, pop, Adx_Coeff_1, Adx_Coeff_2, Adx_Coeff_3) %>% #edit for the number of coeff you found!
  dplyr::arrange(pop) %>%
  melt(., id.vars=c("Sample_ID", "pop"))
head(df)

#define sample size:
n_sample = 16

#define graph
p = ggplot(data=df, mapping = aes(x=factor(Sample_ID),factor(pop),
                                  y= value*100,
                                  fill = factor(variable))) +
  geom_bar(stat="identity", width = 1, size=0.3, colour="black") +
  scale_fill_manual(values= colors_snmf,
                    labels= labels_snmf) +
  theme_minimal() +
  guides(fill=guide_legend(title="Populations"))+
  scale_x_discrete(limits=df$sample_name[1:n_sample], labels =df$pop[1:n_sample], guide = guide_axis(n.dodge=1)) + #guide_axis(n.dodge=1) number of lines in labels
  xlab("") + ylab("")

#check
p


pdf(paste0("./Results/sNMF/sNMF_pipegarph_", geno_name, ".pdf"), onefile =F)
plot(p)
dev.off()

write.csv(geo_data_final, paste0("Results/sNMF/ancestry_coef_snmf_", geno_name, ".csv"), quote = F) #save result as table

#END
