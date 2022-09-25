####################### VALE INSTITUTE OF TECHNOLOGY ##########################

############################     DAPC TUTORIAL     #############################

#By Jeronymo Dalapicolla, 2022


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
if (!require('LEA'))          BiocManager::install("LEA");            library('LEA')
if (!require('qvalue'))       BiocManager::install("qvalue");         library('qvalue')

#From CRAN R:
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require("SNPfiltR"))     install.packages("SNPfiltR");          library("SNPfiltR")
if (!require("dartR"))        install.packages("dartR");             library("dartR")
if (!require('pcadapt'))      install.packages("pcadapt");           library('pcadapt')
if (!require('sf'))           install.packages("sf");                library('sf')


####1. LOAD FILES ----
#A. VCF
snps_raw = read.vcfR("Inputs/vcf/dioclea.vcf", verbose = T)
snps_raw

#B. Remove individuals before filtering SNPs if you need. In this case we have 1 repeated sample: "CC354_rep_sorted" and we do not have geographical data for "CC087_sorted".
colnames(snps_raw@gt)
remove_ind = c("CC354_rep_sorted",  "CC087_sorted")
ind_to_keep = colnames(snps_raw@gt)[!colnames(snps_raw@gt) %in% remove_ind]

snps_raw = snps_raw[samples=ind_to_keep]
snps_raw


#C. Format a PopMAP. A data frame with two column popmap with the same format as stacks, and the columns must be named 'id' and 'pop'
# Load information on samples
geo_data = read.csv("Inputs/metadata/df_dioclea.csv") 
head(geo_data)
tail(geo_data)

#select rows in the same order than vcf
genomic_order = as.data.frame(colnames(snps_raw@gt)[-1]) #always 1st line will be FORMAT. You need to remove
names(genomic_order) = "sample_name"

geo_data_reorder = semi_join(geo_data, genomic_order)
geo_data_reorder = inner_join(genomic_order, geo_data_reorder)
head(geo_data_reorder)

#verify the order:
identical(as.character(colnames(snps_raw@gt)[-1]), as.character(geo_data_reorder$sample_name))
setdiff(as.character(colnames(snps_raw@gt)[-1]), as.character(geo_data_reorder$sample_name))
setdiff(as.character(geo_data_reorder$sample_name), as.character(colnames(snps_raw@gt)[-1]))

popmap = geo_data_reorder[,c(1,4)] #canga as population/location
names(popmap) = c("id", "pop")
class(popmap)
head(popmap)






###2. FILTERING ----
#A. Remove indels
snps_F1 = extract.indels(snps_raw, return.indels = FALSE)


#B. Keep Biallelic SNPs
snps_F2 = filter_biallelic(snps_F1)
snps_F2


#C. By Quality
QUAL = 30 #set a threshold
snps_to_keep = which(snps_F2@fix[,6] > QUAL)
snps_F3 = snps_F2[snps_to_keep, ]
snps_F3


#D. By Minimum Depth by site
#Usually you do not have information on Genotype quality, if you have you can use gq filter
snps_F4 = hard_filter(vcfR = snps_F3, depth = 20, gq = NULL)
snps_F4


#E. By Maximum Depth by site
snps_F5 = max_depth(snps_F4, maxdepth = 200)
snps_F5


#F. By Allele Balance:
snps_F6 = filter_allele_balance(snps_F5)
snps_F6


#G. by MAC or MAF
#estimate MAC by MAF:
MAF_1 = (length(colnames(snps_F6@gt)[-1])*1)/100
MAF_3 = (length(colnames(snps_F6@gt)[-1])*3)/100
MAF_5 = (length(colnames(snps_F6@gt)[-1])*5)/100
MAF_10 = (length(colnames(snps_F6@gt)[-1])*10)/100

snps_F7 = min_mac(snps_F6, min.mac = MAF_5)
snps_F7


#H. By Missing Data per site
missing_by_snp(vcfR=snps_F7) #verify statistics

snps_F8 = missing_by_snp(vcfR=snps_F7, cutoff = .7)
snps_F8


#H. By Missing Data per sample
missing_by_sample(vcfR=snps_F7, popmap = popmap)  #verify statistics
snps_F9 = missing_by_sample(vcfR=snps_F8, cutoff = .4)
snps_F9

popmap = popmap[popmap$id %in% colnames(snps_F9@gt),]
head(popmap)



#I. By Linkage Disequilibrium (LD). If you will use NeEstimate do not use this step, jump to #J.
snps_F10 = distance_thin(vcfR = snps_F9, min.distance = 150) #1 SNP by contig/Chrom
snps_F10

#For NeEstimate
snps_F10 = snps_F9


#J. By Hardy Weinberg Equilibrium (HWE)
#convert to genlight
gl_vcf = vcfR2genlight(snps_F10)
gl_vcf
#add pop/local information
pop(gl_vcf) = popmap$pop
popNames(gl_vcf)

#filter by in the genlight as recommended by Pearman et al 2022:
gl_vcf_fil = gl.filter.hwe(
  gl_vcf,
  subset = "each",
  n.pop.threshold = round(length(popNames(gl_vcf))/2,0), # most of pops
  method_sig = "Exact", # Wigginton et al. (2005)
  multi_comp = FALSE,
  multi_comp_method = "BH", # Benjamini & Hochberg, 1995
  alpha_val = 0.05,
  pvalue_type = "midp", # Graffelman & Moreno, 2013
  cc_val = 0.5,
  min_sample_size = 5,
  verbose = NULL
)


#keep in vcf the same snps in the gl filtered object.
filtered_snps = separate(as.data.frame(gl_vcf_fil@loc.names), col = 1, into = c("CHROM1", "CHROM2", "SNP", "POS"), sep = "_" )
head(filtered_snps)
snps_F11= snps_F10[as.integer(filtered_snps$POS), ]
snps_F11



#J. By outliers SNPs - possibly under selection (Adapt)
#convert VCF to Genotypes
genotypes= t(extract.gt(snps_F11, element = "GT", mask = FALSE, as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = TRUE))
#check:
genotypes[1:5,1:5]
row.names(genotypes)

#Convert Genotypes to LFMM:
dim(genotypes)
((sum(is.na(genotypes)))/(dim(genotypes)[1]*dim(genotypes)[2]))*100
#13.25316%  #Amount of missing data

genotypes[is.na(genotypes)] <- 9 #The missing genotypes have to be encoded with the value 9
genotypes[1:10,1:10]
write.lfmm(genotypes,"Inputs/lfmm/F11.lfmm")

#Load LFMM:
lfmm_input = read.lfmm("Inputs/lfmm/F11.lfmm")
class(lfmm_input)

#Convert LFMM to a PCadapt matrix and run the analyses
pcadapt.test = lfmm_input %>%
  read.pcadapt(., type="lfmm") %>%
  pcadapt(., K=10, ploidy=2)

#"Cattell's Rule" for interpreting the scree plot (PC to the left of the flat line)
plot(pcadapt.test, option="screeplot") #2 PCs
plot(pcadapt.test, option="scores") #PC1 and PC2 #3 clusters
plot(pcadapt.test, option = "scores", i = 3, j = 2) #PC2 and PC3 # 3 clusters
plot(pcadapt.test, option = "scores", i = 3, j = 4) #PC3 and PC4 # 1 cluster
plot(pcadapt.test, option = "scores", i = 5, j = 4) #PC5 and PC4 # 1 cluster
plot(pcadapt.test, option = "scores", i = 5, j = 6) #PC5 and PC6 # 1 cluster
#With PC3 
#So we can use 3 PCs, K = 3.

#K = 3. Run PCadapt again. Test K and see the p-values graphs
gen.pcadapt = read.pcadapt(lfmm_input, type = c("lfmm"))
class(gen.pcadapt)
pcadapt.test = pcadapt(gen.pcadapt, K=3, ploidy=2, min.maf=0.05, method="mahalanobis")
summary(pcadapt.test)
# Use Mahalanobis distance to compute the (rescaled) test statistics (z-scores in this case).
# The robust Mahalanobis distance is a metric that identifies outliers in multidimensional space. "Robust" means the estimation is not sensitive to outliers in the covariance matrix of the z-scores.

#Graphical tools.
plot(pcadapt.test, option = "manhattan")
plot(pcadapt.test, option = "qqplot")
hist(pcadapt.test$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(pcadapt.test, option = "stat.distribution")

#Choosing a cutoff for outlier detection
pcadapt.test$gif #  1.445385
# The GIF indicates how well the test is "calibrated".
# It corrects for inflation of the test score at each locus, which can occur when population
# structure or other confounding factors are not appropriately accounted for in the model.
#
# GIF of 1=well calibrated, >1 =liberal (too many small p-values), <1=conservative (too few small p-values) # Note: GIFs > 2 indicate poor test calibration.
#
#For a given alpha (real valued number between 0 and 1), SNPs with q-values less than alpha will be considered as outliers with an expected false discovery rate bounded by alpha. The false discovery rate is defined as the percentage of false discoveries among the list of candidate SNPs. Here is an example of how to provide a list of candidate SNPs, for an expected false discovery rate lower than 10%.


#q-values
qval = qvalue(pcadapt.test$pvalues)$qvalues
alpha = 0.1
outliers1 = which(qval < alpha)
length(outliers1) #It will be eliminated 9 SNPs

#Benjamini-Hochberg Procedure
padj = p.adjust(pcadapt.test$pvalues, method="BH")
alpha = 0.1
outliers2 = which(padj < alpha)
length(outliers2) #It will be eliminated 9 SNPs

#Bonferroni correction
padj2 = p.adjust(pcadapt.test$pvalues, method="bonferroni")
alpha = 0.1
outliers3 = which(padj2 < alpha)
length(outliers3) #It will be eliminated 5 SNPs


#Choose one approach to eliminate outlier SNPs. In the case Benjamini-Hochberg Procedure the "outliers2"
#Remove SNPs with adaptive signals
snps_to_keep = which(!c(1:dim(snps_F11@gt)[1]) %in% outliers2)
snps_F12 = snps_F11[snps_to_keep, ]
snps_F12




###3. SAVE VCFs ----
#A. For Adaptive analyses save without HWE and PCadapt filters (snps_F10)
write.vcf(snps_F10, "Inputs/vcf/dioclea_Adapt.vcf")

#B. For Neutral SNPs save final VCF (snps_F12)
write.vcf(snps_F12, "Inputs/vcf/dioclea_Neutral.vcf")

#C. For NeEstimate you need save without the filter by linkage disequilibrium (snps_F9), so jump from step F09 to F11
write.vcf(snps_F12, "Inputs/vcf/dioclea_Ne.vcf")



###4. ESTIMATE BASICS METRICS FROM THE VCFs ----
#A. Choose a vcfR
vcf_subset = snps_F12

#B. Estimate average depth by ind and site:
dp = extract.gt(vcf_subset, element = "DP", as.numeric=TRUE)
head(dp)

#Individual
dp_ind = c()
for (i in 1:length(colnames(dp))){dp_ind = c(dp_ind, mean(dp[,i], na.rm = T))}
mean(dp_ind) #137.1958
max(dp_ind) #365.3462
min(dp_ind) #65.49631

#Site
dp_site = c()
for (i in 1:length(rownames(dp))){dp_site = c(dp_site, mean(dp[i,], na.rm = T))}
mean(dp_site) #141.2887
max(dp_site) #211.963
min(dp_site) #41.84746



#C. Estimate average missing data total, by ind, and by site:
MD = extract.gt(vcf_subset, element = "GT", as.numeric=TRUE)
head(MD)

#Total
((sum(is.na(MD)))/(dim(MD)[1]*dim(MD)[2]))*100 #13.9821


#Individual
MD_ind = c()
for (i in 1:length(colnames(MD))){MD_ind = c(MD_ind, sum(is.na(MD[,i]))/length(rownames(MD))*100)}
mean(MD_ind) #13.9821
max(MD_ind) #39.98021
min(MD_ind) #1.583375


#site
MD_site = c()
for (i in 1:length(rownames(MD))){MD_site = c(MD_site, sum(is.na(MD[i,]))/length(colnames(MD))*100)}
mean(MD_site) # 13.9821
max(MD_site) #28.57143
min(MD_site) #0

#END
