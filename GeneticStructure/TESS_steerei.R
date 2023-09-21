#################### UNIVERSIDADE FEDERAL DA PARAÍBA ##########################
################ DEPARTAMENTO DE SISTEMÁTICA E ECOLOGIA #######################
######################     FILTERING SNPs TUTORIAL     ########################

#By Jeronymo Dalapicolla, 2023: Guia de genômica de populações aplicada a mamíferos Neotropicais



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
if (!require('tess3r'))       devtools::install_github("bcm-uga/TESS3_encho_sen");   library('tess3r')

#From CRAN R:
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('ggrepel'))      install.packages("ggrepel");           library('ggrepel')
if (!require('reshape2'))     install.packages("reshape2");          library('reshape2')
if (!require('maps'))         install.packages("maps");              library('maps')
if (!require('raster'))       install.packages("raster");            library('raster')



##IV. AUXILIARY FUNCTIONS ----
create_dir = function(dir_names){
  for (i in 1:length(dir_names)){
    if (dir.exists(file.path(getwd(),dir_names[i])) == FALSE) {dir.create(dir_names[i], recursive = T, showWarnings = T) }
    else (print (paste0(dir_names[i], " has already been created. Be careful you can overwrite files")))}
  
}



#### 2. LOAD FILES ----
# Following this tutorial: https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html
#A. Load VCF
geno_name = "steerei_Neutral"
vcf = read.vcfR("Inputs/vcf/steerei_Neutral.vcf", verbose = T)
vcf

#B. Convert VCF to GENOTYPES 
genotypes= t(extract.gt(vcf, element = "GT", mask = FALSE, as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = TRUE))
#check:
genotypes[1:5,1:5]
row.names(genotypes)

#C. Load information on samples
geo_data = read.csv("Inputs/metadata/coord_steerei.csv") 
head(geo_data)
tail(geo_data)

#D. Filter data according to the genetic samples order
genomic_order = as.data.frame(row.names(genotypes))
names(genomic_order) = "Sample_ID"

geo_data_reorder = semi_join(geo_data, genomic_order)
geo_data_reorder = inner_join(genomic_order, geo_data_reorder)
head(geo_data_reorder)

identical(as.character(row.names(genotypes)), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(row.names(genotypes)), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(geo_data_reorder$Sample_ID), as.character(row.names(genotypes)))


#E. Create a Matrix with long and lat 
coordinates = geo_data_reorder %>%
  dplyr::select(Longitude, Latitude) %>%
  data.matrix(., rownames.force = NA)
#verify the coords
plot(coordinates, pch = 19, cex = .5, xlab = "Longitude", ylab = "Latitude")





#### 3. RUNNING TESS3 ----
#A. Customize values for run TESS3
lambda_values = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5) #test lambda values around 1.
K = c(1:10) # set the number of K to be tested
replications = 10 # number of replication in each K
ploidy = 2 # species ploidy
CPU = 4 #Number of cores for run in parallel
mask = 0.05 #proportion of masked values

set.seed(13)
path = "./Results/TESS/"
if (dir.exists(file.path(getwd(), path)) == FALSE)
{dir.create(path, recursive = T, showWarnings = F)} else (print (paste0(path, " has already been created. Be careful with overwritting")))

for (i in lambda_values){
  tess3.ls = tess3(genotypes, coord = coordinates, K = K, mask = mask, lambda = i,
                   method = "projected.ls", max.iteration = 5000, rep = replications,
                   ploidy = ploidy, openMP.core.num = CPU)
  save(tess3.ls, file = paste0("./Results/TESS/Tess3ls_Lambda_", i,"_", geno_name, ".RData")) 
  
}

#B. Choose the lambda with the minimum cross-validation value:
#create a matrix to save results
cross_value = matrix(NA, max(K)*length(lambda_values), 3)
colnames(cross_value) = c("K", "Lambda", "Crossvalidation")

loop = 0 #Always set loop = 0.
for (i in lambda_values){
  load(file = paste0("./Results/TESS/Tess3ls_Lambda_", i,"_", geno_name, ".RData"))
  for (j in 1:max(K)){
    loop=loop+1
    res = Gettess3res(tess3.ls, K=j)
    cross_value[loop,1] = j
    cross_value[loop,2] = i
    cross_value[loop,3] = min(res$crossvalid.crossentropy)}
}
#save as csv
write.csv(cross_value,paste0("./Results/TESS/Tess3ls_Crossvalidation_values_", geno_name, ".csv"))
#choose best lambda
lambda_tess = as.vector(cross_value[cross_value[,3] == min(cross_value[,3]), ][2])
lambda_tess



#C. Choose best K by graphics:
#load the best lambda project:
load(file = paste0("./Results/TESS/Tess3ls_Lambda_", lambda_tess,"_", geno_name, ".RData"))
#plot results, best K is more distance from dashed line. 4 or 5 clusters.
pdf(paste0("./Results/TESS/TESS3_RE_PlotK_Lambda_", lambda_tess, ".pdf"), onefile =F)
plot.new()
plot(tess3.ls, pch = 19, col = "blue", type="b",lwd=2,
     crossvalid = T, crossentropy = T,
     xlab = "Number of ancestral populations", cex.lab=1.5,
     ylab = "Cross-validation score")
segments(x0 = 1, y0 = median(tess3.ls[[1]]$crossvalid.crossentropy),
        x1 = 10, y1 = median(tess3.ls[[10]]$crossvalid.crossentropy), lty =2, lwd = 2)
dev.off()



#D. Choose best K by assessement greater than 80%, following Bastos et al. 2016 doi: 10.1371/journal.pone.0145165:
#create the dataframe:
results_acess = as.data.frame(matrix(NA, 10, 3))
colnames(results_acess) = c("K", "H_Assessment", "Percentage")
for (l in K){
  res_tess = Gettess3res(tess3.ls, K=l)
  z = res_tess$Q > 0.80
  z = table(z)["TRUE"]
  pc = (z/length(row.names(genotypes)))*100
  results_acess[l,] = c(l, z, pc)
}
#plot results. After K=2 which K has the assessement > 0.85 to   
pdf(paste0("./Results/TESS/TESS3_Assessment_Lambda_", lambda_tess, ".pdf"), onefile =F)
ggplot(data = results_acess, mapping = aes(x= K, y = Percentage)) + 
  geom_line() +
  geom_point(size= 2) +
  xlab("Number of Clusters (K)") + ylab("Assessement (%)") +
  scale_x_continuous(breaks = c(1:10)) +
  theme_bw() +
  geom_label_repel(mapping = aes(x= K, y = Percentage, label=round(Percentage,1)))
dev.off()





#### 4. RESULTS OF TESS3 ----
#A.Choose the Best K and create the Q-matrix
#Best K
optimal_K = 3
res_tess3 = Gettess3res(tess3.ls, K=optimal_K)


#B. Add admixture coefficient and replace the population ID to vcf file
coeff = c()
for (i in 1:optimal_K){coeff[i] = paste0("Adx_Coeff_", i)}
coeff

Qmat =
  res_tess3$Q %>%
  as.data.frame() %>%
  setNames(coeff) %>%
  mutate(PopID_TESS3 = apply(., 1, which.max))
head(Qmat)


#C. define colors and name for populations:
geo_data_final = cbind(geo_data_reorder, Qmat)
geo_data_final[geo_data_final$PopID_TESS3 == 1,] #Madeira
geo_data_final[geo_data_final$PopID_TESS3 == 2,] #Juruá
geo_data_final[geo_data_final$PopID_TESS3 == 3,] #Bolivia


colors_tess = c('red', "yellow", 'blue') #color for the genetic cluster
labels_tess = c("Madeira", "Juruá", "Bolivia") #name for the genetic cluster


#G. Save new vcf file with TESS results and pop ID
write.csv(geo_data_final, paste0("Results/TESS/ancestry_coef_TESS_", geno_name, ".csv"), quote = F) #save result as table


#D. Create a dataframe (df) for ggplot2, all localities
coeff
df = geo_data_final %>%
  dplyr::select(Sample_ID, pop, Adx_Coeff_1, Adx_Coeff_2, Adx_Coeff_3) %>%
  dplyr::arrange(pop) %>%
  melt(., id.vars=c("Sample_ID", "pop"))
head(df)
names(df) = c("Sample_ID", "pop", "variable", "value")
head(df)

#number of samples
nsample = 16

p = ggplot(data=df, mapping = aes(x=factor(Sample_ID),factor(pop),
                                  y= value*100,
                                  fill = factor(variable))) +
  geom_bar(stat="identity", width = 1, size=0.3, colour="black") +
  scale_fill_manual(values= colors_tess,
                    labels= labels_tess) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(limits=df$Sample_ID[1:nsample], labels =df$pop[1:nsample], guide = guide_axis(n.dodge=1)) +
  guides(fill=guide_legend(title="Populations"))+
  xlab("") + ylab("")
#check
p

pdf(paste0("./Results/TESS/Pipegraph_TESS_", geno_name, ".pdf"), onefile =F)
plot(p)
dev.off()




#F. Spatial interpolation of ancestry coefficient
##Maps with TESS
source("Functions_TESS_maps.R")

pdf("./Results/TESS/TESS_MAP2.pdf", onefile =F)
plot.new()
asc.raster="http://membres-timc.imag.fr/Olivier.Francois/RasterMaps/South_America.asc"
grid=createGridFromAsciiRaster(asc.raster)
constraints=getConstraintsFromAsciiRaster(asc.raster, cell_value_min=0)
maps(matrix = res_tess3$Q, coordinates, grid, method = "max",
     onemap=T, onepage=T,  xlim = c(-80, -50), ylim = c(3, -12),
     main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude", cex = .5)
dev.off()


#END
