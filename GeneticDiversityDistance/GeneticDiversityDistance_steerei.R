#################### UNIVERSIDADE FEDERAL DA PARAÍBA ##########################
################ DEPARTAMENTO DE SISTEMÁTICA E ECOLOGIA #######################
######################     FILTERING SNPs TUTORIAL     ########################

#Jeronymo Dalapicolla, 2023: Guia de genômica de populações aplicada a mamíferos Neotropicais


##I. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT ----
##I. REMOVA QUALQUER OBJETO OU FUNÇÃO DO AMBIENTE ----
rm(list=ls())

##II. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! ----
##II. ESCOLHA UMA PASTA PARA EXECUTAR AS ANÁLISES. OS ARQUIVOS DEVEM ESTAR LÁ! ----
# IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H
# NO RStudio VÁ PARA A SESSÃO >> DEFINIR DIRETÓRIO DE TRABALHO >> ESCOLHA O DIRETÓRIO.. NA BARRA DE FERRAMENTAS DO RStudio OU USE O ATALHO CTRL+SHIFT+H

##III. INSTALL AND LOAD THE PACKAGES ----
##III. INSTALAR E CARREGAR OS PACOTES ----
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('adegenet'))     install.packages("adegenet");          library('adegenet')
if (!require('diveRsity'))    devtools::install_github("kkeenan02/diveRsity");         library("diveRsity")
#if (!require('diveRsity'))    install.packages("diveRsity");         library("diveRsity")
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('dartR'))        install.packages("dartR");             library('dartR')
if (!require('raster'))       install.packages("raster");            library('raster')
if (!require('poppr'))        install.packages("poppr");             library('poppr')


####1. PREPARING INPUTS ---- 
####1. PREPARANDO OS ARQUIVOS DE ENTRADA ----
#A. Load neutral .vcf file after filtering SNPs: #A. Carregar o arquivo .vcf após a etapa de filtragem, só com SNPs neutros:
vcf = read.vcfR("Inputs/vcf/steerei_Neutral.vcf")
#check file #Verificar o arquivo
vcf

#B. Convert "VCF" to "GENIND"
genind_data = vcfR2genind(vcf) #CONVERTER VCF PARA GENIND
genind_data

#C. Convert "VCF" to "GENLIGHT" #CONVERTER VCF PARA GENLIGHT
genlight_data = vcfR2genlight(vcf)
genlight_data

#D. Load information on samples #Carregar as informações das amostras
geo_data = read.csv("Inputs/metadata/coord_steerei.csv") 
head(geo_data)
tail(geo_data)

#select rows in the same order than vcf #selecinar as linhas da tabela na mesma ordem que o vcf
genomic_order = as.data.frame(rownames(genind_data@tab))
names(genomic_order) = "Sample_ID"

geo_data_reorder = semi_join(geo_data, genomic_order)
geo_data_reorder = inner_join(genomic_order, geo_data_reorder)
head(geo_data_reorder)

#verify the order: #Verificar a ordem
identical(as.character(rownames(genind_data@tab)), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(rownames(genind_data@tab)), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(geo_data_reorder$Sample_ID), as.character(rownames(genind_data@tab)))



####2. GENETIC DIVERSIY ----
####2. DIVERSIDADE GENÉTICA ----
#A.Ho, He, uHe and FIS:
#A.Ho, He, uHe e FIS:
genlight_data_POP = genlight_data
genlight_data_POP@pop = as.factor(geo_data_reorder$pop)
genlight_data_POP = gl.compliance.check(genlight_data_POP)
div_pop = gl.report.heterozygosity (genlight_data_POP, method = "pop")
write.csv(as.data.frame(div_pop), "./Results/Diversity/metrics_STE_POP.csv")

genlight_data_LOC = genlight_data
genlight_data_LOC@pop = as.factor(geo_data_reorder$State.Province)
genlight_data_LOC = gl.compliance.check(genlight_data_LOC)
div_loc = gl.report.heterozygosity (genlight_data_LOC, method = "pop")
write.csv(as.data.frame(div_loc), "./Results/Diversity/metrics_STE_LOC.csv")



#### 3. GENETIC DISTANCE----
#### 3. DISTÂNCIA GENÉTICA----
#A.BY CLUSTER FST: #FST POR CLUSTER E LOCALIDADE
#Use genlight data
#Usar genlight
genlight_data_POP
genlight_data_LOC

# Run FST #Executar o FST
fst_pop = gl.fst.pop(genlight_data_POP, nboots = 100, percent = 95, nclusters = 1)
fst_loc = gl.fst.pop(genlight_data_LOC, nboots = 100, percent = 95, nclusters = 1)

#save results #Salvar os resultados
write.csv(fst_pop$Fsts, file=paste0("Results/Distance/FSTpop_Fstas_STE.csv"))
write.csv(fst_pop$Pvalues, file=paste0("Results/Distance/FSTpop_pvalues_STE.csv"))
write.csv(fst_pop$Bootstraps, file=paste0("Results/Distance/FSTpop_bootstrap_STE.csv"))

write.csv(fst_loc$Fsts, file=paste0("Results/Distance/FSTloc_Fstas_STE.csv"))
write.csv(fst_loc$Pvalues, file=paste0("Results/Distance/FSTloc_pvalues_STE.csv"))
write.csv(fst_loc$Bootstraps, file=paste0("Results/Distance/FSTloc_bootstrap_STE.csv"))

#END
#FIM
