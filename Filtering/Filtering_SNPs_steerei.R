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
#Basic Packages for installation:
#Pacotes básicos para instalação de outros pacotes
if (!require('remotes'))      install.packages('remotes');           library('remotes')
if (!require('BiocManager'))  install.packages('BiocManager');       library('BiocManager')
if (!require('pacman'))       install.packages('pacman');            library('pacman')
if (!require('devtools'))     install.packages('devtools');          library('devtools')

#From Github or BiocManager:
#Do Github ou BiocManager:
if (!require('LEA'))          BiocManager::install("LEA");            library('LEA')
if (!require('qvalue'))       BiocManager::install("qvalue");         library('qvalue')
if (!require('SNPRelate'))    BiocManager::install("SNPRelate");      library('SNPRelate')

#From CRAN R:
#Do CRAN do R:
if (!require('tidyverse'))     install.packages("tidyverse");         library('tidyverse')
if (!require('vcfR'))          install.packages("vcfR");              library('vcfR')
if (!require("SNPfiltR"))      install.packages("SNPfiltR");          library("SNPfiltR")
if (!require("dartR"))         install.packages("dartR");             library("dartR")
if (!require('pcadapt'))       install.packages("pcadapt");           library('pcadapt')
if (!require('sf'))            install.packages("sf");                library('sf')
if (!require('HardyWeinberg')) install.packages("HardyWeinberg");     library('HardyWeinberg')


#################### FILTERING USING RAD-PARAMETERS ################################## ----
################# FILTRAGEM USANDO PARÂMETROS PARA RAD ############################### ----

####1. LOAD FILES ----
####1. CARREGAR ARQUIVOS ----
#A. VCF
snps_raw = read.vcfR("Inputs/vcf/steerei_raw.vcf", verbose = T)
snps_raw

#B. Remove individuals before filtering SNPs if you need. In this case we will remove two samples as an example:
#B. Remover indivíduos antes de filtrar SNPs, se necessário. Neste caso retiraremos duas amostras como exemplo:
colnames(snps_raw@gt)
remove_ind = c("MVZ194874", "MVZ194879")
ind_to_keep = colnames(snps_raw@gt)[!colnames(snps_raw@gt) %in% remove_ind]

snps_raw = snps_raw[samples=ind_to_keep]
snps_raw


#C. Format a PopMAP. A data frame with two column popmap with the same format as Stacks, and the columns must be named 'id' and 'pop'
# Load information on samples
#C. Formatar um arquivo PopMAP. PopMap é Uma tabela com duas colunas com o mesmo formato do Stacks, e as colunas devem ser nomeadas 'id' e 'pop'
# Carregar informações nas amostras
geo_data = read.csv("Inputs/metadata/coord_steerei.csv") 
head(geo_data)
tail(geo_data)

#select rows in the same order than vcf
#selecionar as linhas da tabela na mesma ordem de indivíduos do arquivo vcf
genomic_order = as.data.frame(colnames(snps_raw@gt)[-1]) #always 1st line will be FORMAT. You need to remove #1ª linha sempre é esse FORMAT e não a primeira amostra
names(genomic_order) = "Sample_ID"

geo_data_reorder = semi_join(geo_data, genomic_order)
geo_data_reorder = inner_join(genomic_order, geo_data_reorder)
head(geo_data_reorder)

#verify the order:
#verificar a ordem:
identical(as.character(colnames(snps_raw@gt)[-1]), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(colnames(snps_raw@gt)[-1]), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(geo_data_reorder$Sample_ID), as.character(colnames(snps_raw@gt)[-1]))

popmap = geo_data_reorder[,c(1,12)] #col "pop" as population/location. #Coluna "pop" vai ser usada como população/localidade nesse exemplo
names(popmap) = c("id", "pop")
class(popmap)
head(popmap)




###2. FILTERING ----
###2. FILTRAGEM ----

#A. Remove indels
#A. Remover indels
snps_F1 = extract.indels(snps_raw, return.indels = FALSE)


#B. Keep Biallelic SNPs
#B. Manter apenas SNPs bialélicos
snps_F2 = filter_biallelic(snps_F1)
snps_F2


#C. Reduce dataset with filter by Missing Data per site
#C. Reduzir o dataset com filtro de dados faltantes por sítio/SNPs
missing_by_snp(snps_F2) #verify statistics  #verificar as estatísticas antes de escolher um limiar. No exemplo, 60%
snps_F3 = missing_by_snp(vcfR=snps_F2, cutoff = .6)
snps_F3


#D. Reduce dataset with filter by Missing Data per sample
#D. Reduzir o dataset com filtro de dados faltantes por indivíduo/Amostra
missing_by_sample(vcfR=snps_F3, popmap = popmap)  #verify statistics #verificar as estatísticas antes de escolher um limiar. No exemplo, 60%
snps_F4 = missing_by_sample(vcfR=snps_F3, cutoff = .6)
snps_F4

popmap = popmap[popmap$id %in% colnames(snps_F1@gt),]
head(popmap)


#C. By Quality
#C. Por Qualidade.
#QUAL = 30 #definir um threshold. O Stacks completa o campo QUAL do vcf com NA, então é preciso fazer essa filtragem dentro da etapa de identificação dos SNPs
#snps_to_keep = which(snps_F4@fix[,6] > QUAL)
#snps_F5 = snps_F4[snps_to_keep, ]
#snps_F5
snps_F5 = snps_F4


#D. By Minimum Depth by site
#D. Por cobertura bruta mínima por sítio. Similar ao parâmetro -m do Stacks.
#Usually you do not have information on Genotype quality, if you have you can use gq filter
#Normalmente você não tem informações sobre a qualidade do genótipo, se tiver pode usar o filtro "gq"
snps_F6 = hard_filter(vcfR = snps_F5, depth = 4, gq = NULL)
snps_F6


#E. By Maximum Depth by site
#E. Por cobertura máxima por sítio.
#dp by site # calcular a cobertura de cada sítio:
dp = extract.gt(snps_F6, element = "DP", as.numeric=TRUE)
head(dp)
dp_site = c()
for (i in 1:length(rownames(dp))){dp_site = c(dp_site, mean(dp[i,], na.rm = T))}
maxdep = median(dp_site, na.rm = T) + (2*sd(dp_site, na.rm = T))
maxdep
#filter by max depth #Filtro pela cobertura máxima
snps_F7 = max_depth(snps_F6, maxdepth = maxdep)
snps_F7


#F. by MAC or MAF
#F. by MAC ou MAF
#estimate MAC by MAF #estimar o MAC a partir do MAF
MAF_1 = (length(colnames(snps_F7@gt)[-1])*1)/100 #1% das amostras
MAF_3 = (length(colnames(snps_F7@gt)[-1])*3)/100 #3% das amostras
MAF_5 = (length(colnames(snps_F7@gt)[-1])*5)/100 #5% das amostras
MAF_10 = (length(colnames(snps_F7@gt)[-1])*10)/100 #10% das amostras

snps_F8 = min_mac(snps_F7, min.mac = MAF_5)
snps_F8



#G. By Allele Balance:
#G. Equilíbrio Alélico:
snps_F9 = filter_allele_balance(snps_F8, min.ratio = 0.2, max.ratio = 0.8)
snps_F9


#H. By Missing Data per site
#H. Por dados faltantes por sítio
missing_by_snp(vcfR=snps_F9) #verify statistics
snps_F10 = missing_by_snp(vcfR=snps_F9, cutoff = .8) #Remover acima de 20%
snps_F10


#I. By Missing Data per sample:
#I. Por dados faltantes por indivíduo/amostra:
missing_by_sample(vcfR=snps_F10, popmap = popmap)  #verify statistics
snps_F11 = missing_by_sample(vcfR=snps_F10, cutoff = .2) # Manter até 20% 
snps_F11

#redefinir o popmap depois de eliminar indivíduos
popmap = popmap[popmap$id %in% colnames(snps_F11@gt),]
head(popmap)



#J. By Linked SNPs.
#J. Remover SNPs ligados. Manter só um SNP por contig ou Chrom
snps_F12 = distance_thin(vcfR = snps_F11, min.distance = 150) #1 SNP by contig/Chrom #min.distance deve ser o tamanho do fragmento
snps_F12



#K. By Hardy Weinberg Equilibrium (HWE)
#K. Por Equilíbrio de Hardy-Weinberg (HWE)
#convert to genlight
#Converter para genlight
gl_vcf = vcfR2genlight(snps_F12)
gl_vcf

#add pop/local information #adicionar as informações de população/localidade
pop(gl_vcf) = popmap$pop
popNames(gl_vcf)

#filter in the genlight as recommended by Pearman et al 2022:
#filtros como os recomendados por Pearman et al. 2022:
gl_vcf_fil = gl.filter.hwe(
  gl_vcf,
  subset = "each",
  n.pop.threshold = round(length(popNames(gl_vcf))/2,0), # most of pops #maioria das populações
  method_sig = "Exact", # Wigginton et al. (2005)
  multi_comp = FALSE,
  multi_comp_method = "BH", # Benjamini & Hochberg, 1995
  alpha_val = 0.05,
  pvalue_type = "midp", # Graffelman & Moreno, 2013
  cc_val = 0.5,
  min_sample_size = 3, #minimum sample size # mínimo de indivíduos por população pra análise rodar
  verbose = NULL
)

head(gl_vcf_fil@loc.names)
head(gl_vcf@loc.names)

#keep in vcf the same snps in the gl filtered object.
#manter no vcf os mesmos snps que ficaram no objeto filtrado gl.
filtered_snps = which(gl_vcf@loc.names %in% gl_vcf_fil@loc.names)
snps_F13= snps_F12[filtered_snps, ]
snps_F13




#L. By outliers SNPs - possibly under selection (Adapt)
#L. Remover SNPs atípicos - possivelmente sob seleção (Adaptativo)
#convert VCF to Genotypes #Converter para Genotypes
genotypes= t(extract.gt(snps_F13, element = "GT", mask = FALSE, as.numeric = TRUE, return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = TRUE))
#check:
genotypes[1:5,1:5]
row.names(genotypes)

#Convert Genotypes to LFMM #Converter para LFMM
dim(genotypes)
((sum(is.na(genotypes)))/(dim(genotypes)[1]*dim(genotypes)[2]))*100
#4.59%  #Amount of missing data #quantidade de dados faltantes

genotypes[is.na(genotypes)] = 9 #The missing genotypes have to be encoded with the value 9 or -9 #Os genótipos faltantes devem ser codificados com o valor 9 ou -9
genotypes[1:10,1:10]
write.lfmm(genotypes,"Inputs/lfmm/F13.lfmm")

#Load LFMM: #carregar LFMM
lfmm_input = read.lfmm("Inputs/lfmm/F13.lfmm")
class(lfmm_input)

#Convert LFMM to a PCadapt matrix and run the analyses. Converta LFMM em uma matriz PCadapt e executar as análises
pcadapt.test = lfmm_input %>%
  read.pcadapt(., type="lfmm") %>%
  pcadapt(., K=10, ploidy=2)

#"Cattell's Rule" for interpreting the scree plot (PC to the left of the flat line)
#"Regra de Cattell" para interpretar o gráfico de scree (PC à esquerda da linha plana)
plot(pcadapt.test, option="screeplot") #3 PCs
plot(pcadapt.test, option="scores") #PC1 and PC2 #3 clusters
plot(pcadapt.test, option = "scores", i = 3, j = 2) #PC2 and PC3 # 3 clusters
plot(pcadapt.test, option = "scores", i = 3, j = 4) #PC3 and PC4 # 1 cluster
plot(pcadapt.test, option = "scores", i = 5, j = 4) #PC5 and PC4 # 1 cluster
plot(pcadapt.test, option = "scores", i = 5, j = 6) #PC5 and PC6 # 1 cluster
#With PC3 we see a clear clustering pattern #Com 3 PCs temos a maior separação das amostras
#So we can use 3 PCs, K = 3. #então usamos 3PCs

#K = 3. Run PCadapt again. Test K and see the p-values graphs #Rodamos novamente fixando o melhor valor de PCs que é 3
gen.pcadapt = read.pcadapt(lfmm_input, type = c("lfmm"))
class(gen.pcadapt)
pcadapt.test = pcadapt(gen.pcadapt, K=3, ploidy=2, min.maf=0.05, method="mahalanobis")
summary(pcadapt.test)
# Use Mahalanobis distance to compute the (rescaled) test statistics (z-scores in this case).
# The robust Mahalanobis distance is a metric that identifies outliers in multidimensional space. "Robust" means the estimation is not sensitive to outliers in the covariance matrix of the z-scores.
# Use a distância de Mahalanobis para calcular as estatísticas de teste (redimensionadas) (pontuações z neste caso).
# A distância de Mahalanobis robusta é uma métrica que identifica valores discrepantes no espaço multidimensional. "Robusto" significa que a estimativa não é sensível a valores discrepantes na matriz de covariância das pontuações z.


#Graphical tools. #gráficos para visualização. Tem muito outliers
plot(pcadapt.test, option = "manhattan")
plot(pcadapt.test, option = "qqplot")
hist(pcadapt.test$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(pcadapt.test, option = "stat.distribution")

#Choosing a cutoff for outlier detection
pcadapt.test$gif #   1.584102
# The GIF indicates how well the test is "calibrated".
# It corrects for inflation of the test score at each locus, which can occur when population
# structure or other confounding factors are not appropriately accounted for in the model.
#
# GIF of 1=well calibrated, >1 =liberal (too many small p-values), <1=conservative (too few small p-values) # Note: GIFs > 2 indicate poor test calibration.
#
#For a given alpha (real valued number between 0 and 1), SNPs with q-values less than alpha will be considered as outliers with an expected false discovery rate bounded by alpha. The false discovery rate is defined as the percentage of false discoveries among the list of candidate SNPs. Here is an example of how to provide a list of candidate SNPs, for an expected false discovery rate lower than 10%.


# O GIF indica quão bem o teste está "calibrado".
# Corrige a inflação da pontuação do teste em cada locus, o que pode ocorrer quando a população é
# estruturada ou outros fatores de intrínsecos não são devidamente contabilizados no modelo.
#
# GIF de 1 = modelo bem calibrado, >1 = liberal (muitos valores p pequenos), <1 = conservador (poucos valores p pequenos) # Observação: GIFs > 2 indicam calibração de teste ruim.
#
#Para um determinado alfa (número real com valor entre 0 e 1), os SNPs com valores q menores que alfa serão considerados como outliers com uma taxa de descoberta falsa esperada limitada por alfa. A taxa de falsas descobertas é definida como a porcentagem de falsas descobertas na lista de SNPs candidatos. Aqui está um exemplo de como fornecer uma lista de SNPs candidatos, para uma taxa esperada de falsas descobertas inferior a 10%.


#q-values
qval = qvalue(pcadapt.test$pvalues)$qvalues
alpha = 0.1
outliers1 = which(qval < alpha)
length(outliers1) #It will be eliminated 3769 SNPs #nº de SNPs que serão eliminados

#Benjamini-Hochberg Procedure
padj = p.adjust(pcadapt.test$pvalues, method="BH")
alpha = 0.1
outliers2 = which(padj < alpha)
length(outliers2) #It will be eliminated 3769 SNPs #nº de SNPs que serão eliminados

#Bonferroni correction
padj2 = p.adjust(pcadapt.test$pvalues, method="bonferroni")
alpha = 0.1
outliers3 = which(padj2 < alpha)
length(outliers3) #It will be eliminated 3008 SNPs #nº de SNPs que serão eliminados


#Choose one approach to eliminate outlier SNPs. In the case Benjamini-Hochberg Procedure the "outliers2"
#Remove SNPs with adaptive signals
#Escolha uma abordagem para eliminar SNPs discrepantes. No caso do Procedimento Benjamini-Hochberg é o arquivo “outliers2”
#Remover SNPs com sinais adaptativos
snps_to_keep = which(!c(1:dim(snps_F13@gt)[1]) %in% outliers2)
snps_F14 = snps_F13[snps_to_keep, ]
snps_F14




#M. Remove by missing data for 20%
#M. Remover dados faltantes superiores a 20%
missing_by_sample(vcfR=snps_F14, popmap = popmap)  #verify statistics #verificar as estatísticas
snps_F15 = missing_by_sample(vcfR=snps_F14, cutoff = .2)
snps_F15

missing_by_snp(vcfR=snps_F15)  #verify statistics #verificar as estatísticas
snps_F16 = missing_by_snp(vcfR=snps_F15, cutoff = .8)
snps_F16


#N. By monomorphic SNPs:
#N. Remover sítios Monomórficos
gl_vcf_final = vcfR2genlight(snps_F16)
gl_vcf_final2 = gl.filter.monomorphs(gl_vcf_final)
gl_vcf_final2

filtered_snps2 = which(gl_vcf_final@loc.names %in% gl_vcf_final2@loc.names)
snps_F17= snps_F16[filtered_snps2, ]
snps_F17



###3. SAVE VCFs ----
###3. SALVAR VCFs ----
#A. For Adaptive analyses save without HWE and PCadapt filters (snps_F12)
#a. Para análises adaptativas, salve sem filtros HWE e PCadapt (snps_F12)
write.vcf(snps_F12, "Inputs/vcf/steerei_Adapt.vcf")

#B. For Neutral SNPs save final VCF (snps_F17)
#B. Para SNPs neutros, salvar o arquivo final de VCF
write.vcf(snps_F17, "Inputs/vcf/steerei_Neutral.vcf") 



###4. ESTIMATE BASICS METRICS FROM THE VCFs ----
###4. ESTIMAR MÉTRICAS BÁSICAS DOS VCFs ----
#A. Choose a vcfR #Escolher o arquivo
vcf_subset = snps_F17

#B. Estimate average depth by ind and site:
#B. Estimar a cobertura média por indivíduos e por SNP 
dp = extract.gt(vcf_subset, element = "DP", as.numeric=TRUE)
head(dp)

#Individual
dp_ind = c()
for (i in 1:length(colnames(dp))){dp_ind = c(dp_ind, mean(dp[,i], na.rm = T))}
mean(dp_ind) #33.80
max(dp_ind) #50.81
min(dp_ind) #17.97

#SNP
dp_site = c()
for (i in 1:length(rownames(dp))){dp_site = c(dp_site, mean(dp[i,], na.rm = T))}
mean(dp_site) #34.106
max(dp_site) #54.786
min(dp_site) #12



#C. Estimate average missing data total, by ind, and by site:
#C. Estimar a média de dados faltantes total, por indivíduo, and by SNP:
MD = extract.gt(vcf_subset, element = "GT", as.numeric=TRUE)
head(MD)

#Total
((sum(is.na(MD)))/(dim(MD)[1]*dim(MD)[2]))*100 #5.17


#Individual
MD_ind = c()
for (i in 1:length(colnames(MD))){MD_ind = c(MD_ind, sum(is.na(MD[,i]))/length(rownames(MD))*100)}
mean(MD_ind) #5.17
max(MD_ind) #14.37
min(MD_ind) #1.31


#site
MD_site = c()
for (i in 1:length(rownames(MD))){MD_site = c(MD_site, sum(is.na(MD[i,]))/length(colnames(MD))*100)}
mean(MD_site) #5.17
max(MD_site) #18.75
min(MD_site) #0

#END
#FIM
