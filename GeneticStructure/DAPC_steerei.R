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
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('dartR'))        install.packages("dartR");             library('dartR')
if (!require('adegenet'))     install.packages("adegenet");          library('adegenet')




#### 1.GENETIC STRUCTURE USING DAPC ----
#### 1.ESTRUTURAÇÃO GENÉTICA USANDO DAPC
#DA Discriminant Analysis focus on between-group variability, while neglecting within-group variation. this method also allows for a probabilistic assignment of individuals to each group, as in Bayesian clustering methods. the method requires the number of variables (alleles) to be less than the number of observations (individuals). This condition is generally not fulfilled in Single Nucleotide Polymorphism (SNP) or re-sequencing datasets. Second, it is hampered by correlations between variables, which necessarily occur in allele frequencies due to the constant-row sum constraint [i.e., compositional data]. Uncorrelated variables will be even more blatant in the presence of linkage disequilibrium
#DAPC relies on data transformation using PCA as a prior step to DA, which ensures that variables submitted to DA are perfectly uncorrelated, and that their number is less than that of analysed individuals (1:10 proportion). Without implying a necessary loss of genetic information, this transformation allows DA to be applied to any genetic data.
#K-means relies on the same model as DA to partition genetic variation into a between-group and a within-group component, and attempts to find groups that minimize the latter. We use Bayesian Information Criterion (BIC) to assess the best supported model, and therefore the number and nature of clusters.

#DA A Análise Discriminante concentra-se na variabilidade entre grupos, enquanto negligencia a variação dentro do grupo. Este método também permite uma atribuição probabilística de indivíduos a cada grupo, como nos métodos de agrupamento bayesiano. o método exige que o número de variáveis (alelos) seja menor que o número de observações (indivíduos). Esta condição geralmente não é atendida no Polimorfismo de Nucleotídeo Único (SNP) ou em conjuntos de dados de re-sequenciamento. Em segundo lugar, é dificultado por correlações entre variáveis, que ocorrem necessariamente em frequências alélicas devido à restrição de soma de linhas constantes [isto é, dados de composição]. Variáveis não correlacionadas serão ainda mais evidentes na presença de desequilíbrio de ligação
#DAPC conta com a transformação de dados utilizando PCA como etapa anterior à AD, o que garante que as variáveis submetidas à AD sejam perfeitamente não correlacionadas e que seu número seja menor que o dos indivíduos analisados. Sem implicar uma perda necessária de informação genética, esta transformação permite que a DA seja aplicada a quaisquer dados genéticos.
#K-means depende do mesmo modelo que DA para particionar a variação genética em um componente entre grupos e dentro do grupo, e tenta encontrar grupos que minimizem o último. Usamos o Critério de Informação Bayesiano (BIC) para avaliar o modelo mais bem suportado e, portanto, o número e a natureza dos clusters.



#A. Load neutral .vcf file after removing outlier SNPs: #A. Carregar o arquivo .vcf após a etapa de filtragem, só com SNPs neutros:
vcf = read.vcfR("Inputs/vcf/steerei_Neutral.vcf", verbose = FALSE)
project = "steerei_Neutral"

#B. Convert "VCF" to "GENIND" #Converter o arquivo vcf para genind
input = vcfR2genind(vcf)
input


#C. Perform a PCA to choose the number of PC in the DAPC: #Executar uma PCA para escolher o número de PCs usados na DAPC
input_scaled = scaleGen (input, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, scale = FALSE, center = TRUE, scannf = FALSE, nf = length(row.names(input@tab))-1)

#D. % of PC variation 
pc_pca = as.data.frame(pca_input$eig)
pc_pca[,2] = (pc_pca/sum(pc_pca))*100

#D1. Rule of 100% of variance: 
index_100 = length(rownames(input@tab))-1
index_100 # number of PC to reach to 100% of explained variance

#D2. Rule of at least 95% of variance:
index_95 = length(which(cumsum(pc_pca[,2]) <= 95))
index_95

#D3. Rule of at least 70% of variance:
index_70 = length(which(cumsum(pc_pca[,2]) <= 70))
index_70 

#D4. Rule of minimum variance per PCs:
variance_pc = 100/(nrow(input@tab)-1)
variance_pc #PCs that increase the explained variance bellow this threshold will be removed
#calculate number of PCs
index_min = length(pc_pca[,2][pc_pca[,2] >= variance_pc])
index_min

#E. Identification of the clusters (We specify that we want to evaluate up to k = 10 groups (max.n.clust = 40)
#If you see a plateau on the graph, you can choose a number of PCs in the beginning of this plateau. In Pilocarpus there's no plateau. We will test different number of PC's
set.seed(13) #set a seed
grp = find.clusters(input, max.n.clust=10, scale = TRUE) # center=T by default.
#Digit the number of PCs to retain. Complete by the results
index_100 #15 PCs - 3 clusters
index_95  #12 PCs - 3 clusters
index_70  #5 PCs - 6 clusters
index_min #2 PCs - 8-10 clusters


#Choose the number of PCs and cluster and verify Group (DPCA) of the first 10 individuals and the Size of Groups
set.seed(13) #set a seed
grp = find.clusters(input, max.n.clust=10, scale = TRUE) # center=T by default.
head(grp$grp, 10)
grp$size



#F. Select the ’best’ BIC is often indicated by an below in the curve of BIC values as a function of k
#save best k graph
dapc_df = as.data.frame(grp$Kstat) %>%
  mutate(K = c(1:10)) %>%
  setNames(c("BIC", "K"))

p = ggplot(data=dapc_df, mapping = aes(x=K,y= BIC)) +
  geom_segment(aes(x=K[1],y= BIC[1], xend =K[10] , yend =BIC[10]), size = 0.5, linetype = 'dashed') +
  geom_line() +
  geom_point(size=4, shape=21, fill="red", color="darkred") +
  theme_bw()+
  scale_x_continuous(breaks = c(1:10)) +
  xlab("Best K") + ylab("Bayesian Information Criterion (BIC)") #set labels

#check the graphics
p
pdf("./Results/DAPC/Bestk_DAPC.pdf", onefile = T)
p
dev.off()


#G. ATTENTION, if your dataset is k = 1 It's finished! Graphs for DAPC need datasets with K >= 2

#H. Choose the best PCs number to recover correctly the clusters. The input should be a Genind object without missing data. Maximum number of PCs is number of individuals -1. Replication by default is 30. Save automatically the graph
pdf("./Results/DAPC/Best_PCs_Number_DAPC.pdf", onefile = T)
number_PCs = xvalDapc(tab(input, NA.method = "mean"), grp$grp, scale = T, n.pca.max = (nrow(input@tab)-1))
dev.off()

#I. Verify the number of PCs and DA used and summary of DAPC, and percentage explained by the three first dapc
number_PCs$DAPC$n.pca #5
number_PCs$DAPC$n.da #2
summary(number_PCs$DAPC)
DAPC1 = paste0("DAPC1 (",round((number_PCs$DAPC$eig[1]/sum(number_PCs$DAPC$eig))*100,2),"%)")
DAPC2 = paste0("DAPC2 (",round((number_PCs$DAPC$eig[2]/sum(number_PCs$DAPC$eig))*100,2),"%)")
DAPC3 = paste0("DAPC3 (",round((number_PCs$DAPC$eig[3]/sum(number_PCs$DAPC$eig))*100,2),"%)")


#J. Verify plots and define the group colors and names:
#plot graph
table(grp$grp)

rownames(input@tab)[which(grp$grp == 1)] #Bolivia
rownames(input@tab)[which(grp$grp == 2)] #Madeira
rownames(input@tab)[which(grp$grp == 3)] #Juruá

#define colors and name for populations
colors_dapc = c('blue', 'red', "yellow" ) #color for the four genetic cluster
labels_dapc = c("Bolivia", "Madeira", "Juruá") #name for the four genetic cluster

#ggplot graph:
df_dapc = number_PCs$DAPC$ind.coord %>%
  as.data.frame() %>%
  rownames_to_column("IND") %>%
  mutate(POP = as.data.frame(grp$grp)[,1])

d_gg = ggplot(df_dapc, aes(x=LD1, y=LD2, fill=POP))+
  geom_point(size=4, pch=21)+
  scale_fill_manual(values= colors_dapc,
                    labels= labels_dapc) +
  theme_bw()+
  guides(fill=guide_legend(title="Populations"))+
  xlab(DAPC1) + ylab(DAPC2) #set labels


#Only One LD
#d_gg = ggplot(df_dapc, aes(x=LD1, fill=POP))+
#  geom_density()+
#  scale_fill_manual(values= colors_dapc,
#                    labels= labels_dapc) +
#  theme_bw()+
#  theme_genetics +
#  xlab("DAPC1 (100%)")# + ylab("DAPC2 (36.51%)") #set labels

d_gg

pdf("./Results/DAPC/Scatter_DAPC_GGPLOT_STE.pdf", onefile = T)
d_gg
dev.off()


#classic plot graph
pdf("./Results/DAPC/Scatter_DAPC_Classical.pdf", onefile = T)
scatter(number_PCs$DAPC, cex = 2, legend = TRUE, col = colors_dapc, txt.leg = labels_dapc,
        clabel = FALSE, posi.leg = "topleft", scree.pca = TRUE, pch=19:20,
        posi.pca = "topright", posi.da = "bottomleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)
dev.off()



#K. Add DAPC clusters ID according the VCF file
identical(row.names(as.data.frame(grp$grp)), as.character(rownames(input@tab)))
setdiff(row.names(as.data.frame(grp$grp)), as.character(rownames(input@tab)))
setdiff(as.character(rownames(input@tab)), row.names(as.data.frame(grp$grp)))

as.data.frame(grp$grp)

#L. Save the dataframe with sample name and pop ID file
write.csv(as.data.frame(grp$grp), paste0("Results/DAPC/DAPC_clusters_", project, ".csv"))

#END
#FIM
