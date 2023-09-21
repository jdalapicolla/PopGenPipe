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
if (!require('adegenet'))     install.packages("adegenet");          library('adegenet')
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('dartR'))        install.packages("dartR");             library('dartR')
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')




####1. GENETIC STRUCTURE USING PCA ----
####1. ESTRUTURAÇÃO GENÉTICA USANDO PCA ---- 
#PCA has the ability to identify genetic structures in very large datasets within negligible computational time, and the absence of any assumption about the underlying population genetic model. PCA aims to summarize the overall variability among individuals, which includes both the divergence between groups (i.e., structured genetic variability), and the variation occurring within groups (‘random’ genetic variability).
#PCA tem a capacidade de identificar estruturas genéticas em conjuntos de dados muito grandes dentro de um tempo computacional curto, tem a ausência de qualquer suposição sobre o modelo genético populacional. O PCA visa resumir a variabilidade geral entre os indivíduos, que inclui tanto a divergência entre grupos (isto é, variabilidade genética estruturada) quanto a variação que ocorre dentro dos grupos (variabilidade genética 'aleatória').

#A. Load neutral .vcf file after filtering SNPs #A. Carregar o arquivo .vcf após a etapa de filtragem, só com SNPs neutros:
vcf = read.vcfR("Inputs/vcf/steerei_Neutral.vcf")
#check file # Verificar o arquivo
vcf

#B. Convert "VCF" to "GENIND" #Converter o arquivo vcf para genind
input = vcfR2genind(vcf)
input
input@tab[1:5,1:5]


#C. Load information on samples #Carregar as informações sobre as amostras
geo_data = read.csv("Inputs/metadata/coord_steerei.csv") 
head(geo_data)
tail(geo_data)

#select rows in the same order than vcf #selecinar as linhas da tabela na mesma ordem que o vcf
genomic_order = as.data.frame(rownames(input@tab))
names(genomic_order) = "Sample_ID"

geo_data_reorder = semi_join(geo_data, genomic_order)
geo_data_reorder = inner_join(genomic_order, geo_data_reorder)
head(geo_data_reorder)

#verify the order:
identical(as.character(rownames(input@tab)), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(rownames(input@tab)), as.character(geo_data_reorder$Sample_ID))
setdiff(as.character(geo_data_reorder$Sample_ID), as.character(rownames(input@tab)))


#D. Perform a PCA: #Executar o PCA
input_scaled = scaleGen (input, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, scale = FALSE, center = FALSE, scannf = FALSE, nf = length(row.names(input@tab))-1)


#E. % of contribution for three first PCs #Contribuição dos 3 primeiros PCs
PC1 = paste0("PC1 (",round(pca_input$eig[1]/sum(pca_input$eig)*100,2),"%)")
PC2 = paste0("PC2 (",round(pca_input$eig[2]/sum(pca_input$eig)*100,2),"%)")
PC3 = paste0("PC3 (",round(pca_input$eig[3]/sum(pca_input$eig)*100,2),"%)")

#F. Min and max values for axes #valores mínimos e máximos por PCs
pc1_range = c(min(pca_input$li[,1]), max(pca_input$li[,1]))
pc2_range = c(min(pca_input$li[,2]), max(pca_input$li[,2]))
pc3_range = c(min(pca_input$li[,3]), max(pca_input$li[,3]))

#G. save results as .CSV #Salvar os resultados em csv
#coordinates in the PCA by axis #as coordenadas das amostras
write.csv(pca_input$li, file = "Results/PCA/Axis_coord_STE.csv")
#% of PC variation #% de variação por PC
perc_pca = as.data.frame(pca_input$eig)
soma = sum (perc_pca)
perc_pca[,2] = (perc_pca/soma)*100
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
write.csv(perc_pca, file = "Results/PCA/contribution_pc_eig_STE.csv")






#H. PCA graphs #Gráficos do PCA
#create a data frame for graphs # criando uma tabela de input para os gráficos comas informações de Localidade/população
df_pca = pca_input$li %>%
  mutate(Local = geo_data_reorder$pop)
head(df_pca)  

pca_graph = ggplot(df_pca, aes(x=Axis1, y=Axis2, fill=Local))+
  geom_point(size=5, pch=21)+
  theme_bw()+
  xlab(PC1) + ylab(PC2)

#check #verificar
pca_graph

table(geo_data_reorder$pop)

#set colors and legends in alphabetical order #definir cores e legendas em ordem alfabética
legends = c("Bolivia", "Juruá River", "Madeira River")
colors = c("blue", "yellow", "red")


pca_graph2 = ggplot(df_pca, aes(x=Axis1, y=Axis2, fill=Local))+
  geom_point(size=5, pch=21)+
  theme_bw()+
  xlab(PC1) + ylab(PC2)+
  scale_fill_manual(values = colors, label = legends)+
  #geom_text(label=rownames(df_pca), nudge_x = 0.25, nudge_y = 0.25, check_overlap = F) + #add label for sample ID #adicionar o nome da amostra no gráfico
  #geom_text(label=df_pca$Local, nudge_x = 0.25, nudge_y = 0.25, check_overlap = F) + #add label for local/pop # adicionar o nome da localidade no gráfico 
  labs(fill = "Populations")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5),
        axis.title.y = element_text(size=12, face="bold", family = "Helvetica", color = "Black"),
        axis.title.x = element_text(size=12, face="bold", family = "Helvetica", color = "Black"),
        axis.text.x = element_text(size=10, family = "Helvetica", color = "Black"),
        axis.text.y = element_text(size=10, family = "Helvetica", color = "Black"),
        plot.tag = element_text(size=16, face="bold", family = "Helvetica", color = "Black"),
        legend.text = element_text(size=10, family = "Helvetica", color = "Black"),
        legend.title= element_text(size=12, face="bold", family = "Helvetica", color = "Black"),
        legend.title.align = 0.5,
        legend.box.background = element_rect(colour = "black", linewidth = 0.5))


#check #verificar
pca_graph2



#save #salvar os resultados
pdf("Results/PCA/PCA_STE.pdf")
pca_graph2
dev.off()


#END
#FIM
