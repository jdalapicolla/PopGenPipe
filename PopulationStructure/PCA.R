####################### VALE INSTITUTE OF TECHNOLOGY ##########################

############################     PCA TUTORIAL     #############################

#By Jeronymo Dalapicolla, 2022


##I. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT ----
rm(list=ls())


##II. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! ----
# IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H


##III. INSTALL AND LOAD THE PACKAGES ----
if (!require('adegenet'))     install.packages("adegenet");          library('adegenet')
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('dartR'))        install.packages("dartR");             library('dartR')
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')




#### GENETIC STRUCTURE USING PCA ---- 
#PCA is its ability to identify genetic structures in very large datasets within negligible computational time, and the absence of any assumption about the underlying population genetic model. PCA aims to summarize the overall variability among individuals, which includes both the divergence between groups (i.e., structured genetic variability), and the variation occurring within groups (‘random’ genetic variability).




####1. LOAD FILES ---- 
#A. Load neutral .vcf file after filtering SNPs:
vcf = read.vcfR("Inputs/vcf/lonchorhina_filtered_75ind.vcf")
#check file
vcf

#B. Convert "VCF" to "GENIND"
input = vcfR2genind(vcf)
input
input@tab[1:5,1:5]


#C. Load information on samples
geo_data = read.csv("Inputs/metadata/Lonchorhina_Sequenciados_BancodeDados.csv") 
head(geo_data)
tail(geo_data)

#select rows in the same order than vcf
genomic_order = as.data.frame(rownames(input@tab))
names(genomic_order) = "sample_name"

geo_data_reorder = semi_join(geo_data, genomic_order)
geo_data_reorder = inner_join(genomic_order, geo_data_reorder)
head(geo_data_reorder)

#verify the order:
identical(as.character(rownames(input@tab)), as.character(geo_data_reorder$sample_name))
setdiff(as.character(rownames(input@tab)), as.character(geo_data_reorder$sample_name))
setdiff(as.character(geo_data_reorder$sample_name), as.character(rownames(input@tab)))


####2. Perform a PCA ----
#D. Scale variables and run PCA 
input_scaled = scaleGen (input, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE, nf = length(row.names(input@tab))-1)

#E. % of contribution for three first PCs
PC1 = paste0("PC1 (",round(pca_input$eig[1]/sum(pca_input$eig)*100,2),"%)")
PC2 = paste0("PC2 (",round(pca_input$eig[2]/sum(pca_input$eig)*100,2),"%)")
PC3 = paste0("PC3 (",round(pca_input$eig[3]/sum(pca_input$eig)*100,2),"%)")

#F. Min and max values for axes
pc1_range = c(min(pca_input$li[,1]), max(pca_input$li[,1]))
pc2_range = c(min(pca_input$li[,2]), max(pca_input$li[,2]))
pc3_range = c(min(pca_input$li[,3]), max(pca_input$li[,3]))

#G. save results as .CSV
#coordinates in the PCA by axis
write.csv(pca_input$li, file = "Results/PCA/Axis_coord_lon_75inv.csv")
#% of PC variation
perc_pca = as.data.frame(pca_input$eig)
soma = sum (perc_pca)
perc_pca[,2] = (perc_pca/soma)*100
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
write.csv(perc_pca, file = "Results/PCA/contribution_pc_eig_lon_75ind.csv")



#H. PCA graphs
#create a data frame for graphs using Serra
df_pca = pca_input$li %>%
  mutate(Local = geo_data_reorder$serra)
head(df_pca)  

pca_graph = ggplot(df_pca, aes(x=Axis1, y=Axis2, fill=Local))+
  geom_point(size=5, pch=21)+
  theme_bw()+
  xlab(PC1) + ylab(PC2)
  
#check
pca_graph


#set colors and legends in alphabetical order
legends = c("Serra Bocaina", "Serra Norte", "Serra Sul")
colors = c("red", "blue", "yellow")


pca_graph2 = ggplot(df_pca, aes(x=Axis1, y=Axis2, fill=Local))+
  geom_point(size=5, pch=21)+
  theme_bw()+
  xlab(PC1) + ylab(PC2)+
  scale_fill_manual(values = colors, label = legends)+
 #geom_text(label=rownames(df_pca), nudge_x = 0.25, nudge_y = 0.25, check_overlap = F) +
  #geom_text(label=df_pca$Local, nudge_x = 0.25, nudge_y = 0.25, check_overlap = F) +
  labs(fill = "Populations")

#check
pca_graph2



#save
pdf("Results/PCA/PCA_Lon_75ind_serra.pdf")
pca_graph2
dev.off()

#END
