## Load Libraries

library(vcfR)
library(adegenet)
library(poppr)
library(ape) ## plot tree
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

library(fields)
library(mapplots)
library(LEA)

## Load some LEA functions
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

library(SNPRelate)
library(dartR)
library(reshape2)
library(ggpubr)

## useful for figure arranging
library(png)
library(grid)
library(gridExtra)

## Keep the default margins before we change them
default_par_mar = par("mar")
default_par_oma = par("oma")

## Create colour blind friendly palette
#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#000000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#56B4E9")


#########################################################################################
## SET UP
#########################################################################################
## Set up the working directory
setwd("/dir/with_data")

#########################################################################################
## LOAD THE SAMPLES FOR THE POPULATION STRUCTURE ANALYSIS
#########################################################################################

## 1.1- Load the samples

OleaPanel = read.vcfR("Filtered_Olea.vcf")

## Generate different data formats and upload the population information
# 1.2- Load the population information

GL_OleaPanel = vcfR2genlight(OleaPanel)
ploidy(GL_OleaPanel_GBS.Pic) = 2
pop.data_GBS.Pic = read.csv("/Pop.data.csv", sep = ",", header = TRUE ,
                    stringsAsFactors=FALSE)
colnames(pop.data_GBS.Pic) = c("name", "origin")
# check sample names in pop file and vcf match
all(colnames(OleaPanel@gt)[-1] == pop.data$name)
pop(GL_OleaPanel) = pop.data$origin


## 1.3- Create the Structure file
## Export the GL Object to a Structure file for further analysis (it may take some time)

gl2faststructure(GL_OleaPanel, outfile = "OleaDiversityPanel.Variants.fs.str", outpath = "/dir/with_data/")


#########################################################################################
## STRUCTURE ANALYSIS WITH FASTSTRUCTURE
#########################################################################################

## 2.1- Create the geno File
struct2geno_GBS.Pic=struct2geno("OleaDiversityPanel.fs.str" , FORMAT = 2,  extra.row = 0, extra.col = 0)

## 2.2- Run the Structure analysis with number of clusters (K) from 1 to 20
obj_str = snmf("genotype.geno", K = 1:20, ploidy = 2, entropy = T, CPU = 5, project = "new")


## 2.3 - Check for the best K
plot(obj_str, col = "blue4", cex = 1.4, pch = 19)
abline(v=2, col="red")
## the cross entropy score for a given K value
cross.entropy(obj_str, K = 2)
cross.entropy(obj_str, K = 3)
plot(obj_str, col = "blue4", cex = 1.4, pch = 19 , main = "TITLE")


## plot admixture

qmatrix_k2 = Q(obj_str, K = 2)
row.names(qmatrix_k2)=GL_OleaPanel$ind.names
qmatrix_k2 <- qmatrix_k2[order(qmatrix_k2[,2],qmatrix_k2[,1] ),]
par(mar=c(2,5,1,1))
barplot(t(as.matrix(qmatrix_k2)), las =2 , horiz=T, col = admix2 , las=1, cex.names = 0.7,
        main = "TITLE")



#########################################################################################
## Neighbor-joining tree
#########################################################################################
## 3.1- Generate a distance tree

Olea_tree = aboot(GL_OleaPanel, tree = "nj", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 0, quiet = T)


rainbow(n = nPop(GL_OleaPanel))
cols = rainbow(n = nPop(GL_OleaPanel))

par(oma=c(3,2,1,1), mar=c(0,0,0,0))
par(bg="white")
plot.phylo(Olea_tree, cex = 0.8, font = 4, adj = 0, tip.color = cols[pop(GL_OleaPanel)] )
nodelabels(Olea_tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend('bottomright', legend = levels(pop(GL_OleaPanel)), fill = cols, border = FALSE, bty = "n", cex = 0.7)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")


ape::write.tree(Olea_tree, file = "./Olea.newick")



#########################################################################################
## PCA ANALYSIS
#########################################################################################

## 4.1- Run PCA analysis
pop(GL_OleaPanel) = pop.data$grp_PCA

Olea_PCA = glPca(GL_OleaPanel, nf=3)
str(Olea_PCA)
## 4.2-Plot the PCA

## Set up new margins
par(oma=c(3,3,3,3), mar=c(1,1,1,1), font=1 )

## Plot the barplot for the eigenvalues
barplot(100*Olea_PCA$eig/sum(Olea_PCA$eig), col = rep(c("red","grey"),c(3,1000)), main="TITLE", xlab="Eigenvalues", ylab="Percent of variance explained")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
Olea_PCA_scores = as.data.frame(Olea_PCA$scores)
Olea_PCA_scores$pop = pop(GL_OleaPanel)


p <- ggplot(Olea_PCA_scores, aes(x=PC1, y=PC2, colour = pop.data$origin))
p <- p + geom_point(size=3)

## To show the labels we will use the geom_text_repel function
p <- p + geom_text_repel(aes(label=gsub("OE._..._", "", rownames(Olea_PCA_scores))), show.legend = FALSE)
p <- p + scale_color_manual(values = Okabe_Ito)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + theme_classic()

## to add the PCA variance percentage use the eignvalues

# sum(Olea_PCA$eig[1]/sum(Olea_PCA$eig))*100
p <- p + xlab("PCA 1 (13.4%)")
p <- p + ylab("PCA 2 (9.7%)")
p <- p + ggtitle("TITLE")
p <- p + theme(
  plot.title = element_text(color="Black", size=12, face="bold"),
  axis.title.x = element_text(color="black", size=12, face="bold"),
  axis.title.y = element_text(color="black", size=12, face="bold"),
  legend.title = element_blank ())
p



#########################################################################################
## DAPC ANALYSIS
#########################################################################################

## 5.0 use find.clusters to identify clusters
## before checking for clusters optim.a.score to see how many pcas to include
dapc2 <- dapc(GL_OleaPanel, n.da=100, n.pca=20)
temp <- optim.a.score(dapc2)

grp <- find.clusters(GL_OleaPanel,n.pca = 6, max.n.clust=40)

names(grp)
head(grp$Kstat, 16) ## BIC values of k
grp$stat           ## selected number of clusters and associated values
head(grp$grp)      ## Group membership
grp$size           ## goup sizes

table(pop(GL_OleaPanel), grp$grp)

## write out the grouping and add to the pop data
write.csv(grp$grp, "grouped.csv")



## to check the representation of each group with a n.pca (x) use:
temp <- summary(dapc(GL_OleaPanel, n.da=100, n.pca=6))$assign.per.pop*100

barplot(temp, xlab= "% of reassignment to group",
        horiz=TRUE, las=1, main = "Group Representation at 6 PCs")

## 5.1- Run the DAPC analysis with n components and groups
Olea_DAPC <- dapc(GL_OleaPanel, n.pca = 6, n.da = 2)


## 5.2- Plot components
par(oma=c(3,1,1,1), mar=c(6,1,1,1), font=2 )

scatter(Olea_DAPC, col = Okabe_Ito, cex = 2, legend = TRUE, clabel = F, size = 4.5 , posi.leg = "bottomleft", 
        scree.pca = TRUE, posi.pca = "topright", cleg = 0.75)

compoplot(Olea_DAPC, col = Okabe_Ito, space = 0.09, posi = 'top', show.lab = TRUE, legend=TRUE )


## 5.3- Generate the barplot with populations assigned to each sector
dapc.results <- as.data.frame(Olea_DAPC$posterior)
dapc.results$pop <- pop(GL_OleaPanel)
dapc.results$indNames <- rownames(dapc.results)
dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity')
p <- p + scale_fill_manual(values = Okabe_Ito)
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p <- p + reorder(dapc.results$Posterior_membership_probability)
p

## The DAPC can be used also to infer the number of groups and then plot them
## in a similar fashion than the Structure. To do it:

## 5.4- Check the number of groups and their variation to find the optimal K

default_par_mar = par("mar")
default_par_oma = par("oma")

## Use as maximum K=40 and create a matrix with the Kstat results of find clusters 
maxK = 24
myMat = matrix(nrow=5, ncol=maxK)
colnames(myMat) = 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp = find.clusters(GL_OleaPanel, n.pca = 6, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

## Generate a data.frame for the plot
my_df = melt(myMat)
colnames(my_df)[1:3] = c("Group", "K", "BIC")
my_df$K = as.factor(my_df$K)
head(my_df)

## Plot the BIC values associated with find.clusters
p1 = ggplot(my_df, aes(x = K, y = BIC))
p1 = p1 + geom_boxplot()
p1 = p1 + theme_bw()
p1 = p1 + xlab("Number of groups (K)")
p1

## Results: It looks like the optimal number of groups is 2 from dapc and 2 from structure

## 5.5- Create the accession group assignments based in the LD 

## We will assay 2 groups
pop(GL_OleaPanel) = pop.data$origin

my_k = c(2)

grp_l = vector(mode = "list", length = length(my_k))
dapc_l = vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] = find.clusters(GL_OleaPanel, n.pca = 6, n.clust = my_k[i])
  dapc_l[[i]] = dapc(GL_OleaPanel, pop = grp_l[[i]]$grp, n.pca = 6, n.da = my_k[i])
}

my_df = as.data.frame(dapc_l[[1 ]]$ind.coord)
my_df$Group = dapc_l[[ 1 ]]$grp
head(my_df)



p2 = ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 = p2 + geom_point(size = 8, shape = 17)
p2 = p2 + theme_bw()
p2 = p2 + scale_color_manual(values=c(Okabe_Ito))
p2 = p2 + scale_fill_manual(values=c(paste(Okabe_Ito, "66", sep = "")))
p2

## 5.6- Create the accession group assignments based in the posterior probability
tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop(GL_OleaPanel)
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- pop(GL_OleaPanel)
  
  my_df <- rbind(my_df, tmp)
}


grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ Region, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
p3 <- p3 + scale_fill_manual(values=c(Okabe_Ito))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3

########################################################################################
## BASIC STATS USING DartR and POPGENOME
#########################################################################################
## Once the different individuals have assigned to different populations, last step will be to redo the population 
## assignment and calculate some population genetic parameters.


my_dfK2 = my_df[my_df$K == 2,]
K2groups = my_dfK2[my_dfK2$Posterior > 0.6,]
K2groups
length(K2groups$Isolate)

## Assign groups from the my_df object or from the pop.data datatable
pop(GL_OleaPanelCp) = pop.data$grp

## Generate the stats
GL_OleaPanelCp_BStats = gl.basic.stats(GL_OleaPanel)
GL_OleaPanelCp_BStats
GL_OleaPanelCp_BStats$Ho
## Calculate the mean Fst per group [,1], [,2]....
mean(GL_OleaPanelCp_BStats$Ho[,1])
mean(GL_OleaPanelCp_BStats$Ho[,2])
mean(GL_OleaPanelCp_BStats$Hs[,1])
mean(GL_OleaPanelCp_BStats$Hs[,2])

### It is necessary to break the VCF into sequences (for readData). To do it it is necessary to use bcftools
### in Linux:
### 1- Compress the file: "bgzip -c myfile.vcf > myfile.vcf.gz"
### 2- Index the compressed file "tabix -p vcf myfile.vcf.gz"
### 3- Create a script to divide the VCF file per sequences
###    grep -v "#" myfile.vcf | cut -f1 | sort -u | awk 'BEGIN{ print "#!/bin/bash\n\nmkdir vcf_by_seq\n"}{ print "bcftools view myfile.vcf.gz "$1" > vcf_by_seq/"$1}' > RunGetVcfBySeq.sh
### 4- Executate the script   

library(PopGenome)
## 6.1- Load the data
PG_OleaPanel = readData("variants_file", format = "VCF")


## Check the total number of SNPs
PG_OleaPanel_sumdata = as.data.frame(get.sum.data(PG_OleaPanel))
sum(PG_OleaPanel_sumdata$n.biallelic.sites)
## It gives 20,000 SNPs lesss and I am not sure why. Nevertheless, the stats should be similar

## Add the populations
K2group1list = K2groups[K2groups$Group == 1,1]
K2group2list = K2groups[K2groups$Group == 2,1]


PG_OleaPanel = set.populations(PG_OleaPanel, list(K2group1list, K2group2list), diploid = TRUE)


## 6.2- Estimate the different parameters
PG_OleaPanel_concatenated = concatenate.regions(PG_OleaPanel)
PG_OleaPanel_concatenated = neutrality.stats(PG_OleaPanel_concatenated, FAST=TRUE)
get.neutrality(PG_OleaPanel_concatenated)
PG_OleaPanel_concatenated = diversity.stats(PG_OleaPanel_concatenated)
get.diversity(PG_OleaPanel_concatenated)
PG_OleaPanel_concatenated = F_ST.stats(PG_OleaPanel_concatenated)
get.F_ST(PG_OleaPanel_concatenated)
PG_OleaPanel_concatenated = detail.stats(PG_OleaPanel_concatenated, site.spectrum=TRUE, site.FST=TRUE) 
PG_OleaPanel_results = get.detail(PG_OleaPanel_concatenated) 

## Get the segregating sites
PG_OleaPanel_concatenated@n.segregating.sites

## Get nucleotide diversity, as Pi, for each population
PG_OleaPanel_concatenated@Pi

## Get the Tajima D for each of the populations
PG_OleaPanel_concatenated@Tajima.D

## Get the Waterson Theta for each of the populations
PG_OleaPanel_concatenated@theta_Watterson


