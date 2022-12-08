library(tidyverse)
library(igraph)
library(reshape2)
library(cowplot)

x <- readRDS("TARA_reference.rds")
n <-dim(x)[1]
K <- dim(x)[2] - 1
genera <- readRDS("TARA_genus_names.rds")

results <- readRDS("TARA_StARS.rds")

genus_1 <- matrix(rep(genera[-(K + 1)], times = K), nrow = K, ncol = K)
genus_2 <- matrix(rep(genera[-(K + 1)], each = K), nrow = K, ncol = K)
row_index <- matrix(rep(seq(1, K), times = K), nrow = K, ncol = K)
col_index <- matrix(rep(seq(1, K), each = K), nrow = K, ncol = K)
genus_1.unique <- genus_1[row_index > col_index]
genus_2.unique <- genus_2[row_index > col_index]

# Literature
pairs_lit <- readRDS("TARA_pairs_literature.rds")
network_lit = data.frame(genus_1 = pairs_lit$genus_org1, genus_2 = pairs_lit$genus_org2)
all_vertices = unique(c(pairs_lit$genus_org1, pairs_lit$genus_org2))

#create graph object
i.graph.lit <- graph.data.frame(network_lit, directed = F)
E(i.graph.lit)$color <- "Purple"
V(i.graph.lit)$size <- 6
#set color of nodes
V(i.graph.lit)$color <- "grey"

pdf(file = "TARA_StARS.pdf", width = 8, height = 8)

#change par to fit labels, might not be necessary for shorter labels
par(mar=c(6,6,6,6), mfrow = c(2, 2))
# par(mar=c(6,6,6,6))

#plot the graph
plot(i.graph.lit,layout=layout.circle,vertex.label="",edge.curved=0.5)

#use the layout info to position node labels
la.lit <- layout.circle(i.graph.lit)
x.lit = la.lit[,1]*1.45
y.lit = la.lit[,2]*1.45

#create vector of angles for text based on number of nodes 
#flipping the orientation of the words half way around so none appear upside down
angle.lit = ifelse(atan(-(la.lit[,1]/la.lit[,2]))*(180/pi) < 0,  
                   90 + atan(- (la.lit[,1]/la.lit[,2]))*(180/pi), 
                   270 + atan(-la.lit[,1]/la.lit[,2])*(180/pi))

#Apply the text labels with a loop with angle as srt
for (i in 1:length(x.lit)) {
  text(x=x.lit[i], y=y.lit[i], labels=V(i.graph.lit)$name[i], adj=NULL, 
       pos=NULL, cex=.7, col="black", srt=angle.lit[i], xpd=T)
}

# Figure S3
title(paste0("Literature","\n","\n","\n","\n")) #the newlines place it at the top

#reset par to default
# par(mar=c(5.1, 4.1, 4.1, 2.1))

# Compo-glasso
Omega.1.ksi <- results$Omegas.1.ksi[, , results$i.rho.1]
Omega.1.ksi.unique <- Omega.1.ksi[row_index > col_index]
Omega.1.selected <- (results$Omega.1 != 0)[row_index > col_index]
Omega.1.unique <- (cov2cor(results$Omega.1))[row_index > col_index]
Omega.1.ksi.all <- data.frame(genus_1 = genus_1.unique, genus_2 = genus_2.unique, instability = Omega.1.ksi.unique, selected = Omega.1.selected, weight = Omega.1.unique)
Omega.1.StARS.selected <- Omega.1.ksi.all %>% filter(selected == TRUE)
Omega.1.StARS.sorted <- Omega.1.StARS.selected %>% arrange(instability, desc(abs(weight))) %>% dplyr::select(-selected)
Omega.1.StARS.sorted[1:100,]
dim(Omega.1.StARS.sorted)[1]
unique(Omega.1.StARS.sorted[, "genus_1"])
dim(unique(Omega.1.StARS.sorted[, c("genus_1", "genus_2")]))

i.graph.1 <- graph.data.frame(Omega.1.StARS.sorted[1:100, ], directed = F, vertices = all_vertices)

#set the edge values to the weights
E(i.graph.1)$value <- abs((Omega.1.StARS.sorted[1:100, ])$weight)
min(E(i.graph.1)$value)
max(E(i.graph.1)$value)

#color the edges based on edge value
library(fields)
col.1 = color.scale(E(i.graph.1)$value, col=two.colors(start="lightblue", end="darkblue", middle = "blue", alpha=1.0),
                    zlim =c(0.009, 0.63), transparent.color="grey", eps= 1e-8)
E(i.graph.1)$color <- col.1

#set color of nodes
V(i.graph.1)$color <- "grey"
V(i.graph.1)$size <- 6

#plot the graph
plot(i.graph.1,layout=layout.circle,vertex.label="",edge.curved=0.5)

#use the layout info to position node labels
la.1 <- layout.circle(i.graph.1)
x.1 = la.1[,1]*1.45
y.1 = la.1[,2]*1.45

#create vector of angles for text based on number of nodes 
#flipping the orientation of the words half way around so none appear upside down
angle.1 = ifelse(atan(-(la.1[,1]/la.1[,2]))*(180/pi) < 0,  
               90 + atan(- (la.1[,1]/la.1[,2]))*(180/pi), 
               270 + atan(-la.1[,1]/la.1[,2])*(180/pi))

#Apply the text labels with a loop with angle as srt
for (i in 1:length(x.1)) {
  text(x=x.1[i], y=y.1[i], labels=V(i.graph.1)$name[i], adj=NULL, 
       pos=NULL, cex=.7, col="black", srt=angle.1[i], xpd=T)
}

#graph title
title(paste0("Comp-gLASSO","\n","\n","\n","\n")) #the newlines place it at the top

#reset par to default
# par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow = c(1, 1))

# Glasso
Omega.2.ksi <- results$Omegas.2.ksi[, , results$i.rho.2]
Omega.2.ksi.unique <- Omega.2.ksi[row_index > col_index]
Omega.2.selected <- (results$Omega.2 != 0)[row_index > col_index]
Omega.2.unique <- (cov2cor(results$Omega.2))[row_index > col_index]
Omega.2.ksi.all <- data.frame(genus_1 = genus_1.unique, genus_2 = genus_2.unique, instability = Omega.2.ksi.unique, selected = Omega.2.selected, weight = Omega.2.unique)
Omega.2.StARS.selected <- Omega.2.ksi.all %>% filter(selected == TRUE)
Omega.2.StARS.sorted <- Omega.2.StARS.selected %>% arrange(instability, desc(abs(weight))) %>% dplyr::select(-selected)
Omega.2.StARS.sorted[1:100,]
dim(Omega.2.StARS.sorted)[1]
unique(Omega.2.StARS.sorted[, "genus_1"])
dim(unique(Omega.2.StARS.sorted[, c("genus_1", "genus_2")]))

#create graph object
i.graph.2 <- graph.data.frame(Omega.2.StARS.sorted[1:100, ], directed = F, vertices = all_vertices)
E(i.graph.2)$value <- abs((Omega.2.StARS.sorted[1:100, ])$weight)
min(E(i.graph.2)$value)
max(E(i.graph.2)$value)
#color the edges based on edge value
col.2 = color.scale(E(i.graph.2)$value, col=two.colors(start="lightblue", end="darkblue", middle = "blue", alpha=1.0),
                    zlim = c(0.009, 0.63), transparent.color="white", eps= 1e-8)
E(i.graph.2)$color <- col.2

#set the size of the nodes based on degree
V(i.graph.2)$color <- "grey"
V(i.graph.2)$size <- 6

#change par to fit labels, might not be necessary for shorter labels
# par(mar=c(6,6,6,6))

#plot the graph
plot(i.graph.2,layout=layout.circle,vertex.label="",edge.curved=0.5)

#use the layout info to position node labels
la.2 <- layout.circle(i.graph.2)
x.2 = la.2[,1]*1.45
y.2 = la.2[,2]*1.45

#create vector of angles for text based on number of nodes 
#flipping the orientation of the words half way around so none appear upside down
angle.2 = ifelse(atan(-(la.2[,1]/la.2[,2]))*(180/pi) < 0,  
                 90 + atan(- (la.2[,1]/la.2[,2]))*(180/pi), 
                 270 + atan(-la.2[,1]/la.2[,2])*(180/pi))

#Apply the text labels with a loop with angle as srt
for (i in 1:length(x.2)) {
  text(x=x.2[i], y=y.2[i], labels=V(i.graph.2)$name[i], adj=NULL, 
       pos=NULL, cex=.7, col="black", srt=angle.2[i], xpd=T)
}

#graph title
title(paste0("gLASSO","\n","\n","\n","\n")) #the newlines place it at the top

#reset par to default
# par(mar=c(5.1, 4.1, 4.1, 2.1))


# MB
Omega.3.ksi <- results$Omegas.3.ksi[, , results$i.rho.3]
Omega.3.ksi.unique <- Omega.3.ksi[row_index > col_index]
Omega.3.selected <- (results$Omega.3 != 0)[row_index > col_index]
Omega.3.unique <- ((results$Omega.3))[row_index > col_index]
Omega.3.ksi.all <- data.frame(genus_1 = genus_1.unique, genus_2 = genus_2.unique, instability = Omega.3.ksi.unique, selected = Omega.3.selected, weight = Omega.3.unique)
Omega.3.StARS.selected <- Omega.3.ksi.all %>% filter(selected == TRUE)
Omega.3.StARS.sorted <- Omega.3.StARS.selected %>% arrange(instability, desc(abs(weight))) %>% dplyr::select(-selected)
Omega.3.StARS.sorted[1:100,]
dim(Omega.3.StARS.sorted)[1]
unique(Omega.3.StARS.sorted[, "genus_1"])
dim(unique(Omega.3.StARS.sorted[, c("genus_1", "genus_2")]))

# MB:
#create graph object
i.graph.3 <- graph.data.frame((Omega.3.StARS.sorted[1:100, ]), directed = F, vertices = all_vertices)
#set the edge values to the weights
E(i.graph.3)$value <- abs((Omega.3.StARS.sorted[1:100, ])$weight)

min(E(i.graph.3)$value)
max(E(i.graph.3)$value)
#color the edges based on edge value
col.3 = color.scale(E(i.graph.3)$value, col=two.colors(start="lightblue", end="darkblue", middle = "blue", alpha=1.0),
                    zlim = c(0.009, 0.63), transparent.color="white", eps= 1e-8)

E(i.graph.3)$color <- col.3


V(i.graph.3)$color <- "grey"
V(i.graph.3)$size <- 6

#plot the graph
plot(i.graph.3, layout=layout.circle,vertex.label="",edge.curved=0.5)

#use the layout info to position node labels
la.3 <- layout.circle(i.graph.3)
x.3 = la.3[,1]*1.45
y.3 = la.3[,2]*1.45

#create vector of angles for text based on number of nodes 
#flipping the orientation of the words half way around so none appear upside down
angle.3 = ifelse(atan(-(la.3[,1]/la.3[,2]))*(180/pi) < 0,  
                 90 + atan(- (la.3[,1]/la.3[,2]))*(180/pi), 
                 270 + atan(-la.3[,1]/la.3[,2])*(180/pi))

#Apply the text labels with a loop with angle as srt
for (i in 1:length(x.3)) {
  text(x=x.3[i], y=y.3[i], labels=V(i.graph.3)$name[i], adj=NULL, 
       pos=NULL, cex=.7, col="black", srt=angle.3[i], xpd=T)
}

#graph title
title(paste0("MB","\n","\n","\n","\n")) #the newlines place it at the top

dev.off()

#reset par to default
par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow = c(1, 1))

# Table S1
pairs_1 <- paste(as.character(Omega.1.StARS.sorted[1:100,]$genus_1), as.character(Omega.1.StARS.sorted[1:100,]$genus_2))
pairs_2 <- paste(as.character(Omega.2.StARS.sorted[1:100,]$genus_1), as.character(Omega.2.StARS.sorted[1:100,]$genus_2))
pairs_3 <- paste(as.character(Omega.3.StARS.sorted[1:100,]$genus_1), as.character(Omega.3.StARS.sorted[1:100,]$genus_2))
pairs_all <- intersect(intersect(pairs_1, pairs_2), pairs_3)

Omega.1.intersect = Omega.1.StARS.sorted[1:100,][pairs_1 %in% pairs_all, ] %>% arrange(instability, desc(abs(weight)))
Omega.2.intersect = Omega.2.StARS.sorted[1:100,][pairs_2 %in% pairs_all, ] %>% arrange(instability, desc(abs(weight)))
Omega.3.intersect = Omega.3.StARS.sorted[1:100,][pairs_3 %in% pairs_all, ] %>% arrange(instability, desc(abs(weight)))

intersect_from_top_100 <- full_join(Omega.1.intersect, Omega.2.intersect, by = c("genus_1", "genus_2"), copy = FALSE, suffix = c("_1", "_2")) %>%
  full_join(Omega.3.intersect, by = c("genus_1", "genus_2"), copy = FALSE) %>% 
  rename("instability_3" = "instability", "weight_3" = "weight") %>%
  mutate(sum_instability = instability_1 + instability_2 + instability_3, 
         sum_abs_weight = abs(weight_1) + abs(weight_2) + abs(weight_3)) %>%
  arrange(sum_instability, desc(sum_abs_weight))

intersect_from_top_100

pairs_lit_paste = paste(pairs_lit$genus_org1, pairs_lit$genus_org2)
paste(intersect_from_top_100$genus_1, intersect_from_top_100$genus_2) %in% pairs_lit_paste #None of the intersect edges from the top-100 lists were in the literature list!

pairs_1 %in% pairs_lit_paste
pairs_2 %in% pairs_lit_paste
pairs_3 %in% pairs_lit_paste

final_list = data.frame("genus pair" = paste(intersect_from_top_100$genus_1, intersect_from_top_100$genus_2))
library(xtable)
xtable(data.frame(final_list))

# Figure 3(b)
i.graph.1.all <- graph.data.frame(Omega.1.StARS.sorted, directed = F, vertices = all_vertices)
i.graph.2.all <- graph.data.frame(Omega.2.StARS.sorted, directed = F, vertices = all_vertices)
i.graph.3.all <- graph.data.frame(Omega.3.StARS.sorted, directed = F, vertices = all_vertices)
data_degree = data.frame(degree = c(as.numeric(degree(i.graph.1.all)), as.numeric(degree(i.graph.2.all)), as.numeric(degree(i.graph.3.all))),
                         method = rep(c("Compo-glasso", "glasso", "MB"), each = 81))
pdf(file = "TARA_degree_half.pdf", width = 6, height = 6)
p <- ggplot(data = data_degree, aes(x = degree, color = method, linetype = method, fill = method)) + 
  geom_density(alpha = 0.1) + 
  # geom_histogram(aes(y = ..density..), alhpa = 0.01, bandwidth = 1) + 
  scale_linetype_manual(values=c("solid", "dashed", "dotted"), name = "Method", labels = c('Compo-glasso', "glasso", "MB")) +
  scale_color_discrete(name = "Method") + 
  scale_fill_discrete(name = "Method") + 
  # scale_color_manual(name = "Method", values = c("blue", "red", "black"), labels = c('Compo-glasso', "glasso", "MB")) + 
  # scale_fill_manual(name = "Method", values = c("blue", "red", "black"), labels = c('Compo-glasso', "glasso", "MB")) + 
  labs(x = "Degree of Vertices", y = "Density") + 
  # theme(legend.position = "bottom")
  theme(legend.position = "none")   #remove the legend
print(p)
dev.off()

# Table 1

order_list = matrix(0, nrow = 7, ncol = 3)
degree_list = matrix(0, nrow = 7, ncol = 3)
for (i in 1:7)
{
  hub = names(sort(degree(i.graph.lit), decreasing = TRUE)[i])
  
  order_1 = which(names(sort(degree(i.graph.1.all), decreasing = TRUE)) == hub)
  order_list[i, 1] = order_1
  degree_list[i, 1] = sort(degree(i.graph.1.all), decreasing = TRUE)[order_1]
  
  order_2 = which(names(sort(degree(i.graph.2.all), decreasing = TRUE)) == hub)
  order_list[i, 2] = order_2
  degree_list[i, 2] = sort(degree(i.graph.2.all), decreasing = TRUE)[order_2]
  
  order_3 = which(names(sort(degree(i.graph.3.all), decreasing = TRUE)) == hub)
  order_list[i, 3] = order_3
  degree_list[i, 3] = sort(degree(i.graph.3.all), decreasing = TRUE)[order_3]
}
colnames(order_list) = c("Compo-glasso", "glasso", "MB")
rownames(order_list) = names(sort(degree(i.graph.lit), decreasing = TRUE)[1:7])
colnames(degree_list) = c("Compo-glasso", "glasso", "MB")
rownames(degree_list) = names(sort(degree(i.graph.lit), decreasing = TRUE)[1:7])
order_list
degree_list

combined_list = matrix("a", nrow = 7, ncol = 3)
for (i in 1:7)
{
  for (j in 1:3)
  {
    combined_list[i, j] = paste(order_list[i, j]," (",degree_list[i,j],")", sep = "")
  }
}
colnames(combined_list) = c("Compo-glasso", "glasso", "MB")
rownames(combined_list) = paste(names(sort(degree(i.graph.lit), decreasing = TRUE)[1:7]), " (",
                                sort(degree(i.graph.lit), decreasing = TRUE)[1:7], ")", sep = "")
sort(degree(i.graph.lit), decreasing = TRUE)[1:7]
xtable(combined_list)
