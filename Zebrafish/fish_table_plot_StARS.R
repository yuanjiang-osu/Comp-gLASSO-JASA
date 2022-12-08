rm(list = ls())

library(tidyverse)
library(igraph)
library(reshape2)
library(cowplot)

genera <- readRDS("genus_names.rds")
K <- length(genera) - 1

# rename the long genus name to a shorter name
genera[which(genera == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium")] <- "A-N-P-R"
genera[which(genera == "Tychonema_CCAP_1459-11B")] <- "Tychonema"

genus_1 <- matrix(rep(genera[-(K + 1)], times = K), nrow = K, ncol = K)
genus_2 <- matrix(rep(genera[-(K + 1)], each = K), nrow = K, ncol = K)
row_index <- matrix(rep(seq(1, K), times = K), nrow = K, ncol = K)
col_index <- matrix(rep(seq(1, K), each = K), nrow = K, ncol = K)
genus_1.unique <- genus_1[row_index > col_index]
genus_2.unique <- genus_2[row_index > col_index]

all_vertices <- genera[-length(genera)]

#change par to fit labels, might not be necessary for shorter labels
pdf(file = "networks.pdf", width = 8, height = 10)
par(mar = c(6, 2, 6, 2), mfcol = c(3, 2))
# par(mar=c(6,6,6,6))

groups <- c("not_infected", "infected")
labels <- c("Uninfected", "Infected")

for(g in 1 : 2)
{
  group <- groups[g]
  results <- readRDS(paste(group, "_StARS.rds", sep = ""))
  
  # Compo-glasso
  Omega.1.ksi <- results$Omegas.1.ksi[, , results$i.rho.1]
  Omega.1.ksi.unique <- Omega.1.ksi[row_index > col_index]
  Omega.1.selected <- (results$Omega.1 != 0)[row_index > col_index]
  # change sign to get partial correlations
  Omega.1.unique <- -(cov2cor(results$Omega.1))[row_index > col_index]
  Omega.1.ksi.all <- data.frame(genus_1 = genus_1.unique, genus_2 = genus_2.unique, instability = Omega.1.ksi.unique, selected = Omega.1.selected, weight = Omega.1.unique)
  Omega.1.StARS.selected <- Omega.1.ksi.all %>% filter(selected == TRUE)
  Omega.1.StARS.sorted <- Omega.1.StARS.selected %>% arrange(instability, desc(abs(weight))) %>% dplyr::select(-selected)
  Omega.1.StARS.sorted[1:100,]
  dim(Omega.1.StARS.sorted)[1]
  unique(Omega.1.StARS.sorted[, "genus_1"])
  dim(unique(Omega.1.StARS.sorted[, c("genus_1", "genus_2")]))
  
  i.graph.1 <- graph.data.frame(Omega.1.StARS.sorted[1:100,], directed = F, vertices = all_vertices)
  
  #set the edge values to the weights
  E(i.graph.1)$value <- abs((Omega.1.StARS.sorted)$weight)
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
  title(paste0(labels[g], ": Comp-gLASSO", "\n","\n","\n","\n")) #the newlines place it at the top
  
  #reset par to default
  # par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow = c(1, 1))
  
  # Glasso
  Omega.2.ksi <- results$Omegas.2.ksi[, , results$i.rho.2]
  Omega.2.ksi.unique <- Omega.2.ksi[row_index > col_index]
  Omega.2.selected <- (results$Omega.2 != 0)[row_index > col_index]
  # change sign to get partial correlations
  Omega.2.unique <- -(cov2cor(results$Omega.2))[row_index > col_index]
  Omega.2.ksi.all <- data.frame(genus_1 = genus_1.unique, genus_2 = genus_2.unique, instability = Omega.2.ksi.unique, selected = Omega.2.selected, weight = Omega.2.unique)
  Omega.2.StARS.selected <- Omega.2.ksi.all %>% filter(selected == TRUE)
  Omega.2.StARS.sorted <- Omega.2.StARS.selected %>% arrange(instability, desc(abs(weight))) %>% dplyr::select(-selected)
  Omega.2.StARS.sorted[1:100,]
  dim(Omega.2.StARS.sorted)[1]
  unique(Omega.2.StARS.sorted[, "genus_1"])
  dim(unique(Omega.2.StARS.sorted[, c("genus_1", "genus_2")]))
  
  #create graph object
  i.graph.2 <- graph.data.frame(Omega.2.StARS.sorted[1:100,], directed = F, vertices = all_vertices)
  E(i.graph.2)$value <- abs((Omega.2.StARS.sorted)$weight)
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
  title(paste0(labels[g],": gLASSO", "\n","\n","\n","\n")) #the newlines place it at the top
  
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
  i.graph.3 <- graph.data.frame((Omega.3.StARS.sorted), directed = F, vertices = all_vertices)
  #set the edge values to the weights
  E(i.graph.3)$value <- abs((Omega.3.StARS.sorted)$weight)
  
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
  title(paste0(labels[g], ": MB","\n","\n","\n","\n")) #the newlines place it at the top

  write.csv(Omega.1.StARS.sorted, paste(group, "_genus_pairs_Comp-gLASSO.csv", sep = ""), row.names = FALSE)
  write.csv(Omega.2.StARS.sorted, paste(group, "_genus_pairs_gLASSO.csv", sep = ""), row.names = FALSE)
  write.csv(Omega.3.StARS.sorted, paste(group, "_genus_pairs_MB.csv", sep = ""), row.names = FALSE)
  
  # Calculate the degree of vertex
  i.graph.1.all <- graph.data.frame(Omega.1.StARS.sorted, directed = F, vertices = all_vertices)
  i.graph.2.all <- graph.data.frame(Omega.2.StARS.sorted, directed = F, vertices = all_vertices)
  i.graph.3.all <- graph.data.frame(Omega.3.StARS.sorted, directed = F, vertices = all_vertices)
  write.csv(sort(degree(i.graph.1.all), decreasing = TRUE), paste(group, "_genus_degrees_Comp-gLASSO.csv", sep = ""), row.names = TRUE)
  write.csv(sort(degree(i.graph.2.all), decreasing = TRUE), paste(group, "_genus_degrees_gLASSO.csv", sep = ""), row.names = TRUE)
  write.csv(sort(degree(i.graph.3.all), decreasing = TRUE), paste(group, "_genus_degrees_MB.csv", sep = ""), row.names = TRUE)
}

dev.off()

# change g to 2 for infected group
for(g in 1 : 2)
{
  group <- groups[g]
  results <- readRDS(paste(group, "_StARS.rds", sep = ""))
  
  # Compo-glasso
  Omega.1.ksi <- results$Omegas.1.ksi[, , results$i.rho.1]
  Omega.1.ksi.unique <- Omega.1.ksi[row_index > col_index]
  Omega.1.selected <- (results$Omega.1 != 0)[row_index > col_index]
  # change sign to get partial correlations
  Omega.1.unique <- -(cov2cor(results$Omega.1))[row_index > col_index]
  Omega.1.ksi.all <- data.frame(genus_1 = genus_1.unique, genus_2 = genus_2.unique, instability = Omega.1.ksi.unique, selected = Omega.1.selected, weight = Omega.1.unique)
  Omega.1.StARS.selected <- Omega.1.ksi.all %>% filter(selected == TRUE)
  Omega.1.StARS.sorted <- Omega.1.StARS.selected %>% arrange(instability, desc(abs(weight))) %>% dplyr::select(-selected)
  
  # Glasso
  Omega.2.ksi <- results$Omegas.2.ksi[, , results$i.rho.2]
  Omega.2.ksi.unique <- Omega.2.ksi[row_index > col_index]
  Omega.2.selected <- (results$Omega.2 != 0)[row_index > col_index]
  # change sign to get partial correlations
  Omega.2.unique <- -(cov2cor(results$Omega.2))[row_index > col_index]
  Omega.2.ksi.all <- data.frame(genus_1 = genus_1.unique, genus_2 = genus_2.unique, instability = Omega.2.ksi.unique, selected = Omega.2.selected, weight = Omega.2.unique)
  Omega.2.StARS.selected <- Omega.2.ksi.all %>% filter(selected == TRUE)
  Omega.2.StARS.sorted <- Omega.2.StARS.selected %>% arrange(instability, desc(abs(weight))) %>% dplyr::select(-selected)
  
  # MB
  Omega.3.ksi <- results$Omegas.3.ksi[, , results$i.rho.3]
  Omega.3.ksi.unique <- Omega.3.ksi[row_index > col_index]
  Omega.3.selected <- (results$Omega.3 != 0)[row_index > col_index]
  Omega.3.unique <- ((results$Omega.3))[row_index > col_index]
  Omega.3.ksi.all <- data.frame(genus_1 = genus_1.unique, genus_2 = genus_2.unique, instability = Omega.3.ksi.unique, selected = Omega.3.selected, weight = Omega.3.unique)
  Omega.3.StARS.selected <- Omega.3.ksi.all %>% filter(selected == TRUE)
  Omega.3.StARS.sorted <- Omega.3.StARS.selected %>% arrange(instability, desc(abs(weight))) %>% dplyr::select(-selected)
  
  # Calculate the degree of vertex
  i.graph.1.all <- graph.data.frame(Omega.1.StARS.sorted, directed = F, vertices = all_vertices)
  i.graph.2.all <- graph.data.frame(Omega.2.StARS.sorted, directed = F, vertices = all_vertices)
  i.graph.3.all <- graph.data.frame(Omega.3.StARS.sorted, directed = F, vertices = all_vertices)
  
  data_degree = data.frame(degree = c(as.numeric(degree(i.graph.1.all)), as.numeric(degree(i.graph.2.all)), as.numeric(degree(i.graph.3.all))),
                           method = rep(c("Compo-glasso", "glasso", "MB"), each = length(degree(i.graph.1.all))))
  
  pdf(file = paste(group, "_degree.pdf", sep = ""), width = 6, height = 6)
  
  p <- ggplot(data = data_degree, aes(x = degree, color = method, linetype = method, fill = method)) + 
    geom_density(alpha = 0.1) + 
    # geom_histogram(aes(y = ..density..), alhpa = 0.01, bandwidth = 1) + 
    scale_linetype_manual(values=c("solid", "longdash", "twodash"), name = "Method", labels = c('Compo-glasso', "glasso", "MB")) +
    scale_color_discrete(name = "Method") + 
    scale_fill_discrete(name = "Method") + 
    # scale_color_manual(name = "Method", values = c("blue", "red", "black"), labels = c('Compo-glasso', "glasso", "MB")) + 
    # scale_fill_manual(name = "Method", values = c("blue", "red", "black"), labels = c('Compo-glasso', "glasso", "MB")) + 
    labs(x = "Degree of Vertices", y = "Density") + 
    # theme(legend.position = "bottom") +
    theme_bw() +
    theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))   #remove the legend
  print(p)

  dev.off()
}