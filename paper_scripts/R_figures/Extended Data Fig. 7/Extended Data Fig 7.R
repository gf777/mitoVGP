setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(corrplot)
library(devtools)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(gridtext)
library(dplyr)

cor.pmat <- function(x, ...) {
  mat <- as.matrix(x)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# Get lower triangle of the matrix
getLower.tri<-function(mat){
  upper<-mat
  upper[upper.tri(mat)]<-""
  mat<-as.data.frame(upper)
  mat
}
# Get upper triangle of the matrix
getUpper.tri<-function(mat){
  lt<-mat
  lt[lower.tri(mat)]<-""
  mat<-as.data.frame(lt)
  mat
}
# Get flatten matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                            "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                            "#4393C3", "#2166AC", "#053061"))(200)
col<-rev(col)

DNAcolor <- function(string){
  colkmer = NULL
  for (char in unlist(strsplit(string,""),use.names=FALSE)) {
    if(char == "A"){
      colkmer = paste(colkmer,"<span style='color:blue'>A</span>", sep="")
    }else if(char == "C"){
      colkmer = paste(colkmer,"<span style='color:red'>C</span>", sep="")
    }else if(char == "T"){
      colkmer = paste(colkmer,"<span style='color:#CCCC00'>T</span>", sep="")
    }else if(char == "G"){
      colkmer = paste(colkmer,"<span style='color:#006400'>G</span>", sep="")
    }else{
      colkmer = paste(colkmer,char, sep="")
    }
  }
  return(colkmer)
}

matrix <- read.csv("matrix7.txt", header = TRUE, row.names = 1, sep = "\t", na.strings="")
kmers <- read.csv("kmers.txt", header = FALSE, sep = "\t", na.strings="")

rowcol<-data.frame(c(rep("a",19),rep("b",9),rep("c",30),rep("d",10)))
colnames(rowcol)<-"type"

# Correlation matrix
cormat<-signif(cor(matrix, use = "complete.obs"),2)
pmat<-signif(cor.pmat(matrix),2)
# Reorder correlation matrix
ord<-corrMatOrder(cormat, order="hclust")
cormat<-cormat[ord, ord]
pmat<-pmat[ord, ord]
# Replace correlation coeff by symbols
sym<-symnum(cormat, abbr.colnames=FALSE)

kmers<-kmers[order(match(kmers$V1,rownames(cormat))),]

kmers$class <- strtrim(kmers$V1, 1)

rownames(kmers)<- NULL

png(height=1200, width=1800, pointsize=15, file="heatmap.png")

ht<-rowAnnotation(
    gap = unit(8, "mm"),
    length = anno_text(gt_render(as.character(kmers$V5), align_widths = TRUE), 
                    just = "right", 
                    gp = gpar(fontsize = 16),
                    location = unit(1, "npc")),
    show_annotation_name = FALSE,
    gc_bar = anno_barplot(round(kmers$V6, digits = 2), border = FALSE, ylim = c(0, 1),
    axis_param = list(
      side = "bottom",
      at = c(0, 0.5, 1),
      labels = c("0", "50", "100"),
      labels_rot = 45,
      gp = gpar(fontsize = 12)
    )),
    #gc = anno_text(gt_render(as.character(round(kmers$V6, digits = 2)), align_widths = TRUE), 
    #                just = "right", 
    #                gp = gpar(fontsize = 16),
    #                location = unit(1, "npc")),
    kmers = anno_text(gt_render(sapply(as.character(kmers$V2), DNAcolor), align_widths = TRUE), 
                    just = "left", 
                    gp = gpar(fontsize = 16, fontfamily="Courier"),
                    location = unit(0, "npc")),
    species = anno_text(gt_render(as.character(kmers$V7), align_widths = TRUE), 
                    just = "right", 
                    gp = gpar(fontsize = 16, fontface="italic"),
                    location = unit(1, "npc"))
    )+
    rowAnnotation(show_annotation_name = FALSE, bar = as.character(kmers$class),
                col = list(bar = c("a" = "#ff0000", "m" = "#ffd700", "f" = "#00ff00", 
                                   "s" = "#1e90ff", "r" = "#ff00ff", "b" = "#66cdaa", 
                                   "k" = "black")),
                annotation_legend_param = list(
    bar = list(
      nrow = 1,
      title = "",
      at = c("a", "m", "f","r","b","k"),
      labels = c("Amphibians", "Mammals", "Fish","Reptiles","Birds","Agnatha"),
      title_gp = gpar(fontsize = 25, fontface="bold"),
      labels_gp = gpar(fontsize = 20),
      grid_height = unit(8, "mm")
    ))) +
    Heatmap(cormat,
          row_labels = rownames(cormat),
          col=col,
          row_dend_side = c("right"), 
          show_column_dend = FALSE, row_names_side = "left", 
          show_column_names = FALSE,
          heatmap_legend_param = list(legend_direction = "horizontal",
                                      title = NULL, labels_gp = gpar(fontsize = 12), grid_height = unit(5, "mm"), 
                                      at = c(-1,0,1),
                                      labels = c(round(min(matrix), digits =2), round(mean(as.matrix(matrix)), digits =2), round(max(matrix), digits =2)))
  )

draw(ht,  heatmap_legend_side = "bottom")

dev.off()

