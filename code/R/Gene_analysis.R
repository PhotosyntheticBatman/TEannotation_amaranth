library(dplyr)
library(ggplot2)
library(chromoMap)
source("D:\\Documents\\Hồ sơ\\Germany\\CEPLAS\\Qualification phase\\Rotation2_Stetter\\Amaranth_TE_annotation\\code\\R\\calcDensity.R")

# Gene annotation file
gene_chr9=Chr9_gen[grep("transcript", Chr9_gen$V3),]

#Chomosome file
Chr_file= data.frame(chr_name = c("chr9"),start_position=c(1),end_position=c(22697164))
# Position file
pos_info = data.frame(chr_name = "chr9",
                      position=(gene_chr9$V4+gene_chr9$V5)/2
)
#Density calculation

density_gene=calcDensity(pos_info,
                       Chr_file,
                       slide_window = 1e5)
density_gene=cbind(No=1:nrow(density_gene),density_gene)

#HEatmap

chromoMap(list(Chr_file),
          list(density_gene),
          n_win.factor = 1,
          win.summary.display = T,
          data_based_color_map = T,
          data_type = "numeric",
          plots = "bar",
          data_colors = list(c("white","orange","red")),
          chr_length = 5,
          chr_width = 20,
          ch_gap = 20,
          chr_color = c("gray"),
          legend = c(T),
          lg_x = 1350,
          lg_y = 275,
          canvas_width = 2000,
          #title = "Helitron along chromosome 9",
          title_font_size = 24,
          top_margin = 100)
