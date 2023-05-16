#Rscript to analyze resonaTE step3 results

#Library
library(dplyr)
library(ggplot2)

# Add the length column

transposon_summarized=cbind(CandidatesB_it2,Length=CandidatesB_it2$V5-CandidatesB_it2$V4)

# Extract the transposon from the data
transposon_summarized_only_Transposon=transposon_summarized[transposon_summarized$V3=="transposon",]

# Create summary table

info_transposon= transposon_summarized_only_Transposon %>% 
  group_by(V2) %>% 
  summarise( number= dplyr::n(),
             avg_length= round(mean(Length)),
             median_length= median(Length),
             min_length= min(Length),
             max_length= max(Length),
             percentage=round(sum(Length)*100/22697164, digit=2),
             )
#Export summary table

write.table(info_transposon, file = "info_transposon_chr9.txt", sep = ",", quote = FALSE, row.names = F)
info
# Creating plots 


#Boxplot for all transposon# Applying ggplot function
ggplot(transposon_summarized_only_Transposon, aes(x = V2, y = Length, fill = V2)) +            
  geom_boxplot()
# Creating columm graph for transposon number
ggplot(info_transposon, aes( y= number, x=V2, fill = V2)) +                       
  geom_col()

# Create sperate graph for each transposon type
plot_all =transposon_summarized_only_Transposon %>% 
  group_by(V2) %>% 
  do(histogram=ggplot(data=.) 
     + aes(Length) 
     + geom_histogram(position = "identity", alpha = 0.2 ) 
     + ggtitle(unique(.$V2)),
     boxplot=ggplot(data=.) 
     + aes(x=V2, y=Length) 
     + geom_boxplot(stat = "boxplot",
                    position = "dodge2",) 
     + ggtitle(unique(.$V2)))

plot_all[[3]]

# Map along chromosome length

#Import gff3 file
library(ape)
read.gff("D:/Documents/CandidatesB_it2.gff3",
         na.strings = c(".", "?"), 
         GFF3 = TRUE)
# Create input dataframe

#Create name tags for each transposon
tmp <- paste(transposon_summarized_only_Transposon$V2, transposon_summarized_only_Transposon$V9, sep='_')
transposon_summarized_only_Transposon <- cbind(transposon_summarized_only_Transposon, tmp)
#transposon_summarized_only_Transposon = subset.data.frame(transposon_summarized_only_Transposon, select=-c(tmp))


Chr_file= data.frame(chr_name = c("chr9"),start_position=c(1),end_position=c(22697164))
Annotation_file= data.frame(element_name=transposon_summarized_only_Transposon$tmp, 
                            chr_name = rep(c("chr9"),each=nrow(transposon_summarized_only_Transposon)),
                            start_position = transposon_summarized_only_Transposon$V4,
                            end_position= transposon_summarized_only_Transposon$V5,
                            software= transposon_summarized_only_Transposon$V2)
#Chromo map generation
library(chromoMap)

chromoMap(list(Chr_file),list(Annotation_file),          
          n_win.factor = 10,
          win.summary.display = T,
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("orange","yellow","blue","green","red","pink")),
          chr_length = 1,
          chr.2D.plot = T,
          plot_filter = list(c("col","byCategory")),
          ch2D.colors = c("orange","yellow","blue","green","red","pink"))
 
# ChromoMap for each transposon software results

#Extract the main data frame into a list of subset dataframe
software_group=split.data.frame(Annotation_file, 
                                Annotation_file$software)

# ChromoMap for each transposon
chromoMap(list(Chr_file,Chr_file,Chr_file,Chr_file,Chr_file,Chr_file),
          software_group,
          ploidy = 6,
          n_win.factor = 10,
          win.summary.display = T,
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("orange","yellow","blue","green","red","pink")),
          chr_length = 1,
          legend = c(T,F,F,F,F,F),
          lg_x = 70,
          lg_y = 150)

# Create bar plot
chromoMap(list(Chr_file),
          list(density_transposon_mod),
          n_win.factor = 1,
          win.summary.display = T,
          data_based_color_map = T,
          data_type = "numeric",
          plots = "bar",
          data_colors = list(c("white","orange","red")),
          #chr_length = 1,
          legend = c(T),
          lg_x = 70,
          lg_y = 150)
ggplot(info_transposon, aes(fill=V2, y=number, x=V2)) + 
  geom_bar(position="stack", stat="identity")

#Other way
Annotation_mod=cbind(select(Annotation_file,-software),
                     copy_no=rep(1,nrow(Annotation_file)),
                     software=Annotation_file$software)
software_mod=split.data.frame(Annotation_mod, 
                                Annotation_mod$software)
chromoMap(list(Chr_file,Chr_file,Chr_file,Chr_file,Chr_file,Chr_file),
          software_mod,
          ploidy = 6,
          n_win.factor = 10,
          win.summary.display = T,
          data_based_color_map = T,
          data_type = "numeric",
          aggregate_func = "sum",
          plots = "bar",
          data_colors = list(c("orange","yellow","blue","green","red","pink")),
          chr_length = 1,
          legend = c(T,F,F,F,F,F),
          lg_x = 70,
          lg_y = 150)

# Create density plot
library(devtools)
devtools::install_github("liuyujie0136/tinyfuncr")
1
