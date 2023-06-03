#Rscript to analyze resonaTE step3 results

#Library
library(dplyr)
library(ggplot2)
# Import 
library(data.table)
TE_summary=fread("D:\\Documents\\Hồ sơ\\Germany\\CEPLAS\\Qualification phase\\Rotation2_Stetter\\Amaranth_TE_annotation\\Raw_data\\CandidatesB_it2(2).gff3", sep="\t",skip="seq1",
                    fill = TRUE)
TE_summary=rbind(TE_summary,fread("D:\\Documents\\Hồ sơ\\Germany\\CEPLAS\\Qualification phase\\Rotation2_Stetter\\Amaranth_TE_annotation\\Raw_data\\CandidatesB_it2(2).gff3", sep="\t",skip="ltrHarvest",
                 fill = TRUE),fill=TRUE)
# Add the length column
TE_summary[TE_summary=="seq1"] <- "chr9"
TE_summary=cbind(TE_summary,
                 Length=TE_summary$V5-TE_summary$V4,
                 density_pos= (TE_summary$V5+TE_summary$V4)/2
                 )

# Extract the transposon from the data
TE_summary_only_Transposon=TE_summary[TE_summary$V3=="transposon",]

# Create summary table

info_transposon= TE_summary_only_Transposon %>% 
  group_by(V2) %>% 
  summarise( number= dplyr::n(),
             avg_length= round(mean(Length)),
             median_length= median(Length),
             min_length= min(Length),
             max_length= max(Length),
             percentage=round(sum(Length)*100/22697164, digit=2),
             )
 info_transposon=cbind(info_transposon,TE=ID_TE)
#Export summary table
setwd("/home/ceplasrotation/TE_annotation/Results")
write.table(info_transposon, file = "info_transposon_chr9.txt", sep = ",", quote = FALSE, row.names = F)

# Creating plots 


#Boxplot for all transposon# Applying ggplot function
ggplot(TE_summary_only_Transposon, aes(x = V2, y = log(Length), fill = V2)) + 
                                   geom_boxplot() + 
                                   theme(text = element_text(size = 16), axis.text.x = element_text(angle = 70, hjust=1)) 
#Dotplot
ggplot(TE_summary_only_Transposon, aes(x = V2, y = Length, fill = V2)) +            
  geom_dotplot(binaxis = "y", stackdir = "center",dotsize = 0.5)
# Creating columm graph for transposon number
a=ggplot(info_transposon, aes( y= number, x=V2, fill = V2)) +                       
  geom_col()
a
# Create sperate graph for each transposon type
plot_all =TE_summary_only_Transposon %>% 
  group_by(V2) %>% 
  do(histogram=ggplot(data=.) 
     + aes(x=Length) 
     + geom_histogram(position = "identity", alpha = 0.7 ) 
     + ggtitle(unique(.$V2))
     + labs( x = "Length(bp)")
     + theme(text = element_text(size = 20),
             axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
             axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
     ),
     boxplot=ggplot(data=.) 
     + aes(x=V2, y=Length) 
     + geom_boxplot(stat = "boxplot",
                    position = "dodge2",) 
     + ggtitle(unique(.$V2))
     + labs(x = "Software", y = "Length(bp)")
     + theme(text = element_text(size = 20),
             axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
             axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
     )

plot_all[[3]] #Indexing from the list to show the plot
plot_all[[2]]
# Map along chromosome length


# Create input dataframe

#Create name tags for each transposon
TE_summary_only_Transposon = cbind(ID = paste(TE_summary_only_Transposon$V2, 
                                              TE_summary_only_Transposon$V9, 
                                              sep='_'),
                                   TE_summary_only_Transposon
                                   )
#transposon_summarized_only_Transposon = subset.data.frame(transposon_summarized_only_Transposon, select=-c(tmp))

#Chomosome file
Chr_file= data.frame(chr_name = c("chr9"),start_position=c(1),end_position=c(22697164))

#Annotation file
Annotation_file= data.frame(element_name=TE_summary_only_Transposon$ID, 
                            chr_name = "chr9",
                            start_position = TE_summary_only_Transposon$V4,
                            end_position= TE_summary_only_Transposon$V5,
                            software= TE_summary_only_Transposon$V2,
                            density_pos=TE_summary_only_Transposon$density_pos,
                            count=1)

#Extract the main data frame into a list of subset dataframe
software_group=split.data.frame(Annotation_file, 
                                Annotation_file$software)

# Density calculation
source("/home/ceplasrotation/TE_annotation/Code/R/calcDensity.R")

pos_info = data.frame(chr_name = "chr9",
                      position=TE_summary_only_Transposon$density_pos
                      )
density_transposon=calcDensity(pos_info,
                               Chr_file,
                               slide_window = 1e5)
density_transposon=cbind(No=1:nrow(density_transposon),density_transposon)

#Indexing from df list to calculate density for each software

software_density=lapply(software_group, select,chr_name,density_pos)

density_software = lapply(software_density, calcDensity,Chr_file,slide_window = 1e05)

ID= c("helitronScanner","ltrHarvest", "miteTracker","must","sineFinder","sineScan","TIRvish")
density_software= mapply(cbind, density_software, "software"=ID, SIMPLIFY=F)

#Merge data frame lists into big data frames
density_software_mod= do.call("rbind", density_software)

#numbering ID
Number=c(1:nrow(density_software_mod))
#create ID column
tmp2=paste(density_software_mod$software,Number,sep = "_")
#Binding
density_software_mod=cbind(ID_win=tmp2,density_software_mod)
#Receate a list
Software_chromo=split.data.frame(density_software_mod,density_software_mod$software)

#Chromo map generation
library(chromoMap)
Annotation_file=Annotation_file[,-6]
chromoMap(list(Chr_file),list(Annotation_file),          
          n_win.factor = 5,
          win.summary.display = T,
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("orange","yellow","blue","green","red","pink","navy")),
          #chr_length = 1,
          ##chr.2D.plot = T,
          #plot_filter = list(c("col","byCategory")),
          #ch2D.colors = c("orange","yellow","blue","green","red","pink","navy")
          )
 
# ChromoMap for each transposon software results



# ChromoMap for each transposon
TE_plot_1=chromoMap(list(Chr_file,Chr_file,Chr_file,Chr_file,Chr_file,Chr_file,Chr_file),
          software_group,
          ploidy = 7,
          n_win.factor = 10,
          win.summary.display = T,
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("orange","yellow","blue","green","red","pink","navy")),
          chr_length = 1,
          chr_width = 20,
          ch_gap = 20,
          chr_color = c("gray"),
          legend = c(T,F,F,F,F,F,F),
          lg_x = 880,
          lg_y = 390,
          canvas_width = 2000,
          title = "TE distribution along chromosome 9",
          title_font_size = 24,
          top_margin = 100,
          )
#pdf(file="/home/ceplasrotation/TE_annotation/Plot/TE_plot.pdf",width=10,height=10)
#TE_plot
#dev.off()


# Draw seperate plot
TE_plot_2=mapply(chromoMap,
                 list(list(Chr_file),list(Chr_file),list(Chr_file),list(Chr_file),list(Chr_file),list(Chr_file),list(Chr_file)),
                 list(list(software_group$helitronScanner),list(software_group$ltrHarvest),list(software_group$miteTracker),list(software_group$must),list(software_group$sineFinder),list(software_group$sineScan),list(software_group$TIRvish)),
                 n_win.factor = 1,
                 win.summary.display = T,
                 data_based_color_map = T,
                 data_type = "numeric",
                 aggregate_func = "sum",
                 data_colors = list(list(c("yellow","orange","red"))),
                 chr_length = 5,
                 legend = c(T,F,F,F,F,F),
                 lg_x = 500,
                 lg_y = 150
                 )
TE_plot_3=as.list(TE_plot_2)
TE_plot
TE_plot_3

deparse(substitute(Software_chromo$helitronScanner))
# Other way
for (i in 1:7){
  assign(paste("TE",ID[i],sep = "_"),
  chromoMap(list(Chr_file),
                         list(Software_chromo[[i]]),
                         ploidy = 1,
                         n_win.factor = 1,
                         win.summary.display = T,
                         data_based_color_map = T,
                         data_type = "numeric",
                         aggregate_func = "sum",
                         data_colors = list(c("yellow","orange","red")) ,
                         chr_length = 5,
                         chr_width = 20,
                         ch_gap = 20,
                         chr_color = c("gray"),
                         legend = c(T),
                         lg_x = 1350,
                         lg_y = 275,
                         canvas_width = 2000,
                         title = "Helitron along chromosome 9",
                         title_font_size = 24,
                         top_margin = 100,))
                         }

#Bar plot for each software based on density

chromoMap(list(Chr_file,Chr_file,Chr_file,Chr_file,Chr_file,Chr_file,Chr_file),
          Software_chromo,
          ploidy = 7,
          n_win.factor = 1,
          win.summary.display = T,
          data_based_color_map = T,
          data_type = "numeric",
          aggregate_func = c("sum"),
          #plots = "bar",
          #data_type = "categorical",
          data_colors = list(c("yellow","orange","red")),
          chr_length = 5,
          legend = c(T,F,F,F,F,F),
          lg_x = 500,
          lg_y = 150
          )

# Create bar plot
library(chromoMap)
chromoMap(list(Chr_file),
          list(density_transposon),
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
##

ggplot(info_transposon, 
       aes(fill=V2, y=number, x=V2)) + 
  geom_bar(position="stack", stat="identity")

print(1)
# barplot for each software
chromoMap(list(Chr_file),
          list(density_transposon),
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

#Other way - fail :(
Annotation_mod=cbind(select(Annotation_file,-software),
                     copy_no=rep(1,nrow(Annotation_file)), 
                     software=Annotation_file$software)
software_mod=split.data.frame(Annotation_mod, 
                                Annotation_mod$software)
chromoMap(list(Chr_file,Chr_file,Chr_file,Chr_file,Chr_file,Chr_file),
          software_group,
          ploidy = 6,
          n_win.factor = 10,
          win.summary.display = T,
          data_type = "categorical",
          #data_type = "numeric",
          #aggregate_func = "sum",
          data_based_color_map = T,
          #plots = "bar",
          data_colors = list(c("orange","yellow","blue","green","red","pink")),
          #chr.2D.plot = T,
          #plot_filter = list(c("col","byCategory")),
          #ch2D.colors = c("orange","yellow","blue","green","red","pink"),
          chr_length = 1,
          legend = c(T,F,F,F,F,F),
          lg_x = 70,
          lg_y = 150)

# Create density plot by ggplot and calcDensity


density_ggplot=data.frame(Chr="chr9", Window=density_transposon$Start,Count=density_transposon$Count)
ggplot(data = density_ggplot, aes(x=Window, y=Count)) + geom_col(position=position_dodge(0.5))  #+ scale_x_continuous(breaks = seq(0,25e6,5e6), labels = c("0 Mb", "15 Mb","25 Mb","25 Mb","25 Mb","25 Mb"))

Dist =density_software_mod %>% 
  group_by(software) %>% 
  do(dist=ggplot(data=.) 
     + aes(x=Start,y=Count) 
     + geom_col() 
     + ggtitle(unique(.$software))
     + labs( x = "Chr9")
     + theme(text = element_text(size = 20),
             axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
             axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
     )
  )
Dist[[2]]
