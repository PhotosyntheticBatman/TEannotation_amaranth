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
# Sort by TE classification
ID_TE= c("TIR","Helitron","LTR", "MITE","MITE","SINE","SINE")
ID_TE1= c("Helitron","LTR", "MITE","MITE","SINE","SINE","TIR")
TE_group=split.data.frame(TE_summary_only_Transposon, 
                          TE_summary_only_Transposon$V2)
new_TE= mapply(cbind,TE_group, "TE"=ID_TE1, SIMPLIFY=F)
new_TE= do.call("rbind", new_TE)
# Create summary table

info_transposon_group= new_TE %>% 
  group_by(TE) %>% 
  summarise( number= dplyr::n(),
             avg_length= round(mean(Length)),
             median_length= median(Length),
             min_length= min(Length),
             max_length= max(Length),
             percentage=round(sum(Length)*100/22697164, digit=2),
  )

ID= c("helitronScanner","ltrHarvest", "miteTracker","must","sineFinder","sineScan","TIRvish")

#Print summary table
#setwd("/home/ceplasrotation/TE_annotation/Results")
write.table(info_transposon, file = "info_transposon_chr9.txt", sep = ",", quote = FALSE, row.names = F)

#Boxplot for all transposon# Applying ggplot function
ggplot(new_TE, aes(x = TE, y = log(Length), fill = TE)) + 
  geom_boxplot() + 
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 70, hjust=1))

# Creating columm graph for transposon number
ggplot(info_transposon, aes( y= number, x=TE, fill = V2)) +                       
  geom_col()

# Length histogramplot
TE_his_box =new_TE %>% 
  group_by(TE) %>% 
  do(histogram=ggplot(data=.) 
     + aes(x=Length, fill=V2) 
     + geom_histogram(position = "identity", alpha = 0.7 ) 
     + ggtitle(unique(.$V2))
     + labs( x = "Length(bp)")
     + theme(text = element_text(size = 20),
             axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
             axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))
     ),
     boxplot=ggplot(data=.) 
     + aes(x=TE, y=Length, fill=TE) 
     + geom_boxplot(stat = "boxplot",
                    position = "dodge2",) 
     + ggtitle(unique(.$V2))
     + labs(x = "Software", y = "Length(bp)")
     + theme(text = element_text(size = 20),
             axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
             axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
  )
TE_his_box[[2]]
TE_his_box[[3]]

#Create name tags for each transposon
new_TE = cbind(ID = paste(new_TE$V2, 
                          new_TE$V9, 
                                              sep='_'),
               new_TE
)
#transposon_summarized_only_Transposon = subset.data.frame(transposon_summarized_only_Transposon, select=-c(tmp))

#Chomosome file
Chr_file= data.frame(chr_name = c("chr9"),start_position=c(1),end_position=c(22697164))

#Annotation file
Annotation_file_gr= data.frame(element_name=new_TE$ID, 
                            chr_name = "chr9",
                            start_position = new_TE$V4,
                            end_position= new_TE$V5,
                            count=1,
                            TE= new_TE$TE,
                            density_pos=new_TE$density_pos
                            )

#Extract the main data frame into a list of subset dataframe
TE_group=split.data.frame(Annotation_file_gr, 
                                Annotation_file_gr$TE)

# Density calculation
source("D:\\Documents\\Hồ sơ\\Germany\\CEPLAS\\Qualification phase\\Rotation2_Stetter\\Amaranth_TE_annotation\\code\\R\\calcDensity.R")

pos_info = data.frame(chr_name = "chr9",
                      position=TE_summary_only_Transposon$density_pos
)
density_gr=calcDensity(pos_info,
                               Chr_file,
                               slide_window = 1e5)
density_transposon=cbind(No=1:nrow(density_transposon),density_transposon)

#Indexing from df list to calculate density for each software

TE_gr_density=lapply(TE_group, select,chr_name,density_pos)

density_TE_gr = lapply(TE_gr_density, calcDensity,Chr_file,slide_window = 1e05)

ID3= c("Helitron","LTR", "MITE","SINE","TIR")
density_TE_gr= mapply(cbind, density_TE_gr, "TE"=ID3, SIMPLIFY=F)

#Merge data frame lists into big data frames
density_TE_gr_mod= do.call("rbind", density_TE_gr)

#numbering ID
Number2=c(1:nrow(density_TE_gr_mod))
#create ID column
tmp3=paste(density_TE_gr_mod$TE,Number2,sep = "_")
#Binding
density_TE_gr_mod=cbind(ID_win=tmp3,density_TE_gr_mod)
#Receate a list
TE_chromo=split.data.frame(density_TE_gr_mod,density_TE_gr_mod$TE)

# Draw seperate plot
library(chromoMap)


for (i in 1:5){
  assign(paste("TE",ID3[i],sep = "_"),
         chromoMap(list(Chr_file),
                   list(TE_chromo[[i]]),
                   ploidy = 1,
                   n_win.factor = 1,
                   win.summary.display = T,
                   data_based_color_map = T,
                   data_type = "numeric",
                   aggregate_func = "sum",
                   data_colors = list(c("white","orange","red")) ,
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
                   top_margin = 100,))
}
TE_Helitron
TE_LTR
TE_MITE
TE_SINE
TE_TIR
chromoMap(list(Chr_file,Chr_file,Chr_file,Chr_file,Chr_file),
          TE_chromo,
          ploidy = 5,
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
