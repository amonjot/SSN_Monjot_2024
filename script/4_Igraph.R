# !/usr/bin/env Rscript
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 15/05/2023
#
# Script Metatrans_Table analyze
Cstack_info()["size"]
# Set input and output -----------------------------------------------------------
# Detect R or Rstudio
se <- Sys.getenv()
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) > 0 ) { inputmode <- TRUE }
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) == 0 ) { inputmode <- FALSE }
#
# Input argument if using Rstudio
if (inputmode == TRUE) {
input <- "80_80_1e-50"
}
#
# Input argument if using R
args = commandArgs(trailingOnly=TRUE)
if ( inputmode == FALSE ) {
  input <- args[1]
}
print(input)
#
# Import package -----------------------------------------------------------
pkg <- c("ggplot2","plyr","dplyr","tidyr","stringr","cowplot","elementalist","svglite","ggupset","ggstatsplot","treemapify","paletteer","fmsb","rstatix","ggpubr","FactoMineR","factoextra","ggrepel","igraph","ggraph","tidygraph")
lapply(pkg, require, character.only = TRUE)
#
# Set directory, create file result -----------------------------------------------------------
if (inputmode == TRUE) {
  current <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current)
  output <- paste(current,"../result/Lagoon_output",input, sep = "/")
}
if (inputmode == FALSE) {
  output <- paste("result/Lagoon_output",input, sep = "/")
}
# Set working directory
if (dir.exists(output) == FALSE) {dir.create(output)}
if (dir.exists(output) == TRUE) { setwd(output) }
# Create result files
  dir.create("Igraph")
  dir.create("Analysis")
#
# Theme unique Dark perso -------------------------------------------------------
  theme_unique_darkcris <- function (base_size = 12, base_family = "") {
    ret <- (theme_bw(base_size = base_size, base_family = base_family) +
              theme(text = element_text(colour = "black"),
                    title = element_text(color = "black", face = "bold"),
                    axis.ticks = element_blank(),
                    line = element_line(color = "black"),
                    rect = element_rect(fill = "white", color = "black"),
                    axis.title = element_blank(),
                    #axis.text.y = element_blank(),
                    #axis.text.x = element_blank(),
                    axis.text.x = element_text(color = "black", size = 10, vjust = 1, hjust=1,angle = 45,face="bold"),
                    #axis.line = element_line(color = "#969696", linetype = 1),
                    legend.background = element_blank(),
                    panel.background = element_rect_round(radius = unit(0.5, "cm"),fill = "white", color = NULL,size = 1),
                    #legend.position = "bottom",
                    legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 11, face="bold"),
                    #strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                    strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5),
                    strip.text.y = element_text(angle=0),
                    #panel.background = element_rect(fill = "white", color = NULL),
                    panel.border = element_blank(),
                    panel.grid = element_line(color = "#252525"),
                    panel.grid.major = element_line(color = "white",linetype="blank"),
                    panel.grid.minor = element_line(color = "#252525",linetype="dotted"),
                    plot.background = element_rect(fill = "white", colour = "#252525", linetype = 0)))
    ret
  } 
  theme_unique_darktris <- function (base_size = 12, base_family = "") {
    ret <- (theme_bw(base_size = base_size, base_family = base_family) +
              theme(text = element_text(colour = "black"),
                    title = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5,size = 10),
                    axis.ticks = element_blank(),
                    line = element_line(color = "black"),
                    rect = element_rect(fill = "white", color = "black"),
                    axis.title = element_text(color = "black", face = "bold"),
                    axis.text.y = element_blank(),
                    axis.text.x = element_blank(),
                    legend.background = element_rect(fill = NULL, color = NULL),
                    legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                    legend.text = element_text(size = 12,face ="bold"),
                    legend.title = element_text(size = 0, face="bold"),
                    strip.background = element_rect(fill=NA,colour=NA,linewidth=NA,linetype = NULL),
                    strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5,size = 14),
                    #strip.text.y.left = element_text(angle=0),
                    panel.background = element_rect(fill = "white", color = NULL),
                    panel.border = element_blank(),
                    panel.grid.major = element_line(color = "white"),
                    panel.grid.minor = element_line(color = "white"),
                    plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
    ret
  } 
  theme_unique_art <- function (base_size = 12, base_family = "") {
    ret <- (theme_minimal(base_size = base_size, base_family = base_family) +
              theme(text = element_text(colour = "black"),
                    title = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
                    axis.ticks = element_line(color = "lightgrey", linetype = 1),
                    line = element_line(color = "black"),
                    rect = element_rect(fill = "white", color = "black"),
                    axis.title = element_text(color = "black", face = "bold"),
                    #axis.text.y = element_blank(),
                    #axis.text.x = element_blank(),
                    #axis.line = element_line(color = "#969696", linetype = 1),
                    legend.background = element_blank(),
                    legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                    legend.text = element_text(size = 8),
                    legend.title = element_text(size = 12, face="bold"),
                    panel.background = element_rect_round(radius = unit(0.5, "cm"),fill = "white", color = NULL,size = 1),
                    #legend.position = "none",
                    #strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                    strip.text = element_text(size = 12,color="black",face="bold",vjust=.5,hjust=.5),
                    #panel.background = element_rect(fill = "white", color = NULL),
                    panel.border = element_blank(),
                    panel.grid = element_line(color = "#252525"),
                    panel.grid.major = element_line(color = "grey86"),
                    panel.grid.minor = element_line(color = "grey86"),
                    plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
    ret
  }
  theme_unique_art_graph <- function (base_size = 12, base_family = "") {
    ret <- (theme_minimal(base_size = base_size, base_family = base_family) +
              theme(text = element_text(colour = "black"),
                    title = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
                    axis.ticks = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    line = element_line(color = "black"),
                    rect = element_rect(fill = "white", color = "black"),
                    legend.background = element_blank(),
                    legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                    legend.text = element_text(size = 8),
                    legend.title = element_text(size = 12, face="bold"),
                    panel.background = element_rect_round(radius = unit(0.5, "cm"),fill = "grey", color = NULL,size = 1),
                    strip.text = element_text(size = 12,color="black",face="bold",vjust=.5,hjust=.5),
                    panel.border = element_blank(),
                    plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
    ret
  }
#
# Graph_input -----------------------------------------------------------------
# Edges "from to + stat"
edges_input <- paste0("diamond_ssn_",input,".edges")
edges_input_df <- read.csv(edges_input,sep = ";")
# All input proteins with metadata
attributes_input <- "../Metadata_Unicellular.attributes"
attributes_input_df <- read.csv(attributes_input,sep = ";")

# Result of clusterisation: CCs + Metadata
stat_input <- paste0("ssn_composantes_connexes_graph_",input,".csv")
stat_input_df <- read.csv(stat_input,sep = ";")
stat_input_df$CC <- stat_input_df$CC+1
#
stat_input_df[,"TRT_1"][stat_input_df[,"TRT_2"]=="PARA"] <- "PARA" # use parasitism as trophic mode
stat_input_df[,"TRT_1"][stat_input_df[,"TRT_1"]=="SAP/HET"] <- "HET" # normalize trophic modes (Uncertainty about sap and het become het)
stat_input_df[stat_input_df==""]<-"Unassigned"
stat_input_df[stat_input_df=="NA"]<-"Unassigned"
stat_input_df[,"TRT_1"][stat_input_df[,"TRT_1"]=="SAP" & stat_input_df[,"TYP_1"]=="Obligate"]<- "PARA"
#
# Input and format kegg database ------------------------------------------
ko_db <- read.csv("../../../rawdata/ko_to_hierarchy.txt",sep="\t")
ko_db <- ko_db %>% filter(lvl_A_val %in% c("Metabolism"))
#
# General_statistic -------------------------------------------------------
for (input_General in c("Total","In_CC")) {
  if (input_General == "Total") {tablearchive <- attributes_input_df}
  if (input_General == "In_CC") {tablearchive <- attributes_input_df %>% filter(ID_Protein %in% stat_input_df$name)}
  # Process -----------------------------------------------------------------
# Count Transcript by Taxref
tablearchive[tablearchive==""] <- NA
tablearchive[,"TRT_1"][tablearchive[,"TRT_2"]=="PARA"] <- "PARA" # use parasitism as trophic mode
tablearchive[,"TRT_1"][tablearchive[,"TRT_1"]=="SAP/HET"] <- "HET" # normalize trophic modes (Uncertainty about sap and het become het)
tablearchive[,"TRT_1"][tablearchive[,"TRT_1"]=="SAP" & tablearchive[,"TYP_1"]=="Obligate"]<- "PARA" # normalize trophic modes (obligate sap become para)
nrow(tablearchive %>% filter(TRT_1=="PARA")) # 24 232 < 104 922
nrow(tablearchive %>% filter(TRT_1=="HET")) # 77 955 < 281 370
nrow(tablearchive %>% filter(TRT_1=="PHOTO")) # 136 946 < 397 999
nrow(tablearchive %>% filter(TRT_1=="SAP")) # 5 952 < 31 362
nrow(tablearchive %>% filter(TRT_1=="MIXO")) # 156 485 < 528 433
nrow(tablearchive %>% filter(is.na(TRT_1)==TRUE)) # 1 763 536 < 6 138 708
nrow(tablearchive %>% filter(is.na(TRT_1)==FALSE)) # 401 570 < 1 344 092
nrow(tablearchive) # 2165106 < 7482800
#
# Affiliation
tablearchivex <- tablearchive %>% select(Taxonomy) %>% distinct(Taxonomy,.keep_all = TRUE)
for (i in row.names(tablearchivex)) { 
  j = 2
  while (str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][j] == "") { j<-j+1
  if (j == 3) {break}}
  k = j+1
  while (str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][k] == "") { k<-k+1
  if (k >= 5) {break}}
  l = k+1
  while (str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][l] == "") { l<-l+1
  if (l >= 7) {break}}
  tablearchivex[i,"Division"] <- str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][j]
  tablearchivex[i,"Phylum"] <- str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][k]
  tablearchivex[i,"Class"] <- str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][l]}
tablearchivex[tablearchivex==""]<-NA
tablearchivex[tablearchivex=="NA"]<-NA
#
#  Process Inconsistencies of taxonomy
summarize_data <- tablearchivex %>% select(Taxonomy,Division,Phylum,Class)
summarize_data[,"Division"][summarize_data[,"Division"]=="unclassified eukaryotes"] <- NA
summarize_data[,"Division"][summarize_data[,"Division"]=="environmental samples"] <- "Environmental samples"
summarize_data[,"Division"][summarize_data[,"Phylum"]=="Jakobida"] <- "Jakobida"
summarize_data[,"Division"][summarize_data[,"Phylum"]=="Malawimonadidae"] <- "Malawimonadidae"
#
# Prepare color for general "Total" stat figure and attribute taxonomy
if (input_General == "Total") {color_definition <- sort(unique(summarize_data$Division))}
summarize_datax <- left_join(tablearchive,summarize_data,"Taxonomy")
#
summarize_datax[,"TYP_1"][summarize_datax[,"TRT_1"]=="MIXO"]<- NA # Remove description of Mixotroph
#
colnames(summarize_datax)[colnames(summarize_datax)=="TRT_1"] <- "Grp_1"
colnames(summarize_datax)[colnames(summarize_datax)=="TYP_1"] <- "Grp_2"
#
  # Prepare Plot ------------------------------------------------------------
    # Prepare Polar data -----------------------------------------------
Polar_data_all <- summarize_datax %>%
  arrange(Division) %>%
  select(Grp_1,Grp_2,Division,Phylum,Class,) %>%
  group_by(Grp_1,Grp_2,Division,Phylum,Class) %>% count(name="count")
#
Polar_data_all[is.na(Polar_data_all)==TRUE] <- "Unassigned"
Polar_data_all[,"Division"][Polar_data_all[,"Division"]=="Unassigned"] <- "Unaffiliated"
Polar_data_all[,"Phylum"][Polar_data_all[,"Phylum"]=="Unassigned"] <- "Unaffiliated"
Polar_data_all[,"Class"][Polar_data_all[,"Class"]=="Unassigned"] <- "Unaffiliated"
# To complete Fig_1
data_summary_para <- as.data.frame(Polar_data_all %>% filter(Grp_1=="PARA"))
data_summary_para$percent <- ""
for (i in row.names(data_summary_para)) {
  Division_i <- data_summary_para[i,"Division"]
  Grp_1_i <- data_summary_para[i,"Grp_1"]
  Grp_2_i <- data_summary_para[i,"Grp_2"]
  data_summary_para[i,"percent"] <- data_summary_para[i,"count"]*100/sum(data_summary_para[,"count"][data_summary_para[,"Division"]==Division_i & data_summary_para[,"Grp_1"]==Grp_1_i & data_summary_para[,"Grp_2"]==Grp_2_i])
}
write.table(data_summary_para,paste0("Analysis/data_summary_para",input_General,".csv"),row.names = FALSE,sep="\t")
#
Polar_data_all <- as.data.frame(Polar_data_all,stringsAsFactors=FALSE)
Proportion_aff_x <- Polar_data_all %>% select(-Phylum,-Class) %>% group_by(Grp_1,Grp_2,Division) %>% summarise_all(sum)
Prop_value <- 100-Proportion_aff_x$count[Proportion_aff_x[,"Division"]=="Unaffiliated"] * 100 / sum(Proportion_aff_x$count)
Proportion_aff <- paste("Affiliated:",round(Prop_value,1),"%")
Polar_data_all <- Polar_data_all %>% filter(Division != "Unaffiliated")
Polar_data <- Polar_data_all %>% select(-Grp_2,-Phylum,-Class) %>% group_by(Grp_1,Division) %>% summarise_all(sum)
Polar_data <- as.data.frame(Polar_data,stringsAsFactors=FALSE)
for (i in rownames(Polar_data)) { 
  Polar_data[i,"Proportion"] <- (Polar_data[i,"count"] * 100) / sum(Polar_data$count)
  Division_i <- Polar_data[i,"Division"]
  Polar_data[i,"Proportion_byDiv"] <- (Polar_data[i,"count"] * 100) / colSums(Polar_data %>% filter(Division == Division_i) %>% select(count))
}
#
    # Polar_data_all (Grp_2) -----------------------------------------------
Polar_data_all <-  Polar_data_all %>% select(-Phylum,-Class) %>% group_by(Grp_1,Grp_2,Division) %>% summarise_all(sum)
for (i in rownames(Polar_data_all)) { 
  Polar_data_all[i,"Proportion"] <- (Polar_data_all[i,"count"] * 100) / sum(Polar_data_all$count)
}
Polar_data_all <- Polar_data_all %>% arrange(Division)
Polar_data_all$color <- NA
Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="Obligate"] <- "lightgreen"
Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="Facultative"] <- "lightblue"
Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="CM"] <- "#CCCCCC"
Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="eSNCM"] <- "darkgrey"
Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="pSNCM"] <- "#646464"
Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="eSNCM/pSNCM"] <- "#333333"
Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="Unassigned"] <- "white"
Polar_data_all <- Polar_data_all %>%
  arrange(Division) %>%
  mutate(lab.ypos = cumsum(Proportion) - 0.5*Proportion)
Polar_data_all$label <- paste(round(Polar_data_all$Proportion,1), "%", sep = "")
for (i in rownames(Polar_data_all)) {
  if (Polar_data_all[i,"label"] == "0%") { Polar_data_all[i,"label"] <- NA}}
for (i in rownames(Polar_data_all)) {
  if (is.na(Polar_data_all[i,"label"]) == FALSE) { Polar_data_all[i,"label"] <- paste0(Polar_data_all[i,"label"])}
}
#
    # Polar_data (Grp_1) -----------------------------------------------
Polar_data <- Polar_data %>%
  arrange(Division) %>%
  mutate(lab.ypos = cumsum(Proportion) - 0.5*Proportion)
Polar_data$label <- paste(round(Polar_data$Proportion_byDiv,1), "%", sep = "")
for (i in rownames(Polar_data)) {
  if (Polar_data[i,"label"] == "0%") { Polar_data[i,"label"] <- NA}}
for (i in rownames(Polar_data)) {
  if (is.na(Polar_data[i,"label"]) == FALSE) { Polar_data[i,"label"] <- paste0(Polar_data[i,"label"])}
}
Polar_data[,"label"][Polar_data[,"Grp_1"]=="Unassigned"]<-NA
#
    # Polar_data_1 (Division) -----------------------------------------------
color <- data.frame(color_definition,c("#c62728cc","#6b1b9acc","white","#273493cc","white","#0277bdcc","#558b30cc","#f9a825cc","white","#fdffcccc","#ef6c00cc","white","white","white","white","#75252595","#880e4fcc","#90ee90cc","#20606490","#add8e6cc","#e53368cc"))
colnames(color) <- c("Division","color")
#
Polar_data_a <- Polar_data %>% select(-Grp_1,-label) %>% 
  arrange(Division) %>%
  group_by(Division) %>% 
  summarise_all(sum)

for (i in row.names(Polar_data_a)) { if (Polar_data_a[i,"Proportion"] > 0.1) {
  Polar_data_a[i,"color"] <- color[,"color"][color[,"Division"]==Polar_data_a[i,"Division"]$Division]}
  else {Polar_data_a[i,"color"] <- "white"}
}
#
Polar_data$color <- NA
Polar_data[,"color"][Polar_data[,"Grp_1"]=="HET"] <- "grey48"
Polar_data[,"color"][Polar_data[,"Grp_1"]=="Unassigned"] <- "white"
Polar_data[,"color"][Polar_data[,"Grp_1"]=="PHOTO"] <- "blue4"
Polar_data[,"color"][Polar_data[,"Grp_1"]=="MIXO"] <- "red4"
Polar_data[,"color"][Polar_data[,"Grp_1"]=="SAP"] <- "grey0"
Polar_data[,"color"][Polar_data[,"Grp_1"]=="PARA"] <- "darkgreen"
#
    # legend ------------------------------------------------------------------
colnames(Polar_data)[colnames(Polar_data)=="Grp_1"] <- "Trophic mode"
colnames(Polar_data_all)[colnames(Polar_data_all)=="Grp_2"] <- "Parasitic lifestyle"
colnames(Polar_data_all)[colnames(Polar_data_all)=="Grp_1"] <- "Trophic mode"
colnames(Polar_data_a)[colnames(Polar_data_a)=="Division"] <- "Taxonomic group"
#
#Legend
Polar_data_b <- Polar_data_a %>% filter(Proportion > 0.1)
Polar_data_b <- Polar_data_b %>% mutate(`Taxonomic group`=paste0(`Taxonomic group`,": ",round(Proportion,1),"%"))
legend_a <- get_legend(ggplot(data = Polar_data_b) + theme_unique_darkcris() +
                         geom_bar(mapping = aes(y= Proportion, fill = `Taxonomic group`)) + guides(fill = guide_legend(override.aes = list(color="black",shape=1))) + 
                         scale_fill_manual(values=noquote(Polar_data_b$color)))
Polar_data_legend <- Polar_data %>% select(`Trophic mode`,Proportion) %>% group_by(`Trophic mode`) %>% summarise_all(sum)
Polar_data_legend <- Polar_data_legend %>% mutate(`Trophic mode` = paste0(`Trophic mode`,": ",round(Proportion,1),"%"))
legend_b <- get_legend(ggplot(data = Polar_data_legend %>% arrange(`Trophic mode`)) + theme_unique_darkcris() + guides(fill="Functional assignment") +
                         geom_bar(mapping = aes(y= Proportion, fill = `Trophic mode`)) + guides(fill = guide_legend(override.aes = list(color="black",shape=1))) + 
                         scale_fill_manual(values=noquote(unique((Polar_data %>% arrange(`Trophic mode`))$color))))
Polar_data_all_legend <- as.data.frame(Polar_data_all) %>% select(`Parasitic lifestyle`,Proportion) %>% group_by(`Parasitic lifestyle`) %>% summarise_all(sum)
Polar_data_all_legend <- Polar_data_all_legend %>% mutate(`Parasitic lifestyle` = paste0(`Parasitic lifestyle`,": ",round(Proportion,1),"%"))
legend_c <- get_legend(ggplot(data = Polar_data_all_legend %>% arrange(`Parasitic lifestyle`)) + theme_unique_darkcris() + guides(fill="Parasitic lifestyle") +
                         geom_bar(mapping = aes(y= Proportion, fill = `Parasitic lifestyle`)) + guides(fill = guide_legend(override.aes = list(color="black",shape=1))) + 
                         scale_fill_manual(values=noquote((Polar_data_all %>% select(`Parasitic lifestyle`, color) %>% distinct() %>% arrange(`Parasitic lifestyle`,.locale="en"))$color)))
#
  # plot --------------------------------------------------------------------
#Plot large graph
Polar_data_all <- as.data.frame(Polar_data_all)
Polar_data_a <- as.data.frame(Polar_data_a)
Polar_data <- as.data.frame(Polar_data)
bxCop <- ggplot(Rowv = NA, col = colMain, scale = "column") +
  theme_unique_darkcris() +
  labs(y = "",x="") + coord_polar("y") + 
  geom_bar(data = Polar_data_all, mapping = aes(y= Proportion, x = 2.85, fill = label,width = 0.1), stat="identity", color = "black", width = 1, linewidth = 0.3, fill=Polar_data_all[,"color"]) + 
  geom_bar(data = Polar_data, mapping = aes(y= Proportion, x = 2.65, fill = label,width = 0.3), stat="identity", color = "black", width = 1, linewidth= 0.3, fill=Polar_data[,"color"]) + 
  scale_y_continuous(limits=c(0,sum(Polar_data %>% select(Proportion))),labels = scales::percent_format(scale=1),breaks = c(25,50,75,100)) +
  geom_bar(data = Polar_data_a,  mapping = aes(y= Proportion, x = 2, fill = rev(Division),width = 1), stat="identity", color = "black", width = 1, linewidth =0.3,fill=Polar_data_a[,"color"]) + 
  xlim(1,3) + theme(axis.text.y=element_blank(),
                    axis.ticks.y=element_blank()) +
  annotate("text",label=paste0(Proportion_aff,"\n(",sum(Proportion_aff_x$count),")"),x=1,y=25, fontface =2)
print(bxCop)
#
svglite(paste0("Analysis/Stat_",input_General,".svg"),width = 10.00,height = 8.00)
legend <- plot_grid(NULL,legend_a, legend_b,legend_c, NULL, ncol = 1, nrow = 5,align = "v",rel_heights = c(0.3,2.3,1.5,0.5,0.5))
b_plot <- plot_grid(bxCop, legend, ncol = 2, nrow = 1 , rel_widths = c(6,2),rel_heights = c(3))
print(b_plot)
dev.off()
#
# plot annot / no annot
axCop_data <- data.frame(Value=c(Prop_value,100-Prop_value),Variable=c("Affiliated","Unaffiliated"),color=c("black","white"))
axCop_data <- axCop_data %>% mutate(Statistics=paste0(Variable,": ",round(Value,1)," %"))
svglite(paste0("Analysis/Stat_bis_",input_General,".svg"),width = 2.00,height = 3.10)
axCop <- ggplot(Rowv = NA, col = colMain, scale = "column") +
  theme_void() + theme(legend.position.inside = c(0.45,-0.1)) +
  labs(y = "",x="") + coord_polar("y",direction = -1) + scale_fill_manual(values=axCop_data$color) +
  geom_bar(data = axCop_data, mapping = aes(y= Value, x = 1, fill = Statistics), stat="identity", color = "black", width = 1, linewidth =0.5) + 
  scale_y_continuous(limits=c(0,sum(axCop_data %>% select(Value)))) +
  theme(axis.text.y=element_blank(),legend.text = element_text(size = 10),axis.ticks.y=element_blank(),legend.title = element_text(size = 11, face="bold"))
print(axCop)
dev.off()
legend_taxo <- get_legend(axCop)
}
#
# Stat on CC --------------------------------------------------------------
  # Prot by CC --------------------------------------------------------------
table_CC <- stat_input_df %>% select(CC) %>% 
    mutate(CC=paste0("CC_",CC)) %>% 
    mutate(length=1) %>% group_by(CC) %>% 
    summarise_all(sum)
  table_CC_mean <- table_CC %>% select(-CC) %>% summarise_all(mean) %>% mutate(length_label=paste("Mean:",round(length,2)))
  table_CC_max <- table_CC %>% select(-CC) %>% summarise_all(max) %>% mutate(length_label=paste("Max:",round(length,2)))
  
  #plot
  c <- ggplot(data = table_CC,aes(x = length,y=1)) +
    geom_boxplot(width = 0.1,fill = NA,color = "black",alpha = 0.2) + 
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 2),breaks = c(0,1,2,4,8,16,32,64,256,512,2048,16368)) +
    scale_x_continuous(breaks = c(3,200,400,800,1600,3200,5000,6000)) +
    geom_point(data=table_CC_mean, aes(y=1, x = length), shape = 21, stroke =2, size = 4, color = "#880e4fcc") +
    geom_label(data=table_CC_mean, aes(y=1.025, x = length, label=length_label), size = 4, color = "#880e4fcc") +
    geom_point(data=table_CC_max, aes(y=1, x = length), shape = 21, stroke =2, size = 4, color = "#880e4fcc") +
    geom_label(data=table_CC_max, aes(y=1.025, x = length, label=length_label), size = 4, color = "#880e4fcc") +
    theme_unique_art() + labs(y="Proteins/CCs",x="Proteins #") + 
    theme(legend.title = element_text(face="bold"),
          strip.text.x = element_blank(),
          axis.text.y = element_blank())
  
  d <- ggplot(table_CC, aes(x = length)) +
    geom_histogram(binwidth = 50,color = "black", fill = "#880e4fcc",linewidth=0.5, center=25) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 2),breaks = c(0,1,2,4,8,16,32,64,256,512,2048,16368)) +
    scale_x_continuous(breaks = c(3,200,400,800,1600,3200,5000,6000)) +
    labs(y="CCs #",x="") + theme_unique_art() +
    theme(legend.title = element_text(face="bold"),
          axis.title = element_text(color = "black", face = "bold"),
          axis.title.x = element_blank(),
          strip.text.x = element_blank(),
          axis.text.x = element_blank())

  svglite("Analysis/Number_Prot_by_CC.svg",width = 12.00,height = 6.00)
  print(plot_grid(d,c, ncol = 1, nrow = 2,rel_heights = c(3,2), align ="v"))
  dev.off()
#
  # Trophic mode by CC -------------------------------------------------------------
  table_Grp <- stat_input_df %>% select(CC,TRT_1,TRT_2,TYP_1)
  # rename column
  colnames(table_Grp)[colnames(table_Grp)=="TRT_1"] <- "Grp_1"
  colnames(table_Grp)[colnames(table_Grp)=="TYP_1"] <- "Grp_2"
  # calcul trophic mode diversity through CC
  table_Grp <- table_Grp %>% select(CC,Grp_1) %>% group_by(CC) %>% summarize(Grp_1=paste(unique(Grp_1),collapse=","))
  table_Grp$Grp_1 <- sapply(str_replace_all(as.character(table_Grp$Grp_1),",,",","), function(x) x)
  table_Grp$Grp_1 <- sapply(str_remove(as.character(table_Grp$Grp_1),"^,"), function(x) x)
  table_Grp$Grp_1 <- sapply(str_remove(as.character(table_Grp$Grp_1),",$"), function(x) x)
  table_Grp$Grp_1 <- sapply(as.list(strsplit(as.character(table_Grp$Grp_1), ",")), function(x) x)
  table_Grp <- table_Grp %>% mutate(count=1)
  # Plot
  svglite("Analysis/Grp_Intersec_in_CC_UPSET_data.svg",width = 11.00,height = 4.0)
  ggplot(table_Grp,aes(x=Grp_1)) +
    geom_bar(color = "black", fill = "#880e4fcc",linewidth=0.5) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 2),breaks = c(0,2,8,32,128,512,2048,16368,131072)) +
    scale_x_upset() +
    #geom_text(stat = "count", aes(label = after_stat(count)), vjust = -1) +
    labs(y="CCs #",x="Trophic mode") + theme_unique_art() +
    theme(legend.title = element_text(face="bold"),
          axis.title = element_text(color = "black", face = "bold"),
          strip.text.x = element_text(color = "black", face = "bold", size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()
  # with data number
  svglite("Analysis/Grp_Intersec_in_CC_UPSET_data_n.svg",width = 22.00,height = 8.0)
  ggplot(table_Grp,aes(x=Grp_1)) +
    geom_bar(color = "black", fill = "#880e4fcc",linewidth=0.5) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 2),breaks = c(0,2,8,32,128,512,2048,16368,131072)) +
    scale_x_upset() +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -1) +
    labs(y="CCs #",x="Trophic mode") + theme_unique_art() +
    theme(legend.title = element_text(face="bold"),
          axis.title = element_text(color = "black", face = "bold"),
          strip.text.x = element_text(color = "black", face = "bold", size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()
#
# Network statistics for each trophic modes -----------------------------------------------------------------
  # Create CCs subset for each trophic modes ------------------------------------------------------------
  for (GRP in c("PARA","SAP","HET","PHOTO","MIXO")) {
    table <- paste("interest",GRP,"id",sep ="_")
    assign(table,stat_input_df %>% filter(TRT_1==GRP) %>% select(CC) %>% distinct())
    TRT_Stat_data <- stat_input_df %>% select(CC,TRT_1) %>% filter(CC %in% get(table)$CC)
    TRT_Stat_data <- TRT_Stat_data %>% group_by(CC, TRT_1) %>% 
    summarise(total_count=n(),
              .groups = 'drop')
  assign(paste0('interest_CC_75percent_major',GRP),vector())
  assign(paste0('interest_CC_50percent_major',GRP),vector())
  assign(paste0('interest_CC_100percent_major',GRP),vector())
  assign(paste0('interest_CC_100percent_only',GRP),vector())
  for (i in unique(TRT_Stat_data$CC)) {
    TRT_df <- TRT_Stat_data %>% filter(CC == i) %>% select(CC,TRT_1,total_count) %>% filter(TRT_1 != "Unassigned")
    unknown_i <- as.numeric(paste(TRT_Stat_data %>% filter(CC == i) %>% select(CC,TRT_1,total_count) %>% filter(TRT_1 == "Unassigned") %>% select(total_count)))
    GRP_i <- as.numeric(paste(TRT_df %>% filter(TRT_1 == GRP) %>% select(total_count)))
    other_i <- as.numeric(paste(TRT_df %>% filter(TRT_1 != GRP) %>% select(-TRT_1) %>% group_by(CC) %>% summarise_all(sum) %>% select(total_count)))
    if (is.na(unknown_i)==TRUE) {unknown_i <- 0}
    if (is.na(GRP_i)==TRUE) {GRP_i <- 0}
    if (is.na(other_i)==TRUE) {other_i <- 0}
    if (GRP_i/3 >= other_i) { assign(paste0('interest_CC_75percent_major',GRP),c(get(paste0('interest_CC_75percent_major',GRP)),i)) } # at least 75 % of the associated proteins are PARA
    if (GRP_i >= other_i) { assign(paste0('interest_CC_50percent_major',GRP),c(get(paste0('interest_CC_50percent_major',GRP)),i)) } # at least 50 % of the associated proteins are PARA
    if (GRP_i > 0 & other_i == 0) { assign(paste0('interest_CC_100percent_major',GRP),c(get(paste0('interest_CC_100percent_major',GRP)),i)) } # 100% of the associated proteins are PARA
    if (GRP_i > 0 & other_i == 0 & unknown_i == 0) { assign(paste0('interest_CC_100percent_only',GRP),c(get(paste0('interest_CC_100percent_only',GRP)),i)) } # 100% of proteins are PARA
  }
#
  assign(paste("data",GRP,"stat_i",sep="_"),tibble())
  for (i in list(paste0('interest_CC_75percent_major',GRP),
                 paste0('interest_CC_50percent_major',GRP),
                 paste0('interest_CC_100percent_major',GRP),
                 paste0('interest_CC_100percent_only',GRP))) {
  CC_withGRP <- stat_input_df %>% filter(CC %in% get(i))
  CC_withGRP$sum <- 1
  assign(paste("data",GRP,"stat_i",sep="_"),rbind(get(paste("data",GRP,"stat_i",sep="_")),CC_withGRP %>% select(CC,sum) %>% group_by(CC) %>% summarise_all(sum) %>% mutate(Category=paste(i))))
  }
  assign("input_i",get(paste("data",GRP,"stat_i",sep="_")))
  input_i[input_i==paste0('interest_CC_100percent_only',GRP)] <- "100%*"
  input_i[input_i==paste0('interest_CC_100percent_major',GRP)] <- "100%"
  input_i[input_i==paste0('interest_CC_75percent_major',GRP)] <- "≥75%"
  input_i[input_i==paste0('interest_CC_50percent_major',GRP)] <- "≥50%"
  
# Plot
  svglite(paste0("Analysis/CC_containing",GRP,"_Stat.svg"),width = 5.00,height = 4.50)
  input_i$Category <- factor(input_i$Category,levels=c("≥50%","≥75%","100%","100%*"))
  plt <- ggbetweenstats(
    data = input_i,
    x = Category,
    y = sum,
    results.subtitle =FALSE,
    bf.message = FALSE,
    centrality.label.args = list(size  = 4,nudge_x=0.1,nudge_y=-1.5)
    ) + scale_y_continuous(trans = scales::pseudo_log_trans(base = 2),breaks = c(0,2,4,8,16,32,64,128,256,512,1024),limits =c(0,1024)) +
    #stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    scale_color_manual(values = c("#ea716e","red2","red3","red4")) + theme_unique_art() + 
    labs(x=paste(GRP,"proportion across CC"),y="CC Abundance") + guides(color = "none")
  print(plt)
  dev.off()
}
#
  # Functional and taxonomic diversity through trophic strategies -----------------
stat_specific <- data.frame()
ko_rate_table <- tibble()
stat_taxo <- data.frame()
CC_count_tax  <- data.frame()
stat_ko <- data.frame()
CC_count_ko  <- data.frame()
protein_count <- data.frame()
for (Tmode in c("PARA","SAP","HET","PHOTO","MIXO")) {
  assign("input_75",get(paste("data",Tmode,"stat_i",sep="_")))
  # Select 75% threshold
  input_75 <- input_75 %>% filter(Category == paste0("interest_CC_75percent_major",Tmode))
  subset_cat <- (input_75 %>% select(CC))$CC
  # Format table
  stat_input_df_plus <- stat_input_df %>% select(-TRT_2,-TYP_1,-TYP_2,-ZONE_PA,-FRACTION_PA,-PERIODS_PA,-Short_name_Genes,-GOterms,-PfamDom)
  stat_compo <- stat_input_df_plus %>% filter(CC %in% subset_cat)
  stat_compo[,"ko"][grepl(stat_compo[,"ko"],pattern=",") == TRUE] <- "Unassigned"
  stat_compo[,"ko"][stat_compo[,"ko"] == ""] <- "Unassigned"

#
    # Taxonomy ----------------------------------------------------------------
  stat_compo_tax <- left_join(stat_compo,summarize_data,"Taxonomy")
  stat_compo_tax$count <- 1
  CC_count_tax <- rbind(CC_count_tax,data.frame(CC_count_tax=length(unique(stat_compo_tax$CC)),GRP=Tmode))
  stat_compo_tax$Phylum[is.na(stat_compo_tax[,"Phylum"])==TRUE & is.na(stat_compo_tax[,"Class"])==FALSE] <- stat_compo_tax$Class[is.na(stat_compo_tax[,"Phylum"])==TRUE & is.na(stat_compo_tax[,"Class"])==FALSE]
  stat_compo_tax$Class[is.na(stat_compo_tax[,"Phylum"])==FALSE & is.na(stat_compo_tax[,"Class"])==TRUE & is.na(stat_compo_tax[,"Division"])==FALSE] <- paste0("Unaffiliated_",stat_compo_tax$Phylum[is.na(stat_compo_tax[,"Phylum"])==FALSE & is.na(stat_compo_tax[,"Class"])==TRUE & is.na(stat_compo_tax[,"Division"])==FALSE])
  stat_compo_tax$Phylum[is.na(stat_compo_tax[,"Phylum"])==TRUE & is.na(stat_compo_tax[,"Class"])==TRUE & is.na(stat_compo_tax[,"Division"])==FALSE] <- paste0("Unaffiliated_",stat_compo_tax$Division[is.na(stat_compo_tax[,"Phylum"])==TRUE & is.na(stat_compo_tax[,"Class"])==TRUE & is.na(stat_compo_tax[,"Division"])==FALSE])
  stat_compo_tax$Class[is.na(stat_compo_tax[,"Class"])==TRUE & is.na(stat_compo_tax[,"Phylum"])==FALSE & is.na(stat_compo_tax[,"Division"])==FALSE] <- stat_compo_tax$Phylum[is.na(stat_compo_tax[,"Phylum"])==FALSE & is.na(stat_compo_tax[,"Class"])==TRUE & is.na(stat_compo_tax[,"Division"])==FALSE]
  stat_compo_tax[is.na(stat_compo_tax)==TRUE] <- "Unaffiliated"
  stat_compo_tax <- stat_compo_tax %>% select(Division,Phylum,Class,count) %>% group_by(Division,Phylum,Class,) %>% summarise_all(sum)
  stat_compo_tax$GRP <- Tmode
  stat_compo_tax <- stat_compo_tax %>% mutate(percent = count*100/sum(stat_compo_tax$count))
  stat_compo_tax <- as.data.frame(stat_compo_tax)
  tempx <- c(t(stat_compo_tax %>% filter(Division != "Unaffiliated") %>% select(-Division,-Phylum,-Class) %>% group_by(GRP) %>% summarise_all(sum) %>% select(count,GRP,percent)))
  stat_compo_tax <- rbind(stat_compo_tax,c("Affiliated","Affiliated","Affiliated",tempx))
  stat_compo_tax$percent <- as.numeric(stat_compo_tax$percent)
  stat_compo_tax$count <- as.numeric(stat_compo_tax$count)
  stat_compo_tax <- stat_compo_tax %>% mutate(label = paste0(Division,": ",round(percent,1),"%"))
  stat_taxo <- rbind(stat_taxo,stat_compo_tax)
#
    # Function (KO) -----------------------------------------------------------
  ko_rate_table_i <- tibble(GRP=Tmode,Nprot=nrow(stat_compo),NAnnot=nrow(stat_compo %>% filter(ko != "Unassigned")))
  ko_rate_table <- rbind(ko_rate_table,ko_rate_table_i)
  stat_compo_ko <- left_join(stat_compo,ko_db %>% select(KO_id,lvl_B_val,lvl_C_val) %>% distinct(),by = c("ko"="KO_id"),relationship = "many-to-many")
  stat_compo_ko[,"lvl_B_val"][is.na(stat_compo_ko[,"lvl_B_val"]) == TRUE] <- "Unassigned"
  stat_compo_ko[,"lvl_C_val"][is.na(stat_compo_ko[,"lvl_C_val"]) == TRUE] <- "Unassigned"
  stat_compo_ko$count <- 1
  protein_count <- rbind(protein_count,data.frame(protein_count=length(unique(stat_compo_ko$name)),GRP=Tmode))
  CC_count_ko <- rbind(CC_count_ko,data.frame(CC_count_ko=length(unique(stat_compo_ko$CC)),GRP=Tmode))
  stat_compo_ko <- stat_compo_ko %>% select(lvl_B_val,lvl_C_val,count) %>% group_by(lvl_B_val,lvl_C_val) %>% summarise_all(sum)
  stat_compo_ko$GRP <- Tmode
  stat_compo_ko <- stat_compo_ko %>% mutate(percent = count*100/sum(stat_compo_ko$count))
  stat_compo_ko <- as.data.frame(stat_compo_ko)
  stat_compo_ko$lvl_B_val <- sapply(str_replace_all(as.character(stat_compo_ko$lvl_B_val),"_"," "), function(x) x)
  stat_compo_ko$lvl_C_val <- sapply(str_replace_all(as.character(stat_compo_ko$lvl_C_val),"_"," "), function(x) x)
  stat_compo_ko$percent <- as.numeric(stat_compo_ko$percent)
  stat_compo_ko$count <- as.numeric(stat_compo_ko$count)
  tempx <- c(t(stat_compo_ko %>% filter(lvl_B_val != "Unassigned") %>% select(-lvl_B_val,-lvl_C_val) %>% group_by(GRP) %>% summarise_all(sum) %>% select(count,GRP,percent)))
  stat_compo_ko <- rbind(stat_compo_ko,c("Assigned","Assigned",tempx))
  stat_compo_ko$percent <- as.numeric(stat_compo_ko$percent)
  stat_compo_ko$count <- as.numeric(stat_compo_ko$count)
  stat_compo_ko <- stat_compo_ko %>% mutate(label = paste0(lvl_C_val,": ",round(percent,1),"%"))
  stat_ko <- rbind(stat_ko,stat_compo_ko)
}

    # Plot Taxonomy -----------------------------------------------------------
  val_tab <- as.data.frame(stat_taxo %>% filter(Division %in% c("Unaffiliated","Affiliated")) %>% select(count,GRP) %>% group_by(GRP) %>% summarise_all(sum))
  val_tab <- merge(val_tab,as.data.frame(stat_taxo %>% filter(Division %in% c("Affiliated")) %>% select(GRP,label)),by="GRP")
  val_tab <- merge(val_tab,CC_count_tax,by="GRP")
  # Pre-plot
  u <- ggplot() +
    geom_treemap(data = stat_taxo %>% filter(Division %in% c("Unaffiliated","Affiliated")), aes(area = percent, fill = Division,subgroup = Division, label = label),color="white",size=2) +
    scale_fill_manual(values=c("black","white")) +
    facet_grid(factor(GRP,levels = c("PHOTO","HET","MIXO","SAP","PARA"), labels = c("Phototrophs","Heterotrophs","Mixotrophs","Saprotrophs","Parasites"))~., switch="both") +
    geom_treemap_subgroup_border(size=2,color="black") + theme_unique_darktris() +
    theme(panel.border = element_rect(fill=NA,linewidth=2),legend.position="none") +
    geom_label(data = val_tab, mapping=aes(label = paste0(count," proteins","\n",label,"\n",CC_count_tax," CCs")),fontface = "bold", x = 0.5 ,y = 0.5,color = "white",fill="black",size=4,alpha=0.9)
  print(u)
  # Color
  stat_taxo_i <- stat_taxo %>% filter(!Division %in% c("Unaffiliated","Affiliated"))
  color <- data.frame(color_definition,c("#c62728cc","#6b1b9acc","white","#273493cc","white","#0277bdcc","#558b30cc","#f9a825cc","white","#fdffcccc","#ef6c00cc","white","white","white","white","#75252595","#880e4fcc","#90ee90cc","#20606490","#add8e6cc","#e53368cc"))
  colnames(color) <- c("Division","color")
  stat_taxo_i <- merge(stat_taxo_i,color,by="Division")
  # Legend
  order <- stat_taxo_i %>% select(Division,percent) %>% group_by(Division) %>% summarise_all(sum) %>% arrange(-percent)
  color <- left_join(order,color,by="Division")
  color$Division <- factor(color$Division,levels=order)
  legend_a <- get_legend(ggplot(data = stat_taxo_i %>% filter(percent > 0.02) %>% arrange(GRP,percent),mapping = aes(y= percent, x = GRP,fill = factor(Division,levels= order$Division),group=factor(Division,levels= order$Division)),color="black", Rowv = NA, scale = "column",alpha=0.8) + 
                           theme_unique_darkcris() + facet_grid(.~GRP,scale="free") + labs(fill ="Taxonomy") +
                           geom_bar(stat="identity") + guides(fill = guide_legend(ncol=6,override.aes = list(color="black",shape=1))) + 
                           scale_fill_manual(values=noquote(color$color)))
  # Final plot
  v <- ggplot(stat_taxo_i, aes(area = percent, fill = Division,subgroup = Division, label = Class)) + 
    geom_treemap(alpha=0.8,color="white",size=2, fill=stat_taxo_i[,"color"]) + 
    #scale_fill_paletteer_d("MetBrewer::Signac",direction = -1) +
    facet_grid(factor(GRP,levels = c("PHOTO","HET","MIXO","SAP","PARA"), labels = c("Phototrophs","Heterotrophs","Mixotrophs","Saprotrophs","Parasites"))~., switch="both") + 
    guides(fill=guide_legend(ncol=6,title.position = "top")) +
    geom_treemap_subgroup_border(size=2,color="black") + theme_unique_darktris() +
    geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=5) +
    theme(panel.border = element_rect(fill=NA,linewidth=2),strip.text.y.left = element_blank()) + theme(legend.position="none")
  #
  svglite(paste0("Analysis/taxo_compo_75.svg"),width = 14.00,height = 8.00)
  v_plot <- plot_grid(u,v,NA,legend_a, ncol = 2, nrow = 2, align = "h",rel_widths = c(1,6),rel_heights = c(7,1))
  print(v_plot)
  dev.off()
  # Save table
  write.table(stat_taxo_i %>% arrange(-percent,Division), file = paste0("Analysis/taxo_compo_75.csv"),row.names = FALSE,col.names = TRUE,quote = FALSE,sep= ";")
#
    # Plot Function ----------------------------------------------------------------------
  val_tab <- as.data.frame(stat_ko %>% filter(lvl_B_val %in% c("Unassigned","Assigned")) %>% select(count,GRP) %>% group_by(GRP) %>% summarise_all(sum))
  val_tab <- merge(val_tab,as.data.frame(stat_ko %>% filter(lvl_B_val %in% c("Assigned")) %>% select(GRP,label)),by="GRP")
  val_tab <- merge(val_tab,CC_count_ko,by="GRP")
  val_tab <- merge(val_tab,protein_count,by="GRP")
  # Pre-plot
  t <- ggplot() + 
    geom_treemap(data = stat_ko %>% filter(lvl_B_val %in% c("Unassigned","Assigned")), aes(area = percent, fill = lvl_B_val,subgroup = lvl_B_val, label = label),color="white",size=2) +
    scale_fill_manual(values=c("black","white")) +
    facet_grid(factor(GRP,levels = c("PHOTO","HET","MIXO","SAP","PARA"), labels = c("Phototrophs","Heterotrophs","Mixotrophs","Saprotrophs","Parasites"))~., switch="both") + 
    geom_treemap_subgroup_border(size=2,color="black") + theme_unique_darktris() +
    geom_text(data = val_tab, mapping=aes(label = paste0(protein_count," proteins","\n",label,"\n",CC_count_ko," CCs")),fontface = "bold", x = 0,hjust =0 ,y = 0.5,color="black",size=4) +
    theme(panel.border = element_rect(fill=NA,linewidth=2),legend.position="none")
  print(t)
  # Color
  stat_ko_i <- stat_ko %>% filter(!lvl_B_val %in% c("Unassigned","Assigned"))
  color_ko <- data.frame(unique(stat_ko$lvl_B_val),c(paletteer_d("MetBrewer::Signac",direction = -1,n=length(unique(stat_ko$lvl_B_val)))))
  colnames(color_ko) <- c("lvl_B_val","color")
  stat_ko_i <- merge(stat_ko_i,color_ko,by="lvl_B_val")
  # Legend
  order <- stat_ko_i %>% select(lvl_B_val,percent) %>% group_by(lvl_B_val) %>% summarise_all(sum) %>% arrange(-percent)
  color_ko <- left_join(order,color_ko,by="lvl_B_val")
  color_ko$lvl_B_val <- factor(color_ko$lvl_B_val,levels=order)
  legend_b <- get_legend(ggplot(data = stat_ko_i %>% filter(percent > 0.02) %>% arrange(GRP,percent),mapping = aes(y= percent, x = GRP,fill = factor(lvl_B_val,levels= order$lvl_B_val) ,group=factor(lvl_B_val,levels= order$lvl_B_val)),color="black", Rowv = NA, scale = "column",alpha=0.8) + 
                           theme_unique_darkcris() + facet_grid(.~GRP,scale="free") + labs(fill ="KEGG category") +
                           geom_bar(stat="identity") + guides(fill = guide_legend(ncol=4,override.aes = list(color="black",shape=1))) + 
                           scale_fill_manual(values=noquote(color_ko$color)))
  # Final plot
  s <- ggplot(stat_ko_i, aes(area = percent, fill = lvl_B_val,subgroup = lvl_B_val, label = lvl_C_val)) + 
    geom_treemap(alpha=1,color="white",size=2, fill=stat_ko_i[,"color"]) + 
    facet_grid(factor(GRP,levels = c("PHOTO","HET","MIXO","SAP","PARA"), labels = c("Phototrophs","Heterotrophs","Mixotrophs","Saprotrophs","Parasites"))~., switch="both") + 
    geom_treemap_subgroup_border(size=2,color="black") + theme_unique_darktris() +
    geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=5) +
    theme(panel.border = element_rect(fill=NA,linewidth=2),strip.text.y.left = element_blank()) + theme(legend.position="none")
  #
  svglite(paste0("Analysis/KO_compo_75.svg"),width = 15.00,height = 12.00)
  b_plot <- plot_grid(t,s,NA,legend_b, ncol = 2, nrow = 2, align = "h",rel_widths = c(1.2,6),rel_heights = c(6,1))
  print(b_plot)
  dev.off()
  # Save table
  write.table(stat_ko_i %>% arrange(-percent,lvl_B_val), file = paste0("Analysis/KO_compo_75.csv"),row.names = FALSE,col.names = TRUE,quote = FALSE,sep= ";")
#
    # Specific metabolisms ------------------------------------------------------------
  if (unique(ko_db$lvl_A_val)=="Metabolism") {
    stat_specific <- stat_ko_i %>% filter(lvl_C_val %in% c("Fatty acid degradation",
                                                                                         "Biosynthesis of unsaturated fatty acids",
                                                                                         "Glycerophospholipid metabolism",
                                                                                         "Glycerolipid metabolism",
                                                                                         "Sphingolipid metabolism",
                                                                                         "Sulfur metabolism",
                                                                                         "Photosynthesis - antenna proteins",
                                                                                         "Prodigiosin biosynthesis",
                                                                                         "Carotenoid biosynthesis")) %>% select(lvl_C_val,count,percent,GRP)
  }
  #
  # Hist
  svglite(paste0("Analysis/KO_compo_special_category.svg"),width = 10.50,height = 6.00)
  ggplot(stat_specific, aes(y = GRP, x = percent, fill = lvl_C_val)) + 
    geom_bar(color="black",linewidth=0.5,stat="identity") + 
    facet_wrap(lvl_C_val~., ncol = 3) + scale_fill_manual(values=color$color) +
    theme_unique_art() + theme(legend.position="none") + labs(x="Proportion",y="Functional Groups")
  dev.off()
  #
  # Radar
  ndata <- reshape2::dcast(data = stat_specific,formula = GRP~lvl_C_val,fun.aggregate = sum,value.var = "percent")
  coul <- c("#c62728cc","#6b1b9acc","#0277bdcc","#558b30cc","#f9a825cc")
  svglite(paste0("Analysis/KO_compo_special_radar.svg"),width = 10.50,height = 12.00)
  row.names(ndata) <- ndata$GRP ; ndata <- ndata %>% select(-GRP)
  colors_border <- coul
  library(scales)
  colors_in <- alpha(coul,0.2)
  radarchart( ndata  , axistype=1 , maxmin=F,
              #custom polygon
              pcol=colors_border , pfcol=colors_in , plwd=2 , 
              caxislabels=seq(0,1.2,0.3),
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
              #custom labels
              vlcex=0.8 
  )
  legend(x=-1.2, y=1.2, legend = rownames(ndata[,]), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=0.8, pt.cex=3, ncol = 3)
  dev.off()
  # ko_rate_table
  ko_rate_table <- ko_rate_table %>% mutate(rate=NAnnot*100/Nprot)
  write.table(ko_rate_table, file = paste0("Analysis/ko_rate_table.csv"),row.names = FALSE,col.names = TRUE,quote = FALSE,sep= ";")
#
  # Save R Image for tree stat ------------------------------------------------------------
  save.image("R_Tree_stat.RData")
#
# Statistic test -----------------------------------------------------------
  stat_input_df_pcoa <- stat_input_df %>% select(CC,TRT_1)
  table_pcoa <- as.data.frame.matrix(table(stat_input_df_pcoa))
  x <- table_pcoa %>% select(-Unassigned)
  # Transform table
  x <- as.data.frame(x)
  x <- reshape2::melt(x)
  colnames(x)[colnames(x)=="variable"] <- "group"
  # Bartlett
    bartlett.test(data = x, value ~ group) # p-value < 0.05 : Bartlett test is not validate => necessity to make non-parametric test
  # Wilcox 
  wilcox_x <- x %>% wilcox_test(value ~ group, paired = TRUE, p.adjust.method = "bonferroni")
  wilcox_x
  # graph
    wilcox_x_data <- wilcox_x %>% add_xy_position(x = "group", fun = "max")
    ggboxplot(x, x = "group", y = "value") +
      stat_pvalue_manual(wilcox_x_data, hide.ns = TRUE)
#
# Multidimensional analyse GRP ----------------------------------------------------------------
  HET <- (data_HET_stat_i %>% filter(Category == "interest_CC_75percent_majorHET"))$CC
  MIXO <- (data_MIXO_stat_i %>% filter(Category == "interest_CC_75percent_majorMIXO"))$CC
  PHOTO <- (data_PHOTO_stat_i %>% filter(Category == "interest_CC_75percent_majorPHOTO"))$CC
  PARA <- (data_PARA_stat_i %>% filter(Category == "interest_CC_75percent_majorPARA"))$CC
  SAP <- (data_SAP_stat_i %>% filter(Category == "interest_CC_75percent_majorSAP"))$CC
  CC_int <- unique(c(HET,MIXO,PHOTO,PARA,SAP))
#
    # CA ggplot2 TRT_1 ----------------------------------------------------------------
  stat_input_df_pcoa <- stat_input_df %>% select(CC,TRT_1)
  table_pcoa <- as.data.frame.matrix(table(stat_input_df_pcoa))
  #
  # Distance
  res.ca <- CA(table_pcoa,graph = FALSE, ncp = 3)
  #
  #row
  p <- get_ca_row(res.ca)
  coord.y <- p$coord
  coord.y <- as.data.frame(coord.y)
  coord.y$CC <- row.names(coord.y)
  coord.y$TYPE[coord.y[,"CC"] %in% PARA] <- "Parasites"
  coord.y$TYPE[coord.y[,"CC"] %in% MIXO] <- "Mixotrophs"
  coord.y$TYPE[coord.y[,"CC"] %in% HET] <- "Heterotrophs"
  coord.y$TYPE[coord.y[,"CC"] %in% PHOTO] <- "Phototrophs"
  coord.y$TYPE[coord.y[,"CC"] %in% SAP] <- "Saprotrophs"
  coord.y$TYPE[is.na(coord.y[,"TYPE"]) == TRUE] <- "Unassigned"
  #
  #col
  q <- get_ca_col(res.ca)
  coord.x <- q$coord
  coord.x <- as.data.frame(coord.x)
  coord.x$Functional_group <- row.names(coord.x)
  coord.x$Functional_group[coord.x[,"Functional_group"] == "PARA"] <- "Parasites"
  coord.x$Functional_group[coord.x[,"Functional_group"] == "MIXO"] <- "Mixotrophs"
  coord.x$Functional_group[coord.x[,"Functional_group"] == "HET"] <- "Heterotrophs"
  coord.x$Functional_group[coord.x[,"Functional_group"] == "PHOTO"] <- "Phototrophs"
  coord.x$Functional_group[coord.x[,"Functional_group"] == "SAP"] <- "Saprotrophs"
  Xseq <- fviz_screeplot(res.ca,addlabels = TRUE, ylim = c(0,50), barfill = "maroon", barcolor = "black") + theme_unique_art()
  svglite(paste0("Analysis/CA_TRT_eig.svg"),width = 6.00,height = 4.00)
  print(Xseq)
  dev.off()
  Xseq <- Xseq$data
  Dim1 <- paste0("Dim 1 [",round(Xseq %>% filter(dim == 1) %>% select(eig),1),"%]")
  Dim2 <- paste0("Dim 2 [",round(Xseq %>% filter(dim == 2) %>% select(eig),1),"%]")
  Dim3 <- paste0("Dim 3 [",round(Xseq %>% filter(dim == 3) %>% select(eig),1),"%]")
  #
  # palet
  palet_tree_grp <- c("#92C051FF","#2B9B81FF","#1F6E9CFF","#633372FF","#E87B89FF","#D8443CFF")
  #
  #Dim1 x Dim 2
  df <- coord.y %>% select(-CC,-`Dim 3`)
  find_hull <- function(df) df[chull(df$`Dim 1`, df$`Dim 2`), ]
  hulls <- ddply(df, "TYPE", find_hull)
  #
  LEGEND_a <- get_legend(pcoa_plot_c <- ggplot() + theme_unique_art() + labs(title="",color="Trophic mode",fill = "CCs", x= Dim1, y=Dim2) + 
                           theme(legend.text=element_text(size=12),legend.position = "right") + scale_fill_manual(values = palet_tree_grp) +
                           geom_point(data=coord.y,aes(x=`Dim 1`,y=`Dim 2`,fill = TYPE),  shape = 21, color =c("black")) + guides(fill = guide_legend(override.aes = list(size = 3))))
  pcoa_plot_c <- ggplot() + theme_unique_art() + 
    theme(legend.text=element_text(size=12),legend.position = "right") + 
    geom_point(data=coord.y,aes(x=`Dim 1`,y=`Dim 2`,fill = TYPE),  shape = 21, color =c("black")) + 
    geom_point(data=coord.x,aes(x=`Dim 1`,y=`Dim 2`,color = Functional_group), size = 3, shape = 21, stroke =2, fill ='white') + 
    scale_color_manual(values = palet_tree_grp) +
    scale_fill_manual(values = palet_tree_grp) +
    geom_polygon(data = hulls %>% filter(TYPE != "Unassigned"), aes(x=`Dim 1`,y=`Dim 2`,color=TYPE), alpha = 0.05, linetype = 4) +
    geom_label_repel(data=coord.x,aes(x=`Dim 1`,y=`Dim 2`,label = Functional_group),fontface = "bold",color = "black",fill = "white",size = 4,show.legend = FALSE, nudge_x = -1, nudge_y = 0.5) +
    labs(title="",color="Trophic mode",fill = "CCs", x= Dim1, y=Dim2) + 
    scale_y_continuous(position = 'left')
  #LEGEND_a <- get_legend(pcoa_plot_c)
  pcoa_plot_c <- pcoa_plot_c + theme(legend.position = "none")
  #Dim3 x Dim 2
  df <- coord.y %>% select(-CC,-`Dim 1`)
  find_hull <- function(df) df[chull(df$`Dim 3`, df$`Dim 2`), ]
  hulls <- ddply(df, "TYPE", find_hull)
  #
  pcoa_plot_d <- ggplot() + theme_unique_art() + 
    theme(legend.text=element_text(size=12),legend.position = "right") + 
    geom_point(data=coord.y,aes(x=`Dim 3`,y=`Dim 2`,fill = TYPE), size = 2, shape = 21, color =c("black")) + 
    geom_point(data=coord.x,aes(x=`Dim 3`,y=`Dim 2`,color = Functional_group), size = 3, shape = 21, stroke =2, fill ='white') + 
    scale_color_manual(values = palet_tree_grp) + 
    scale_fill_manual(values = palet_tree_grp) +
    geom_polygon(data = hulls %>% filter(TYPE != "Unassigned"), aes(x=`Dim 3`,y=`Dim 2`,color=TYPE), alpha = 0.05, linetype = 4) +
    geom_label_repel(data=coord.x,aes(x=`Dim 3`,y=`Dim 2`,label = Functional_group),fontface = "bold",color = "black",fill = "white",size = 4,show.legend = FALSE, nudge_x = -1, nudge_y = -0.5) +
    labs(title="",color="Trophic mode", x= Dim3, y=Dim2) + 
    scale_y_continuous(position = 'right') + theme(legend.position = "none")
  
  b_plot <- plot_grid(pcoa_plot_c, pcoa_plot_d, LEGEND_a,ncol = 3, nrow = 1, rel_widths = c(3,3,1),rel_heights = c(3))
  svglite(paste0("Analysis/CA_TRT.svg"),width = 14.00,height = 8.00)
  plot(b_plot)
  dev.off()
    # CA barycentre KO -------------------------------------------------------
  #Filter obligate
  Obligate_PARA_CC <- unique(stat_input_df[,"CC"][stat_input_df[,"TRT_1"]=="PARA" & stat_input_df[,"TYP_1"]=="Obligate"])
  interest_CC_majorPARA <- PARA[PARA %in% Obligate_PARA_CC]
  CC_bary <- stat_input_df
  CC_bary$sum <- 1
  summary(CC_bary %>% select(CC,sum) %>% group_by(CC) %>% summarise_all(sum) %>% select(sum))
  length(unique(CC_bary$CC))
  #
  #ko
  df_bary <- CC_bary %>% select(CC, ko) %>% distinct(CC, ko)
  df_bary <- separate_rows(df_bary,"ko",sep=",") %>% mutate(sum = 1) %>% group_by(ko,sum) %>% summarize(CC=paste(unique(CC),collapse=","),sum(sum))
  df_bary_subset <- df_bary %>% mutate(CDT = length(intersect(unlist(str_split(CC,pattern=",")),interest_CC_majorPARA))) %>% filter(CDT > 0) %>% filter(ko != "Unassigned")
  #
  #Barycentre
  for (i in row.names(df_bary_subset)) { 
    Dim1j <- vector()
    Dim2j <- vector()
    Dim3j <- vector()
    print(df_bary_subset[i,"CC"])
    for (j in df_bary_subset[i,"CC"] %>% str_split(pattern = ",") %>% unlist()) {
      Dim1j <- c(Dim1j,coord.y[,"Dim 1"][coord.y[,"CC"]==j])
      Dim2j <- c(Dim2j,coord.y[,"Dim 2"][coord.y[,"CC"]==j])
      Dim3j <- c(Dim3j,coord.y[,"Dim 3"][coord.y[,"CC"]==j])
    }
    Dim1j <- mean(Dim1j)
    df_bary_subset[i,"Dim1j"] <- Dim1j
    Dim2j <- mean(Dim2j)
    df_bary_subset[i,"Dim2j"] <- Dim2j
    Dim3j <- mean(Dim3j)
    df_bary_subset[i,"Dim3j"] <- Dim3j
  }
  df_bary_subsetx <- df_bary_subset %>% filter(Dim1j >= 2)
  #Dim1 x Dim 2
  df <- coord.y %>% select(-CC,-`Dim 3`)
  find_hull <- function(df) df[chull(df$`Dim 1`, df$`Dim 2`), ]
  hulls <- ddply(df, "TYPE", find_hull)
  #
  df_bary_subset_legend <- df_bary_subset %>% mutate(uptwo = ifelse(ko %in% df_bary_subsetx$ko,"≥ 2","< 2"))
  LEGEND_a <- get_legend(pcoa_plot_c <- ggplot() + theme_unique_art() + labs(title="",color="Trophic mode",fill = "CCs",alpha = "Barycenter coordinates\non the dimension 1", x= Dim1, y=Dim2) + 
                           geom_point(data=df_bary_subset_legend,aes(x=`Dim1j`,y=`Dim2j`,alpha = uptwo), size = 2, shape = 8, stroke =1,color ='black') +
                           theme(legend.text=element_text(size=12),legend.position = "right") + scale_fill_manual(values = palet_tree_grp) + scale_alpha_manual(values = c(0.5,1)) +
                           geom_point(data=coord.y,aes(x=`Dim 1`,y=`Dim 2`,fill = TYPE),  shape = 21, color =c("black")) + guides(fill = guide_legend(override.aes = list(size = 3)),
                                                                                                                                  alpha = guide_legend(override.aes = list(size = 3))))
  #
  pcoa_plot_c <- ggplot() + theme_unique_art() + 
    theme(legend.text=element_text(size=12),legend.position = "right") + 
    geom_point(data=coord.y,aes(x=`Dim 1`,y=`Dim 2`,fill = TYPE), size = 2, shape = 21, color =c("black")) + 
    geom_point(data=df_bary_subset %>% filter(Dim1j < 2),aes(x=`Dim1j`,y=`Dim2j`), size = 2, shape = 8, stroke =1,color ='black',alpha=0.5) +
    geom_point(data=df_bary_subset %>% filter(Dim1j >= 2),aes(x=`Dim1j`,y=`Dim2j`), size = 2, shape = 8, stroke =1,color ='black') +
    geom_point(data=coord.x,aes(x=`Dim 1`,y=`Dim 2`,color = Functional_group), size = 3, shape = 21, stroke =2, fill ='white') +
    scale_color_manual(values = palet_tree_grp) +
    scale_fill_manual(values = palet_tree_grp) +
    geom_vline(xintercept = 2, linetype = 2, color = "#383e42") +
    geom_polygon(data = hulls %>% filter(TYPE != "Unassigned"), aes(x=`Dim 1`,y=`Dim 2`,color=TYPE), alpha = 0.05, linetype = 4) +
    geom_label_repel(data=df_bary_subsetx,aes(x=`Dim1j`,y=`Dim2j`,label = ko),size = 3,show.legend = FALSE, nudge_x = 0.5, nudge_y = 0.5,max.overlaps = 500,color="black",fill = "white") +
    geom_label_repel(data=coord.x,aes(x=`Dim 1`,y=`Dim 2`,label = Functional_group),color = "black",size = 4,fontface = "bold",show.legend = FALSE, nudge_x = -0.8, nudge_y = -0.5) +
    labs(title="",color="Trophic mode",fill = "CCs", x= Dim1, y=Dim2) +
    scale_y_continuous(position = 'left') + theme(legend.position = "none")
  #
  b_plot <- plot_grid(pcoa_plot_c, LEGEND_a,ncol = 2, nrow = 1, rel_widths = c(4,0.8),rel_heights = c(3))
  svglite(paste0("Analysis/CA_TRT_Barycentre_obligatePARA.svg"),width = 16.00,height = 10.00)
  plot(b_plot)
  dev.off()
#
  # Save R Image for barycentre ------------------------------------------------------------
  save.image("R_Barycentre_checkpoint.RData")
#
# Processes_network -------------------------------------------------------
  graph_input <- graph_from_data_frame(edges_input_df, directed = FALSE, vertices = attributes_input_df)
  graph_decomposed <- decompose(graph_input, mode = c("weak", "strong"), max.comps = NA,
                                      min.vertices = 3)
  #
# Display statistics graph of network ---------------------------------------------------------
  # Select PARA CCs which contains sequences labeled with obligate parasite --------------------------------------------------
  Obligate_PARA_CC <- unique(stat_input_df[,"CC"][stat_input_df[,"TRT_1"]=="PARA" & stat_input_df[,"TYP_1"]=="Obligate"])
  interest_CC_majorPARA <- PARA[PARA %in% Obligate_PARA_CC]
  CC_withPARA <- stat_input_df %>% filter(CC %in% interest_CC_majorPARA)
  CC_withPARA$sum <- 1
  summary(CC_withPARA %>% select(CC,sum) %>% group_by(CC) %>% summarise_all(sum) %>% select(sum))
  length(unique(CC_withPARA$CC))
  stat_CC_withPARA <- CC_withPARA %>% select(CC,sum) %>% group_by(CC) %>% summarise_all(sum) %>% filter(sum <= 226 & sum >= 3) %>% select(CC)
  nrow(stat_CC_withPARA)
  summarize_datax <- summarize_datax %>% filter(ID_Protein %in% CC_withPARA$name)
#
  # Add functional information ---------------------------------
  df_PARA <- CC_withPARA %>% select(CC, ko) %>% group_by(CC) %>% summarize(ko=paste(unique(ko),collapse=","))
  df_PARA$ko <- sapply(str_remove(as.character(df_PARA$ko),"^,"), function(x) x)
  df_PARA$ko <- sapply(str_remove(as.character(df_PARA$ko),",$"), function(x) x)
  df_PARA$ko <- sapply(str_replace_all(as.character(df_PARA$ko),",,",","), function(x) x)
  df_PARA$ko <- sapply(str_split(df_PARA$ko,","), function(x) paste(unique(x),collapse = ","))
  df_PARA <- as.data.frame(df_PARA)
  for (i in row.names(df_PARA)) { 
    if (grepl(",",df_PARA[i,"ko"])==TRUE){ df_PARA[i,"ko"] <- paste(unique(unlist(str_split(df_PARA[i,"ko"],","))),collapse=",")}
    if (str_count(df_PARA[i,"ko"],",") > 2) { df_PARA[i,"ko"] <- "multiple"}
    if (df_PARA[i,"ko"] == "") { df_PARA[i,"ko"] <- "unknown"}
  }
  #
  # Process subset network --------------------------------------------------
  stat_CC_withPARA <- as.data.frame(stat_CC_withPARA)
  first_graph <- stat_CC_withPARA[1,1]
  graph_ix <- graph_decomposed[[first_graph]]
  V(graph_ix)$CC <- first_graph
  V(graph_ix)$ko_compo <- paste(first_graph,df_PARA %>% filter(CC==first_graph) %>% select(ko), sep = " - ")
  for (i in stat_CC_withPARA[2:1679,"CC"]) { 
    V(graph_decomposed[[i]])$CC <- i
    V(graph_decomposed[[i]])$ko_compo <- paste(i, df_PARA %>% filter(CC == i) %>% select(ko), sep = " - ")
    graph_ix <- as_tbl_graph(graph_ix) %>% tidygraph::bind_graphs(graph_decomposed[[i]])}
  # Modification of attribute  ----------------------------------------------
  # Trophic role annotations
  V(graph_ix)$TRT_1[V(graph_ix)$TRT_1=="SAP/HET"] <- "HET"
  V(graph_ix)$TRT_1[V(graph_ix)$TRT_1=="SAP" & V(graph_ix)$TYPE_1=="Obligate"] <- "PARA"
  V(graph_ix)$TRT_1[V(graph_ix)$TRT_2=="PARA"] <- "PARA"
  V(graph_ix)$TRT_2[V(graph_ix)$TRT_1=="PARA"] <- "yes"
  V(graph_ix)$TRT_2[is.na(V(graph_ix)$TRT_2)==TRUE] <- "no"
  V(graph_ix)$TRT_1[is.na(V(graph_ix)$TRT_1)==TRUE] <- "unknown"
  #
  # Taxonomy & Function
  for ( i in unique(V(graph_ix)$Taxonomy)) {
    temp <- tablearchivex[,"Phylum"][tablearchivex[,"Taxonomy"]==i]
    print(temp)
    V(graph_ix)$Class[V(graph_ix)$Taxonomy==i] <- temp
  }
  V(graph_ix)$Class[is.na(V(graph_ix)$Class)==TRUE] <- "unknown"
  V(graph_ix)$ko_yn <- "known"
  V(graph_ix)$ko_yn[is.na(V(graph_ix)$ko)==TRUE] <- "unknown"
  #
  # Cross taxonomy plot -----------------------------------------------------
  tax_filter <- cbind(V(graph_ix)$CC,V(graph_ix)$Class)
  tax_filter <- as.data.frame(tax_filter) %>% group_by(V1) %>% summarize(V2=paste(unique(V2),collapse=","))
  tax_filter$V2 <- str_split(tax_filter$V2,",")
  svglite("Igraph/CC_containing75PARA_CC_PARA_obligate.svg",width = 10.00,height = 10.00)
  ggplot(tax_filter,aes(x=V2)) +
    geom_bar(color = "black", fill = "#880e4fcc",linewidth=0.5) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 2),breaks = c(0,2,8,32,128,512,1024)) +
    scale_x_upset() +
    #geom_text(stat = "count", aes(label = after_stat(count)), vjust = -1) +
    labs(y="CC #",x="Phylum") + theme_unique_art() +
    theme(legend.title = element_text(face="bold"),
          axis.title = element_text(color = "black", face = "bold"),
          strip.text.x = element_text(color = "black", face = "bold", size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()
  #
  # Statistics plot --------------------------------------------------------------------
  stat_subset <- CC_withPARA
  stat_subset$sum <- 1
  stat_subset <- stat_subset %>% select(CC,sum) %>% group_by(CC) %>% summarise_all(sum)
  stat_subset$Category <- "obligate PARA"
  #
  # Plot
  svglite(paste0("Analysis/CC_containing_obligate_PARA_Stat.svg"),width = 3.00,height = 4.50)
  plt <- ggbetweenstats(
    data = stat_subset,
    x = Category,
    y = sum,
    results.subtitle =FALSE,
    bf.message = FALSE,
    centrality.label.args = list(size  = 4,nudge_x=0.1,nudge_y=-1.5)
  ) + scale_y_continuous(trans = scales::pseudo_log_trans(base = 2),breaks = c(0,2,4,8,16,32,64,128,256),limits =c(2,256)) +
    theme_unique_art() + 
    labs(x="",y="CC Abundance") + guides(color = "none")
  print(plt)
  dev.off()
#
  # Function plot --------------------------------------------------------------------
  summarize_function  <- separate_rows(summarize_datax,"ko",sep=",")
  summarize_function <- left_join(summarize_function,ko_db %>% select(KO_id,lvl_B_val,lvl_C_val) %>% distinct(),by = c("ko"="KO_id"),relationship = "many-to-many")
  summarize_function$lvl_B_val <- sapply(str_replace_all(as.character(summarize_function$lvl_B_val),"_"," "), function(x) x)
  summarize_function$lvl_C_val <- sapply(str_replace_all(as.character(summarize_function$lvl_C_val),"_"," "), function(x) x)
  summarize_function[,"lvl_B_val"][is.na(summarize_function[,"lvl_B_val"]) == TRUE] <- "Unassigned"
  summarize_function[,"lvl_C_val"][is.na(summarize_function[,"lvl_C_val"]) == TRUE] <- "Unassigned"
  summarize_function$count <- 1
  summarize_function <- summarize_function %>% select(lvl_B_val,lvl_C_val,count) %>% group_by(lvl_B_val,lvl_C_val) %>% summarise_all(sum)
  summarize_function <- summarize_function %>% mutate(percent = count*100/sum(summarize_function$count))
  summarize_function <- as.data.frame(summarize_function)
  summarize_function$percent <- as.numeric(summarize_function$percent)
  summarize_function$count <- as.numeric(summarize_function$count)
  tempx <- c(t(summarize_function %>% filter(lvl_B_val != "Unassigned") %>% select(-lvl_B_val,-lvl_C_val)  %>% summarise_all(sum) %>% select(count,percent)))
  summarize_function <- rbind(summarize_function,c("Assigned","Assigned",tempx))
  summarize_function$percent <- as.numeric(summarize_function$percent)
  summarize_function$count <- as.numeric(summarize_function$count)
  summarize_function <- summarize_function %>% mutate(label = paste0(lvl_C_val,": ",round(percent,1),"%"))
  #
  val_tab <- as.data.frame(summarize_function %>% filter(lvl_B_val %in% c("Unassigned","Assigned")) %>% select(count) %>% group_by() %>% summarise_all(sum))
  val_tab <- cbind(val_tab,as.data.frame(summarize_function %>% filter(lvl_B_val %in% c("Assigned")) %>% select(label)))
  # Pre-plot
  t <- ggplot() + 
    geom_treemap(data = summarize_function %>% filter(lvl_B_val %in% c("Unassigned","Assigned")), aes(area = percent, fill = lvl_B_val,subgroup = lvl_B_val, label = label),color="white",size=2) +
    scale_fill_manual(values=c("black","white")) +
    geom_treemap_subgroup_border(size=2,color="black") + theme_unique_darktris() +
    geom_text(data = val_tab, mapping=aes(label = paste0(count," proteins","\n",label)),fontface = "bold", x = 0,hjust =0 ,y = 0.5,color="black",size=4) +
    theme(panel.border = element_rect(fill=NA,linewidth=2),legend.position="none")
  print(t)
  # Color
  summarize_function_i <- summarize_function %>% filter(!lvl_B_val %in% c("Unassigned","Assigned"))
  color_ko <- data.frame(unique(summarize_function$lvl_B_val),c(paletteer_d("MetBrewer::Signac",direction = -1))[1:length(unique(summarize_function$lvl_B_val))])
  colnames(color_ko) <- c("lvl_B_val","color")
  summarize_function_i <- merge(summarize_function_i,color_ko,by="lvl_B_val")
  # Legend
  order <- summarize_function_i %>% select(lvl_B_val,percent) %>% group_by(lvl_B_val) %>% summarise_all(sum) %>% arrange(-percent)
  color_ko <- left_join(order,color_ko,by="lvl_B_val")
  color_ko$lvl_B_val <- factor(color_ko$lvl_B_val,levels=order)
  legend_b <- get_legend(ggplot(data = summarize_function_i %>% filter(percent > 0.02) %>% arrange(percent),mapping = aes(y= percent, x = 1,fill = factor(lvl_B_val,levels= order$lvl_B_val) ,group=factor(lvl_B_val,levels= order$lvl_B_val)),color="black", Rowv = NA, scale = "column",alpha=0.8) + 
                           theme_unique_darkcris() + labs(fill ="KEGG category") +
                           geom_bar(stat="identity") + guides(fill = guide_legend(ncol=4,override.aes = list(color="black",shape=1))) + 
                           scale_fill_manual(values=noquote(color_ko$color)))
  # Final plot
  s <- ggplot(summarize_function_i, aes(area = percent, fill = lvl_B_val,subgroup = lvl_B_val, label = lvl_C_val)) + 
    geom_treemap(alpha=1,color="white",size=2, fill=summarize_function_i[,"color"]) + 
    geom_treemap_subgroup_border(size=2,color="black") + theme_unique_darktris() +
    geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=5) +
    theme(panel.border = element_rect(fill=NA,linewidth=2),strip.text.y.left = element_blank()) + theme(legend.position="none")
  svglite(paste0("Analysis/KO_compo_subset75obligate_category.svg"),width = 15.00,height = 6.00)
  b_plot <- plot_grid(t,s,NA,legend_b, ncol = 2, nrow = 2, align = "h",rel_widths = c(1.2,6),rel_heights = c(6,2))
  print(b_plot)
  dev.off()
#
  # Taxonomy plot -----------------------------------------------------------
  summarise_tax <- summarize_datax
  summarise_tax$count <- 1
  summarise_tax$Phylum[is.na(summarise_tax[,"Phylum"])==TRUE & is.na(summarise_tax[,"Class"])==FALSE] <- summarise_tax$Class[is.na(summarise_tax[,"Phylum"])==TRUE & is.na(summarise_tax[,"Class"])==FALSE]
  summarise_tax$Class[is.na(summarise_tax[,"Phylum"])==FALSE & is.na(summarise_tax[,"Class"])==TRUE & is.na(summarise_tax[,"Division"])==FALSE] <- paste0("Unaffiliated_",summarise_tax$Phylum[is.na(summarise_tax[,"Phylum"])==FALSE & is.na(summarise_tax[,"Class"])==TRUE & is.na(summarise_tax[,"Division"])==FALSE])
  summarise_tax$Phylum[is.na(summarise_tax[,"Phylum"])==TRUE & is.na(summarise_tax[,"Class"])==TRUE & is.na(summarise_tax[,"Division"])==FALSE] <- paste0("Unaffiliated_",summarise_tax$Division[is.na(summarise_tax[,"Phylum"])==TRUE & is.na(summarise_tax[,"Class"])==TRUE & is.na(summarise_tax[,"Division"])==FALSE])
  summarise_tax$Class[is.na(summarise_tax[,"Class"])==TRUE & is.na(summarise_tax[,"Phylum"])==FALSE & is.na(summarise_tax[,"Division"])==FALSE] <- summarise_tax$Phylum[is.na(summarise_tax[,"Phylum"])==FALSE & is.na(summarise_tax[,"Class"])==TRUE & is.na(summarise_tax[,"Division"])==FALSE]
  summarise_tax[is.na(summarise_tax)==TRUE] <- "Unaffiliated"
  summarise_tax <- summarise_tax %>% select(Division,Phylum,Class,count) %>% group_by(Division,Phylum,Class,) %>% summarise_all(sum)
  summarise_tax <- summarise_tax %>% mutate(percent = count*100/sum(summarise_tax$count))
  summarise_tax <- as.data.frame(summarise_tax)
  tempx <- c(t(summarise_tax %>% filter(Division != "Unaffiliated") %>% select(-Division,-Phylum,-Class) %>% summarise_all(sum) %>% select(count,percent)))
  summarise_tax <- rbind(summarise_tax,c("Affiliated","Affiliated","Affiliated",tempx))
  summarise_tax$percent <- as.numeric(summarise_tax$percent)
  summarise_tax$count <- as.numeric(summarise_tax$count)
  summarise_tax <- summarise_tax %>% mutate(label = paste0(Division,": ",round(percent,1),"%"))
  #
  val_tab <- as.data.frame(summarise_tax %>% filter(Division %in% c("Unaffiliated","Affiliated")) %>% select(count) %>% summarise_all(sum))
  val_tab <- cbind(val_tab,as.data.frame(summarise_tax %>% filter(Division %in% c("Affiliated")) %>% select(label)))
  # Pre-plot
  u <- ggplot() + 
    geom_treemap(data = summarise_tax %>% filter(Division %in% c("Unaffiliated","Affiliated")), aes(area = percent, fill = Division,subgroup = Division, label = label),color="white",size=2) +
    scale_fill_manual(values=c("black","white")) +
    geom_treemap_subgroup_border(size=2,color="black") + theme_unique_darktris() +
    theme(panel.border = element_rect(fill=NA,linewidth=2),legend.position="none") +
    geom_label(data = val_tab, mapping=aes(label = paste0(count," proteins","\n",label)),fontface = "bold", x = 0.5 ,y = 0.5,color = "white",fill="black",size=4,alpha=0.9)
  print(u)
  # Color
  summarise_tax_i <- summarise_tax %>% filter(!Division %in% c("Unaffiliated","Affiliated"))
  color <- data.frame(color_definition,c("#c62728cc","#6b1b9acc","white","#273493cc","white","#0277bdcc","#558b30cc","#f9a825cc","white","#fdffcccc","#ef6c00cc","white","white","white","white","#75252595","#880e4fcc","#90ee90cc","#20606490","#add8e6cc","#e53368cc"))
  colnames(color) <- c("Division","color")
  summarise_tax_i <- merge(summarise_tax_i,color,by="Division")
  # Legend
  order <- summarise_tax_i %>% select(Division,percent) %>% group_by(Division) %>% summarise_all(sum) %>% arrange(-percent)
  color <- left_join(order,color,by="Division")
  color$Division <- factor(color$Division,levels=order)
  legend_a <- get_legend(ggplot(data = summarise_tax_i %>% filter(percent > 0.02) %>% arrange(percent),mapping = aes(y= percent, x = 1,fill = factor(Division,levels= order$Division),group=factor(Division,levels= order$Division)),color="black", Rowv = NA, scale = "column",alpha=0.8) + 
                           theme_unique_darkcris() +  labs(fill ="Taxonomy") +
                           geom_bar(stat="identity") + guides(fill = guide_legend(ncol=6,override.aes = list(color="black",shape=1))) + 
                           scale_fill_manual(values=noquote(color$color)))
  # Final plot
  v <- ggplot(summarise_tax_i, aes(area = percent, fill = Division,subgroup = Division, label = Class)) + 
    geom_treemap(alpha=0.8,color="white",size=2, fill=summarise_tax_i[,"color"]) +
    guides(fill=guide_legend(ncol=6,title.position = "top")) +
    geom_treemap_subgroup_border(size=2,color="black") + theme_unique_darktris() +
    geom_treemap_text(fontface = "bold", colour = "black", place = "centre",grow = FALSE,size=12,min.size=5) +
    theme(panel.border = element_rect(fill=NA,linewidth=2),strip.text.y.left = element_blank()) + theme(legend.position="none")
  svglite(paste0("Analysis/taxo_compo_subset75obligate_category.svg"),width = 14.00,height = 6.00)
  v_plot <- plot_grid(u,v,NA,legend_a, ncol = 2, nrow = 2, align = "h",rel_widths = c(1,6),rel_heights = c(7,1))
  print(v_plot)
  dev.off()
#
  # Polar Taxonomy Plot ------------------------------------------------------------
    # Prepare Polar data -----------------------------------------------
    Polar_data_all <- summarize_datax %>%
      arrange(Division) %>%
      select(Grp_1,Grp_2,Division,Phylum,Class,) %>%
      group_by(Grp_1,Grp_2,Division,Phylum,Class) %>% count(name="count")
    #
    Polar_data_all[is.na(Polar_data_all)==TRUE] <- "Unassigned"
    Polar_data_all[,"Division"][Polar_data_all[,"Division"]=="Unassigned"] <- "Unaffiliated"
    Polar_data_all[,"Phylum"][Polar_data_all[,"Phylum"]=="Unassigned"] <- "Unaffiliated"
    Polar_data_all[,"Class"][Polar_data_all[,"Class"]=="Unassigned"] <- "Unaffiliated"
    #
    Polar_data_all <- as.data.frame(Polar_data_all,stringsAsFactors=FALSE)
    Proportion_aff_x <- Polar_data_all %>% select(-Phylum,-Class) %>% group_by(Grp_1,Grp_2,Division) %>% summarise_all(sum)
    Prop_value <- 100-Proportion_aff_x$count[Proportion_aff_x[,"Division"]=="Unaffiliated"] * 100 / sum(Proportion_aff_x$count)
    Proportion_aff <- paste("Affiliated:",round(Prop_value,1),"%")
    Polar_data_all <- Polar_data_all %>% filter(Division != "Unaffiliated")
    Polar_data <- Polar_data_all %>% select(-Grp_2,-Phylum,-Class) %>% group_by(Grp_1,Division) %>% summarise_all(sum)
    Polar_data <- as.data.frame(Polar_data,stringsAsFactors=FALSE)
    for (i in rownames(Polar_data)) { 
      Polar_data[i,"Proportion"] <- (Polar_data[i,"count"] * 100) / sum(Polar_data$count)
      Division_i <- Polar_data[i,"Division"]
      Polar_data[i,"Proportion_byDiv"] <- (Polar_data[i,"count"] * 100) / colSums(Polar_data %>% filter(Division == Division_i) %>% select(count))
    }
    #
    # Polar_data_all (Grp_2) -----------------------------------------------
    Polar_data_all <-  Polar_data_all %>% select(-Phylum,-Class) %>% group_by(Grp_1,Grp_2,Division) %>% summarise_all(sum)
    for (i in rownames(Polar_data_all)) { 
      Polar_data_all[i,"Proportion"] <- (Polar_data_all[i,"count"] * 100) / sum(Polar_data_all$count)
    }
    Polar_data_all <- Polar_data_all %>% arrange(Division)
    Polar_data_all$color <- NA
    Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="Obligate"] <- "lightgreen"
    Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="Facultative"] <- "lightblue"
    Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="CM"] <- "#CCCCCC"
    Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="eSNCM"] <- "darkgrey"
    Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="pSNCM"] <- "#646464"
    Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="eSNCM/pSNCM"] <- "#333333"
    Polar_data_all[,"color"][Polar_data_all[,"Grp_2"]=="Unassigned"] <- "white"
    Polar_data_all <- Polar_data_all %>%
      arrange(Division) %>%
      mutate(lab.ypos = cumsum(Proportion) - 0.5*Proportion)
    Polar_data_all$label <- paste(round(Polar_data_all$Proportion,1), "%", sep = "")
    for (i in rownames(Polar_data_all)) {
      if (Polar_data_all[i,"label"] == "0%") { Polar_data_all[i,"label"] <- NA}}
    for (i in rownames(Polar_data_all)) {
      if (is.na(Polar_data_all[i,"label"]) == FALSE) { Polar_data_all[i,"label"] <- paste0(Polar_data_all[i,"label"])}
    }
    #
    # Polar_data (Grp_1) -----------------------------------------------
    Polar_data <- Polar_data %>%
      arrange(Division) %>%
      mutate(lab.ypos = cumsum(Proportion) - 0.5*Proportion)
    Polar_data$label <- paste(round(Polar_data$Proportion_byDiv,1), "%", sep = "")
    for (i in rownames(Polar_data)) {
      if (Polar_data[i,"label"] == "0%") { Polar_data[i,"label"] <- NA}}
    for (i in rownames(Polar_data)) {
      if (is.na(Polar_data[i,"label"]) == FALSE) { Polar_data[i,"label"] <- paste0(Polar_data[i,"label"])}
    }
    Polar_data[,"label"][Polar_data[,"Grp_1"]=="Unassigned"]<-NA
    #
    # Polar_data_1 (Division) -----------------------------------------------
    color <- data.frame(color_definition,c("#c62728cc","#6b1b9acc","white","#273493cc","white","#0277bdcc","#558b30cc","#f9a825cc","white","#fdffcccc","#ef6c00cc","white","white","white","white","#75252595","#880e4fcc","#90ee90cc","#20606490","#add8e6cc","#e53368cc"))
    colnames(color) <- c("Division","color")
    #
    Polar_data_a <- Polar_data %>% select(-Grp_1,-label) %>% 
      arrange(Division) %>%
      group_by(Division) %>% 
      summarise_all(sum)
    
    for (i in row.names(Polar_data_a)) { if (Polar_data_a[i,"Proportion"] > 0.1) {
      Polar_data_a[i,"color"] <- color[,"color"][color[,"Division"]==Polar_data_a[i,"Division"]$Division]}
      else {Polar_data_a[i,"color"] <- "white"}
    }
    #
    Polar_data$color <- NA
    Polar_data[,"color"][Polar_data[,"Grp_1"]=="HET"] <- "grey48"
    Polar_data[,"color"][Polar_data[,"Grp_1"]=="Unassigned"] <- "white"
    Polar_data[,"color"][Polar_data[,"Grp_1"]=="PHOTO"] <- "blue4"
    Polar_data[,"color"][Polar_data[,"Grp_1"]=="MIXO"] <- "red4"
    Polar_data[,"color"][Polar_data[,"Grp_1"]=="SAP"] <- "grey0"
    Polar_data[,"color"][Polar_data[,"Grp_1"]=="PARA"] <- "darkgreen"
    #
    # legend ------------------------------------------------------------------
    colnames(Polar_data)[colnames(Polar_data)=="Grp_1"] <- "Trophic mode"
    colnames(Polar_data_all)[colnames(Polar_data_all)=="Grp_2"] <- "Parasitic lifestyle"
    colnames(Polar_data_all)[colnames(Polar_data_all)=="Grp_1"] <- "Trophic mode"
    colnames(Polar_data_a)[colnames(Polar_data_a)=="Division"] <- "Taxonomic group"
    #
    #Legend
    Polar_data_b <- Polar_data_a %>% filter(Proportion > 0.1)
    Polar_data_b <- Polar_data_b %>% mutate(`Taxonomic group`=paste0(`Taxonomic group`,": ",round(Proportion,1),"%"))
    legend_a <- get_legend(ggplot(data = Polar_data_b) + theme_unique_darkcris() +
                             geom_bar(mapping = aes(y= Proportion, fill = `Taxonomic group`)) + guides(fill = guide_legend(override.aes = list(color="black",shape=1))) + 
                             scale_fill_manual(values=noquote(Polar_data_b$color)))
    Polar_data_legend <- Polar_data %>% select(`Trophic mode`,Proportion) %>% group_by(`Trophic mode`) %>% summarise_all(sum)
    Polar_data_legend <- Polar_data_legend %>% mutate(`Trophic mode` = paste0(`Trophic mode`,": ",round(Proportion,1),"%"))
    legend_b <- get_legend(ggplot(data = Polar_data_legend %>% arrange(`Trophic mode`)) + theme_unique_darkcris() + guides(fill="Functional assignment") +
                             geom_bar(mapping = aes(y= Proportion, fill = `Trophic mode`)) + guides(fill = guide_legend(override.aes = list(color="black",shape=1))) + 
                             scale_fill_manual(values=noquote(unique((Polar_data %>% arrange(`Trophic mode`))$color))))
    Polar_data_all_legend <- as.data.frame(Polar_data_all) %>% select(`Parasitic lifestyle`,Proportion) %>% group_by(`Parasitic lifestyle`) %>% summarise_all(sum)
    Polar_data_all_legend <- Polar_data_all_legend %>% mutate(`Parasitic lifestyle` = paste0(`Parasitic lifestyle`,": ",round(Proportion,1),"%"))
    legend_c <- get_legend(ggplot(data = Polar_data_all_legend %>% arrange(`Parasitic lifestyle`)) + theme_unique_darkcris() + guides(fill="Parasitic lifestyle") +
                             geom_bar(mapping = aes(y= Proportion, fill = `Parasitic lifestyle`)) + guides(fill = guide_legend(override.aes = list(color="black",shape=1))) + 
                             scale_fill_manual(values=noquote((Polar_data_all %>% select(`Parasitic lifestyle`, color) %>% distinct() %>% arrange(`Parasitic lifestyle`,.locale="en"))$color)))
    #
    # plot --------------------------------------------------------------------
    Polar_data_all <- as.data.frame(Polar_data_all)
    Polar_data_a <- as.data.frame(Polar_data_a)
  #Plot large graph
  bxCop <- ggplot(Rowv = NA, col = colMain, scale = "column") +
    theme_unique_darkcris() +
    labs(y = "",x="") + coord_polar("y") + 
    geom_bar(data = Polar_data_all, mapping = aes(y= Proportion, x = 2.85, fill = label,width = 0.1), stat="identity", color = "black", width = 1, linewidth = 0.3, fill=Polar_data_all[,"color"]) + 
    geom_bar(data = Polar_data, mapping = aes(y= Proportion, x = 2.65, fill = label,width = 0.3), stat="identity", color = "black", width = 1, linewidth= 0.3, fill=Polar_data[,"color"]) + 
    scale_y_continuous(limits=c(0,sum(Polar_data %>% select(Proportion))),labels = scales::percent_format(scale=1),breaks = c(25,50,75,100)) +
    geom_bar(data = Polar_data_a,  mapping = aes(y= Proportion, x = 2, fill = rev(Division),width = 1), stat="identity", color = "black", width = 1, linewidth =0.3,fill=Polar_data_a[,"color"]) + 
    xlim(1,3) + theme(axis.text.y=element_blank(),
                      axis.ticks.y=element_blank()) +
    annotate("text",label=paste0(Proportion_aff,"\n(",sum(Proportion_aff_x$count),")"),x=1,y=25, fontface =2)
  print(bxCop)
  #
  svglite("Analysis/Stat_CC_75percent_obligate.svg",width = 10.00,height = 8.00)
  legend <- plot_grid(NULL,legend_a, legend_b,legend_c, NULL, ncol = 1, nrow = 5,align = "v",rel_heights = c(0.3,2.3,1.5,0.5,0.5))
  b_plot <- plot_grid(bxCop, legend, ncol = 2, nrow = 1 , rel_widths = c(6,2),rel_heights = c(3))
  print(b_plot)
  dev.off()
  #
# Display Network graph --------------------------------------------------------------
  # Select PARA CCs which contains sequences labeled with obligate parasite ----------
    Obligate_PARA_CC <- unique(stat_input_df[,"CC"][stat_input_df[,"TRT_1"]=="PARA" & stat_input_df[,"TYP_1"]=="Obligate"])
    interest_CC_majorPARA <- PARA[PARA %in% Obligate_PARA_CC]
    CC_withPARA <- stat_input_df %>% filter(CC %in% interest_CC_majorPARA)
    CC_withPARA$sum <- 1
    summary(CC_withPARA %>% select(CC,sum) %>% group_by(CC) %>% summarise_all(sum) %>% select(sum))
    length(unique(CC_withPARA$CC))
    stat_CC_withPARA <- CC_withPARA %>% select(CC,sum) %>% group_by(CC) %>% summarise_all(sum) %>% filter(sum <= 226 & sum >= 36) %>% select(CC)
    nrow(stat_CC_withPARA)
    summarize_datax <- summarize_datax %>% filter(ID_Protein %in% CC_withPARA$name)
  # Process subset network --------------------------------------------------
  stat_CC_withPARA <- as.data.frame(stat_CC_withPARA)
  first_graph <- stat_CC_withPARA[1,1]
  graph_ix <- graph_decomposed[[first_graph]]
  V(graph_ix)$CC <- first_graph
  V(graph_ix)$ko_compo <- paste(first_graph,df_PARA %>% filter(CC==first_graph) %>% select(ko), sep = " - ")
  for (i in stat_CC_withPARA[2:50,"CC"]) { 
    V(graph_decomposed[[i]])$CC <- i
    V(graph_decomposed[[i]])$ko_compo <- paste(i, df_PARA %>% filter(CC == i) %>% select(ko), sep = " - ")
    graph_ix <- as_tbl_graph(graph_ix) %>% tidygraph::bind_graphs(graph_decomposed[[i]])}
#
  # Modification of attribute  ----------------------------------------------
    # SAP/HET -> HET
  V(graph_ix)$TRT_1[V(graph_ix)$TRT_1=="SAP/HET"] <- "HET"
  V(graph_ix)$TRT_1[V(graph_ix)$TRT_1=="SAP" & V(graph_ix)$TYPE_1=="Obligate"] <- "PARA"
  V(graph_ix)$TRT_1[V(graph_ix)$TRT_2=="PARA"] <- "PARA"
  V(graph_ix)$TRT_2[V(graph_ix)$TRT_1=="PARA"] <- "yes"
  V(graph_ix)$TRT_2[is.na(V(graph_ix)$TRT_2)==TRUE] <- "no"
  V(graph_ix)$TRT_1[is.na(V(graph_ix)$TRT_1)==TRUE] <- "unknown"
#
    # Taxonomy
  for ( i in unique(V(graph_ix)$Taxonomy)) {
    temp <- tablearchivex[,"Class"][tablearchivex[,"Taxonomy"]==i]
    print(temp)
    V(graph_ix)$Class[V(graph_ix)$Taxonomy==i] <- temp
  }
  V(graph_ix)$Class[is.na(V(graph_ix)$Class)==TRUE] <- "unknown"
  V(graph_ix)$ko_yn <- "known"
  V(graph_ix)$ko_yn[is.na(V(graph_ix)$ko)==TRUE] <- "unknown"
#
  # plot --------------------------------------------------------------------
    graph_ix <- graph_ix %>%
      mutate(centrality = centrality_authority())
  svglite("Igraph/CC_containing75PARA_CC_PARA_obligate_KO_50most_abundant.svg",width = 26.00,height = 14.00)
    ggraph(graph_ix, layout = 'stress') + 
      geom_edge_link(color="#7b7b7b7b") + 
      geom_node_point(aes(shape=TRT_1,fill=Class,color=ko_yn),size = 3) +
      #geom_node_text(aes(label = ko), colour = 'black', vjust = 0.4) + 
      scale_fill_manual(values = c("#303f9f","#ef5350","#fbc02c","maroon","grey","#ee7620","red4","limegreen","#4cc6da","black","#98a5ff","lightgreen","blue4","#2B9B81FF","#1F6E9CFF","#41acc1","#E87B89FF","#424242","#E69F00","white","green","#92C051FF","green4","#633372FF","#56B4E9","#009E73","pink","green","#F0E442","#0072B2")) + 
      scale_shape_manual(values = c(25,22,23,21,24,20)) + 
      scale_color_manual(values = c("black","white")) + 
      scale_size_manual(values = c(3, 5)) + 
      facet_nodes(.~ko_compo, scales = "free",nrow = 7) + 
      labs(title="",fill="Taxonomy", shape = "Trophic modes", size = "Parasite Y/N", color = "Kegg ID") + 
      guides(fill = guide_legend(ncol=1,override.aes = list(shape = 21, size = 4)), color = guide_legend(override.aes = list(shape = 21, size = 4, fill = "grey")),size = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(size = 4))) +
      theme_unique_art_graph()
  dev.off()
  #
# Display CCs Network clustering interested KO (40 barycenters) --------------------------------------
  # Select PARA CCs which contains sequences labeled with parasite and with interested KO (40 barycenters) --------
  Obligate_PARA_CC <- unique(stat_input_df[,"CC"][stat_input_df[,"TRT_1"]=="PARA"])# & stat_input_df[,"TYP_1"]=="Obligate"])
  interest_CC_majorPARA <- PARA[PARA %in% Obligate_PARA_CC]
  # Select CCs appearing in Barycenters
  CC_Barycenter <- c(unlist(str_split(paste(df_bary_subsetx$CC,collapse=","),pattern=",")))
  interest_CC_majorPARA <- interest_CC_majorPARA[interest_CC_majorPARA %in% CC_Barycenter]
  #
  CC_withPARA <- stat_input_df %>% filter(CC %in% interest_CC_majorPARA)
  CC_withPARA$sum <- 1
  summary(CC_withPARA %>% select(CC,sum) %>% group_by(CC) %>% summarise_all(sum) %>% select(sum))
  length(unique(CC_withPARA$CC))
  stat_CC_withPARA <- CC_withPARA %>% select(CC,sum) %>% group_by(CC) %>% summarise_all(sum) %>% filter(sum <= 226 & sum >= 3) %>% select(CC)
  nrow(stat_CC_withPARA)
  summarize_datax <- summarize_datax %>% filter(ID_Protein %in% CC_withPARA$name)
#
  # Add functional information ---------------------------------
  df_PARA <- CC_withPARA %>% select(CC, ko) %>% group_by(CC) %>% summarize(ko=paste(unique(ko),collapse=","))
  df_PARA$ko <- sapply(str_remove(as.character(df_PARA$ko),"^,"), function(x) x)
  df_PARA$ko <- sapply(str_remove(as.character(df_PARA$ko),",$"), function(x) x)
  df_PARA$ko <- sapply(str_replace_all(as.character(df_PARA$ko),",,",","), function(x) x)
  df_PARA$ko <- sapply(str_split(df_PARA$ko,","), function(x) paste(unique(x),collapse = ","))
  df_PARA <- as.data.frame(df_PARA)
  for (i in row.names(df_PARA)) { 
    if (grepl(",",df_PARA[i,"ko"])==TRUE){ df_PARA[i,"ko"] <- paste(unique(unlist(str_split(df_PARA[i,"ko"],","))),collapse=",")}
    if (str_count(df_PARA[i,"ko"],",") > 2) { df_PARA[i,"ko"] <- "multiple"}
    if (df_PARA[i,"ko"] == "") { df_PARA[i,"ko"] <- "unknown"}
  }
  # Process subset network --------------------------------------------------
  stat_CC_withPARA <- as.data.frame(stat_CC_withPARA)
  first_graph <- stat_CC_withPARA[1,1]
  graph_ix <- graph_decomposed[[first_graph]]
  V(graph_ix)$CC <- first_graph
  V(graph_ix)$ko_compo <- paste(first_graph,df_PARA %>% filter(CC==first_graph) %>% select(ko), sep = " - ")
  for (i in stat_CC_withPARA[2:74,"CC"]) { 
    V(graph_decomposed[[i]])$CC <- i
    V(graph_decomposed[[i]])$ko_compo <- paste(i, df_PARA %>% filter(CC == i) %>% select(ko), sep = " - ")
    graph_ix <- as_tbl_graph(graph_ix) %>% tidygraph::bind_graphs(graph_decomposed[[i]])}
  #
  # Modification of attribute  ----------------------------------------------
  # SAP/HET -> HET
  V(graph_ix)$TRT_1[V(graph_ix)$TRT_1=="SAP/HET"] <- "HET"
  V(graph_ix)$TRT_1[V(graph_ix)$TRT_1=="SAP" & V(graph_ix)$TYPE_1=="Obligate"] <- "PARA"
  V(graph_ix)$TRT_1[V(graph_ix)$TRT_2=="PARA"] <- "PARA"
  V(graph_ix)$TRT_2[V(graph_ix)$TRT_1=="PARA"] <- "yes"
  V(graph_ix)$TRT_2[is.na(V(graph_ix)$TRT_2)==TRUE] <- "no"
  V(graph_ix)$TRT_1[is.na(V(graph_ix)$TRT_1)==TRUE] <- "unknown"
  #
  # Taxonomy
  for ( i in unique(V(graph_ix)$Taxonomy)) {
    temp <- tablearchivex[,"Class"][tablearchivex[,"Taxonomy"]==i]
    print(temp)
    V(graph_ix)$Class[V(graph_ix)$Taxonomy==i] <- temp
  }
  V(graph_ix)$Class[is.na(V(graph_ix)$Class)==TRUE] <- "unknown"
  V(graph_ix)$ko_yn <- "known"
  V(graph_ix)$ko_yn[is.na(V(graph_ix)$ko)==TRUE] <- "unknown"
  #
  # plot --------------------------------------------------------------------
  graph_ix <- graph_ix %>%
    mutate(centrality = centrality_authority())
  svglite("Igraph/CC_containing75PARA_CC_PARA_obligate_KO_correlated.svg",width = 26.00,height = 14.00)
  ggraph(graph_ix, layout = 'stress') + 
    geom_edge_link(color="#7b7b7b7b") + 
    geom_node_point(aes(shape=TRT_1,fill=Class, size=TRT_2,color=ko_yn)) +
    #geom_node_text(aes(label = ko), colour = 'black', vjust = 0.4) + 
    scale_fill_manual(values = c("#303f9f","#ef5350","#fbc02c","#0072B2","grey","#ee7620","red4","limegreen","#4cc6da","black","#98a5ff","lightgreen","maroon","blue4","#2B9B81FF","#1F6E9CFF","#41acc1","#E69F00","white","pink","#424242","green","#92C051FF","green4","#633372FF","#56B4E9","#009E73","pink","green","#F0E442","#0072B2")) + 
    scale_shape_manual(values = c(25,22,23,24,21,20)) + 
    scale_color_manual(values = c("black","white")) + 
    scale_size_manual(values = c(3, 5)) + facet_nodes(.~ko_compo, scales = "free",ncol = 7) + 
    labs(title="",fill="Taxonomy (fill)", shape = "Trophic modes", size = "Parasite Y/N", color = "Kegg ID") + 
    guides(fill = guide_legend(ncol=1,override.aes = list(shape = 21, size = 4)), color = guide_legend(override.aes = list(shape = 21, size = 4, fill = "grey")),size = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(size = 4))) +
    theme_unique_art_graph()
  dev.off()
  #
# Save R Image ------------------------------------------------------------
  save.image("R_Igraph.RData")
  #
# Annotation estimation -------------------------------------------------
  # Function -----------------------------------------------------------
  #ko
  df_network_function <- stat_input_df %>% select(CC, ko) %>% group_by(CC) %>% summarize(ko=paste(unique(ko),collapse=","))
  df_network_function$ko <- sapply(str_remove(as.character(df_network_function$ko),"^,"), function(x) x)
  df_network_function$ko <- sapply(str_remove(as.character(df_network_function$ko),",$"), function(x) x)
  df_network_function$ko <- sapply(str_replace_all(as.character(df_network_function$ko),",,",","), function(x) x)
  df_network_function$ko <- sapply(str_split(df_network_function$ko,","), function(x) paste(unique(x),collapse = ","))
  df_network_function <- as.data.frame(df_network_function)
  #
  df_network_function$count <- sapply(str_count(df_network_function$ko,","), function(x) x+1)
  df_network_function$count_u <- sapply(str_count(df_network_function$ko,"Unassigned"), function(x) x)
  df_network_function$count <- df_network_function$count - df_network_function$count_u
#
  # Taxonomy ----------------------------------------------------------------
  tablearchive <- stat_input_df
  tablearchive[tablearchive==""] <- NA
  # Affiliation
  tablearchivex <- tablearchive %>% select(Taxonomy) %>% distinct(Taxonomy,.keep_all = TRUE)
  for (i in row.names(tablearchivex)) { 
    j = 2
    while (str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][j] == "") { j<-j+1
    if (j == 3) {break}}
    k = j+1
    while (str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][k] == "") { k<-k+1
    if (k >= 5) {break}}
    l = k+1
    while (str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][l] == "") { l<-l+1
    if (l >= 7) {break}}
    tablearchivex[i,"Division"] <- str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][j]
    tablearchivex[i,"Phylum"] <- str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][k]
    tablearchivex[i,"Class"] <- str_split(tablearchivex[i,"Taxonomy"],pattern = "/[dkpcofg]_")[[1]][l]}
  tablearchivex[tablearchivex==""]<-NA
  tablearchivex[tablearchivex=="NA"]<-NA
  #  Process Inconsistencies of taxonomy
  summarize_data <- tablearchivex %>% select(Taxonomy,Division,Phylum,Class)
  summarize_data[,"Division"][summarize_data[,"Division"]=="unclassified eukaryotes"] <- NA
  summarize_data[,"Division"][summarize_data[,"Division"]=="environmental samples"] <- "Environmental samples"
  summarize_data[,"Division"][summarize_data[,"Phylum"]=="Jakobida"] <- "Jakobida"
  summarize_data[,"Division"][summarize_data[,"Phylum"]=="Malawimonadidae"] <- "Malawimonadidae"
  summarize_datax <- left_join(tablearchive,summarize_data,"Taxonomy")
  #
  for ( i in unique(tablearchive$Taxonomy)) {
    temp_Phylum <- summarize_datax[,"Phylum"][summarize_datax[,"Taxonomy"]==i]
    temp_Class <- summarize_datax[,"Class"][summarize_datax[,"Taxonomy"]==i]
    temp_Division <- summarize_datax[,"Division"][summarize_datax[,"Taxonomy"]==i]
    tablearchive$Division[tablearchive$Taxonomy==i] <- temp_Division
    tablearchive$Phylum[tablearchive$Taxonomy==i] <- temp_Phylum
    tablearchive$Class[tablearchive$Taxonomy==i] <- temp_Class
  }
  #
  tablearchive[is.na(tablearchive)==TRUE] <- "unknown"
  tablearchive[tablearchive==""] <- "unknown"
  df_network_tax <- tablearchive %>% select(CC, Division, Phylum, Class) %>% group_by(CC) %>% summarize(Class=paste(unique(Class),collapse=","),Phylum=paste(unique(Phylum),collapse=","),Division=paste(unique(Division),collapse=","))
  for (i in c("Class","Phylum","Division")) {
    df_network_tax[[i]] <- sapply(str_remove(as.character(df_network_tax[[i]]),"^,"), function(x) x)
    df_network_tax[[i]] <- sapply(str_remove(as.character(df_network_tax[[i]]),",$"), function(x) x)
    df_network_tax[[i]] <- sapply(str_replace_all(as.character(df_network_tax[[i]]),",,",","), function(x) x)
    df_network_tax[[i]] <- sapply(str_split(df_network_tax[[i]],","), function(x) paste(unique(x),collapse = ","))
  }
  df_network_tax <- as.data.frame(df_network_tax)
  #Class
  df_network_tax$count_Class <- sapply(str_count(df_network_tax$Class,","), function(x) x+1)
  df_network_tax$count_Class_u <- sapply(str_count(df_network_tax$Class,"unknow"), function(x) x)
  df_network_tax$count_Class <- df_network_tax$count_Class - df_network_tax$count_Class_u
  #Phylum
  df_network_tax$count_Phylum <- sapply(str_count(df_network_tax$Phylum,","), function(x) x+1)
  df_network_tax$count_Phylum_u <- sapply(str_count(df_network_tax$Phylum,"unknow"), function(x) x)
  df_network_tax$count_Phylum <- df_network_tax$count_Phylum - df_network_tax$count_Phylum_u
  #Division
  df_network_tax$count_Division <- sapply(str_count(df_network_tax$Division,","), function(x) x+1)
  df_network_tax$count_Division_u <- sapply(str_count(df_network_tax$Division,"unknow"), function(x) x)
  df_network_tax$count_Division <- df_network_tax$count_Division - df_network_tax$count_Division_u
#
  # Trophic mode -----------------------------------------------------------
  #TRT_1
  df_network_mode <- stat_input_df %>% select(CC, TRT_1) %>% group_by(CC) %>% summarize(TRT_1=paste(unique(TRT_1),collapse=","))
  df_network_mode$TRT_1 <- sapply(str_remove(as.character(df_network_mode$TRT_1),"^,"), function(x) x)
  df_network_mode$TRT_1 <- sapply(str_remove(as.character(df_network_mode$TRT_1),",$"), function(x) x)
  df_network_mode$TRT_1 <- sapply(str_replace_all(as.character(df_network_mode$TRT_1),",,",","), function(x) x)
  df_network_mode$TRT_1 <- sapply(str_split(df_network_mode$TRT_1,","), function(x) paste(unique(x),collapse = ","))
  df_network_mode <- as.data.frame(df_network_mode)
  #
  df_network_mode$count.grp <- sapply(str_count(df_network_mode$TRT_1,","), function(x) x+1)
  df_network_mode$count.grp_u <- sapply(str_count(df_network_mode$TRT_1,"Unassigned"), function(x) x)
  df_network_mode$count.grp <- df_network_mode$count.grp - df_network_mode$count.grp_u
  #
  # Merging and re-assignation -----------------------------------------------------------------
  # Merging
  stat_CC <- merge(df_network_function, df_network_tax, by="CC")
  stat_CC <- merge(stat_CC, df_network_mode, by="CC")
  # Count number of CCs and proteins in each conformation ((un)knwon taxonomy / (un)known function / (un)knwon trophic mode)
  # Taking into account SSN
  cat("Statistics using SSN:",file="Statistics_of_annotations_using_or_not_SSN.txt",sep="\n")
  cat(paste("Taxonomically affiliated & Functionnally annotated CCs:",nrow(stat_CC %>% filter(Division != "unknown") %>% filter(ko != "Unassigned"))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # 124455 full annotated OK
  cat(paste("Taxonomically affiliated & Functionnally annotated proteins:",nrow(stat_input_df %>% filter(CC %in% (stat_CC %>% filter(Division != "unknown") %>% filter(ko != "Unassigned"))$CC))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # representing 1252112 proteins
  cat(paste("Taxonomically unaffiliated & Functionnally unannotated CCs:",nrow(stat_CC %>% filter(Division == "unknown") %>% filter(ko == "Unassigned"))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # 48773 CCs unnanotated OK
  cat(paste("Taxonomically unaffiliated & Functionnally unannotated proteins:",nrow(stat_input_df %>% filter(CC %in% (stat_CC %>% filter(Division == "unknown") %>% filter(ko == "Unassigned"))$CC))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # representing 207330 proteins OK
  cat(paste("Taxonomically affiliated & Functionnally unannotated CCs:",nrow(stat_CC %>% filter(Division != "unknown") %>% filter(ko == "Unassigned"))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # 124455 full annotated OK
  cat(paste("Taxonomically affiliated & Functionnally unannotated proteins:",nrow(stat_input_df %>% filter(CC %in% (stat_CC %>% filter(Division != "unknown") %>% filter(ko == "Unassigned"))$CC))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # representing 1252112 proteins
  cat(paste("Taxonomically unaffiliated & Functionnally annotated CCs:",nrow(stat_CC %>% filter(Division == "unknown") %>% filter(ko != "Unassigned"))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # 48773 CCs unnanotated OK
  cat(paste("Taxonomically unaffiliated & Functionnally annotated proteins:",nrow(stat_input_df %>% filter(CC %in% (stat_CC %>% filter(Division == "unknown") %>% filter(ko != "Unassigned"))$CC))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # representing 207330 proteins OK
  cat(paste("Trophically assigned CCs:",nrow(stat_CC %>% filter(TRT_1 != "Unassigned"))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # 140303 CCS annotated and 162001 CCs unannotated
  cat(paste("Trophically assigned proteins:",nrow(stat_input_df %>% filter(CC %in% (stat_CC %>% filter(TRT_1 != "Unassigned"))$CC))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # representing 207330 proteins
  # Without SSN
  cat(print("Statistics without SSN:"),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n")
  cat(paste("Taxonomically affiliated & Functionnally annotated proteins:",nrow(summarize_datax %>% filter(is.na(Division) == FALSE) %>% filter(ko != "Unassigned"))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # 340482 full annotated
  cat(paste("Taxonomically unaffiliated & Functionnally annotated proteins:",nrow(summarize_datax %>% filter(is.na(Division) == TRUE) %>% filter(ko != "Unassigned"))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # 826681 with unknown taxonomy and known function
  cat(paste("Taxonomically affiliated & Functionnally unannotated proteins:",nrow(summarize_datax %>% filter(is.na(Division) == FALSE) %>% filter(ko == "Unassigned"))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # 363448 with unknown function and known taxonomy
  cat(paste("Taxonomically unaffiliated & Functionnally unannotated proteins:",nrow(summarize_datax %>% filter(is.na(Division) == TRUE) %>% filter(ko == "Unassigned"))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n") # 634495 without any annotations
  cat(paste("Trophically assigned proteins:",nrow(summarize_datax %>% filter(TRT_1 != "Unassigned"))),file="Statistics_of_annotations_using_or_not_SSN.txt",append=TRUE,sep="\n")
  #
# Re-annotation estimation using CCs with the same type of each annotations to annot unknown proteins -------------------------------------------------------------
  stat_CC <- stat_CC %>% filter(count == 1 & count_Class == 1 & count_Phylum == 1 & count_Division == 1 & count.grp == 1)
  CC_reassignable <- stat_CC$CC
  # correct_annotation
  CC_correct <- as.data.frame(stat_CC)
  CC_correct$newko <- sapply(unlist(str_split(CC_correct$ko,","))[!unlist(str_split(CC_correct$ko,",")) %in% "Unassigned"], function(x) paste(unique(x),collapse=","))
  CC_correct$newTRT_1 <- sapply(unlist(str_split(CC_correct$TRT_1,","))[!unlist(str_split(CC_correct$TRT_1,",")) %in% "Unassigned"], function(x) paste(unique(x),collapse=","))
  CC_correct$newClass <- sapply(unlist(str_split(CC_correct$Class,","))[!unlist(str_split(CC_correct$Class,",")) %in% "unknown"], function(x) paste(unique(x),collapse=","))
  CC_correct$newPhylum <- sapply(unlist(str_split(CC_correct$Phylum,","))[!unlist(str_split(CC_correct$Phylum,",")) %in% "unknown"], function(x) paste(unique(x),collapse=","))
  CC_correct$newDivision <- sapply(unlist(str_split(CC_correct$Division,","))[!unlist(str_split(CC_correct$Division,",")) %in% "unknown"], function(x) paste(unique(x),collapse=","))
  CC_correct <- as.data.frame(CC_correct %>% select(CC,newko,newDivision,newPhylum,newClass,newTRT_1)) %>% mutate(reannotation = "Yes")
  #
  stat_CCx <- tablearchive %>% select(CC, name,Division,Phylum,Class,TRT_1,ko) %>% mutate(reannotation = "No")
  stat_CCx$newko <- stat_CCx$ko
  stat_CCx$newTRT_1 <- stat_CCx$TRT_1
  stat_CCx$newDivision <- stat_CCx$Division
  stat_CCx$newPhylum <- stat_CCx$Phylum
  stat_CCx$newClass <- stat_CCx$Class
  # Statistics
  stat_CCx <- stat_CCx %>% rows_update(CC_correct, by = "CC")
  # Write file
  cat("Statistics of re-annotation using CCs with the same type of each annotations to annot unknown proteins:",file="Statistics_of_re-annotations.txt",sep="\n")
  cat(paste0("Functionnally annotated proteins with re-annotation: ",nrow(stat_CCx %>% filter(newko != "Unassigned")) *100 /nrow(stat_CCx),"% and without: ",nrow(stat_CCx %>% filter(ko != "Unassigned")) *100 /nrow(stat_CCx),"%"),file="Statistics_of_re-annotations.txt",append=TRUE,sep="\n")
  cat(paste0("Trophically assigned proteins with re-annotation: ",nrow(stat_CCx %>% filter(newTRT_1 != "Unassigned")) *100 /nrow(stat_CCx),"% and without: ",nrow(stat_CCx %>% filter(TRT_1 != "Unassigned")) *100 /nrow(stat_CCx),"%"),file="Statistics_of_re-annotations.txt",append=TRUE,sep="\n")
  cat(paste0("Taxonomically (Class) affiliated proteins with re-annotation: ",nrow(stat_CCx %>% filter(newClass != "unknown")) *100 /nrow(stat_CCx),"% and without: ",nrow(stat_CCx %>% filter(Class != "unknown")) *100 /nrow(stat_CCx),"%"),file="Statistics_of_re-annotations.txt",append=TRUE,sep="\n")
  cat(paste0("Taxonomically (Phylum) affiliated proteins with re-annotation: ",nrow(stat_CCx %>% filter(newPhylum != "unknown")) *100 /nrow(stat_CCx),"% and without: ",nrow(stat_CCx %>% filter(Phylum != "unknown")) *100 /nrow(stat_CCx),"%"),file="Statistics_of_re-annotations.txt",append=TRUE,sep="\n")
  cat(paste0("Taxonomically (Division) affiliated proteins with re-annotation: ",nrow(stat_CCx %>% filter(newDivision != "unknown")) *100 /nrow(stat_CCx),"% and without: ",nrow(stat_CCx %>% filter(Division != "unknown")) *100 /nrow(stat_CCx),"%"),file="Statistics_of_re-annotations.txt",append=TRUE,sep="\n")
#
# Save R Image ------------------------------------------------------------
save.image("R_reannot.RData")
#