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
# Import package -----------------------------------------------------------
pkg <- c("ggplot2","dplyr","tidyr","cowplot","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","svglite","stringr","elementalist","gplots","grid","ggh4x")
lapply(pkg, require, character.only = TRUE)
#
# Set directory, create file result -----------------------------------------------------------
if (inputmode == TRUE) {
  current <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current)
  output <- paste(current,"../result/Lagoon", sep = "/")
}
if (inputmode == FALSE) {
  output <- paste("result/Lagoon", sep = "/")
}
# Set working directory
if (dir.exists(output) == FALSE) {dir.create(output)}
if (dir.exists(output) == TRUE) { setwd(output) }
#
# Theme unique Dark perso -------------------------------------------------------
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
#
# Lists inputS ------------------------------------------------------------
input_list <- list.dirs(".", recursive = FALSE)

# Graph_input -----------------------------------------------------------------
list_CC_byparam <- vector()
for (input in input_list) { 
  temp_file <- str_split(input,pattern = "/")[[1]][2]
  temp_tab <- read.csv(paste0(input,"/ssn_composantes_connexes_graph_",temp_file,".csv"),sep = ";")
  temp_tab$CC <- temp_tab$CC+1
  list_CC_byparam <- c(list_CC_byparam,paste0("input_df_",temp_file))
  assign(paste0("input_df_",temp_file),temp_tab)
}
rm(temp_tab)
#
# Processing Homogeneity scores -------------------------------------------
table_list_Total <- data.frame()
#for (annot_type in c("ko")) { print(annot_type)
for (annot_type in c("ko","TRT_1")) { print(annot_type)
  #input <- "input_df_80_80_1e-50"
  table_list_byparam <- as.data.frame(list_CC_byparam)
  table_list_byparam$annot_type <- annot_type
  table_list_byparam$parameters <- sapply(str_split(as.character(table_list_byparam$list_CC_byparam),"df_"), function(x) x[[2]])
  row.names(table_list_byparam) <- paste0(table_list_byparam$list_CC_byparam,"_",annot_type)
  for (input in list_CC_byparam) { print(input)
    CC_list_temp <- get(input) %>% mutate(Nprot = 1) %>% select(CC,Nprot) %>% group_by(CC) %>% summarise_all(sum)
    CC_list <- merge(get(input) %>% select(CC, all_of(annot_type)) %>% group_by(CC) %>% summarize("{annot_type}":=paste(unique(get(annot_type)),collapse=",")),CC_list_temp,by="CC")
    CC_list[,annot_type] <- sapply(str_remove(as.character(CC_list[,annot_type]),"^,"), function(x) x)
    CC_list[,annot_type] <- sapply(str_remove(as.character(CC_list[,annot_type]),",$"), function(x) x)
    CC_list[,annot_type] <- sapply(str_replace_all(as.character(CC_list[,annot_type]),",,",","), function(x) x)
    CC_list[,annot_type] <- sapply(str_split(CC_list[,annot_type],","), function(x) paste(unique(x),collapse = ","))
    CC_list[CC_list==""] <- NA
    CC_list <- CC_list %>% mutate(Nannot=str_count(get(annot_type),",")+1)
    CC_list <- CC_list %>% mutate(Fhom=ifelse(Nannot>1,1-(Nannot/Nprot),1))
    # append stat
    #table_list_byparam[paste0(input,"_",annot_type),"resume"] <- sum((1-CC_list[,"Fhom"][is.na(CC_list[,"Fhom"])==FALSE])^2)/length(CC_list[,"Fhom"][is.na(CC_list[,"Fhom"])==FALSE])
    table_list_byparam[paste0(input,"_",annot_type),"resume"] <- mean(CC_list[,"Fhom"],na.rm =TRUE)/length(CC_list[,"Fhom"][is.na(CC_list[,"Fhom"])==FALSE])
    table_list_byparam[paste0(input,"_",annot_type),"mean Fhom"] <- mean(CC_list[,"Fhom"],na.rm =TRUE)
    table_list_byparam[paste0(input,"_",annot_type),"mean Nprot"] <- mean(CC_list$Nprot,na.rm =TRUE)
    table_list_byparam[paste0(input,"_",annot_type),"mean Nannot"] <- mean(CC_list$Nannot,na.rm =TRUE)
    table_list_byparam[paste0(input,"_",annot_type),"nclust"] <- nrow(CC_list)
    table_list_byparam[paste0(input,"_",annot_type),"#NA"] <- nrow(CC_list %>% filter(is.na(Fhom)==TRUE))*100/nrow(CC_list)
  }
  table_list_Total <- rbind(table_list_Total,table_list_byparam)
}
#
# graph stat --------------------------------------------------------------
input_data_ggplot <- reshape2::melt(table_list_Total,id.vars = c("list_CC_byparam","annot_type","parameters","mean Nprot","nclust")) #%>% select(list_CC_byparam,annot_type,parameters,`mean Fhom`))
order <- c("80_65_1e-50","80_70_1e-50","80_75_1e-50","80_80_1e-50","80_85_1e-50","80_90_1e-50","80_95_1e-50","80_100_1e-50")
input_data_ggplot[,"value"][input_data_ggplot[,"annot_type"]=="TRT_1"&input_data_ggplot[,"variable"]=="mean Fhom"] <- NA
input_data_ggplot[,"value"][input_data_ggplot[,"annot_type"]=="TRT_1"&input_data_ggplot[,"variable"]=="resume"] <- NA
input_data_ggplot$variable <- factor(input_data_ggplot$variable,
                                     levels=c("mean Nannot","mean Fhom","resume","#NA"),
                                     labels=c('bold(bar( "X" )[Nannot]/CC)','bold(bar( "X" )[Fhom])','bold(bar( "X" )[Fhom]/N[CCa])','"NA %"'))
#scale
df_scales <- data.frame(
  Panel = c("mean Nannot", "mean Fhom","resume", "#NA"),
  ymin = c(1, 0.875, 1e-6, 25),
  trans =c("identity","identity","identity","identity"),
  ymax = c(1.5, 0.890, 1.5e-5, 65),
  n = c(4, 4, 4, 4)
)
df_scales$Panel <- factor(df_scales$Panel,
                          levels=c("mean Nannot","mean Fhom","resume","#NA"),
                          labels=c('bold(bar( "X" )[Nannot]/CC)','bold(bar( "X" )[Fhom])','bold(bar( "X" )[Fhom]/N[CCa])','"NA %"'))
df_scales <- split(df_scales, df_scales$Panel)
scales <- lapply(df_scales, function(x) {
  scale_y_continuous(limits = c(x$ymin, x$ymax), trans = x$trans,n.breaks = x$n)
})
#plot
q<-ggplot(input_data_ggplot, aes(x = factor(parameters,levels = order) , y = value)) + 
  geom_line(aes(group = annot_type, color = annot_type)) + 
  geom_text(aes(label=paste("max:",round(value,3))),
            data=input_data_ggplot[which.max(input_data_ggplot[,"value"][input_data_ggplot[,"variable"]=='bold(bar( "X" )[Fhom])'&input_data_ggplot[,"annot_type"]=='ko'])+16,], 
            hjust=-0.15, vjust=2,color="#878787",size=3) +
  geom_text(aes(label=paste("min:",signif(value,3))),
            data=input_data_ggplot[which.min(input_data_ggplot[,"value"][input_data_ggplot[,"variable"]=='bold(bar( "X" )[Fhom]/N[CCa])'&input_data_ggplot[,"annot_type"]=='ko']),], 
            hjust=-0.15, vjust=2,color="#878787",size=3) +
  scale_color_manual(values=c("#e53368","#90ee90","#c1dbff"),labels=c("ko"="Function (KO id)","TRT_1"="Trophic mode")) +
  scale_x_discrete(labels=c("80_65_1e-50"="65%","80_70_1e-50"="70%","80_75_1e-50"="75%","80_80_1e-50"="80%",
                            "80_85_1e-50"="85%","80_90_1e-50"="90%","80_95_1e-50"="95%","80_100_1e-50"="100%")) +
  geom_point(aes(color = annot_type)) + facet_grid(variable~.,scales = "free",labeller = label_parsed) + 
  ggh4x::facetted_pos_scales(y = scales) + 
  geom_vline(xintercept= input_data_ggplot[which.max(input_data_ggplot[,"value"][input_data_ggplot[,"variable"]=='bold(bar( "X" )[Fhom])'&input_data_ggplot[,"annot_type"]=='ko']),"parameters"], 
             linetype="dashed", color = "#878787", linewidth=0.4) + 
  geom_vline(xintercept= input_data_ggplot[which.min(input_data_ggplot[,"value"][input_data_ggplot[,"variable"]=='bold(bar( "X" )[Fhom]/N[CCa])'&input_data_ggplot[,"annot_type"]=='ko']),"parameters"], 
             linetype="dashed", color = "#878787", linewidth=0.4) + 
  labs(x="ID%",y="",color="Annotations") + theme_unique_art() 
print(q)

qustom <- q +
  theme(axis.title.x = element_text(margin=margin(50,0,0,0))) +
  coord_cartesian(clip='off')
for (i in 1:length(unique(input_data_ggplot$parameters))) {
  qustom <- qustom + annotation_custom(
    textGrob(label=paste0('N = ',input_data_ggplot$nclust[i],'\n(x̄ = ',round(input_data_ggplot$`mean Nprot`[i],1),")"), rot=0, gp=gpar(fontsize=9,col="#878787")),
    xmin=input_data_ggplot$parameters[i], xmax=input_data_ggplot$parameters[i], ymin=25, ymax=-1
  )
}
sub.div <- 1/length(unique(input_data_ggplot$parameters))
for (i in 1:(length(unique(input_data_ggplot$parameters))-1)) {
  qustom <- qustom + annotation_custom(
    linesGrob(
      x=unit(c(sub.div*i,sub.div*i), 'npc'),
      y=unit(c(0,-0.4), 'npc')   # length of the line below axis
    )
  )
}
svglite("Parameter_analyses.svg",width = 9.50,height = 8.00)
print(qustom)
dev.off()
#
# Save R Image ------------------------------------------------------------
save.image("R_Treshold.RData")

    