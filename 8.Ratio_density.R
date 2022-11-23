library(ggplot2)
library(ggpubr)
library(ggsci)
library(stringr)
library(cowplot)
args=commandArgs(T)
# args <- c("/home/zhl/projects/batch_effect_correction/new_results/Ratio_allspecies_7methods/Protein_Ratio_Combat.csv",
# "/home/zhl/projects/batch_effect_correction/new_results/Ratio_allspecies_7methods/")
df <- read.csv(args[1],row.names = 1,check.names = F)
df_cols = colnames(df)
df_cols <- gsub("swissprot_homo_201604","Human",df_cols)
df_cols <- gsub("Saccharomyces cerevisiae_201602","Yeast",df_cols)
df_cols <- gsub("uniprot-Ecoli_K12_201709","Ecoli",df_cols)
colnames(df) <- df_cols
df_t <- reshape2::melt(df,measure.vars = colnames(df)[!grepl("median",colnames(df))],variable.name = "Species",value.name="Ratio")
df_t <- reshape2::melt(df_t,measure.vars = colnames(df)[grepl("median",colnames(df))],variable.name = "median_name",value.name="median")
df_t <- df_t[abs(df_t$Ratio)<=4,]
df_t <- na.omit(df_t)
if (grepl("Peptide_Ratio_",args[1])){
    title = str_split_fixed(args[1],"Peptide_Ratio_",2)[,2]
    title <- str_split_fixed(title,".csv",2)[,1]
} else if (grepl("Protein_Ratio_",args[1])) {
    title = str_split_fixed(args[1],"Protein_Ratio_",2)[,2]
    title <- str_split_fixed(title,".csv",2)[,1]
} else {title=""}
df_t <- df_t[order(df_t$Species,decreasing = T),]
df_t$Species <- factor(df_t$Species,levels = c("Ecoli","Human","Yeast"))
plot=ggscatter(df_t,x="median",y="Ratio",
               color="Species",font.label=c(10,"plain","black"),
               palette = 'Dark2',size=1.5, alpha=0.75, legend = "bottom",
                xlab = 'Log2 Intensity', ylab = 'Log2 (A:B)')+ggtitle(title)+
    geom_hline(yintercept = c(-1,0,2),lty=2,col="black")+
    theme(plot.title = element_text(hjust=0,vjust = -1),
        axis.text=element_text(color='black',size = 18,face = "bold"),
        axis.title=element_text(size=20,color='black'),
        panel.grid=element_blank(),
        legend.key.size = unit(12, "pt"),
        panel.border=element_blank(),
        panel.grid.major=element_line(colour=NA),
        panel.background=element_rect(fill="white",colour="black"))+
        scale_y_continuous(limits = c(-4,4),breaks = seq(-4,4,1))

yplot <- ggdensity(df_t, x="Ratio", alpha = 0.8,
          color = "Species", fill = "Species",
          palette = 'Dark2')+scale_x_continuous(limits = c(-4,4),breaks = seq(-4,4,1) )+rotate()+clean_theme()

xplot <- gghistogram(df_t,x="median",
            bins=40,
            fill="median_name",position = "stack",size=0.1,
            palette = 'Dark2',alpha = 0.8)+clean_theme()
p1=insert_yaxis_grob(plot,grob = yplot,position="right")
p2=insert_xaxis_grob(p1,grob = xplot,position="top")
p=cowplot::ggdraw(p2)
ggsave(plot=p,file=args[2],width=8.5,height = 7.5)
ggsave(plot=p,file=args[3],width=8.5,height = 7.5)
