setwd("C:/Sauvegarde/fmollandin (JJ-1313-SRVDATAuser)/Résultats")

library(ggplot2)
library(tidyverse)
library(reshape2)
library(pROC)
library(xlsx)
library(knitr)
library(circlize)
library(plotly)
library(wesanderson)
library(ComplexHeatmap)
library(ggnewscale)
library(gridExtra)
library(cowplot)
source("functions_aux.R")

## Read in Data -------------------------------------------------------------------------
scenarios=read.xlsx("Copie de Sim_Param.xlsx",sheetIndex = 1,header=F,encoding = "UTF-8")
scenarios=scenarios[-22,]
rownames(scenarios)=scenarios[,1]
scenarios=scenarios[,-1]
scenarios=scenarios[,17:29]
colnames(scenarios)=c("R","S","T","U","V","W","X","Y","Z","AA","AB","AC","Ref")
liste_names=colnames(scenarios)
nscen=length(liste_names)
liste_names_wo=paste(liste_names,"wo",sep="")
df_pi3<-data.frame(Scenario=colnames(scenarios),pi3=t(scenarios["pi3",]))
df_var<-data.frame(Scenario=colnames(scenarios),Var=t(scenarios["pi3",])*80)
colnames(df_var)=c("Scenario","Var")

## Figure 1 -----------------------------------------------------------------------------

## Notes:
## => Not sure that (A) should be vertical (harder to compare between the two panels)
## => I inversed the order of h2 labels in (A) and plots in (B) to match the order of plot (A)
## => Colors of heritability now sequential to reflect that this is a quantitative value
## => Made scales equal for y-axis in (B)
## => Should we express part of genetic varaicne as per mille (instead of per cent)? Or sum over 5 QTLs?
## => Diagonal labels are really hard to read, so I made them perpendicular

# (A) Corrélations
load("df_cor_all.RData")
load("df_cor_mean.RData")
levels(df_cor_mean$Data) <- c("50k custom", "50k")
df_cor_mean$Data <- relevel(df_cor_mean$Data, ref = "50k")
df_cor_mean$h2 <- factor(df_cor_mean$h2)
levels(df_cor_mean$h2) <- paste("h\u00b2=", c(0.1, 0.3,0.5,0.8))
fig1a <- ggplot(df_cor_mean,
                aes(x=pi3,y=Correlation,color=as.factor(h2)))+
    geom_line(size=1.25)+
    geom_point()+
    geom_errorbar(data=df_cor_mean,
                  aes(x=pi3,ymin=Correlation-SD, ymax=Correlation+SD),
                  position=position_dodge(0.00001),
                  alpha=0.5,size=1)+
    facet_wrap(~Data, ncol = 1)+
    labs(color="Heritability")+
    guides(fill=F)+
    ## Reverse legend order to match plots,
    scale_color_brewer(type="seq",
                       palette=1,
                       limits=c(1:2, levels(df_cor_mean$h2)),
                       breaks=rev(levels(df_cor_mean$h2)))+
    xlab("% genetic variance per QTL")+
    ylab("Correlation (\u03c1)") +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom",
          legend.direction = "vertical"); fig1a

load("df_cor_all.RData")
load("df_cor_mean.RData")
levels(df_cor_mean$Data) <- c("50k custom", "50k")
df_cor_mean$h2 <- factor(df_cor_mean$h2)
levels(df_cor_mean$h2) <- paste("h\u00b2=", c(0.1, 0.3,0.5,0.8))

fig1a_v2 <- ggplot(df_cor_mean,
                aes(x=pi3,y=Correlation,linetype=Data,color=as.factor(h2),group=interaction(as.factor(h2),Data)))+
    geom_line(size=1.25)+
    geom_point()+
    labs(color="Heritability")+
    guides(fill=F)+
    ## Reverse legend order to match plots,
    scale_color_brewer(type="seq",
                       palette=1,
                       limits=c(1:2, levels(df_cor_mean$h2)),
                       breaks=rev(levels(df_cor_mean$h2)))+
    xlab("% genetic variance per QTL")+
    ylab("Correlation (\u03c1)") +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    theme_bw(base_size = 18) +
    theme(legend.direction = "vertical"); fig1a_v2

ggsave(fig1a_v2, filename="Figure_1A_v2.pdf",
       width = 12, height = 9,
       device = cairo_pdf)

# (B) Différence de corrélations
load("df_rapport_all.RData")
df_rapport_all$h2 <- factor(df_rapport_all$h2, levels = c(0.1, 0.3,0.5,0.8),
                            ordered=TRUE,
                            labels = c(expression(paste(h^2, "= 0.1")),
                                       expression(paste(h^2, "= 0.3")),
                                       expression(paste(h^2, "= 0.5")),
                                       expression(paste(h^2, "= 0.8"))))
levels(df_rapport_all$h2) <- rev(levels(df_rapport_all$h2))
fig1b <- ggplot(df_rapport_all,aes(x=pi3,y=Rapport, group = pi3))+
    geom_boxplot()+
    xlab("% genetic variance per QTL")+
    ylab("(50k custom \u03c1) - (50k \u03c1)")+
    # scale_x_discrete(labels=paste(sort(df_pi3$pi3),sep=""))+
    facet_wrap(~h2,
               #scales = "free",
               labeller = label_parsed, ncol = 1)+
    geom_jitter(height = 0, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = 2, col = "red")+
    theme_bw(base_size = 18); fig1b


# Full figure combined using cowplot and exported using ggsave
fig1 <- plot_grid(fig1a, fig1b, labels = c("A", "B"),
          rel_widths = c(1, 1.75), label_size = 20)
ggsave(fig1, filename="Figure_1.pdf",
       width = 16, height = 10,
       device = cairo_pdf)


## Figure 1 alternative 1-----------------------------------------------------------------------------

# (A) Corrélations
load("df_cor_all.RData")
load("df_cor_mean.RData")
levels(df_cor_mean$Data) <- c("50k custom", "50k")
df_cor_mean$Data <- relevel(df_cor_mean$Data, ref = "50k")
df_cor_mean$h2 <- factor(df_cor_mean$h2)
levels(df_cor_mean$h2) <- paste("h\u00b2=", c(0.1, 0.3,0.5,0.8))
fig1a <- ggplot(df_cor_mean,
                aes(x=pi3,y=Correlation,color=as.factor(h2)))+
    geom_line(size=1.25)+
    geom_point()+
    geom_errorbar(data=df_cor_mean,
                  aes(x=pi3,ymin=Correlation-SD, ymax=Correlation+SD),
                  position=position_dodge(0.00001),
                  alpha=0.5,size=1)+
    facet_wrap(~Data, ncol = 2)+
    labs(color="Heritability")+
    guides(fill=F)+
    ## Reverse legend order to match plots,
    scale_color_brewer(type="seq",
                       palette=1,
                       limits=c(1:2, levels(df_cor_mean$h2)),
                       breaks=rev(levels(df_cor_mean$h2)))+
    xlab("% genetic variance per QTL")+
    ylab("Correlation (\u03c1)") +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    theme_bw(base_size = 18); fig1a

# (B) Différence de corrélations
load("df_rapport_all.RData")
df_rapport_all$h2 <- factor(df_rapport_all$h2, levels = c(0.1, 0.3,0.5,0.8),
                            ordered=TRUE,
                            labels = c(expression(paste(h^2, "= 0.8")),
                                       expression(paste(h^2, "= 0.5")),
                                       expression(paste(h^2, "= 0.3")),
                                       expression(paste(h^2, "= 0.1"))))
levels(df_rapport_all$h2) <- rev(levels(df_rapport_all$h2))
fig1b <- ggplot(df_rapport_all,aes(x=pi3,y=Rapport,group=pi3))+
    geom_boxplot()+
    xlab("% genetic variance per QTL")+
    ylab("(50k custom \u03c1) - (50k \u03c1)")+
#    scale_x_discrete(labels=paste(sort(df_pi3$pi3),sep=""))+
    facet_wrap(~h2,
               #scales = "free",
               labeller = label_parsed, ncol = 4)+
    geom_jitter(height = 0, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = 2, col = "red")+
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle=90)); fig1b


# Full figure combined using cowplot and exported using ggsave
fig1 <- plot_grid(fig1a, fig1b, labels = c("A", "B"), ncol = 1,
                   label_size = 20)
ggsave(fig1, filename="Figure_1_alt1.pdf",
       width = 16, height = 10,
       device = cairo_pdf)

## Figure 1 alternative 2 -----------------------------------------------------------------------------

# (A) Corrélations
load("df_cor_all.RData")
load("df_cor_mean.RData")
levels(df_cor_mean$Data) <- c("50k custom", "50k")
df_cor_mean$Data <- relevel(df_cor_mean$Data, ref = "50k")
df_cor_mean$h2 <- factor(df_cor_mean$h2)
levels(df_cor_mean$h2) <- paste("h\u00b2=", c(0.1, 0.3,0.5,0.8))
fig1a <- ggplot(df_cor_mean,
                aes(x=pi3,y=Correlation,color=as.factor(h2)))+
    geom_line(size=1.25)+
    geom_point()+
    geom_errorbar(data=df_cor_mean,
                  aes(x=pi3,ymin=Correlation-SD, ymax=Correlation+SD),
                  position=position_dodge(0.00001),
                  alpha=0.5,size=1)+
    facet_wrap(~Data, ncol = 1)+
    labs(color="Heritability")+
    guides(fill=F)+
    ## Reverse legend order to match plots,
    scale_color_brewer(type="seq",
                       palette=1,
                       limits=c(1:2, levels(df_cor_mean$h2)),
                       breaks=rev(levels(df_cor_mean$h2)))+
    xlab("% genetic variance per QTL")+
    ylab("Correlation (\u03c1)") +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom",
          legend.direction = "vertical"); fig1a

# (B) Différence de corrélations
load("df_rapport_all.RData")
df_rapport_all$h2 <- factor(df_rapport_all$h2, levels = c(0.1, 0.3,0.5,0.8),
                            ordered=TRUE,
                            labels = c(expression(paste(h^2, "= 0.8")),
                                       expression(paste(h^2, "= 0.5")),
                                       expression(paste(h^2, "= 0.3")),
                                       expression(paste(h^2, "= 0.1"))))
levels(df_rapport_all$h2) <- rev(levels(df_rapport_all$h2))
fig1b <- ggplot(df_rapport_all,aes(x=factor(pi3),y=Rapport))+
    geom_boxplot()+
    xlab("% genetic variance per QTL")+
    ylab("(50k custom \u03c1) - (50k \u03c1)")+
    scale_x_discrete(labels=paste(sort(df_pi3$pi3),sep=""))+
    facet_wrap(~h2,
               #scales = "free",
               labeller = label_parsed, ncol = 1)+
    geom_jitter(height = 0, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = 2, col = "red")+
    theme_bw(base_size = 18); fig1b

ggsave(fig1b, filename="Figure_1B.pdf",
       width = 10, height = 9,
       device = cairo_pdf)

# Full figure combined using cowplot and exported using ggsave
fig1 <- plot_grid(fig1a, fig1b, labels = c("A", "B"),
                  rel_widths = c(1, 1.75), label_size = 20)
ggsave(fig1, filename="Figure_1_alt2.pdf",
       width = 16, height = 10,
       device = cairo_pdf)


## Figure 1 alternative 3-----------------------------------------------------------------------------

# (A) Corrélations
load("df_cor_all.RData")
load("df_cor_mean.RData")
levels(df_cor_mean$Data) <- c("50k custom", "50k")
df_cor_mean$Data <- relevel(df_cor_mean$Data, ref = "50k")
df_cor_mean$h2 <- factor(df_cor_mean$h2)
levels(df_cor_mean$h2) <- paste("h\u00b2=", c(0.1, 0.3,0.5,0.8))
fig1a <- ggplot(df_cor_mean,
                aes(x=pi3,y=Correlation,color=as.factor(h2)))+
    geom_line(size=1.25)+
    geom_point()+
    geom_errorbar(data=df_cor_mean,
                  aes(x=pi3,ymin=Correlation-SD, ymax=Correlation+SD),
                  position=position_dodge(0.00001),
                  alpha=0.5,size=1)+
    facet_wrap(~Data, ncol = 2)+
    labs(color="Heritability")+
    guides(fill=F)+
    ## Reverse legend order to match plots,
    scale_color_brewer(type="seq",
                       palette=1,
                       limits=c(1:2, levels(df_cor_mean$h2)),
                       breaks=rev(levels(df_cor_mean$h2)))+
    xlab("% genetic variance per QTL")+
    ylab("Correlation (\u03c1)") +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    theme_bw(base_size = 18); fig1a

# (B) Différence de corrélations
load("df_rapport_all.RData")
df_rapport_all$h2 <- factor(df_rapport_all$h2, levels = c(0.1, 0.3,0.5,0.8),
                            ordered=TRUE,
                            labels = c(expression(paste(h^2, "= 0.8")),
                                       expression(paste(h^2, "= 0.5")),
                                       expression(paste(h^2, "= 0.3")),
                                       expression(paste(h^2, "= 0.1"))))
levels(df_rapport_all$h2) <- rev(levels(df_rapport_all$h2))
fig1b <- ggplot(df_rapport_all,aes(x=factor(pi3),y=Rapport))+
    geom_boxplot()+
    xlab("% genetic variance per QTL")+
    ylab("(50k custom \u03c1) - (50k \u03c1)")+
    scale_x_discrete(labels=paste(sort(df_pi3$pi3),sep=""))+
    facet_wrap(~h2,
               #scales = "free",
               labeller = label_parsed, ncol = 4)+
    geom_jitter(height = 0, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = 2, col = "red")+
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle=90)); fig1b


# Full figure combined using cowplot and exported using ggsave
fig1 <- plot_grid(fig1a, fig1b, labels = c("A", "B"), ncol = 1,
                  label_size = 20)
ggsave(fig1, filename="Figure_1_alt3.pdf",
       width = 16, height = 10,
       device = cairo_pdf)


## Figure 2 -----------------------------------------------------------------------------

## Notes:
## => I left the heritability in decreasing order to match Figure 1

# PIP
load("df_PIP_all.RData")
levels(df_PIP_all$Data) <- c("50k custom", "50k")
df_PIP_all$Data <- relevel(df_PIP_all$Data, ref = "50k")
levels(df_PIP_all$h2) <- paste("h\u00b2=", c(0.8, 0.5,0.3,0.1))
levels(df_PIP_all$gpin) <- c("Null",
                             "Small",
                             "Medium",
                             "Large")
wescol <- wes_palette("Zissou1",5,type="discrete")
fig2 <- ggplot(df_PIP_all,aes(x=as.factor(pi3),y=Moyenne,fill=gpin))+
    geom_bar(stat="Identity",width = 1)+
    xlab("% genetic variance per QTL")+
    scale_fill_manual(values=c("grey",wescol[-c(1,4)])) +
    ylab("Mean inclusion probability")+
    facet_grid(Data~as.factor(h2))+
    # theme(axis.text.x = element_text(face="bold",  size=4, angle=45))+
    labs(fill="SNP effect\nclass")+
    # scale_x_discrete(labels=paste(sort(df_pi3$pi3)*100,"%",sep="")) +
    theme_minimal(base_size = 18) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))

ggsave(fig2, filename="Figure_2.pdf",
       width = 12, height = 8, device = cairo_pdf)

## Figure 3 -----------------------------------------------------------------------------

## Note:
## => Made these into stacked barplots, matched colors for groups to Figure 2,
##  heritability still in same order as other plots

# Effectifs MAP
load("MAP_all.RData")
load("QTL_all.RData")
QTL_all <- QTL_all %>% filter(h2 != 0.1, Classe != "Low") %>% droplevels()
levels(QTL_all$Data) <- c("50kcustom", "50k")
QTL_all$Data <- relevel(QTL_all$Data, ref = "50k")
QTL_all$Data <- factor(QTL_all$Data, labels = c(expression(paste("50k")),
                                                expression(paste("50k custom"))))
QTL_all$h2 <- factor(QTL_all$h2, levels = c(0.8, 0.5,0.3,0.1),
                            ordered=TRUE,
                            labels = c(expression(paste(h^2, "= 0.8")),
                                       expression(paste(h^2, "= 0.5")),
                                       expression(paste(h^2, "= 0.3")),
                                       expression(paste(h^2, "= 0.1"))))
levels(QTL_all$Classe) <- c("Large", "Medium")
QTL_all$Classe <- factor(QTL_all$Classe, levels = c("Large", "Medium"))
wescol <- wes_palette("Zissou1",5,type="discrete")
fig3 <- ggplot(QTL_all)+
    # geom_line(data=QTL_all,aes(y=MAP,x=pi3,color=Classe),alpha=1,size=1)+
    geom_col(data=QTL_all,
             aes(y=MAP,x=factor(pi3),fill=Classe), position = "stack")+
    facet_grid(h2~Data, labeller = label_parsed)+
    scale_fill_manual(values=wescol[c(5,3)], name= "SNP effect class") +
    xlab("% genetic variance per QTL") +
    ylab("Number of assigned SNPs (MAP)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
    guides(fill = guide_legend(reverse = TRUE))

ggsave(fig3, filename="Figure_3.pdf",
       width = 12, height = 9, device = cairo_pdf)

## Figure 4 -----------------------------------------------------------------------------
load("matrix_QTL.RData")
load("matrix_QTLwo.RData")
load("row_annotation_matrix.RData")

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
col_fun2 = colorRamp2(c(0, 1), c("white", "red"))
tmp <- as.character(as.numeric(substr(colnames(matrix_QTL), 1, nchar(colnames(matrix_QTL))-1)) / 100)
colnames(matrix_QTL) <- colnames(matrix_QTLwo) <- tmp

ht3 <- Heatmap(row_annotation_matrix,
               cluster_columns = FALSE,
               col = col_fun2,
               cluster_rows = TRUE,
               row_title = NULL,
               name = "D'")
ht2 <- Heatmap(matrix_QTLwo,
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               col=col_fun,
               show_heatmap_legend = FALSE,
               column_names_gp = grid::gpar(fontsize = 14),
               column_title = "50k",
               show_row_names = FALSE)
ht1 <- Heatmap(matrix_QTL,
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               col=col_fun,
               column_names_gp = grid::gpar(fontsize = 14),
               name="Inclusion\nprobability",
               column_title = "50k custom",
               show_row_names = FALSE)

pdf("Figure_4.pdf", width = 8, height = 6)
draw(ht3 + ht2 + ht1, column_title = "% genetic variance per QTL",
     column_title_side = "bottom" )
dev.off()




## Figure 5 -----------------------------------------------------------------------------

## Notes:
## => Match category colors to Figures 2 and 3 (add a lighter shade of the red for the QTL window)
## => Big plot with lots of overplotted points -- can we omit points with estimated vi < some threshold?

# Exemple variance (0.1)
load("param_app.RData")
levels(param_app$Data) <- c("50k custom", "50k")
param_app$Data <- relevel(param_app$Data, ref = "50k")
levels(param_app$True_Cat) <- c("Null", "Polygenic", "QTL window", "QTL")
wescol <- wes_palette("Zissou1",5,type="discrete")
# https://www.color-hex.com/color/f21a00
fig5 <- ggplot(param_app,aes(x=Position,y=Part_Vi,color=True_Cat))+
    geom_point()+
    geom_point(data = subset(param_app, True_Cat == 'QTL window'),
               aes(x=Position,y=Part_Vi,color=True_Cat), size = 3)+
    geom_point(data = subset(param_app, True_Cat == 'QTL'),
               aes(x=Position,y=Part_Vi,color=True_Cat), size = 3)+
    facet_grid(Data ~ as.factor(k))+
    scale_color_manual(values=c("grey", wescol[2], "#f9a399", wescol[5]),
                       name = "Simulated\neffect class")+
    xlab("SNP index") +
    ylab("Estimated proportion of genetic variance") +
    theme_bw(base_size = 18)

ggsave(fig5, filename="Figure_5.pdf",
       width = 14, height = 8, device = cairo_pdf)

## Figure 6 -----------------------------------------------------------------------------

## Note:
##  => heritability still in same order as other plots
##  => categories are different here than in previous plots, perhaps use a different palette here?
##  => Need to come up with better names for the criteria names

# Comparaison méthodes
load("df_comp_QTL.RData")
levels(df_comp_QTL$Data) <- c("50k custom", "50k")
df_comp_QTL$Data <- relevel(df_comp_QTL$Data, ref = "50k")
df_comp_QTL$Data <- factor(df_comp_QTL$Data, labels = c(expression(paste("50k")),
                                                expression(paste("50k custom"))))
df_comp_QTL$h2 <- factor(df_comp_QTL$h2, levels = c(0.1, 0.3,0.5,0.8),
                            ordered=TRUE,
                            labels = c(expression(paste(h^2, "= 0.8")),
                                       expression(paste(h^2, "= 0.5")),
                                       expression(paste(h^2, "= 0.3")),
                                       expression(paste(h^2, "= 0.1"))))
levels(df_comp_QTL$h2) <- rev(levels(df_comp_QTL$h2))
levels(df_comp_QTL$Method) <- c("CIPi", "CIPi + Vi", "CIPi + Vi + Non-null MAPi")
fig6 <- ggplot(data=df_comp_QTL, aes(x=as.factor(pi3), y=Detected, fill=Method)) +
    geom_bar(stat="identity",width=1) +
    facet_grid(Data~as.factor(h2),labeller = label_parsed)+
    scale_fill_manual(values=wes_palette("Zissou1",10,type="continuous")[c(1,4,7)],
                      name = "Detection criterion")+
    ylab("Number of detected QTL")+
    xlab("% genetic variance per QTL")+
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90),
          legend.direction = "vertical")

ggsave(fig6, filename="Figure_6.pdf",
       width = 12, height = 8, device = cairo_pdf)

# Eventual supplementary plot? -----------------------------------------------------------
## Vi
# load("df_Vi_all.RData")
# ggplot(df_Vi_all,aes(x=pi3,y=Vi,color=h2))+
#     facet_grid(~Data)+
#     geom_jitter(size=0.8)+
#     geom_crossbar(aes(ymin = Vi, ymax = Vi))




