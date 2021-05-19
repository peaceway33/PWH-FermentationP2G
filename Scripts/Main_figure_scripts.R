#### Figure 1  ####
phenotypes_combined_subset <- #Supplementary data Figure 1#

## PLOT ALL TRAITS
theme_custom <- theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.placement = "outside",
        strip.text.x = element_text(angle = 0, size = 18),
        strip.text.y.left = element_text(angle = 0, size = 18, hjust = 1),
        axis.text = element_text(size = 10),
        legend.position = "top"
  )
clean_labels <- clean_names %>% 
  mutate(fieldnames_R = as.character(fieldnames_R)) %>%
  filter(fieldnames_R %in% colnames(phenotypes_combined_subset %>% select(-segregant))) 
clean_labels <- clean_labels[match(phenotypes_combined_subset %>% select(-segregant) %>% colnames(), clean_labels$fieldnames_R),] %>%
  mutate(cleanname = gsub(" \\(","\n(",as.character(cleanname))) %>%
  mutate(suffix = ifelse(str_detect(cleanname, ".*(?<!\\))$"), "\n","")) %>% 
  mutate(cleanname = gsub(" tolerance"," tol.", cleanname)) %>%
  mutate(cleanname = paste0(cleanname,suffix)) %>%
  pull(cleanname)
num_columns <- (length(phenotypes_combined_subset) - 1)

scale_manual_nodiff <- function(g){
  for(i in 1:g$nrow){
    for(j in 1:g$ncol){
      if(i == j){
        g[i,j] <- g[i,j] + 
          scale_fill_manual(values = "#ff99c8") + 
          scale_color_manual(values = "#ff99c8") + 
          scale_alpha_manual(values = 0.4)
      }
    }
  }
  return(g)
}
my_diag_nodiff <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
  geom_histogram(..., position = "identity", 
                   color = "#ff99c8", fill = "#ff99c8", alpha = 0.5)
}
g_combined_subset_nodiff <- ggpairs(phenotypes_combined_subset %>% select(-segregant),
                                    columns = 1:(length(phenotypes_combined_subset)-1),
                                    lower = list(continuous = wrap("points", 
                                    color = "#AAAAAA", fill = "#AAAAAA", alpha = 0.3)),
                                    diag = list(continuous = my_diag_nodiff), 
                                    columnLabels = clean_labels,
                                    switch = "both")  +
  theme_custom


#### Figure 2A ####
variant_freq <- #Supplementary data Figure 2A#

variant_freq_distinct <- distinct(variant_freq, Variant, .keep_all = TRUE)

variant_freq_distinct %>%
  group_by(annotation_fin) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

Figure_2A <-
  ggplot(variant_freq_distinct %>% filter(!Class %in% c("Other")& 
                                            Peak_width==1))+ 
  geom_bar( aes(x = forcats::fct_rev(fct_infreq(Condition_new)), fill=Class),
            width = 0.8,
            color = "white",
            alpha = 0.8) +
  theme_cowplot()+
  theme(
    legend.position = c(0.025,0.85),
    legend.title = element_blank(),
    plot.title=element_text(hjust=0.5),
    axis.line =  element_line(size = 0.1),
    axis.ticks = element_line(size = 0.5),
    axis.text.x = element_text(angle=45,hjust = 1, vjust=1),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", size=0.5),
    strip.background = element_blank(),
    axis.title= element_text(face = "bold"))+
  scale_fill_brewer(direction = -1, palette = "Set2")+
  ylab("# of QTN")

#### Figure 2B ####
variant_compare_all <- #Supplementary data Figure 2B#

Figure_2B <- variant_compare_all %>%
  filter(!Class %in% c("Other")) %>%
  group_by(Group, Class) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(aes(x= Class, y= freq, fill= fct_rev(Group))) +
  geom_col(color= "black", width=0.6, alpha=0.8, position = position_dodge(width=0.75))+
  theme_cowplot()+
  theme(
    legend.position = c(0.05,0.9),
    legend.title=element_blank(),
    axis.line =  element_line(size = 0.1),
    axis.ticks = element_line(size = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45,hjust = 1, vjust=1),
    axis.title= element_text(face="bold"),
    panel.border = element_rect(colour = "black", size=0.5),
    strip.background = element_blank())+
  scale_fill_manual(values = c("grey90", "grey45", "black"))+
  ylab("Fraction of QTN")+
  guides(fill = guide_legend(reverse=TRUE))

#### Figure 2C ####
variant_freq <- #Supplementary data Figure 2C#

variants_fin <- variant_freq %>%
  drop_na()%>%
  mutate(R_minus_Y=R_mean-Y_mean)%>% 
  arrange(Condition_new,R_mean,abs(R_minus_Y)) %>%
  mutate(rank = row_number())

Figure_2C <- ggplot(variants_fin) +
  facet_wrap(Condition_new~., scale="free")+
  geom_segment( aes(x=rank, xend=rank, y=R_mean, yend=Y_mean), color="grey70", size=0.25) +
  geom_hline(yintercept = 0, size=0.25, color="grey70") +
  geom_point( aes(x=rank, y=R_mean, size=as.numeric(freq1)), fill="#0F579A", color="black", shape=21) +
  geom_point( aes(x=rank, y=Y_mean, size=as.numeric(freq2)), fill="#FDD36F", color="black", shape=21) +
  scale_size(range=c(0.5,2))+
  theme_cowplot()+
  theme(
    legend.position = "none",
    legend.box = "vertical",
    axis.text.x = element_blank(),
    axis.line =  element_line(size = 0.5),
    axis.ticks.y = element_line(size = 0.5),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", size=0.5),
    strip.background = element_blank())+
  labs(x="Candidate variant", y="Effect (z-score)")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.25))

#### Figure 2D ####

shared <- #Supplementary data Figure 2D#

ggplot() +
  geom_histogram(data=shared, aes(x=Nr..Strains.carrying.at.least.1.copy.RM, y=stat(count)/sum(stat(count))),binwidth = 10, alpha=0.3, fill="red") +
  geom_histogram(data=shared, aes(x=Nr..strains.arrying.at.least.1.copy.YJM, y=stat(count)/sum(stat(count))), binwidth = 10,, alpha=0.3, fill="blue")

shared.m<-melt(shared, id.vars = c("variantPos.FW."))

facets <- c("Nr..Strains.carrying.at.least.1.copy.RM" =	"RM11-1a" ,
            "Nr..strains.arrying.at.least.1.copy.YJM" =	"YJM975" )

Figure_2D <- ggplot() +
  geom_histogram(data=shared.m, aes(x=value, y=stat(count)/sum(stat(count)), fill=variable),binwidth = 10, alpha=1) +
  theme_bw(base_size = 14) +
  facet_wrap(~variable, ncol = 1, nrow=2, labeller = labeller(variable = facets,label_wrap_gen(multi_line = TRUE) ),scales="free_x") +
  ylab("Proportion of QTLs/QTNs") +
  xlab("Nr. of strains carrying at least 1 allele (binwidth=10)") +
  scale_fill_manual(values=c("#0F579A", "#FDD36F")) +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    axis.line =  element_line(size = 0.5),
    axis.ticks.y = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", size=0.5),
    strip.background = element_blank(), panel.grid = element_blank())

#### Figure 2E ####
top <- #Supplementary data Figure 2E#
top_melt<-melt(top, id.vars = c("Original.order" , "Condition","Condition.simplified",              
                                "R_mean","Condition_new","Group1_medium" , 
                                "Group2_phenotyping", "Group3_sample_type","Group4_final","Y_mean", "Greater_effect_allele",             
                                "combined_effect" , "Rank", "variantPos.FW." ,
                                "classifier_fixed" , "chrom..formal.name.","chromosome..ping.", 
                                "location..ping.", "peak" , "RM_genotype",                       
                                "YJM_genotype"))


lineage.toplot <-top_melt %>% separate(variable, c("key", "lineage"), sep="([.])")
lineage.toplot$key<-factor(lineage.toplot$key)
lineage.toplot$lineage<-factor(lineage.toplot$lineage)

#lineage.m$variable<-as.character(lineage.m$variable)
#lineage.toplot$group<-as.character(lineage.toplot$group)

lineage.toplot$parent <- sapply(strsplit(as.character(lineage.toplot$lineage), "\\_"), tail, 1)
lineage.toplot$key<-factor(lineage.toplot$key)
lineage.toplot$lineage<-factor(lineage.toplot$lineage)
lineage.toplot$parent<-factor(lineage.toplot$parent)


lineage.toplot<-subset(lineage.toplot, parent != "others" & key != "M1" & key != "M2" & key != "M3" & key != "X7")
lineage.toplot$key<-factor(lineage.toplot$key)
lineage.toplot$parent<-factor(lineage.toplot$parent)
new_labels <- c("X1A" =	"Wine European" ,
                "X1B" =	"Wine European subclade1" ,
                "X1C" =	"Wine European subclade2" ,
                "X1D" =	"Wine European subclade3" ,
                "X1E" =	"Wine European subclade4" ,
                "X2" =	"Alpechin" ,
                "X3" =	"Brazilian bioethanol" ,
                "X4" =	"Mediterranean oak" ,
                "X5" =	"French dairy" ,
                "X6" =	"African beer" ,
                "X7" =	"Mosaic beer" ,
                "X8" =	"Mixed origin" ,
                "X9" =	"Mexican agave" ,
                "X10" =	"French Guiana human" ,
                "X11" =	"Ale beer" ,
                "X12" =	"West African cocoa" ,
                "X13" =	"African palm wine" ,
                "X14" =	"CHNIII" ,
                "X15" =	"CHNII" ,
                "X16" =	"CHNI" ,
                "X17" =	"Taiwanese" ,
                "X18" =	"Far East Asia" ,
                "X19" =	"Malaysian" ,
                "X20" =	"CHNV" ,
                "X21" =	"Ecuadorean" ,
                "X22" =	"Far East Russian" ,
                "X23" =	"North American oak" ,
                "X24" =	"Asian islands" ,
                "X25" =	"Sake" ,
                "X26" =	"Asian fermentation" ,
                "M1" =	"Mosaic region1" ,
                "M2" =	"Mosaic region2" ,
                "M3" =	"Mosaic region3")
lineage.toplot$key<-factor(lineage.toplot$key, levels=c("X1A", "X1B", "X1C", "X1D", "X2", "X1E",
                                                        "X3",  "X4" , "X6" , "X5",  "X8" , "X9",
                                                        "X10", "X11" ,"X12", "X13", "X14", "X15", "X16", "X17" ,"X18", "X19",
                                                        "X20" ,"X21" ,"X22", "X23", "X24","X26", "X25" ))

lineage.toplot$Condition<-factor(lineage.toplot$Condition, levels=c("1-propanol (ppm)", "ethyl acetate (ppm)", "ethyl hexanoate (ppm)", "ethyl octanoate (ppm)", "isoamyl alcohol (ppm)", "isobutanol (ppm)", "isopentyl acetate (ppm)", "phenethyl acetate (ppm)", "maltose (2% w/v)", "abv (% v/v)", "acetic acid (g/l)", "glycerol (g/l)", "SO2 (mg/l)", "acetic acid tolerance (50nM)", "ethanol tolerance (6% v/v)", "halotolerance (0.5M NaCl)", "Osmotolerance (40% w/v sorbitol)", "SO3 tolerance (1.0 mM)"))
lineage.toplot$variantPos.FW.<-factor(lineage.toplot$variantPos.FW., levels = lineage.toplot$variantPos.FW.)

ggplot(data=lineage.toplot,  
       aes(y=key, x=variantPos.FW., size=value, fill=parent, 
           group = interaction(factor(parent),factor(variantPos.FW.)))) + 
  geom_point(color="#000000", shape=21, position=position_dodge(.5)) +
  scale_size(range=c(0,3),breaks=c(0,25,75,100),labels=c("0","25%","75%","100%"),guide="legend") +
  scale_fill_manual(values=c("#0F579A", "#FDD36F")) +
  theme_bw(base_size = 10) + 
  theme(axis.text.x=element_text(face="bold", size=8), legend.position='bottom') +
  scale_y_discrete(labels = new_labels, limits = rev(levels(lineage.toplot$key))) +
  xlab("QTN") + facet_wrap(~Condition, scales = "free_x", nrow = 3, ncol = 8)

#save as 16 x 11 or 15 x 11

###plot symbols for peaks
lineage.toplot$peak_symbol<-ifelse(lineage.toplot$peak == 1, "QTN", "Non-QTN")

ggplot() + geom_point(data=subset(lineage.toplot, key == "X1A" & parent == "RM"), aes(x=variantPos.FW., y=1, shape=peak_symbol, color=peak_symbol, fill=peak_symbol)) +
  facet_wrap(~Condition, scales = "free_x",  nrow = 3, ncol = 8) +
  scale_shape_manual(values=c(22, 25)) +
  scale_fill_manual(values=c("#1745FC","#C11B17")) +
  scale_color_manual(values=c("#1745FC","#C11B17")) +
  theme_minimal()+
  scale_y_continuous(limits=c(1), expand = c(0, 0)) +
  theme(panel.grid = element_blank()) +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside")
#### Figure 3  ####
PVAL_Data <- #Supplementary data Figure 3#

QTL_to_change_1 <- "SUC2"
QTL_to_change_2 <- "IMA1"
QTL_to_change_3 <- "ALD6"
QTL_to_change_4 <- "URA5"
QTL_to_change_5 <- "URK1"

trait_to_change_1 <- c("1-propanol conc.", "ethyl acetate conc.", 
                       "isoamyl alcohol conc.",
                       "glycerol conc.", "acetic acid conc.", 
                       "ethanol conc.")
trait_to_change_2 <- c("ethyl acetate conc.",
                       "phenethyl acetate conc.",
                       "isobutanol conc.",
                       "glycerol conc.", 
                       "acetic acid conc." , 
                       "maltose growth")
trait_to_change_3 <- c("acetic acid conc.", 
                       "ethyl octanoate conc.",
                       "isoamyl alcohol conc.")
trait_to_change_4 <- c("halotolerance",
                       "ethyl octanoate conc.",
                       "acetic acid conc.", 
                       "ethanol tolerance",
                       "acetic acid tolerance")
trait_to_change_5 <- c("isobutanol conc.")

chromosome_1 <- "Chr IX"
chromosome_2 <- "Chr VII"
chromosome_3 <- "Chr XVI"
chromosome_4 <- "Chr XIII"
chromosome_5 <- "Chr XIV"

qtl_coord_1 <- PVAL_Data$pos[5461:5463]
qtl_coord_2 <- PVAL_Data$pos[4837:4840]
qtl_coord_3 <- PVAL_Data$pos[11572:11573]
qtl_coord_4 <- PVAL_Data$pos[8458:8459]
qtl_coord_5 <- PVAL_Data$pos[9801:9802]

dat_chr <- data.frame(
  label = c("Chr IX", "Chr VII", "Chr XVI", "Chr XIII", "Chr XIV"),
  QTL   = c("QTL 1","QTL 2", "QTL 3", "QTL 4", "QTL 5")
)

dat_seg <- data.frame(
  pos = c("37379", "1067985", "432058", "56990", "647843"),
  xmax = c("37779", "1068695", "432771", "57813", "648789"),
  PVAL = c(Inf, Inf, Inf, Inf, Inf),
  QTL   = c("QTL 1","QTL 2", "QTL 3", "QTL 4", "QTL 5")
)


subset_data_1 <- PVAL_Data %>% slice(5445:5470) %>% gather(trait,PVAL, 6:26, factor_key = TRUE) %>% dplyr::select(trait, chr, pos, position,PVAL)%>%
  filter(trait %in% trait_to_change_1)%>% mutate(QTL="QTL 1")
subset_data_2 <- PVAL_Data %>% slice(4820:4841) %>% gather(trait,PVAL, 6:26, factor_key = TRUE) %>% dplyr::select(trait, chr, pos, position,PVAL)%>%
  filter(trait %in% trait_to_change_2)%>% mutate(QTL="QTL 2")
subset_data_3 <- PVAL_Data %>% slice(11559:11580) %>% gather(trait,PVAL, 6:26, factor_key = TRUE) %>% dplyr::select(trait, chr, pos, position,PVAL)%>% 
  filter(trait %in% trait_to_change_3)%>% mutate(QTL="QTL 3")
subset_data_4 <- PVAL_Data %>% slice(8453:8475) %>% gather(trait,PVAL, 6:26, factor_key = TRUE) %>% dplyr::select(trait, chr, pos, position,PVAL)%>%
  filter(trait %in% trait_to_change_4)%>% mutate(QTL="QTL 4")
subset_data_5 <- PVAL_Data %>% slice(9786:9812) %>% gather(trait,PVAL, 6:26, factor_key = TRUE) %>% dplyr::select(trait, chr, pos, position,PVAL)%>%
  filter(trait %in% trait_to_change_5)%>% mutate(QTL="QTL 5")

PVAL_subset_all <- rbind(subset_data_1, 
                         subset_data_2,
                         subset_data_3,
                         subset_data_4,
                         subset_data_5)

qtl_coord_1 <- PVAL_Data$pos[5461:5463]
qtl_coord_2 <- PVAL_Data$pos[4837:4840]
qtl_coord_3 <- PVAL_Data$pos[11572:11573]
qtl_coord_4 <- PVAL_Data$pos[8458:8459]
qtl_coord_5 <- PVAL_Data$pos[9801:9802]

Figure_3 <- ggplot(PVAL_subset_all, aes(x = pos, y = PVAL)) +
  facet_wrap(~QTL, scales = "free")+
  geom_rect(xmin = 37385, xmax = 38983, ymin = -Inf, ymax = Inf,
            fill = "grey90", alpha=0.1,inherit.aes=FALSE)+
  geom_rect(xmin = 1067222, xmax = 1068991, ymin = -Inf, ymax = Inf,
            fill = "grey90", alpha=0.1,inherit.aes=FALSE)+
  geom_rect(xmin = 432588, xmax = 434090, ymin = -Inf, ymax = Inf,
            fill = "grey90", alpha=0.1,inherit.aes=FALSE)+
  geom_rect(xmin = 56773, xmax = 57453, ymin = -Inf, ymax = Inf,
            fill = "grey90", alpha=0.1,inherit.aes=FALSE)+
  geom_rect(xmin = 647432, xmax = 648937, ymin = -Inf, ymax = Inf,
            fill = "grey90", alpha=0.1,inherit.aes=FALSE)+
  geom_line(size=0.5)+
  geom_point(aes(fill=trait), size=3, shape= 21)+
  geom_hline(yintercept = 5.2, color="grey40", size=0.25, linetype = "longdash")+
  ylab("PVAL")+
  theme_classic(base_size = 16)+
  theme(legend.position = c(0.85,0.23),
        legend.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_blank(),
        axis.text=element_text(colour = "black"),
        axis.line =  element_line(size=0.5),
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text.x = element_text(hjust=0,size=14, angle = 0),
        strip.background = element_blank(),
        plot.title = element_text(size=14, face="bold"),
        plot.margin = unit(c(1,1,1,1), "lines"))+
  scale_fill_brewer(palette = "Set3")+
  scale_x_continuous(breaks = breaks_pretty(n = 3))+
  guides(fill = guide_legend(override.aes = list(size=3.5), ncol=2))+
  coord_cartesian(clip = "off")+
  geom_text(
    data    = dat_chr,
    mapping = aes(x = Inf, y = Inf, label = label),
    hjust   = 1.1,
    vjust   = 1.5
  )

ggdraw(Figure_3) + 
  draw_label("SUC2", fontface=3, x = 0.20, y = 0.93, size=11.5, color = "grey45") +
  draw_label("IMA1", fontface=3, x = 0.60, y = 0.93, size=11.5, color = "grey45") +
  draw_label("ALD6", fontface=3, x = 0.84, y = 0.93, size=11.5, color = "grey45") +
  draw_label("URA5", fontface=3, x = 0.18, y = 0.45, size=11.5, color = "grey45") +
  draw_label("URK1", fontface=3, x = 0.50, y = 0.45, size=11.5, color = "grey45") 


#### Figure 4  ####
QTL_ver_all <- #Supplementary data Figure 4#
  
KO_compare_all <- QTL_ver_all %>%
  group_by(QTL, date, trait, class) %>%
  mutate(control_mean = mean(value[background == "WT"], na.rm = T), 
         value_norm = value/control_mean) 

# calculate sd
sd_all <- KO_compare_all %>% 
  group_by(QTL,trait, class, background) %>%
  summarise(mean.val = mean(value_norm, na.rm = TRUE),
            sd=sd(value_norm))%>%
  mutate(value_norm=1)

stat.test_all <- KO_compare_all %>% 
  group_by(QTL, trait, class) %>% 
  t_test(value_norm~background)%>%
  add_significance() %>%
  mutate(x_pos=c(1,1,2,2,3,3,4,4,5,5,6,6,1,1,2,2,3,3,4,4,5,5,6,6,
                 1,1,2,2,3,3,1,1,2,2,3,3,4,4,5,5,1,1))%>%
  add_xy_position( x = "x_pos")

Figure_4 <- ggplot(KO_compare_all, 
               aes(trait,value_norm))+
  facet_grid(class~QTL, space= "free_x", 
             scales="free_x")+
  geom_hline(yintercept = 1, color="grey", size=0.5, linetype="dashed")+
  geom_errorbar(data = sd_all, 
                aes(ymin = mean.val-sd, ymax = mean.val+sd, colour = background), 
                width=0, lwd=0.5, position=position_nudge(x=0))+
  stat_summary(fun.y=mean, geom="point", aes(color=background), size=1.5, position=position_nudge(x=0))+ 
  theme_classic(base_size = 16)+
  theme(axis.title.y.right = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.line.y = element_blank(),
        axis.text.y.right = element_blank(),                  
        axis.line =  element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        panel.spacing = unit(1, "mm"),
        panel.border = element_rect(size=0.5, fill="transparent"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(face = "italic"),
        strip.text.y.left = element_text(hjust=0,size=12, color = "grey45", angle = 0),
        legend.position = "none",
        plot.margin=unit(c(1,1,1,1.2),"cm"))+
  scale_color_manual(values = c("tomato3",NA))+
  scale_fill_manual(values = c(NA, "grey30"))+
  scale_y_continuous( breaks = seq(0.5, 2, by = 0.5))+
  ylim(0, 2)+
  stat_pvalue_manual(
    stat.test_all_new, 
    label = "p.signif", tip.length = 0,
    coord.flip = FALSE, y.position = c(rep(1.8)), hide.ns = TRUE)


#### Figure 5B ####
SNV_seg <- #Supplementary data Figure 5B#

violin_SNV_new <-
  mutate(SNV_seg, group  =
           ifelse(V5461 == "R" & V5462 == "R", "RR",
           ifelse(V5461 == "R" & V5462 == "Y", "RY",
           ifelse(V5461 == "Y" & V5462 == "R", "YR",
           ifelse(V5461 == "Y" & V5462 == "Y", "YY", "-"))))) %>%
  filter(group!="-")

my_comparisons <- list( c("RR", "RY"), c("RR", "YY"), c("RY", "YY") )

Figure_5B <- ggplot(violin_SNV_new %>% filter(group!= "-"), aes(group, abv)) +
  geom_violin(aes(fill = group), color= "black", draw_quantiles = c(0.5), alpha=0.8, size=0.5) +
  geom_point(aes(fill = group), shape=21, draw_quantiles = c(0.5),alpha=0.2) +
  stat_compare_means(label = "p.signif", 
                     comparisons = my_comparisons, 
                     method = "t.test", 
                     size=5,
                     tip.length = 0.1,
                     label.y = c(6.75, 6.95, 7.2))+
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(),
    axis.line =  element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    plot.title = element_text(size=14)) +
  scale_fill_brewer(palette = "Set2")+
  ylim(c(5.25,7.5))+
  ylab("ethanol conc. (% v/v)")+
  scale_x_discrete(labels = c("T/A", "T/\u0394", "C/A", "C/\u0394"))+
  xlab("")
#### Figure 5C ####
SUC2_swap <- #Supplementary data Figure 5C#

my_SNV_comparisons_new <- list(c("WT", "SNV5461"),
                               c("WT", "SNV5462"))

SUC2_abv <- SUC2_swap %>%
  mutate(sample = factor(sample, levels = c("SNV5461",
                                            "SNV5462",
                                            "WT",
                                            "Blank")))

p_abv_sample <- ggplot(SUC2_abv, aes(sample,value))+
  facet_grid(~class, scales = "free_x", space = "free_x") +
  stat_summary(fun.y=mean, geom="bar",aes(fill=name),
               width = 0.8,
               color = "black",
               alpha = 0.8,
               size = 0.5) +
  stat_summary(fun.data=mean_sdl, geom="errorbar", 
               width=0.25, color="black", lwd=0.5, 
               position=position_dodge(.7),
               fun.args = list(mult = 1))+
  stat_compare_means(
    aes(group=sample),
    method = "t.test",
    comparisons = my_SNV_comparisons_new,
    label = "p.signif",
    size = 5,
    tip.length = 0.1,
    label.y = c(6.75, 7.5))+
  theme_classic(base_size = 16)+
  theme(legend.position = "none",
        axis.line =  element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        plot.title = element_text(size=14),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(colour = "black", fill=NA, size=0.5))+
  scale_fill_manual(values=c("#8DA0CB","#FC8D62","#66C2A5","#FC8D62","#8DA0CB" ,"#E78AC3"))+
  ylab("ethanol conc. (% v/v)")+
  ylim(c(0,8))+
  xlab("")
#### Figure 5D ####
SNV_barplot <- #Supplementary data Figure 5C#
  
stat.test <- compare_means(
  value ~ background, data = SNV_barplot, group.by = c("class"),
  method = "t.test", ref.group = "WT"
)

SNV5462 <- ggbarplot(SNV_barplot, x = "class", y = "value",
                     fill = "sequence",
                     color="background",
                     add = "mean_sd", 
                     add.params = list(group = "background"),
                     position = position_dodge(0.8),
                     alpha=0.5,
                     size = 0.5)

p_extra_strain_SNV5462 <- SNV5462 + stat_pvalue_manual(
  stat.test, x = "class", y.position = 8,
  label = "p.signif",
  position = position_dodge(0.8))+
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title=element_text(size=12),
    axis.text.x = element_text(size = 14),
    axis.line =  element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, color = "grey45"),
    panel.border = element_rect(colour = "black", fill=NA, size=0.6)) +
    scale_fill_manual(name= "Variant (394)\ngenotype ",
                    labels=c("A","\u0394"),
                    values = c("grey45","tomato3"),
                    guide = guide_legend(reverse=TRUE))+
  scale_color_manual(values = c("black", "black"))+
  scale_y_continuous(limits = c(0, max(SNV_barplot$value) * 1.1))+
  xlab("")+
  ylab("ethanol conc. (% v/v)")+
  guides(color = FALSE)
#### Figure 6B ####
ALD6_df <- #Supplementary data Figure 6B#

my_comparisons <- list(c("R", "Y"))

Figure_6B <- ggplot(ALD6_df %>% filter(V11573!= "-"), aes(V11573, trait)) +
  geom_violin(aes(fill = V11573), color= "black", draw_quantiles = c(0.5), alpha=0.5, size=0.5) +
  geom_point(aes(fill = V11573), shape=21, draw_quantiles = c(0.5),alpha=0.2) +
  stat_compare_means(label = "p.signif", 
                     comparisons = my_comparisons, 
                     method = "t.test", 
                     size=4)+
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.x = element_text(),
    axis.line =  element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    plot.title = element_text(size=14)
  ) +
  ggtitle("acetic acid (g/l)") +
  scale_fill_manual(name="Variant genotype",values = c("tomato3", "grey45"))+
  scale_x_discrete(labels = c("C", "A"))+
  xlab("Variant genotype (184)")+
  ylim(c(0, 0.8))
#### Figure 6C ####
SNV_barplot <- #Supplementary data Figure 6C#

my_SNV_comparisons <- list(c("WT", "SNV11573"))

Figure_6C <- 
  ggplot(SNV_barplot, aes(sample, value)) +
  facet_wrap(class ~ ., scales = "free_x") +
  stat_summary(
    fun.y = mean,
    geom = "bar",
    aes(fill = protein),
    size = 0.5,
    width = 0.7,
    color = "black",
    alpha = 0.5)+
  stat_summary(
    fun.data = mean_sdl,
    geom = "errorbar",
    width = 0.2,
    color = "black",
    lwd = 0.5,
    position = position_dodge(.7),
    fun.args = list(mult = 1)) +
  stat_compare_means(
    method = "t.test",
    comparisons = my_SNV_comparisons,
    label = "p.signif",
    size = 5,
    label.y = c(max(SNV_barplot$value) * 1.1,max(SNV_barplot$value) * 1.2) ,
    tip.length = 0.02) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 14),
    axis.line =  element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, color = "grey45"),
    plot.title = element_text(size=14, hjust=0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_fill_manual(name="Variant genotype",
                    labels=c("C", "A"),
                    values = c("tomato3", "grey45"))+
  ylab("acetic acid (g/l)") +
  scale_y_continuous(limits = c(0, max(SNV_barplot$value) * 1.2))+
  scale_x_discrete(labels = c("Mutant","WT"))+
  xlab("")
#### Figure 6D ####
SNV_barplot_other <- #Supplementary data Figure 6D#
  
my_SNV_comparisons <- list(c("WT", "SNV11573"))

Figure_6D <- 
  ggplot(SNV_barplot, aes(sample, value)) +
  facet_wrap(class ~ ., scales = "free_x") +
  stat_summary(
    fun.y = mean,
    geom = "bar",
    aes(fill = protein),
    size = 0.5,
    width = 0.7,
    color = "black",
    alpha = 0.5)+
  stat_summary(
    fun.data = mean_sdl,
    geom = "errorbar",
    width = 0.2,
    color = "black",
    lwd = 0.5,
    position = position_dodge(.7),
    fun.args = list(mult = 1)) +
  stat_compare_means(
    method = "t.test",
    comparisons = my_SNV_comparisons,
    label = "p.signif",
    size = 5,
    label.y = c(max(SNV_barplot$value) * 1.1,max(SNV_barplot$value) * 1.2) ,
    tip.length = 0.02) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 14),
    axis.line =  element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, color = "grey45"),
    plot.title = element_text(size=14, hjust=0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_fill_manual(name="Variant genotype",
                    labels=c("C", "A"),
                    values = c("tomato3", "grey45"))+
  ylab("acetic acid (g/l)") +
  scale_y_continuous(limits = c(0, max(SNV_barplot$value) * 1.2))+
  scale_x_discrete(labels = c("Mutant","WT"))+
  xlab("")