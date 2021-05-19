#### Figure S1 ####
phenotypes_combined_subset <-  #Supplementary data Figure S1#

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
  #mutate(cleanname = gsub(" tolerance","\ntolerance", cleanname)) %>%
  mutate(cleanname = paste0(cleanname,suffix)) %>%
  pull(cleanname)
num_columns <- (length(phenotypes_combined_subset) - 1)
Figure_S1 <- ggpairs(phenotypes_combined_subset,
                             columns = 1:(length(phenotypes_combined_subset)-1),
                             ggplot2::aes(colour = as.factor(segregant), fill = as.factor(segregant),
                                          alpha = as.factor(segregant)),
                             lower = list(continuous = wrap("points"), mapping = ggplot2::aes(size = segregant)),
                             diag = list(continuous = my_diag),
                             #upper = list(continuous = "blank"), 
                             columnLabels = clean_labels,
                             switch = "both",
                             legend = clean_legend)  +
  theme_custom
Figure_S1 <- scale_manual(Figure_S1)
#### Figure S2 ####
SUC2_swap <- #Supplementary data Figure S2#

my_sucorse_comparisons <- list(c("suc2", "SNV5462"),
                               c("suc2", "WT"),
                               c("SNV5462","WT"))


Figure_S2 <- ggplot(SUC2_swap, aes(sample,value))+
  facet_grid(~class, scales = "free_x", space = "free_x") +
  stat_summary(fun.y=mean, geom="bar",aes(fill=sample),
               width = 0.8,
               color = "white",
               alpha = 0.5,
               size = 0.00) +
  stat_summary(fun.data=mean_sdl, geom="errorbar", 
               width=0.25, color="black", lwd=0.5, 
               position=position_dodge(.7),
               fun.args = list(mult = 1))+
  stat_compare_means(
    method = "t.test",
    comparisons = my_sucorse_comparisons,
    label = "p.signif",
    size = 5,
    tip.length = 0.005,
    label.y = c(4, 4.3, 4.6))+
  theme_classic(base_size = 16)+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.line =  element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        plot.title = element_text(hjust=0.5, size=14),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_blank())+
  scale_fill_manual(values=c("royalblue3", "royalblue3", "grey45", "04151F"))+
  ggtitle("Sucrose (g/l)")+
  coord_cartesian(ylim=c(0, 5))+
  xlab("")
#### Figure S3 ####
SNV_SWAP <- #Supplementary data Figure S3#


trait_to_change_SUC2 <- c("glycerol conc.", 
                          "acetic acid conc.", 
                          "ethanol conc.", 
                          "1-propanol conc.", 
                          "ethyl acetate conc.")
trait_to_change_ALD6 <- c("1-propanol conc.", 
                          "isoamyl alcohol conc.", 
                          "acetic acid conc.", 
                          "ethyl octanoate conc.")

trait_dotplot_all_SUC2 <- SNV_SWAP %>% filter(trait %in% trait_to_change_SUC2, QTL == "SUC2", sample != "SNV5461") %>%
  group_by(date, trait, class) %>%
  mutate(control_mean = mean(value[background == "WT"], na.rm = T), 
         value_norm = value/control_mean) %>%
  mutate(trait = factor(trait, levels = c("ethanol conc.", 
                                          "acetic acid conc.",
                                          "glycerol conc.",
                                          "1-propanol conc.", 
                                          "ethyl acetate conc.")))
trait_dotplot_all_ALD6 <- SNV_SWAP %>% filter(trait %in% trait_to_change_ALD6, QTL == "ALD6", class != "Sb") %>%
  group_by(date, trait, class) %>%
  mutate(control_mean = mean(value[background == "WT"], na.rm = T), 
         value_norm = value/control_mean) %>%
  mutate(trait = factor(trait, levels = c("acetic acid conc.",
                                          "1-propanol conc.", 
                                          "isoamyl alcohol conc.", 
                                          "ethyl octanoate conc.")))

trait_dotplot_all_SUC2$trait <- as.factor(trait_dotplot_all_SUC2$trait)
trait_dotplot_all_ALD6$trait <- as.factor(trait_dotplot_all_ALD6$trait)

stat.test_SUC2 <- trait_dotplot_all_SUC2 %>% group_by(trait, class) %>% 
  t_test(value_norm~background)%>% 
  add_significance("p.adj")
stat.test_ALD6 <- trait_dotplot_all_ALD6 %>% group_by(trait, class) %>% 
  t_test(value_norm~background)%>% 
  add_significance("p.adj")

stat.test_SUC2 <- stat.test_SUC2 %>% add_xy_position(fun = "mean", x = "background")
stat.test_SUC2

stat.test_ALD6 <- stat.test_ALD6 %>% add_xy_position(fun = "mean", x = "background")
stat.test_ALD6


# Figure_S3A 

Figure_S3A <- ggplot(trait_dotplot_all_SUC2, aes(background,value_norm))+
  facet_grid(trait~class, scales="free_x", switch = "y", labeller = label_wrap_gen(width=10))+
  geom_hline(yintercept = 1, color="grey", size=0.5, linetype="dashed")+
  stat_summary(fun.data=mean_sdl, geom="errorbar", aes(color=background), width=0, lwd=0.5, position=position_dodge(.7))+
  stat_summary(fun.y=mean, geom="point", aes(color=background), size=1.5) +
  theme_classic(base_size = 16)+
  theme(axis.title.y.right = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y.right = element_blank(),                  
        axis.ticks.y = element_blank(),
        axis.line =  element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        panel.spacing = unit(1, "mm"),
        panel.border = element_rect(size=0.5, fill="transparent"), 
        strip.background = element_blank(),
        strip.text = element_text(size=14, color = "grey45"),
        legend.position = "none",
        plot.title = element_text(size=12, face = "bold.italic")
  )+ 
  scale_color_manual(values = c("tomato3", "grey40"))+
  scale_y_continuous(limits = c(min(trait_dotplot_all_SUC2$value_norm), 
                                max(trait_dotplot_all_SUC2$value_norm)*1.01),
                     breaks = seq(0, max(trait_dotplot_all_SUC2$value_norm)*1.3, by = 0.2))+
  stat_pvalue_manual(
    stat.test_SUC2, label = "p.adj.signif", tip.length = 0, vjust = -0.05,
    coord.flip = TRUE, 
    y.position = max(trait_dotplot_all_SUC2$value_norm)*1,
    hide.ns = TRUE
  ) +
  coord_flip()


# Figure_S3B 

Figure_S3B <- ggplot(trait_dotplot_all_ALD6, aes(background,value_norm))+
  facet_grid(trait~class, scales="free_x", switch = "y", labeller = label_wrap_gen(width=10))+
  geom_hline(yintercept = 1, color="grey", size=0.5, linetype="dashed")+
  stat_summary(fun.data=mean_sdl, geom="errorbar", aes(color=background), width=0, lwd=0.5, position=position_dodge(.7))+
  stat_summary(fun.y=mean, geom="point", aes(color=background), size=1.5) +  #directlabels::geom_dl(data= trait_dotplot_all_ALD6 %>% filter(trait=="acetic acid" & class == "RM") %>% 
  theme_classic(base_size = 16)+
  theme(axis.title.y.right = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y.right = element_blank(),                  
        axis.ticks.y = element_blank(),
        axis.line =  element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        panel.spacing = unit(1, "mm"),
        panel.border = element_rect(size=0.5, fill="transparent"), 
        strip.background = element_blank(),
        strip.text = element_text(size=14, color = "grey45"),
        legend.position = "none",
        plot.title = element_text(size=12, face = "bold.italic")
  )+
  scale_color_manual(values = c("tomato3", "grey40"))+
  scale_y_continuous(limits = c(min(trait_dotplot_all_ALD6$value_norm), 
                                max(trait_dotplot_all_ALD6$value_norm)*1.01),
                     breaks = seq(0, max(trait_dotplot_all_ALD6$value_norm)*1.2, by = 0.2))+
  stat_pvalue_manual(
    stat.test_ALD6, label = "p.adj.signif", tip.length = 0, vjust = -0.05,
    coord.flip = TRUE, 
    y.position = max(trait_dotplot_all_ALD6$value_norm)*1.001,
    hide.ns = TRUE
  ) +
  coord_flip()

#### Figure S4 ####
SNV_SWAP_sup <- #Supplementary data Figure S4#


trait_dotplot_ima1 <- SNV_SWAP_sup %>% 
  group_by(date, trait, class) %>%
  mutate(control_mean = mean(value[background == "WT"], na.rm = T), 
         value_norm = value/control_mean) 

trait_dotplot_ima1$trait <- as.factor(trait_dotplot_ima1$trait)

my_SNV_comparisons_new <- list(c("WT", "SNV4837"),
                               c("WT", "SNV4835"))

stat.test_ima1 <- trait_dotplot_ima1 %>% 
  group_by(trait, class) %>% 
  t_test(value_norm~sample) %>%
  adjust_pvalue(method = 'bonferroni')%>%
  adjust_pvalue() 

stat.test_ima1 <- stat.test_ima1 %>% add_xy_position(fun = "mean", x = "background")

Figure_S4 <- ggplot(trait_dotplot_ima1, aes(sample,value_norm))+
  facet_grid(trait~class, scales="free_x", switch = "y", labeller = label_wrap_gen(width=15))+
  geom_hline(yintercept = 1, color="grey", size=0.5, linetype="dashed")+
  stat_summary(fun.data=mean_sdl, geom="errorbar", aes(color=sample), width=0, lwd=0.5, position=position_dodge(.7))+
  stat_summary(fun.y=mean, geom="point", aes(color=sample), size=1.5) +
  theme_cowplot()+
  theme_classic(base_size = 16)+
  theme(axis.title.y.right = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y.right = element_blank(),                  
        axis.ticks.y = element_blank(),
        axis.line =  element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        panel.spacing = unit(1, "mm"),
        panel.border = element_rect(size=0.5, fill="transparent"), 
        strip.background = element_blank(),
        strip.text = element_text(size=14, color = "grey45"),
        legend.position = "none"
  )+ 
  scale_color_manual(values = c("tomato3","steelblue", "grey40"))+
  stat_pvalue_manual(
    stat.test_ima1, label = "p.adj.signif", tip.length = 0,
    coord.flip = TRUE, y.position=c(1.2, 1.25))+
  coord_flip()
#### Figure S5 ####
# read data
cov_file = #Supplementary data Figure S5A#
MALgenes = #Supplementary data Figure S5B#

cov_file$V1 = factor(cov_file$V1,
                     levels = c("SRR5634755", "SRR5634774", "SRR5630056", "SRR5634576", "SRR5630345", "SRR5634504", "SRR5630053", "SRR5630190", "SRR5630087", "SRR5634745", "SRR5629932", "SRR5634542", "SRR5630013", "SRR5634491", "SRR5634513", "SRR5634546", "SRR5634530", "SRR5630171", "SRR5634612", "SRR5630175", "SRR5634823", "SRR5630349", "SRR5634650", "SRR5634752", "SRR5629824", "SRR5634416", "SRR5634482", "SRR5629863", "SRR5629972", "SRR5634730", "SRR5630447", "SRR5634498", "SRR5630208", "SRR5634441", "SRR5630164", "SRR5634349", "SRR5634368", "SRR5634725", "SRR5634767", "SRR5629868", "SRR5634391", "SRR5634638",
                                "SRR5629859", "SRR5630263", "SRR5630419", "SRR5634667", "SRR5634507", "SRR5634666", "SRR5630353", "SRR5634460", "SRR5629879", "SRR5629872", "SRR5630424", "SRR5634547", "SRR5634657", "SRR5630200", "SRR5634759", "SRR5634661", "SRR5630362", "SRR5634676", "SRR5629810", "SRR5630147", "SRR5630299", "SRR5629845", "SRR5630331", "SRR5630356", "SRR5630399", "SRR5634466", "SRR5630052", "SRR5630103", "SRR5630291", "SRR5634408", "SRR5630093", "SRR5629950", "SRR5634532", "SRR5629838", "SRR5630368", "SRR5630076", "SRR5634463", "SRR5630211", "SRR5630040", "SRR5629851", "SRR5634421", "SRR5634645",
                                "SRR5634686", "SRR5629801", "SRR5634586", "SRR5630121", "SRR5630034", "SRR5629835", "SRR5630351", "SRR5630023", "SRR5629830", "SRR5630032", "SRR5630444", "SRR5634431", "SRR5630438", "SRR5634596", "SRR5630320", "SRR5634364", "SRR5634459", "SRR5634540", "SRR5630048", "SRR5634786", "SRR5634558", "SRR5629805", "SRR5630192", "SRR5630074", "SRR5630148", "SRR5630105", "SRR5629891", "SRR5634362", "SRR5629920", "SRR5634438", "SRR5634406", "SRR5629864", "SRR5634422", "SRR5630425", "SRR5630028", "SRR5630326", "SRR5634499", "SRR5634706", "SRR5630170", "SRR5630179", "SRR5630233", "SRR5634810",
                                "SRR5630128", "SRR5630225", "SRR5629794", "SRR5630243", "SRR5630155", "SRR5630384", "SRR5630264", "SRR5634570", "SRR5630045", "SRR5629959", "SRR5630374", "SRR5630395", "SRR5634743", "SRR5630246", "SRR5629934", "SRR5634601", "SRR5634402", "SRR5629973", "SRR5634414", "SRR5630129", "SRR5634673", "SRR5634519", "SRR5630276", "SRR5629981", "SRR5634696", "SRR5634740", "SRR5630254", "SRR5630337", "SRR5634369", "SRR5634618", "SRR5629829", "SRR5630304", "SRR5634575", "SRR5634735", "SRR5634640", "SRR5629987", "SRR5634516", "SRR5634736", "SRR5629825", "SRR5629978", "SRR5630184", "SRR5634409",
                                "SRR5634792", "SRR5629791", "SRR5629871", "SRR5634581", "SRR5634616", "SRR5634452", "SRR5634723", "SRR5634512", "SRR5634642", "SRR5630266", "SRR5630448", "SRR5634588", "SRR5634430", "SRR5629874", "SRR5634813", "SRR5634500", "SRR5630280", "SRR5630379", "SRR5634753", "SRR5630270", "SRR5629849", "SRR5630008", "SRR5630375", "SRR5629913", "SRR5634824", "SRR5634623", "SRR5630086", "SRR5630396", "SRR5630365", "SRR5634388", "SRR5630075", "SRR5634567", "SRR5634802", "SRR5634423", "SRR5629941", "SRR5634660", "SRR5634350", "SRR5634474", "SRR5630146", "SRR5630389", "SRR5630415", "SRR5634481",
                                "SRR5634634", "SRR5630242", "SRR5630406", "SRR5629837", "SRR5630257", "SRR5630031", "SRR5630065", "SRR5630322", "SRR5634505", "SRR5630080", "SRR5634539", "SRR5634456", "SRR5630376", "SRR5634415", "SRR5634665", "SRR5630113", "SRR5634476", "SRR5630067", "SRR5629806", "SRR5629886", "SRR5630361", "SRR5630039", "SRR5634494", "SRR5629847", "SRR5629788", "SRR5634424", "SRR5630324", "SRR5630369", "SRR5630259", "SRR5630329", "SRR5634502", "SRR5630012", "SRR5634508", "SRR5634681", "SRR5634777", "SRR5634744", "SRR5630152", "SRR5634687", "SRR5630081", "SRR5634708", "SRR5634535", "SRR5634428",
                                "SRR5634348", "SRR5634716", "SRR5630255", "SRR5630385", "SRR5634822", "SRR5629903", "SRR5634464", "SRR5634403", "SRR5630068", "SRR5634399", "SRR5630109", "SRR5630408", "SRR5634677", "SRR5629968", "SRR5634798", "SRR5630106", "SRR5630073", "SRR5634515", "SRR5630290", "SRR5629971", "SRR5629974", "SRR5629858", "SRR5629979", "SRR5630104", "SRR5634603", "SRR5634492", "SRR5634707", "SRR5630413", "SRR5634384", "SRR5634746", "SRR5629841", "SRR5634779", "SRR5629797", "SRR5629834", "SRR5634347", "SRR5629820", "SRR5630046", "SRR5630295", "SRR5629840", "SRR5630241", "SRR5630205", "SRR5630265",
                                "SRR5634655", "SRR5630357", "SRR5630077", "SRR5634375", "SRR5629815", "SRR5630178", "SRR5629792", "SRR5630101", "SRR5630433", "SRR5630401", "SRR5630403", "SRR5630286", "SRR5630303", "SRR5634719", "SRR5634549", "SRR5634790", "SRR5634554", "SRR5629969", "SRR5634527", "SRR5630102", "SRR5630364", "SRR5634437", "SRR5634580", "SRR5630199", "SRR5634478", "SRR5634643", "SRR5634592", "SRR5634607", "SRR5630273", "SRR5629857", "SRR5630227", "SRR5629782", "SRR5630333", "SRR5629921", "SRR5630334", "SRR5629867", "SRR5630082", "SRR5630232", "SRR5630397", "SRR5629962", "SRR5634713", "SRR5634383",
                                "SRR5630250", "SRR5629882", "SRR5630318", "SRR5634602", "SRR5634633", "SRR5630026", "SRR5630226", "SRR5630132", "SRR5630383", "SRR5634649", "SRR5629823", "SRR5634521", "SRR5634386", "SRR5630177", "SRR5629802", "SRR5634714", "SRR5629951", "SRR5634405", "SRR5629915", "SRR5634511", "SRR5634652", "SRR5634754", "SRR5634747", "SRR5630301", "SRR5630244", "SRR5629964", "SRR5634738", "SRR5630358", "SRR5634390", "SRR5630180", "SRR5629832", "SRR5634750", "SRR5634709", "SRR5630214", "SRR5630278", "SRR5634670", "SRR5629929", "SRR5630268", "SRR5630439", "SRR5634731", "SRR5629911", "SRR5634662",
                                "SRR5634401", "SRR5634587", "SRR5630251", "SRR5630445", "SRR5629984", "SRR5634366", "SRR5630112", "SRR5630138", "SRR5634598", "SRR5630411", "SRR5634632", "SRR5634473", "SRR5634704", "SRR5630123", "SRR5634694", "SRR5630169", "SRR5630154", "SRR5634804", "SRR5630159", "SRR5634659", "SRR5629908", "SRR5630118", "SRR5630017", "SRR5634577", "SRR5629963", "SRR5629822", "SRR5629997", "SRR5634639", "SRR5630441", "SRR5630094", "SRR5630410", "SRR5634664", "SRR5630248", "SRR5630370", "SRR5634674", "SRR5634668", "SRR5630328", "SRR5634720", "SRR5634726", "SRR5630224", "SRR5629831", "SRR5634520",
                                "SRR5629783", "SRR5634626", "SRR5630186", "SRR5630089", "SRR5634454", "SRR5634648", "SRR5634465", "SRR5634372", "SRR5630409", "SRR5630137", "SRR5634757", "SRR5634772", "SRR5630041", "SRR5630124", "SRR5634352", "SRR5629875", "SRR5634742", "SRR5634784", "SRR5630059", "SRR5630209", "SRR5634449", "SRR5634410", "SRR5630422", "SRR5630258", "SRR5634543", "SRR5634773", "SRR5629839", "SRR5634453", "SRR5630236", "SRR5634443", "SRR5629796", "SRR5630191", "SRR5634710", "SRR5629976", "SRR5634551", "SRR5630204", "SRR5630058", "SRR5630267", "SRR5629860", "SRR5634712", "SRR5630157", "SRR5634404",
                                "SRR5630231", "SRR5634815", "SRR5630294", "SRR5634379", "SRR5634518", "SRR5634669", "SRR5634761", "SRR5630055", "SRR5630256", "SRR5629938", "SRR5634729", "SRR5634688", "SRR5630011", "SRR5629826", "SRR5634495", "SRR5630437", "SRR5630153", "SRR5630033", "SRR5634658", "SRR5634411", "SRR5634387", "SRR5634426", "SRR5630070", "SRR5630293", "SRR5634796", "SRR5630330", "SRR5630079", "SRR5629989", "SRR5630402", "SRR5629878", "SRR5629936", "SRR5630090", "SRR5634764", "SRR5630078", "SRR5634357", "SRR5630125", "SRR5630335", "SRR5634485", "SRR5630063", "SRR5634628", "SRR5630271", "SRR5634826",
                                "SRR5629967", "SRR5634819", "SRR5629862", "SRR5630141", "SRR5634733", "SRR5629970", "SRR5629827", "SRR5630061", "SRR5629998", "SRR5630366", "SRR5630016", "SRR5630161", "SRR5630150", "SRR5634442", "SRR5630085", "SRR5634395", "SRR5634617", "SRR5630194", "SRR5630450", "SRR5634544", "SRR5634563", "SRR5634351", "SRR5629896", "SRR5629861", "SRR5629999", "SRR5629870", "SRR5634425", "SRR5634717", "SRR5634522", "SRR5630115", "SRR5630359", "SRR5634553", "SRR5634470", "SRR5630163", "SRR5634806", "SRR5630219", "SRR5630066", "SRR5630037", "SRR5630202", "SRR5634739", "SRR5634446", "SRR5634613",
                                "SRR5630285", "SRR5630350", "SRR5630183", "SRR5634365", "SRR5630160", "SRR5634705", "SRR5630394", "SRR5630310", "SRR5630088", "SRR5629866", "SRR5634417", "SRR5634809", "SRR5634797", "SRR5634615", "SRR5634698", "SRR5629994", "SRR5630313", "SRR5630292", "SRR5634487", "SRR5634367", "SRR5629892", "SRR5634557", "SRR5630393", "SRR5629808", "SRR5634461", "SRR5630340", "SRR5629865", "SRR5634600", "SRR5630346", "SRR5630317", "SRR5634629", "SRR5629807", "SRR5629965", "SRR5629986", "SRR5630275", "SRR5634778", "SRR5629947", "SRR5629996", "SRR5630314", "SRR5629977", "SRR5629955", "SRR5634702",
                                "SRR5630309", "SRR5630440", "SRR5629894", "SRR5629884", "SRR5634548", "SRR5630381", "SRR5634545", "SRR5634419", "SRR5634793", "SRR5634541", "SRR5630151", "SRR5630064", "SRR5630230", "SRR5634765", "SRR5634671", "SRR5634489", "SRR5634560", "SRR5634477", "SRR5630047", "SRR5634497", "SRR5634534", "SRR5630083", "SRR5630127", "SRR5634801", "SRR5634398", "SRR5634825", "SRR5630174", "SRR5634469", "SRR5634816", "SRR5634644", "SRR5634611", "SRR5630238", "SRR5634814", "SRR5630261", "SRR5634631", "SRR5629990", "SRR5629924", "SRR5630418", "SRR5629786", "SRR5630050", "SRR5629933", "SRR5629926",
                                "SRR5630036", "SRR5629833", "SRR5634355", "SRR5634536", "SRR5634400", "SRR5634766", "SRR5629985", "SRR5630018", "SRR5630100", "SRR5634396", "SRR5630193", "SRR5634363", "SRR5629873", "SRR5634599", "SRR5634562", "SRR5634808", "SRR5634697", "SRR5634732", "SRR5634389", "SRR5634734", "SRR5630071", "SRR5630404", "SRR5630042", "SRR5630414", "SRR5630057", "SRR5634353", "SRR5630020", "SRR5634589", "SRR5630315", "SRR5629905", "SRR5634486", "SRR5634371", "SRR5634762", "SRR5630300", "SRR5630114", "SRR5629890", "SRR5634699", "SRR5630006", "SRR5634715", "SRR5629930", "SRR5630185", "SRR5629904",
                                "SRR5630228", "SRR5630136", "SRR5634385", "SRR5634354", "SRR5630038", "SRR5630237", "SRR5630043", "SRR5634471", "SRR5630382", "SRR5630223", "SRR5629927", "SRR5629819", "SRR5634584", "SRR5634721", "SRR5634693", "SRR5630284", "SRR5634376", "SRR5634566", "SRR5629975", "SRR5630387", "SRR5629854", "SRR5630289", "SRR5629993", "SRR5629961", "SRR5634654", "SRR5630072", "SRR5630172", "SRR5630452", "SRR5630378", "SRR5629785", "SRR5630156", "SRR5634680", "SRR5630390", "SRR5630165", "SRR5630434", "SRR5634606", "SRR5630360", "SRR5629945", "SRR5634695", "SRR5630197", "SRR5630176", "SRR5630407",
                                "SRR5634724", "SRR5630427", "SRR5630221", "SRR5630014", "SRR5630312", "SRR5630069", "SRR5629799", "SRR5630116", "SRR5630142", "SRR5634651", "SRR5634625", "SRR5630417", "SRR5629846", "SRR5634768", "SRR5629917", "SRR5629949", "SRR5634433", "SRR5630216", "SRR5634737", "SRR5634683", "SRR5630117", "SRR5634610", "SRR5629781", "SRR5630372", "SRR5634571", "SRR5634653", "SRR5630240", "SRR5629914", "SRR5634451", "SRR5629888", "SRR5630386", "SRR5630166", "SRR5630297", "SRR5630222", "SRR5630262", "SRR5629937", "SRR5634462", "SRR5630099", "SRR5630206", "SRR5629940", "SRR5630054", "SRR5630027",
                                "SRR5630435", "SRR5634514", "SRR5634722", "SRR5634529", "SRR5634555", "SRR5630380", "SRR5630239", "SRR5629943", "SRR5629966", "SRR5630143", "SRR5634475", "SRR5629946", "SRR5629919", "SRR5634356", "SRR5630229", "SRR5634583", "SRR5634427", "SRR5630182", "SRR5634614", "SRR5629880", "SRR5630283", "SRR5630252", "SRR5630188", "SRR5634805", "SRR5634609", "SRR5630311", "SRR5630269", "SRR5629848", "SRR5630051", "SRR5634756", "SRR5630341", "SRR5634800", "SRR5629906", "SRR5634450", "SRR5630388", "SRR5634537", "SRR5634550", "SRR5634561", "SRR5629855", "SRR5630449", "SRR5629885", "SRR5634630",
                                "SRR5629916", "SRR5629821", "SRR5630144", "SRR5630323", "SRR5634552", "SRR5629988", "SRR5630416", "SRR5634435", "SRR5630272", "SRR5629958", "SRR5630287", "SRR5629922", "SRR5630245", "SRR5634807", "SRR5630044", "SRR5634413", "SRR5630135", "SRR5630203", "SRR5629836", "SRR5634803", "SRR5634748", "SRR5630003", "SRR5634727", "SRR5630247", "SRR5634817", "SRR5630338", "SRR5634374", "SRR5629803", "SRR5630426", "SRR5634635", "SRR5630327", "SRR5634509", "SRR5630062", "SRR5629881", "SRR5634622", "SRR5630091", "SRR5634524", "SRR5630139", "SRR5634812", "SRR5634517", "SRR5630095", "SRR5629925",
                                "SRR5630431", "SRR5634760", "SRR5634434", "SRR5630421", "SRR5630354", "SRR5634637", "SRR5634675", "SRR5630428", "SRR5630004", "SRR5629814", "SRR5630371", "SRR5634506", "SRR5634597", "SRR5629843", "SRR5634458", "SRR5634820", "SRR5634627", "SRR5629923", "SRR5630029", "SRR5630325", "SRR5630405", "SRR5634397", "SRR5629877", "SRR5630149", "SRR5634620", "SRR5630344", "SRR5629983", "SRR5629844", "SRR5634501", "SRR5630021", "SRR5630000", "SRR5630281", "SRR5634811", "SRR5629828", "SRR5630084", "SRR5634608", "SRR5634392", "SRR5629991", "SRR5634763", "SRR5634711", "SRR5634393", "SRR5630432",
                                "SRR5630212", "SRR5634559", "SRR5630308", "SRR5634703", "SRR5634582", "SRR5630010", "SRR5634691", "SRR5629850", "SRR5629956", "SRR5634679", "SRR5630126", "SRR5634769", "SRR5629942", "SRR5634782", "SRR5630196", "SRR5630347", "SRR5629842", "SRR5630352", "SRR5630210", "SRR5630145", "SRR5630296", "SRR5630022", "SRR5630363", "SRR5634496", "SRR5630302", "SRR5634569", "SRR5629812", "SRR5634429", "SRR5630195", "SRR5630007", "SRR5630097", "SRR5630181", "SRR5630412", "SRR5630220", "SRR5629995", "SRR5629899", "SRR5630201", "SRR5629817", "SRR5634432", "SRR5630430", "SRR5634556", "SRR5634771",
                                "SRR5630398", "SRR5629931", "SRR5630446", "SRR5634663", "SRR5634821", "SRR5629809", "SRR5634758", "SRR5630282", "SRR5634565", "SRR5634528", "SRR5629960", "SRR5634741", "SRR5630373", "SRR5630215", "SRR5630305", "SRR5630249", "SRR5630158", "SRR5630274", "SRR5629790", "SRR5634573", "SRR5634564", "SRR5629856", "SRR5629869", "SRR5629816", "SRR5634690", "SRR5629907", "SRR5630049", "SRR5629909", "SRR5629928", "SRR5634358", "SRR5634523", "SRR5634636", "SRR5630019", "SRR5634595", "SRR5629954", "SRR5629795", "SRR5629813", "SRR5634420", "SRR5629957", "SRR5634770", "SRR5630015", "SRR5630108",
                                "SRR5630343", "SRR5634799", "SRR5630140", "SRR5630307", "SRR5629887", "SRR5634361", "SRR5634440", "SRR5630107", "SRR5630092", "SRR5630392", "SRR5629910", "SRR5634407", "SRR5634531", "SRR5630420", "SRR5634692", "SRR5634641", "SRR5630024", "SRR5629793", "SRR5634672", "SRR5634525", "SRR5630442", "SRR5629918", "SRR5634467", "SRR5630001", "SRR5634488", "SRR5634377", "SRR5634776", "SRR5634689", "SRR5630096", "SRR5630133", "SRR5630400", "SRR5629883", "SRR5630298", "SRR5630110", "SRR5629953", "SRR5630189", "SRR5634685", "SRR5630332", "SRR5630253", "SRR5630213", "SRR5634503", "SRR5630187",
                                "SRR5634479", "SRR5629980", "SRR5630391", "SRR5634480", "SRR5630451", "SRR5634439", "SRR5629804", "SRR5629889", "SRR5630234", "SRR5630217", "SRR5630235", "SRR5634700", "SRR5634794", "SRR5634455", "SRR5634619", "SRR5630339", "SRR5630377", "SRR5629818", "SRR5634568", "SRR5630306", "SRR5630005", "SRR5634718", "SRR5634780", "SRR5634394", "SRR5630277", "SRR5634510", "SRR5634381", "SRR5629853", "SRR5629798", "SRR5634621", "SRR5634382", "SRR5630218", "SRR5630288", "SRR5634604", "SRR5634728", "SRR5630279", "SRR5634359", "SRR5630002", "SRR5630443", "SRR5634538", "SRR5629948", "SRR5630134",
                                "SRR5629811", "SRR5629992", "SRR5634472", "SRR5634447", "SRR5634682", "SRR5634370", "SRR5634647", "SRR5634701", "SRR5634493", "SRR5634468", "SRR5634605", "SRR5634572", "SRR5630098", "SRR5634749", "SRR5630319", "SRR5629939", "SRR5630316", "SRR5634624", "SRR5630367"))

# prepare heatmap
p1 = ggplot(cov_file) +
  geom_tile(data = subset(cov_file, cov_file$V7 == "A"),
            aes(x = V3, y = V1, fill = V6)) +
  scale_fill_gradientn(na.value = "grey85", limits = c(0, 2),
                       colours = c("#053061", "#0F437B", "#195696", "#2369AD", "#2F79B5", "#3B89BE", "#4E9AC6", "#6AACD0", "#86BDDA", "#9FCBE1", "#B6D7E8", "#CCE2EE", "#DBEAF2", "#E9F0F4", "#F7F7F7", "#F9ECE5", "#FBE3D4", "#FCD7C2", "#F9C3A9", "#F5B090", "#EF9B7A", "#E58267", "#DA6954", "#CE5045", "#C13639", "#B41D2D", "#9C1127", "#810823", "#67001F"),
                       breaks=c(0, 1, 2)) +
  scale_x_continuous(labels = comma) +
  coord_cartesian(expand = FALSE) +
  labs(fill = "Maltose growth",
       x = " ",
       y = "Segregants") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.position = "left",
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.75)) +
  guides(fill = guide_colourbar(ticks.colour = "black",
                                frame.colour = "black",
                                frame.linewidth = 1))

p2 = ggplot(cov_file) +
  geom_tile(data =  subset(cov_file, cov_file$V7 == "B" & cov_file$V8 == "R"),
            aes(x = V3, y = V1, fill = V6, alpha = V9 )) +
  scale_fill_gradientn(na.value = "red", limits = c(0, 2),
                       colours = c("#FFFFCC", "#FFF8BC", "#FFF1AC", "#FEEB9C", "#FEE38C", "#FEDC7D", "#FED16E", "#FEC35F", "#FEB54F", "#FDA747", "#FD9A41", "#FD8D3C", "#FC7635", "#FC5F2E", "#F94928", "#F03623", "#E7231E", "#DC151D", "#CE0B21", "#C00225", "#AC0026", "#960026", "#800026"),
                       values = c(0, 0.00000001, 1, 2),
                       breaks = c(0, 1, 2)) +
  new_scale_fill() +
  geom_tile(data =  subset(cov_file, cov_file$V7 == "B" & cov_file$V8 == "Y"),
            aes(x = V3, y = V1, fill = V6, alpha = V9 )) +
  scale_fill_gradientn(na.value = "black", limits = c(0, 2),
                       colours = c("white", "grey95", "grey90", "grey85", "grey80", "grey75", "grey70", "grey65", "grey60", "grey55", "grey50", "grey45", "grey40", "grey35", "grey30", "grey25", "grey20", "grey15", "grey10", "black"),
                       values = c(0, 0.00000001, 1, 2),
                       breaks=c(0, 1, 2)) +
  scale_x_continuous(labels = comma) +
  coord_cartesian(expand = FALSE) +
  labs(fill = "Coverage",
       x = " ") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.75)) +
  guides(fill = guide_colourbar(ticks.colour = "black",
                                frame.colour = "black",
                                frame.linewidth = 1))

# Figure_S5A
Figure_S5A <- grid.arrange(p1, p2, nrow = 1,
             layout_matrix = rbind(c(1, 2)),
             widths = c(0.12, 0.90))

# Figure_S5B
MALgenes$IMA1_ALLELE = factor(MALgenes$IMA1_ALLELE, levels=c("R", "Y"))
mu = ddply(MALgenes, "IMA1_ALLELE", summarise,
           grp.mean=mean(Maltose_growth))
Figure_S5B <- 
ggplot(MALgenes, aes(x = Maltose_growth, fill = IMA1_ALLELE, color =
                       IMA1_ALLELE)) +
  geom_histogram(bins = 20, position = "identity", alpha = 0.5) +
  geom_vline(data = mu, aes(xintercept = grp.mean), linetype = "dashed",
             colour = c("#FDA747", "grey50"), size = 1.5) +
  scale_fill_manual(values = c("#FDA747", "grey50")) +
  scale_color_manual(values = c("#FDA747", "grey50")) +
  guides(colour = FALSE) +
  ylim(0, 127) +
  coord_cartesian(expand = FALSE) +
  labs(x = "Maltose growth",
       y = "Number of segregants",
       fill = "Genotype") +
  theme(axis.ticks = element_line(size = 0.75),
        axis.ticks.length = unit(10, "points"),
        axis.title = element_text(size = 48),
        axis.text = element_text(size = 36),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 24),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0.95),
        legend.background = element_rect(color = "black", fill = "white"),
        panel.border = element_rect(size = 1.25, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_blank())



#### Figure S6 ####
SNV_SWAP_sup <- #Supplementary data Figure S6B#

trait_dotplot_urk1 <- 
  SNV_SWAP_sup %>% 
  group_by(date, trait, class) %>%
  mutate(control_mean = mean(value[background == "WT"], na.rm = T), 
         value_norm = value/control_mean) 

trait_dotplot_urk1$trait <- as.factor(trait_dotplot_urk1$trait)

stat.test_urk1 <- trait_dotplot_urk1 %>% 
  group_by(trait, class) %>% 
  t_test(value_norm~background)%>% 
  adjust_pvalue(method = 'bonferroni')%>%
  adjust_pvalue() %>%
  add_significance("p.adj")

stat.test_urk1 <- stat.test_urk1 %>% add_xy_position(fun = "mean", x = "background")

#Figure_S6B

Figure_S6B <- ggplot(trait_dotplot_urk1, aes(sample,value_norm))+
  facet_grid(trait~class, scales="free_x", switch = "y")+
  geom_hline(yintercept = 1, color="grey", size=0.5, linetype="dashed")+
  stat_summary(fun.data=mean_sdl, geom="errorbar", aes(color=background), width=0, lwd=0.5, position=position_dodge(.7))+
  stat_summary(fun.y=mean, geom="point", aes(color=background), size=1.5) +
  theme_cowplot()+
  theme_classic(base_size = 16)+
  theme(axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank(),                  
        axis.line =  element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5),
        panel.spacing = unit(1, "mm"),
        panel.border = element_rect(size=0.5, fill="transparent"), 
        strip.background = element_blank(),
        strip.text = element_text(size=14, color = "grey45"),
        strip.text.y.left = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size=12),
        legend.position = c(0.8,0.85))+ 
  scale_color_manual(labels=c("SNV swap","WT"),
                     values = c("tomato3","grey40"))+
  xlab("Variant genotype")+
  ylab("Relative isobutanol conc. (STD)")+
  scale_x_discrete(labels=c("SNV9801" = " ",
                            "SNV9802" = " ",
                            "WT" = " "
  ))

urk1_del <- #Supplementary data Figure S6C#

urk1_del$class <- factor(urk1_del$class, levels = c("S.boulardii", "Ethanol Red", "CEN.PK"))

stat.test <- compare_means(
  value ~ background, data = urk1_del%>%filter(!class %in% c("RM", "YJM")), group.by = "class",
  method = "t.test", ref.group = "WT"
)

urk1 <- ggbarplot(urk1_del, x = "class", y = "value",
                  fill = "background",
                  add = "mean_sd", add.params = list(group = "background"),
                  position = position_dodge(0.8),
                  alpha=0.5)

p_urk1_KO <- urk1 + stat_pvalue_manual(
  stat.test, x = "class", y.position = 22,
  label = "p.signif",
  position = position_dodge(0.8), 
  hide.ns = T
)+
  theme_classic(base_size = 16) +
  theme(
    legend.position = c(0.8,0.8),
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.line =  element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, color = "grey45"),
    plot.title = element_text(size=16),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  ) +
  scale_fill_manual(
    labels=c("urk1 KO", "WT"),
    values = c("tomato3", "grey40"))+
  ylab("isobutanol (ppm)") +
  scale_y_continuous(limits = c(0, max(urk1_del$value) * 0.8))+
  xlab("")