### The BEGINNING ~~~~~
##
# ~ Plots TreeSparrowGenomics--Stats | By George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(scales, extrafont, dplyr, grid, lubridate, cowplot, egg, tidyverse, stringr, reshape, plotly)
library(lemon)


# Loads datasets ~
Stats <- read.csv("TreeSparrow--Stats.txt", sep = "\t", header = TRUE,
                    stringsAsFactors = TRUE)
Adaptors <- read.table("TreeSparrow--Stats_Adaptors.txt", sep = "\t", header = FALSE,
                       stringsAsFactors = FALSE); colnames(Adaptors) <- c("Sample_ID", "reads_adaptors")


# Combines DFs ~
fulldf <- merge(Stats, Adaptors, by = "Sample_ID")


# Gets total_reads ~
fulldf$total_reads <-
       fulldf$seq_reads_pairs * 2


# Gets percentages ~
fulldf$percentage_retained_reads <-
       fulldf$seq_retained_reads * 100 / fulldf$total_reads
fulldf$reads_adaptors <-
  fulldf$reads_adaptors * 100 / fulldf$total_reads


# Fixes percentages ~
fulldf$hits_raw_frac_HouseSparrow <- fulldf$hits_raw_frac_HouseSparrow * 100
fulldf$hits_clonality_HouseSparrow <- fulldf$hits_clonality_HouseSparrow  * 100
fulldf$hits_unique_frac_HouseSparrow <- fulldf$hits_unique_frac_HouseSparrow * 100
fulldf$hits_raw_frac_TreeSparrow <- fulldf$hits_raw_frac_TreeSparrow * 100
fulldf$hits_clonality_TreeSparrow <- fulldf$hits_clonality_TreeSparrow  * 100
fulldf$hits_unique_frac_TreeSparrow <- fulldf$hits_unique_frac_TreeSparrow * 100


# Expands PCA_Annot by adding Population ~
fulldf$Population <- ifelse(grepl("Aizawl", fulldf$Sample_ID), "Aizawl",
                     ifelse(grepl("Eastermar", fulldf$Sample_ID), "Eastermar",
                     ifelse(grepl("PMON2002KOR0002U", fulldf$Sample_ID), "Gyeonggi",
                     ifelse(grepl("TreeSparrow_01", fulldf$Sample_ID), "Corsica",
                     ifelse(grepl("PHL", fulldf$Sample_ID), "Bulacan",
                     ifelse(grepl("TreeSparrow", fulldf$Sample_ID), "Malta",
                     ifelse(grepl("Zhabagly", fulldf$Sample_ID), "Zhabagly", "Japan")))))))


# Fixes Japanese population ~
fulldf$Population <- ifelse(fulldf$Sample_ID %in% c("PMON2011JPN0001U", "PMON2011JPN0002U", "PMON2011JPN0003U", "PMON2011JPN0004U", "PMON2011JPN0005U",
                                                    "PMON2011JPN0006U", "PMON2011JPN0007U"), "Fukushima",
                     ifelse(fulldf$Sample_ID %in% c("PMON1981JPN0001M", "PMON2000JPN0001U", "PMON2001JPN0001F", "PMON2004JPN0001M", "PMON2006JPN0002F",
                                                    "PMON2007JPN0001M", "PMONXXXXJPN0001M"), "Hokkaido",
                     ifelse(fulldf$Sample_ID %in% c("PMON2014JPN0002U", "PMON2016JPN0001U", "PMON2016JPN0002M", "PMON2016JPN0003F", "PMON2016JPN0004U"), "Ibaraki",
                     ifelse(fulldf$Sample_ID %in% c("PMON2008JPN0002U", "PMON2008JPN0003M"), "Kagoshima",
                     ifelse(fulldf$Sample_ID %in% c("PMON2009JPN0001M", "PMON2009JPN0002M"), "Okinawa",
                     ifelse(fulldf$Sample_ID %in% c("PMON1995JPN0001U", "PMON1995JPN0002U", "PMON2005JPN0001U"), "Osaka", fulldf$Population))))))


# Fixes Tokyo population ~
fulldf$Population <- ifelse(grepl("Japan", fulldf$Population), "Tokyo", fulldf$Population)


# Reorders Population ~
fulldf$Population <- factor(fulldf$Population, ordered = T,
                                   levels = c("Eastermar",
                                              "Malta",
                                              "Corsica",
                                              "Zhabagly",
                                              "Bulacan",
                                              "Kagoshima",
                                              "Okinawa",
                                              "Gyeonggi",
                                              "Hokkaido",
                                              "Fukushima",
                                              "Ibaraki",
                                              "Tokyo",
                                              "Osaka",
                                              "Aizawl"))


# Expands PCA_Annot by adding Country ~
fulldf$Country <- ifelse(fulldf$Population %in% c("Eastermar"), "Netherlands",
                         ifelse(fulldf$Population %in% c("Corsica"), "France",
                         ifelse(fulldf$Population %in% c("Malta"), "Malta",
                         ifelse(fulldf$Population %in% c("Zhabagly"), "Kazakhstan",
                         ifelse(fulldf$Population %in% c("Aizawl"), "India",
                         ifelse(fulldf$Population %in% c("Bulacan"), "Philippines",
                         ifelse(fulldf$Population %in% c("Kagoshima", "Okinawa", "Hokkaido", "Fukushima",  "Ibaraki", "Tokyo", "Osaka"), "Japan",
                         ifelse(fulldf$Population %in% c("Gyeonggi"), "South Korea", fulldf$Population))))))))


# Reorders Population ~
fulldf$Country <- factor(fulldf$Country, ordered = T,
                         levels = c("Netherlands",
                                    "France",
                                    "Malta",
                                    "Kazakhstan",
                                    "India",
                                    "Philippines",
                                    "South Korea",
                                    "Japan"))

# Cleans DF ~
fulldf <- fulldf %>%
          select(Sample_ID, Country, Population, total_reads, reads_adaptors, percentage_retained_reads,
                 hits_raw_frac_HouseSparrow, hits_clonality_HouseSparrow, hits_unique_frac_HouseSparrow, hits_coverage_HouseSparrow,
                 hits_raw_frac_TreeSparrow, hits_clonality_TreeSparrow, hits_unique_frac_TreeSparrow, hits_coverage_TreeSparrow)


# Converts DF from wide into long ~
fulldfUp <- gather(fulldf, Stat, Value, "total_reads", "reads_adaptors", "percentage_retained_reads",
                                        "hits_raw_frac_HouseSparrow", "hits_raw_frac_TreeSparrow",
                                        "hits_clonality_HouseSparrow", "hits_clonality_TreeSparrow",
                                        "hits_unique_frac_HouseSparrow", "hits_coverage_HouseSparrow",
                                        "hits_unique_frac_TreeSparrow", "hits_coverage_TreeSparrow")


# Expands fulldf by adding REF ~
fulldfUp$REF <- ifelse(grepl("HouseSparrow", fulldfUp$Stat), "House",
                ifelse(grepl("TreeSparrow", fulldfUp$Stat), "Tree", NA))


# Corrects Afastamento Motivos ~
levels(fulldfUp$Stat <- sub("_HouseSparrow", "", fulldfUp$Stat))
levels(fulldfUp$Stat <- sub("_TreeSparrow", "", fulldfUp$Stat))


# Reorders Stat ~
fulldfUp$Stat <- factor(fulldfUp$Stat, ordered = T,
                        levels = c("hits_coverage",
                                   "hits_unique_frac",
                                   "hits_clonality",
                                   "hits_raw_frac",
                                   "percentage_retained_reads",
                                   "reads_adaptors",
                                   "total_reads"))


# Corrects facet labels ~
ylabels <- c("total_reads" = "# of Reads",
             "reads_adaptors" = "% of Reads With Adaptors",
             "percentage_retained_reads" = "% of Reads Retained",
             "hits_raw_frac" = "% of Mapped Reads",
             "hits_clonality" = "% of Clonality",
             "hits_unique_frac" = "% of Uniquely Mapped Reads",
             "hits_coverage" = "Mean Depth")


# Custom y-axis labels ~
plot_index_labels <- 0
labels_fun <- function(z) {
  plot_index_labels <<- plot_index_labels + 1L
  switch(plot_index_labels,
         scales::label_number(accuracy = 1, suffix = "X")(z),
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(z),
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(z),
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(z),
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(z),
         scales::label_percent(accuracy = 1, scale = 1, big.mark = "")(z),
         scales::label_number(accuracy = 1, scale = 1/1000000, big.mark = "", suffix = "M")(z))}


# Creates the panel ~
TreeSparrowGenomics_Stat <- 
 ggplot(fulldfUp, aes(x = Country, y = Value, fill = REF)) +
  #geom_violin(data = fulldfUp, aes(x = Population, y = Value),
  #            fill = "#ffffff", colour = "#000000", show.legend = FALSE, alpha = .9, size = .45, width = 1) +
  geom_boxplot(position = "dodge", outlier.shape = NA, width = .5, lwd = .25, colour = "#000000", alpha = .7) +
  #stat_summary(data = fulldfUp, aes(x = Population, y = Value),  
  #             fun = mean, geom = "point", shape = 21, size = 2.5, alpha = 1, colour = "#000000", fill = "#df65b0") +
  scale_fill_manual(values = c("#1E90FF", "#856D46"), na.translate = FALSE) + 
  facet_rep_grid(Stat ~. , scales = "free", labeller = labeller(Stat = ylabels)) +
  scale_y_continuous(#limits = limits_fun,
                     #breaks = breaks_fun,
                     labels = labels_fun) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major = element_line(color = "#E5E7E9", linetype = "dashed", linewidth = .005),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.spacing.y = unit(1, "cm"),
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "#000000", size = 10, face = "bold", angle = 45, vjust = 1, hjust = 1),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(color = "#000000", size = 8, face = "bold"),
        axis.ticks.x = element_line(color = "#000000", linewidth = .3),
        axis.ticks.y = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 8, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(fill = guide_legend(title = "Reference", title.theme = element_text(size = 17, face = "bold"),
                             label.theme = element_text(size = 16),
                             override.aes = list(starshape = 15, size = 4.25, alpha = .7, starstroke = .15), nrow = 1, order = 1),
         starshape = "none",
         colour = "none")


# Saves the panel ~
ggsave(TreeSparrowGenomics_Stat, file = "TreeSparrowGenomics--Stats_Country.pdf",
       device = cairo_pdf, width = 12, height = 16, scale = 1, dpi = 600)
ggsave(TreeSparrowGenomics_Stat, file = "TreeSparrowGenomics--Stats_Country.jpeg",
       width = 12, height = 16, scale = 1, dpi = 600)


#
## 
### The END ~~~~~


# Custom y-axis limits ~
limits_fun <- function(x){
  limitVal <- max(x)
  print(x)
  if (limitVal < 60){
    c(0, 60)}
  else if (limitVal < 100){
    c(90, 100)}
  else { 
    c(70000000, 227746600)}}


# Custom y-axis breaks ~
breaks_fun <- function(y){
  caseVal <- min(y)
  print(y)
  if (caseVal < 100 & caseVal > 50){
    seq(20, 60, by = 10)}
  else if (caseVal < 100){
    seq(92, 100, by = 2)}}