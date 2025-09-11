library(rtracklayer)

library(tidyverse)
library(dplyr)

library(ggplot2)
library(patchwork)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
ucsc_id <- args[1]
work_dir <- args[2]

# Upload alignment
alig <- read.csv(paste(work_dir, '/', ucsc_id, '_genome_dir/', ucsc_id, '_mapping.txt', sep = ''), 
                 sep = '\t', col.names = c('ID', 'count', 'chr', 'place', c(5:21)))

# Choose the genome and connection
session <- browserSession("UCSC")
genome(session) <- ucsc_id

# Input chromosome info
query <- ucscTableQuery(session, 
                        table = "ucscToRefSeq",
                        range = ucsc_id)
chrs <- getTable(query) %>% 
  filter(startsWith(name, 'NC_')) %>%     # choose chromosome only
  filter(chrom != 'chrM') %>%             # excluding mitochondrial chromosome
  select(c('chrom', 'chromEnd', 'name')) %>% 
  rename(lens = chromEnd)

# Input centromere info
if (ucsc_id == 'hg38') {
  query <- ucscTableQuery(session, 
                          table = "centromeres",
                          range = ucsc_id)
  centromeres <- getTable(query) %>%
    group_by(chrom) %>%
    summarize(
      centStart = min(chromStart),
      centEnd = max(chromEnd),
      size = max(chromEnd) - min(chromStart)
    ) %>%
    ungroup()
} else {
  query <- ucscTableQuery(session, 
                          table = "gap")
                          #range = ucsc_id)
  centromeres <- getTable(query) %>% 
    filter(type == 'centromere') %>%               # choose centromeres
    select(c('chrom', 'chromStart', 'chromEnd', 'size')) %>% 
    rename(centStart = chromStart, centEnd = chromEnd)
}

# Join chromosome and centromeres info tables
chromosome_info <- left_join(chrs, centromeres, by = 'chrom') %>% 
  replace(is.na(.), 0)


# Pull normalize coefficient
norm_dens <- c()
for (i in c(1:nrow(chromosome_info))){
  alig_chr <- alig %>% 
    filter(chr == chromosome_info$chrom[i]) %>% 
    select(chr, place) %>% 
    mutate(cnt = 1)
  
  seps <- chromosome_info$lens[i]/50
  tics <- seq(0, chromosome_info$lens[i], by = seps)
  dens <- c()
  for (k in tics){
    dens <- append(dens, sum(alig_chr %>% 
                               filter(place >= k & place < k+seps) %>% 
                               pull(cnt)))}
  #print(max(dens))
  norm_dens <- append(norm_dens, max(dens))}

norm_coef <- max(norm_dens)


# Visualize chromosome FISH-colouring
for (i in c(1:nrow(chromosome_info))){
  alig_chr <- alig %>% 
    filter(chr == chromosome_info$chrom[i]) %>% 
    select(chr, place) %>% 
    mutate(cnt = 1)
  
  seps = chromosome_info$lens[i]/50
  tics <- seq(0, chromosome_info$lens[i], by = seps)
  dens <- c()
  for (k in tics){
    dens <- append(dens, sum(alig_chr %>% 
                               filter(place >= k & place < k+seps) %>% 
                               pull(cnt)/norm_coef))}
  
  chr_density <- as_tibble(tics) %>% 
    mutate(dens = dens)
  
  chr_vis <- ggplot(chr_density) +
    geom_bar(aes(x = value, y = 1, fill = dens),
             stat = "identity", position = "dodge") +
    scale_x_continuous(name = 'position', breaks = seq(0, chromosome_info$lens[i], by = 10000000)) +
    scale_fill_gradientn(
      colours = c("darkblue", 'springgreen3', "green1"),
      rescaler = ~ scales::rescale_mid(.x, mid = 0.3),
      limits = c(0, 1)) +
    annotate("rect", xmin = chromosome_info$centStart[i], xmax = chromosome_info$centEnd[i], ymin = 0, ymax = 1,
             alpha = 1, fill = "black") +
    labs(title = paste('chromosome', str_split_i(chromosome_info$chrom[i], 'chr', -1), '\n', 'N =', nrow(alig_chr))) +
    theme_bw()
  
  ggsave(paste(work_dir, '/plots/', ucsc_id, "_coloring_", chromosome_info$chrom[i], ".pdf", sep = ''),
         chr_vis, width=6, height=2, units="in", scale=3)
}


# Creating main plot
chr_info <- chromosome_info %>%
  mutate(chrEnd = lens, .before = lens) %>% 
  mutate(chrStart = 0, .before = chrEnd) %>% 
  mutate(id = as.numeric(str_split_i(chrom, 'chr', -1))) %>% 
  arrange(id, chrom)

chr_info[c('chrStart', 'chrEnd', 'centStart', 'centEnd')] <- 
  chr_info[c('chrStart', 'chrEnd', 'centStart', 'centEnd')] - chr_info$centStart

plot_list <- list()
for (i in c(1:nrow(chr_info))){ 

  print(chr_info$chrom[i])
  alig_chr <- alig %>% 
    filter(chr == chr_info$chrom[i]) %>% 
    select(chr, place) %>% 
    mutate(cnt = 1)
  
  seps = chr_info$lens[i]/50
  tics <- seq(0, chr_info$lens[i], by = seps)
  dens <- c()
  for (k in tics){
    dens <- append(dens, sum(alig_chr %>% 
                               filter(place >= k & place < k+seps) %>% 
                               pull(cnt)/norm_coef))}
  chr_density <- as_tibble(tics) %>% 
    mutate(dens = dens)
  
  chr_vis <- ggplot(chr_density) +
    geom_bar(aes(x = value + chr_info$chrStart[i], y = 1, fill = dens),
             stat = "identity", position = "dodge") +
    scale_fill_gradientn(
      colours = c("darkblue", 'springgreen3', "green1"),
      rescaler = ~ scales::rescale_mid(.x, mid = 0.3),
      limits = c(0, 1)) +
    annotate("rect", alpha = 1, fill = "black",
             xmin = chr_info$centStart[i], 
             xmax = chr_info$centEnd[i], 
             ymin = -0.3, ymax = 1.3) +
    coord_flip() +
    scale_x_reverse() +
    xlim(max(chr_info$chrEnd) +10000000, min(chr_info$chrStart) -10000000) +
    theme_bw() +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          panel.background=element_blank(), panel.border=element_blank(),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          plot.background=element_blank(), legend.position="none")
  chr_vis <- grid.arrange(chr_vis, bottom = str_split_i(chr_info$chrom[i], 'chr', -1))
  plot_list[[i]] <- chr_vis}

main_plot <- wrap_plots(plot_list, ncol = nrow(chromosome_info)/2)
ggsave(paste(work_dir, '/plots/', ucsc_id, "_main.pdf", sep = ''),
       main_plot, width=5, height=4, units="in", scale=3)
