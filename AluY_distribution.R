library(tidyverse)
library(dplyr)
library(jsonlite)

library(ggplot2)
library(RColorBrewer)

# Upload alignment
alig <- read.csv('/home/vera_kanaeva/Vera/lab/rna_seq/bogomazova_task/mmul_view_sam.txt', 
                 sep = '\t', col.names = c('ID', 'count', 'chr', 'place', c(5:21)))

# Input chromosome info
data <- stream_in(file("/home/kanaevavera/Vera/lab/rna_seq/bogomazova_task/sequence_report.jsonl"))
df <- as.data.frame(data) %>% 
  filter(role == 'assembled-molecule') %>% 
  select(c('chrName', 'refseqAccession', 'length'))

# Input centromere info
...

chr_info <- as_tibble(chrs) %>% 
  mutate(ids = ids, lens = lens, gap_start = gap_start, gap_end = gap_end)


# Pull normalize coefficient
norm_dens <- c()
for (i in c(1:24)){
  alig_chr18 <- alig %>% 
    filter(chr == chr_info$value[i]) %>% 
    select(chr, place) %>% 
    mutate(cnt = 1)
  
  print(chr_info$ids[i])
  seps = chr_info$lens[i]/50
  tics <- seq(0, chr_info$lens[i], by = seps)
  dens <- c()
  for (k in tics){
    #print(i)
    dens <- append(dens, sum(alig_chr18 %>% 
                               filter(place >= k & place < k+seps) %>% 
                               pull(cnt)))
  }
  print(max(dens))
  norm_dens <- append(norm_dens, max(dens))
}
max(norm_dens)

# Visualize chromosome FISH-colouring
for (i in c(1:24)){
  print(chr_info$ids[i])
  alig_chr18 <- alig %>% 
    filter(chr == chr_info$value[i]) %>% 
    select(chr, place) %>% 
    mutate(cnt = 1)
  
  seps = chr_info$lens[i]/50
  tics <- seq(0, chr_info$lens[i], by = seps)
  dens <- c()
  for (k in tics){
    #print(i)
    dens <- append(dens, sum(alig_chr18 %>% 
                               filter(place >= k & place < k+seps) %>% 
                               pull(cnt)/15))
  }
  #print(dens)
  
  chr18_density <- as_tibble(tics) %>% 
    mutate(dens = dens)
  
  chr_vis <- ggplot(chr18_density) +
    geom_bar(aes(x = value, y = 1, fill = dens),
             stat = "identity", position = "dodge") +
    scale_x_continuous(name = 'position', breaks = seq(0, chr_info$lens[i], by = 10000000)) +
    scale_fill_gradientn(
      colours = c("darkblue", 'springgreen3', "green1"),
      rescaler = ~ scales::rescale_mid(.x, mid = 0.3),
      limits = c(0, 1)) +
    annotate("rect", xmin = chr_info$gap_start[i], xmax = chr_info$gap_end[i], ymin = 0, ymax = 1,
             alpha = 1, fill = "black") +
    labs(title = paste('chromosome', ids[i], '\n', 'N =', nrow(alig_chr18))) +
    theme_bw()
  
  ggsave(paste("/home/vera_kanaeva/Vera/lab/rna_seq/bogomazova_task/alignments/grch38_coloring_chr", ids[i], ".pdf", sep = ''),
         chr_vis, width=6, height=2, units="in", scale=3)
}