library(tidyverse)

dir <- "~/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/G4PromFinder_outputs/"

list.files(dir, pattern = "coordinates.txt", full.names = T)%>%
  data.table::fread(header = F, fill = TRUE) -> df

df %>%
  mutate_at(c("V4","V5"), as.numeric) %>%
  as_tibble() %>%
  drop_na() -> df
  # readr::parse_number()


df %>%
  select(V1, V4, V5) %>%
  write.table(., sep = "\t",file = paste0(dir, "coordinates.bed"), quote = F, 
              col.names = F, row.names = F)

df %>%
  filter(V1 %in% "Chromosome")


#
mtd_file <- "/Users/cigom/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/genome_metadata.txt"

mtd <- read.delim(mtd_file, sep = "\t", header = F) %>% as_tibble() %>% distinct()
names(mtd) <- c("genome", "Index")


mtd %>% filter(Index %in% "Chromosome") %>% pull(genome) -> rm_genomes



dir <- "/Users/cigom/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/genomes/gtf/"

file <- list.files(path = dir, pattern = "intersect_bed_f50.txt", full.names = T)

# filter(grepl("ASM148338v1", V9))

read.delim(file, sep = "\t", header = F) %>%
  as_tibble() %>%
  rename(V1 = "Index") %>% 
  left_join(., mtd,by = "Index") %>%
  distinct() -> df_viz

df_viz %>%
  filter(!genome %in% rm_genomes) %>%
  group_by(Index, V3, genome) %>%
  tally(sort = T) %>%
  mutate(sp = genome) %>%
  separate(sp, c("Genus", "sp"), "_") %>%
  arrange(sp) %>% 
  group_by(genome) %>%
  mutate(pct = n / sum(n) * 100, sp = as.integer(as.factor(sp))) %>%
  mutate(genome = forcats::fct_reorder(genome, sp)) %>%
  # mutate(genome = forcats::fct_reorder(genome, n)) %>%
  ggplot() +
  geom_col(aes(y = pct, x = genome, fill = V3)) +
  coord_flip() +
  ggsci::scale_fill_aaas(name = "") +
  labs(y = "Prediction (Normalized)", x = "Genome") +
  theme_bw(base_size = 16, base_family = "GillSans") -> p1
  
dir <- "~/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/G4PromFinder_outputs"

ggsave(p1, filename = "S2_figure.png", path = dir, 
       width = 10, height = 8)  
  
