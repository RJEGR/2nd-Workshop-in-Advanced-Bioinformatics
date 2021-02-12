library(tidyverse)

dir <- "~/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/G4PromFinder_outputs/coordinates_files/"

about_files_1 <- list.files(dir, pattern = "about.txt", full.names = T)


about_summary <- function(about_txt) {
  id <- unlist(strsplit(basename(about_txt),'[.]'))[1]
  
  about <- readLines(about_txt)
  about <- about[grepl('Total number of prediction', about)] 
  Predicted <- unlist(strsplit(as.character(unique(about)), ":"))[2]
  data.frame(id, Predicted = as.double(Predicted))
  
}

lapply(about_files_1, about_summary) %>%
  do.call(rbind, .) -> df1

# random genome results

dir <- "~/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/G4PromFinder_outputs/coordinates_files_shuffled_kmer2/"


about_files_2 <- list.files(dir, pattern = "about.txt", full.names = T)


about_summary(about_files_1[1])
about_summary(about_files_2[1])

lapply(about_files_2, about_summary) %>%
  do.call(rbind, .) %>%
  arrange(match(id, df1$id)) -> df2

# Sanity check


identical(df1$id, df2$id)

df <- data.frame(id = df1$id, 
                 genome = df1$Predicted, 
                 shuffled = df2$Predicted) %>%
  pivot_longer(-id, names_to = "group", values_to = "n_predict") %>%
  group_by(group) %>%
  mutate(Rank = rank(n_predict)) %>%
  mutate(id = forcats::fct_reorder(id, Rank))

library(ggplot2)

df %>%
  ggplot(aes(x = id, y = n_predict)) +
  geom_bar(aes(fill = group), position="dodge", stat="identity") +
  labs(y = "# Predictions with G4PromFinder", x = "") +
  coord_flip() + 
  # facet_grid(~ group) +
  scale_fill_manual("", values = c("#3182bd", "#de2d26"))

# dataViz coords ----

dir <- "~/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/G4PromFinder_outputs/coordinates_files/"

coor_files_1 <- list.files(dir, pattern = "coordinates.txt", full.names = T)

load_coord <- function(coor_files) {
  x <- data.table::fread(coor_files, fill = T)
  id <- unlist(strsplit(basename(coor_files),'[.]'))[1]
  data.frame(x, id)
}

lapply(coor_files_1, load_coord) %>%
  do.call(rbind, .) %>% as_tibble() -> df1

df1 %>%
  drop_na() %>%
  group_by(Strand, id) %>%
  tally(sort = TRUE)  %>%
  mutate(group = "Genome") -> df1_summ

dir <- "~/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/G4PromFinder_outputs/coordinates_files_shuffled_kmer2/"

coor_files_2 <- list.files(dir, pattern = "coordinates.txt", full.names = T)

lapply(coor_files_2, load_coord) %>%
  do.call(rbind, .) %>% as_tibble() -> df2

df2 %>%
  drop_na() %>%
  group_by(Strand, id) %>%
  tally(sort = TRUE) %>%
  mutate(group = "Genome (shuffled)") %>%
  arrange(match(id, df1_summ$id)) -> df2_summ

# sanity check
identical(df1_summ$id, df1_summ$id)

rbind(df1_summ, df2_summ) %>%
  group_by(Strand, group) %>%
  mutate(Rank = rank(n)) %>%
  mutate(id = forcats::fct_reorder(id, Rank)) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = id, y = n)) +
  geom_bar(aes(fill = group), position="dodge", stat="identity") +
  labs(y = "# Predictions with G4PromFinder", x = "") +
  coord_flip() + 
  facet_grid(~ Strand) +
  scale_fill_manual("", values = c("#3182bd", "#de2d26")) +
  theme_bw(base_size = 16, base_family = "GillSans") -> p1

# dir <- "~/raw_data/G4PromFinder_outputs/"

ggsave(p1, filename = "G4PromFinder.png", path = dir, 
       width = 10, height = 14)

raf <- function(x) x/sum(x) * 100

rbind(df1_summ, df2_summ) %>%
  group_by(Strand, group) %>%
  mutate(Rank = rank(n)) %>%
  mutate(id = forcats::fct_reorder(id, Rank)) %>%
  group_by(id, group) %>%
  mutate(pct = round(raf(n), digits = 2)) %>%
  # mutate(sep = id) %>%
  # separate(sep, c("Genus", "sp"), "_") %>%
  # mutate(lPos = cumsum(pct) - pct/2) %>%
  ggplot(aes(x = id, y = n)) +
  geom_col(aes(fill = Strand)) +
  ggrepel::geom_text_repel(aes(label= paste0(pct, "%")),
                   position=position_stack(vjust=0.5, reverse = T), 
                   size = 3) +
  # geom_bar(aes(fill = group), position="dodge", stat="identity") +
  labs(y = "# Predictions with G4PromFinder", x = "", caption = "2-kmer shuffling used") +
  coord_flip() + 
  facet_grid(~ group, scales = "free_y") +
  scale_fill_manual("", values = c("#3182bd", "#de2d26")) +
  theme_bw(base_family = "GillSans") -> p2

dir <- "~/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/G4PromFinder_outputs"


ggsave(p2, filename = "G4PromFinder_kmer2.png", path = dir, 
       width = 10, height = 12)


# by GC content

y <- read.csv("~/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/genomes/GC_content.csv") %>% as_tibble()

rm_genomes <- c("Tenacibaculum_dicentrarchi_gca_001483385", "Tenacibaculum_sp_dsm_106434_gca_003867015", "Tenacibaculum_sp_sz_18_gca_002813915", "Tenacibaculum_todarodis_gca_001889045")

rbind(df1_summ, df2_summ) %>%
  ungroup() %>%
  filter(!id %in% rm_genomes) %>%
  group_by(Strand, group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  mutate(sep = id) %>%
  separate(sep, c("Genus", "sp"), "_") %>%
  inner_join(y) %>% 
  filter(group %in% "Genome") %>%
  filter(Strand %in% "positivo") %>%
  mutate(d = n/Size * 1E6, 
         Size = Size * 1E6) %>%
  mutate(so_sort = as.integer(as.factor(sp))) %>%
  arrange(so_sort) %>%
  mutate(sp = forcats::fct_reorder(sp, so_sort)) %>%
  # mutate(sp = forcats::fct_reorder(sp, Size)) %>%
  ggplot() +
  geom_point(aes(y = d, 
                 size = n,
                 x = log(Size), color = GC)) +
  facet_wrap(~ sp) +
  # facet_grid( ~ Strand) +
  labs(y = "N/Mbp", x = "log(Mbp)") +
  theme_bw(base_size = 16, base_family = "GillSans") +
  ggsci::scale_color_gsea(name = "GC %", 
                             na.value = 'white',
                             limits = c(30,34),
                             labels = scales::percent_format(scale = 1)) -> p4

dir <- "~/Documents/GitHub/2nd-Workshop-in-Advanced-Bioinformatics/G4PromFinder_outputs"

ggsave(p4, filename = "S1_figure.png", path = dir, 
       width = 10, height = 10)

rbind(df1_summ, df2_summ) %>%
  group_by(Strand, group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  inner_join(y) %>% 
  filter(group %in% "Genome") %>%
  filter(Strand %in% "positivo") %>%
  separate(id, c("Genus", "sp"), "_") %>%
  mutate(d = n/Size * 1E6, 
         Size = Size * 1E6) %>%
  ggplot() +
  geom_point(aes(y = d, 
                 x = GC, color = sp), 
             size = 3)



