# ==============
## Checking and Load packages ----
# ==============

.bioc_packages <- c("data.table", "tidyr", "igraph", "dplyr", "magrittr", "ggplot2", "NetSci")

.inst <- .bioc_packages %in% installed.packages()

if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(.bioc_packages, require, character.only = TRUE)


url <- "https://raw.githubusercontent.com/deisygysi/NetMed_Workshop/master/data/PPI_Symbol_Entrez.csv"

PPI <- fread(url)

PPI %>%
  select(starts_with("Symbol")) %>%
  # drop_na() %>%
  filter_all(any_vars((. != ""))) %>%
  graph_from_data_frame(., directed = FALSE, v = NULL) %>%
  simplify() -> gPPI

dd = degree(gPPI) %>% table() %>% as.data.frame()
names(dd) = c('Degree', "Nodes")
dd$Degree %<>% as.character %>% as.numeric()
dd$Nodes  %<>% as.character %>% as.numeric()

# esquisse::esquisser()

ggplot(dd) +
  aes(x = Degree, y = Nodes) +
  geom_point(colour = "#1d3557") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  theme_minimal()
  
degree(gPPI) %>% as.data.frame() %>% View()

url <- "https://raw.githubusercontent.com/deisygysi/NetMed_Workshop/master/data/curated_gene_disease_associations.tsv"

GDA <- fread(url)

Cleaned_GDA = GDA %>% filter(diseaseType == 'disease') %>%
  mutate(diseaseName = tolower(diseaseName)) %>%
  select(geneSymbol, diseaseName, diseaseSemanticType) %>%
  unique() 

numGenes = Cleaned_GDA %>% 
  group_by(diseaseName) %>%
  summarise(numGenes = n()) %>%
  ungroup() %>%
  group_by(numGenes) %>%
  summarise(numDiseases = n())

ggplot(numGenes) +
  aes(x = numGenes, y = numDiseases, size = numGenes) +
  geom_point(colour = "#1d3557") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(x = "Genes", y = "Diseases")+
  theme_minimal()

Cleaned_GDA %<>% 
  group_by(diseaseName) %>%
  mutate(numGenes = n()) %>%
  filter(numGenes > 10)

Cleaned_GDA$diseaseName %>%
  unique() %>%  
  length()

url <- "https://raw.githubusercontent.com/deisygysi/NetMed_Workshop/master/data/DB_DrugTargets_1201.csv"

DT = fread(url)

Cleaned_DT = DT %>% 
  filter(organism == 'Humans') %>%
  select(Gene_Target, Name,ID, Type, known_action) %>%
  unique() 

TargetDist = Cleaned_DT %>% 
  group_by(Gene_Target) %>%
  summarise(numDrugs = n()) 

DrugDist = Cleaned_DT %>% 
  group_by(ID) %>%
  summarise(numTargets = n()) 

ggplot(TargetDist) +
  aes(x = numDrugs) +
  geom_histogram(colour = "#1d3557", fill = "#a8dadc" ) +
  labs(x = "Targets", y = "Drugs")+
  theme_minimal()

TargetDist %>%
  arrange(desc(numDrugs)) %>%
  filter(numDrugs > 400)

ggplot(DrugDist) +
  aes(x = numTargets) +
  geom_histogram(colour = "#1d3557", fill = "#a8dadc" ) +
  labs(y = "Targets", x = "Drugs")+
  theme_minimal()

# Chapter 3 Methods for Disease Module Identification and Disease Similarity

#First, let's select genes that are associated with Schizophrenia.

SCZ_Genes = 
  Cleaned_GDA %>% 
  filter(diseaseName %in% 'schizophrenia') %>%
  pull(geneSymbol) %>% 
  unique()

# Next, let's see how they are localized in the PPI.
# Fist, we have to make sure all genes are in the PPI.
# Later, we calculate the LCC.
# And lastly, let's visualize it.

SCZ_PPI = SCZ_Genes[SCZ_Genes %in% V(gPPI)$name]

gScz = gPPI %>%
  induced.subgraph(., SCZ_PPI)

components(gScz)

# Calculates the Largest Connected Component (LCC)
# it takes a time!!
LCC_scz = LCC_Significance(N = 1000, Targets = SCZ_PPI,
                           G = gPPI)
Histogram_LCC(LCC_scz)


V(gScz)$size = degree(gScz) %>% 
  CoDiNA::normalize()
V(gScz)$size = (V(gScz)$size + 0.1)*5
V(gScz)$color = '#83c5be'
V(gScz)$frame.color = '#006d77'
V(gScz)$label = ifelse(V(gScz)$size  > 4, V(gScz)$name, NA )
V(gScz)$label.color = "black" # '#e29578'

E(gScz)$width = edge.betweenness(gScz, directed = F) %>% CoDiNA::normalize()
E(gScz)$width = E(gScz)$width + 0.01
E(gScz)$weight = E(gScz)$width
par(mar = c(0,0,0,0))
plot(gScz)

# calculate edge betweness

gScz %<>% delete.vertices(., degree(.) == 0)

V(gScz)$size = degree(gScz) %>% 
  CoDiNA::normalize()
V(gScz)$size = (V(gScz)$size + 0.1)*5
V(gScz)$color = '#83c5be'
V(gScz)$frame.color = '#006d77'
V(gScz)$label = ifelse(V(gScz)$size  > 4, V(gScz)$name, NA )
V(gScz)$label.color = 'black'

E(gScz)$width = edge.betweenness(gScz, directed = F) %>% CoDiNA::normalize()
E(gScz)$width = E(gScz)$width + 0.01
E(gScz)$weight = E(gScz)$width
par(mar = c(0,0,0,0))
plot(gScz)


Dis_Ex1 = c('schizophrenia',
            "autistic disorder", 
            'obesity',
            'hyperlipidemia',
            'rheumatoid arthritis')
GDA_Interest = Cleaned_GDA %>% 
  filter(diseaseName %in% Dis_Ex1) %>%
  select(diseaseName, geneSymbol) %>%
  unique()

# calculate distance 
# Note the diameter here is the average distance instead of the maximal distance

Jaccard_Ex2 = Jaccard(GDA_Interest)


# Calculates the separation of two set of targets on a network.
# the more negative the value are, it is better
sab = separation(gPPI, GDA_Interest)

Sep_ex2 = sab$Sab %>% as.matrix()

Sep_ex2[lower.tri(Sep_ex2)] = t(Sep_ex2)[lower.tri(Sep_ex2)]

Sep_ex2 %>% heatmap(., symm = T)

require(eulerr)
# install.packages("eulerr")

Euler_List = list (
  SCZ = GDA_Interest$geneSymbol[GDA_Interest$diseaseName == 'schizophrenia'],
  
  ASD = GDA_Interest$geneSymbol[GDA_Interest$diseaseName == 'autistic disorder'],
  
  OB = GDA_Interest$geneSymbol[GDA_Interest$diseaseName == 'obesity'],
  
  HD = GDA_Interest$geneSymbol[GDA_Interest$diseaseName == 'hyperlipidemia'],
  
  RA = GDA_Interest$geneSymbol[GDA_Interest$diseaseName == 'rheumatoid arthritis'])

EULER = euler(Euler_List)
plot(EULER, quantities = TRUE)
