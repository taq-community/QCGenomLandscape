# Exploratory search of NCBI SRA for Quebec eDNA/metabarcoding datasets.
# Not part of the reproducible _targets.R pipeline -- nothing downstream
# consumes this yet.
library(rentrez)

# Recherche avec variantes géographiques
search_terms <- c(
    "(Quebec OR Québec OR QC) AND (eDNA OR environmental DNA OR metabarcoding)",
    "Canada[Geography] AND (Quebec OR Québec) AND amplicon"
)

sra_results <- rentrez::entrez_search(
    db = "sra",
    term = search_terms[2],
    retmax = 99000
)
