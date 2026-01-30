# Avec rentrez
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

