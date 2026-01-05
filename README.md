# Portrait Génomique du Québec

Un projet pour cartographier la disponibilité des données génomiques pour toutes les espèces observées au Québec.

## Vue d'ensemble

Ce projet interroge les principales bases de données génomiques publiques (NCBI et BOLD Systems) pour établir un portrait de la couverture génomique des ~24,800 espèces documentées au Québec.

## Structure du projet

```
QCGenomLandscape/
├── data/
│   └── bdqc_list_01122025.csv       # Liste des espèces du Québec
├── scripts/
│   ├── query_genomic_data.R         # Script principal d'interrogation
│   ├── create_visualizations.R      # Génération des figures
│   ├── test_api.R                   # Test rapide des APIs
│   ├── utils.R                      # Fonctions utilitaires
│   └── README.md                    # Guide d'utilisation
├── results/                         # Résultats des analyses (généré)
├── figures/                         # Visualisations (généré)
├── _quarto.yml                      # Configuration du site web
├── index.qmd                        # Page d'accueil
├── intro.qmd                        # Introduction
├── summary.qmd                      # Résumé
├── references.qmd                   # Références
└── APPROCHE.md                      # Documentation méthodologique
```

## Installation

### Prérequis

- R (≥ 4.0.0)
- Quarto (pour générer le site web)

### Packages R nécessaires

```r
# Installer tous les packages nécessaires
install.packages(c(
  "tidyverse",
  "httr",
  "jsonlite",
  "rentrez",
  "bold",
  "ggplot2",
  "viridis",
  "patchwork",
  "scales"
))
```

## Utilisation

### 1. Tester les APIs

Avant de lancer l'analyse complète, testez que tout fonctionne:

```bash
Rscript scripts/test_api.R
```

### 2. Lancer l'analyse (mode test)

Par défaut, le script traite 100 espèces en mode test:

```bash
Rscript scripts/query_genomic_data.R
```

### 3. Lancer l'analyse complète

Pour traiter toutes les espèces (~6-8 heures):

```r
# Dans query_genomic_data.R, modifier la ligne 270:
TEST_MODE <- FALSE

# Puis lancer:
Rscript scripts/query_genomic_data.R
```

### 4. Créer les visualisations

Après avoir obtenu les résultats:

```bash
Rscript scripts/create_visualizations.R
```

### 5. Générer le site web

```bash
quarto render
quarto preview
```

## Bases de données interrogées

### NCBI (National Center for Biotechnology Information)
- **nucleotide**: Toutes les séquences d'ADN/ARN
- **genome**: Génomes complets assemblés
- **API**: rentrez (package R)

### BOLD Systems (Barcode of Life Data System)
- **Focus**: Barcoding ADN (COI)
- **API**: HTTP REST
- **Particulièrement riche pour**: Arthropodes et vertébrés

## Outputs

### Fichiers de données

Le dossier `results/` contient:

- `genomic_data_test.csv` - Résultats du test (100 espèces)
- `genomic_data_complete.csv` - Résultats complets
- `genomic_data_with_taxonomy.csv` - Résultats avec taxonomie complète
- `summary_statistics.csv` - Statistiques générales
- `statistics_by_taxonomy.csv` - Statistiques par groupe taxonomique
- `species_without_genomic_data.csv` - Espèces sans données
- `quick_report.md` - Rapport rapide

### Visualisations

Le dossier `figures/` contient:

- `00_composite.png` - Figure composite
- `01_overview.png` - Vue d'ensemble de la couverture
- `02_coverage_by_kingdom.png` - Couverture par royaume
- `03_top_classes.png` - Top 20 des classes
- `04_heatmap_coverage.png` - Carte de chaleur phylum × classe
- `05_sequence_distribution.png` - Distribution des séquences
- `06_genomes_by_class.png` - Génomes complets par classe

## Fonctions utilitaires

Le fichier `scripts/utils.R` contient des fonctions pratiques:

```r
source("scripts/utils.R")

# Charger et afficher un résumé
results <- load_results("results/genomic_data_complete.csv")
print_summary(results)

# Identifier les erreurs
errors <- find_errors(results)

# Exporter les espèces sans données
no_data <- export_species_without_data(results)

# Créer un rapport rapide
create_quick_report(results)
```

## Exemples d'analyse

### Charger et explorer les résultats

```r
library(tidyverse)

# Charger les résultats
results <- read_csv("results/genomic_data_with_taxonomy.csv")

# Espèces les mieux couvertes
top_species <- results %>%
  filter(!is.na(ncbi_sequences)) %>%
  arrange(desc(ncbi_sequences)) %>%
  select(species, kingdom, ncbi_sequences, bold_records, genome_available) %>%
  head(20)

# Groupes taxonomiques les mieux couverts
coverage_by_class <- results %>%
  group_by(class) %>%
  summarise(
    n_species = n(),
    n_with_data = sum(ncbi_sequences > 0 | bold_records > 0, na.rm = TRUE),
    pct_coverage = round(n_with_data / n_species * 100, 2)
  ) %>%
  filter(n_species >= 10) %>%
  arrange(desc(pct_coverage))
```

### Identifier les lacunes

```r
# Espèces communes sans données génomiques
common_no_data <- results %>%
  filter(
    obs_count > 100,
    (is.na(ncbi_sequences) | ncbi_sequences == 0),
    (is.na(bold_records) | bold_records == 0)
  ) %>%
  arrange(desc(obs_count)) %>%
  select(species, kingdom, class, obs_count, vernacular_fr)
```

## Limitations

### Techniques
- Dépendance aux APIs publiques (rate limiting)
- Nomenclature taxonomique non standardisée
- Qualité variable des données NCBI

### Biologiques
- Biais taxonomiques (vertébrés sur-représentés)
- Données mondiales (pas spécifiques au Québec)
- Certains groupes peu séquencés (fungi, invertébrés)

Voir [APPROCHE.md](APPROCHE.md) pour plus de détails.

## Contribution

### Auteurs
- Steve Vissault
- Marie Pier Brochu
- Valérie Langlois

### Citation

```
Vissault, S., Brochu, M.P., & Langlois, V. (2025).
Portrait génomique du Québec: Cartographie de la disponibilité
des données génomiques pour les espèces observées au Québec.
```

## Licence

À définir

## Contact

Pour questions ou suggestions, créer une issue sur le dépôt GitHub.

## Références

- NCBI: https://www.ncbi.nlm.nih.gov/
- BOLD Systems: http://www.boldsystems.org/
- Biodiversité Québec: https://biodiversite-quebec.ca/
- Quarto: https://quarto.org/
