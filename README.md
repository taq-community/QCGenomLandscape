# Portrait Génomique du Québec

Cartographie de la disponibilité des données génomiques pour les espèces documentées au Québec. Le projet interroge les bases de données NCBI et BOLD Systems afin d'établir un portrait de la couverture génomique des ~24 800 espèces de la liste BDQC, avec un focus sur les marqueurs moléculaires (COI, cytb, gènes mitochondriaux) et les espèces à statut de conservation.

## Vue d'ensemble

Le pipeline:

1. Télécharge les séquences barcode disponibles sur **BOLD Systems** pour le Québec
2. Interroge **NCBI** pour les séquences nucléotidiques (avec et sans voucher) ainsi que les génomes complets
3. Filtre géographiquement les enregistrements avec coordonnées à l'intérieur du Québec
4. Intègre les statuts de conservation (COSEPAC / LEMV) et produit des figures de couverture

## Structure du projet

```
QCGenomLandscape/
├── data/
│   ├── bdqc_list_01122025.csv                    # Liste des espèces BDQC (~24 800 sp.)
│   ├── CA_especes_en_peril.csv                   # Espèces à statut COSEPAC (Canada)
│   ├── QC_especes_en_peril.csv                   # Espèces à statut LEMV (Québec)
│   ├── primers_keck_2022.csv                     # Amorces de référence (Keck 2022)
│   ├── primers_map_group_bdqc_list_01122025.csv  # Correspondance marqueurs × groupes taxonomiques
│   └── canvec_1M_CA_Admin.gdb/                  # Limites administratives (CanVec 1M)
├── scripts/
│   ├── get_bold_seqs.R              # Télécharge les séquences BOLD pour le Québec
│   ├── get_ncbi_seqs.R              # Requêtes NCBI nucléotide (séquences avec voucher)
│   ├── get_ncbi_seqs_non_voucher.R  # Requêtes NCBI nucléotide (sans filtre voucher) + filtrage QC
│   ├── get_ncbi_full_genome.R       # Requêtes NCBI genome/nucléotide pour génomes complets
│   ├── get_sra.R                    # Récupération des données SRA
│   ├── create_dataframe.R           # Intègre tous les résultats en un tableau récapitulatif
│   ├── taxon_representation.R       # Figures de prévalence génique et couverture espèces à risque
│   ├── map_bold.R                   # Cartographie des données BOLD
│   ├── map_ncbi.R                   # Cartographie des données NCBI
│   ├── map_combined.R               # Application Shiny: carte hexagonale BOLD + NCBI
│   ├── quality_seq_check.R          # Contrôle qualité des séquences
│   └── utils.R                      # Fonctions utilitaires (ex. parse_latlon)
├── results/
│   ├── bold_qc_data.tsv                          # Séquences BOLD brutes (Québec)
│   ├── ncbi_results.rds                          # Séquences NCBI avec voucher
│   ├── ncbi_non_voucher_results.rds              # Séquences NCBI sans voucher
│   ├── deficient_queries.rds                     # Requêtes en erreur (voucher)
│   ├── deficient_queries_non_voucher_results.rds # Requêtes en erreur (non-voucher)
│   ├── high_id_queries.rds                       # Requêtes avec > 500 résultats (voucher)
│   ├── high_id_queries_non_voucher_results.rds   # Requêtes avec > 500 résultats (non-voucher)
│   ├── genes_subsamp_50_df.rds                   # Gènes annotés (sous-échantillon 50 acc./sp.)
│   ├── genes_prevalence.png / .svg               # Figure: prévalence des gènes par groupe
│   └── risk_status_coverage.png / .svg           # Figure: couverture génomique × statut de risque
└── logs/                            # Journaux horodatés des requêtes API
```

## Prérequis

- R ≥ 4.1
- Clé API NCBI (gratuite sur https://www.ncbi.nlm.nih.gov/account/)

### Packages R

```r
install.packages(c(
  "tidyverse", "httr2", "rentrez", "logger", "glue",
  "xml2", "sf", "ggplot2", "purrr",
  "shiny", "bslib", "leaflet", "mapview"
))
```

### Clé API NCBI

Ajouter dans `.Renviron` à la racine du projet:

```
NCBI_API_KEY=votre_clé_ici
```

## Utilisation

Les scripts doivent être exécutés dans l'ordre suivant:

### 1. Télécharger les séquences BOLD

```bash
Rscript scripts/get_bold_seqs.R
```

Interroge le portail BOLD Systems pour toutes les séquences du Québec et enregistre le résultat dans `results/bold_qc_data.tsv`.

### 2. Récupérer les séquences NCBI (avec voucher)

```bash
Rscript scripts/get_ncbi_seqs.R
```

Pour chaque espèce × marqueur moléculaire, interroge la base `nucleotide` de NCBI en filtrant sur les séquences associées à un spécimen voucher. Produit `results/ncbi_results.rds`.

### 3. Récupérer les séquences NCBI (sans voucher)

```bash
Rscript scripts/get_ncbi_seqs_non_voucher.R
```

Même requêtes sans le filtre `voucher[Title]`. Filtre ensuite les enregistrements géoréférencés à l'intérieur du Québec à l'aide des limites CanVec. Produit `results/ncbi_non_voucher_results.rds`.

### 4. Récupérer les génomes complets

```bash
Rscript scripts/get_ncbi_full_genome.R
```

Interroge les bases `genome` (génomes nucléaires) et `nucleotide` (génomes mitochondriaux complets) pour chaque espèce de la liste BDQC. Produit `results/ncbi_genome_results.rds`.

### 5. Analyser et produire les figures

```bash
Rscript scripts/taxon_representation.R
```

Génère les figures de prévalence génique par groupe taxonomique et de couverture des espèces à statut de conservation. Produit `results/genes_prevalence.svg` et `results/risk_status_coverage.svg`.

### 6. Construire le tableau récapitulatif

```bash
Rscript scripts/create_dataframe.R
```

Intègre les résultats NCBI, les gènes annotés, la taxonomie BDQC et les statuts de conservation (COSEPAC/LEMV) en un seul tableau par espèce.

### 7. Visualiser sur la carte interactive

```r
shiny::runApp("scripts/map_combined.R")
```

Lance une application Shiny avec une carte hexagonale (10 km) des occurrences BOLD + NCBI géoréférencées au Québec, filtrable par famille, rang taxonomique et période de collecte.

## Bases de données

| Source | Contenu | Accès |
|--------|---------|-------|
| BOLD Systems | Barcodes ADN (COI principalement), enregistrements géoréférencés | API REST (`httr2`) |
| NCBI Nucleotide | Séquences nucléotidiques par espèce × marqueur | `rentrez` |
| NCBI Genome | Génomes nucléaires assemblés | `rentrez` |

## Outputs principaux

| Fichier | Description |
|---------|-------------|
| `results/bold_qc_data.tsv` | Données BOLD brutes pour le Québec |
| `results/ncbi_results.rds` | Séquences NCBI avec voucher (métadonnées complètes) |
| `results/ncbi_non_voucher_results.rds` | Séquences NCBI sans filtre voucher |
| `results/genes_subsamp_50_df.rds` | Gènes annotés (sous-échantillon 50 acc./sp.) |
| `results/genes_prevalence.svg` | Prévalence des marqueurs géniques par groupe taxonomique |
| `results/risk_status_coverage.svg` | Couverture génomique par statut de risque (QC) |

## Limitations

- Les requêtes NCBI sont soumises au rate limiting de l'API (une clé API augmente le quota)
- Les séquences NCBI ne sont pas spécifiques au Québec; le filtrage géographique ne s'applique qu'aux enregistrements avec coordonnées
- Biais taxonomiques importants: vertébrés et arthropodes mieux représentés que les champignons, invertébrés, et microorganismes
- La nomenclature taxonomique entre BDQC, NCBI et BOLD n'est pas toujours concordante

Voir [APPROCHE.md](APPROCHE.md) pour la justification méthodologique détaillée.

## Auteurs

- Steve Vissault
- Marie Pier Brochu
- Valérie Langlois

## Citation

```
Vissault, S., Brochu, M.P., & Langlois, V. (2025).
Portrait génomique du Québec: Cartographie de la disponibilité
des données génomiques pour les espèces observées au Québec.
```

## Références

- NCBI: https://www.ncbi.nlm.nih.gov/
- BOLD Systems: https://boldsystems.org/
- Biodiversité Québec: https://biodiversite-quebec.ca/
- CanVec: https://open.canada.ca/data/en/dataset/306e5004-534b-4110-9feb-58e3a5c3fd97
