# QCGenomLandscape

Cartographie de la disponibilité des données génomiques pour les espèces
documentées au Québec. Le package interroge les bases de données NCBI et
BOLD Systems afin d’établir un portrait de la couverture génomique des
espèces de la liste BDQC, avec un focus sur les marqueurs moléculaires
(COI, cytb, gènes mitochondriaux) et les espèces à statut de
conservation.

## Vue d’ensemble

Le pipeline:

1.  Télécharge les séquences barcode disponibles sur **BOLD Systems**
    pour le Québec
2.  Interroge **NCBI** pour les séquences nucléotidiques (avec et sans
    voucher) ainsi que les génomes complets
3.  Filtre géographiquement les enregistrements avec coordonnées à
    l’intérieur du Québec
4.  Intègre les statuts de conservation (COSEPAC / LEMV), contrôle la
    qualité des séquences (détection de codons stop) et produit des
    figures de couverture
5.  Exporte optionnellement les livrables vers un conteneur Arbutus
    (Swift)

## Structure du projet

    QCGenomLandscape/
    ├── R/                              # Fonctions exportées du package
    │   ├── bold.R                      # fetch_bold_sequences()
    │   ├── ncbi_sequences.R            # build_ncbi_queries(), fetch_ncbi_sequences()
    │   ├── ncbi_genome.R               # query_full_genome(), fetch_ncbi_genomes()
    │   ├── genbank.R                   # fetch_gene_annotations(), fetch_gb_records(), parse_gb_records()
    │   ├── gene_classification.R       # assign_gene_group()
    │   ├── taxonomy.R                  # classify_taxon_group()
    │   ├── risk_status.R               # load_risk_status()
    │   ├── summary_table.R             # build_summary_dataframe()
    │   ├── plots.R                     # plot_gene_prevalence(), plot_risk_status_coverage()
    │   ├── sequence_qc.R               # has_stop_codon_coi(), score_sequence_quality(), align_sequences()
    │   ├── spatial.R                   # load_canvec_boundary(), flag_within_boundary(), make_hex_grid()
    │   ├── swift.R                     # swift_auth(), upload_to_swift()
    │   └── utils.R                     # parse_latlon()
    ├── tests/testthat/                 # Tests unitaires (aucun n'appelle une vraie API)
    ├── vignettes/
    │   └── approche-methodologique.Rmd # Justification méthodologique détaillée
    ├── _targets.R                      # Définition du pipeline (package `targets`)
    ├── exploration/                    # Scripts exploratoires, hors pipeline reproductible
    │   ├── get_sra.R
    │   ├── map_bold.R
    │   └── map_ncbi.R
    ├── internal/                       # Scripts liés à une feuille Google Sheets privée, hors package
    │   ├── fill_arisque_sheet.R
    │   └── get_ncbi_full_genome_arisque_sheet.R
    ├── data/                            # Données sources (gitignored, voir ci-dessous)
    ├── results/                         # Sorties du pipeline (gitignored)
    └── logs/                            # Journaux horodatés des requêtes API (gitignored)

`data/` et `results/` ne sont pas versionnés (trop volumineux —
`results/` seul dépasse 200 Mo) ; les fonctions du package prennent les
chemins en argument plutôt que de les supposer fixes.

## Installation

### Prérequis

- R ≥ 4.1
- Clé API NCBI (gratuite sur <https://www.ncbi.nlm.nih.gov/account/>)

### Packages Bioconductor

`Biostrings` et `DECIPHER` viennent de Bioconductor, pas de CRAN — à
installer avant le reste :

``` r

install.packages("BiocManager")
BiocManager::install(c("Biostrings", "DECIPHER"))
```

### Le package

``` r

# install.packages("remotes")
remotes::install_local(".")

# Pour lancer le pipeline (targets) et régénérer les figures/vignettes
remotes::install_local(".", dependencies = TRUE)
```

### Clé API NCBI

Ajouter dans `.Renviron` à la racine du projet:

    NCBI_API_KEY=votre_clé_ici

### Export vers Arbutus (optionnel)

L’authentification se fait par **identifiant applicatif** (application
credential), pas par mot de passe — un compte protégé par MFA (comme la
plupart des comptes Alliance Canada) ne peut pas s’authentifier par mot
de passe simple via l’API.

Pour en générer un : tableau de bord Horizon d’Arbutus → *Identity* →
*Application Credentials* → *Create Application Credential* (ou
`openstack application credential create <nom>` en ligne de commande,
après une connexion interactive qui satisfait le MFA). L’ID et le secret
ne sont affichés qu’une seule fois à la création.

Ajouter ensuite dans `.Renviron` :

    SWIFT_APP_CRED_ID=votre_id
    SWIFT_APP_CRED_SECRET=votre_secret

L’URL d’authentification a déjà la valeur par défaut du projet
QCGenomicLandscape sur Arbutus — voir
[`?swift_auth`](https://taq-community.github.io/QCGenomLandscape/reference/swift_auth.md)
pour la surcharger (`SWIFT_AUTH_URL`) si vous pointez vers un autre
projet OpenStack. Sans `SWIFT_APP_CRED_ID`/`SWIFT_APP_CRED_SECRET`, le
pipeline s’exécute normalement — les étapes d’export sont simplement
ignorées (avec un avertissement).

## Utilisation

Le pipeline est orchestré avec
[`targets`](https://books.ropensci.org/targets/) :

``` r

# Aperçu du graphe de dépendances, sans rien exécuter
targets::tar_visnetwork()

# Exécute (ou reprend) le pipeline au complet
targets::tar_make()

# Lire un résultat intermédiaire ou final
targets::tar_read(summary_table)
```

Voir
[`vignette("approche-methodologique")`](https://taq-community.github.io/QCGenomLandscape/vignettes/approche-methodologique.Rmd)
pour le détail de chaque étape et la justification méthodologique.

Les fonctions du package peuvent aussi être utilisées directement, hors
du pipeline :

``` r

library(QCGenomLandscape)

bold_data <- fetch_bold_sequences()

qc_species <- read.csv("data/bdqc_list_01122025.csv") |>
  dplyr::filter(rank == "species")
query_primers <- read.csv2("data/primers_map_group_bdqc_list_01122025.csv")

queries <- build_ncbi_queries(qc_species, query_primers, voucher = TRUE)
ncbi <- fetch_ncbi_sequences(queries)
```

## Bases de données

| Source | Contenu | Accès |
|----|----|----|
| BOLD Systems | Barcodes ADN (COI principalement), enregistrements géoréférencés | API REST (`httr2`) |
| NCBI Nucleotide | Séquences nucléotidiques par espèce × marqueur | `rentrez` |
| NCBI Genome | Génomes nucléaires assemblés | `rentrez` |

## Outputs principaux

| Fichier | Description |
|----|----|
| `results/bold_qc_data.tsv` | Données BOLD brutes pour le Québec |
| `results/ncbi_results.rds` | Séquences NCBI avec voucher (métadonnées complètes) |
| `results/ncbi_non_voucher_results.rds` | Séquences NCBI sans filtre voucher |
| `results/ncbi_genome_results.rds` | Génomes complets (nucléaires/mitochondriaux) par espèce |
| `results/genes_subsamp_50_df.rds` | Gènes annotés (sous-échantillon) |
| `results/sequence_qc.rds` | Métriques de qualité de séquence + détection de codons stop |
| `results/genes_prevalence.svg` | Prévalence des marqueurs géniques par groupe taxonomique |
| `results/risk_status_coverage.svg` | Couverture génomique par statut de risque (QC) |

## Limitations

- Les requêtes NCBI sont soumises au rate limiting de l’API (une clé API
  augmente le quota)
- Les séquences NCBI ne sont pas spécifiques au Québec; le filtrage
  géographique ne s’applique qu’aux enregistrements avec coordonnées
- Biais taxonomiques importants: vertébrés et arthropodes mieux
  représentés que les champignons, invertébrés, et microorganismes
- La nomenclature taxonomique entre BDQC, NCBI et BOLD n’est pas
  toujours concordante

Voir la vignette `approche-methodologique` pour la justification
méthodologique détaillée.

## Auteurs

- Steve Vissault
- Marie Pier Brochu
- Valérie Langlois

## Citation

    Vissault, S., Brochu, M.P., & Langlois, V. (2025).
    Portrait génomique du Québec: Cartographie de la disponibilité
    des données génomiques pour les espèces observées au Québec.

## Références

- NCBI: <https://www.ncbi.nlm.nih.gov/>
- BOLD Systems: <https://boldsystems.org/>
- Biodiversité Québec: <https://biodiversite-quebec.ca/>
- CanVec:
  <https://open.canada.ca/data/en/dataset/306e5004-534b-4110-9feb-58e3a5c3fd97>
- targets: <https://books.ropensci.org/targets/>
