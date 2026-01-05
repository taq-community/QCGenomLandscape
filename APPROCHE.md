# Approche méthodologique - Portrait génomique du Québec

## Objectif

Établir un portrait de la disponibilité des données génomiques pour les ~24,800 espèces observées au Québec, en interrogeant les principales bases de données publiques.

## Sources de données

### 1. Liste des espèces
- **Source**: BDQC (Biodiversité du Québec)
- **Fichier**: `data/bdqc_list_01122025.csv`
- **Contenu**: 24,815 observations d'espèces avec informations taxonomiques complètes
- **Groupes couverts**:
  - Arthropodes (insectes, araignées, crustacés)
  - Plantes (angiospermes, conifères, bryophytes)
  - Fungi (mycètes)
  - Oiseaux
  - Autres vertébrés et invertébrés
  - Bactéries

### 2. Bases de données génomiques interrogées

#### NCBI (National Center for Biotechnology Information)
- **API**: rentrez (package R)
- **Bases interrogées**:
  - **nucleotide**: Toutes les séquences d'ADN/ARN
  - **genome**: Génomes complets assemblés
- **Type de données**:
  - Gènes individuels
  - Régions génomiques
  - Génomes mitochondriaux
  - Génomes nucléaires complets
- **Avantages**: La plus grande base de données publique, tous types d'organismes

#### BOLD Systems (Barcode of Life Data System)
- **API**: HTTP REST API
- **Focus**: Barcoding ADN (principalement COI pour animaux)
- **Type de données**:
  - Séquences COI standardisées
  - Métadonnées d'échantillonnage
  - Photos de spécimens
- **Avantages**: Données standardisées, particulièrement riche pour arthropodes et vertébrés

## Workflow d'analyse

### Phase 1: Interrogation des APIs

```
Pour chaque espèce:
├── Requête NCBI nucleotide
│   ├── Comptage des séquences
│   └── Note: tous types confondus
├── Requête NCBI genome
│   └── Vérification génome complet
└── Requête BOLD
    ├── Comptage des barcodes
    └── Note: principalement COI
```

**Paramètres techniques**:
- Délai entre requêtes: 0.5 secondes (respect rate limiting)
- Sauvegarde intermédiaire: tous les 50 espèces
- Temps estimé pour 24,800 espèces: 6-8 heures

### Phase 2: Traitement des données

**Nettoyage**:
- Utilisation des noms scientifiques valides (`valid_scientific_name`)
- Gestion des synonymes taxonomiques
- Filtrage des doublons

**Enrichissement**:
- Jointure avec données taxonomiques complètes
- Classification en catégories de couverture
- Calcul de métriques par groupe taxonomique

### Phase 3: Analyse et visualisation

**Métriques calculées**:
1. Nombre d'espèces avec séquences NCBI
2. Nombre d'espèces avec barcodes BOLD
3. Nombre d'espèces avec génomes complets
4. Pourcentage de couverture global et par groupe
5. Distribution du nombre de séquences

**Visualisations générées**:
1. Vue d'ensemble de la couverture
2. Couverture par royaume taxonomique
3. Top 20 des classes les mieux représentées
4. Carte de chaleur phylum × classe
5. Distribution du nombre de séquences
6. Espèces avec génomes complets

## Structure des outputs

### Fichiers de données

```
results/
├── genomic_data_test.csv              # Résultats mode test (100 espèces)
├── genomic_data_complete.csv          # Résultats complets
├── genomic_data_with_taxonomy.csv     # Résultats + taxonomie
├── summary_statistics.csv             # Statistiques générales
└── statistics_by_taxonomy.csv         # Stats par groupe
```

### Fichiers de visualisation

```
figures/
├── 00_composite.png                   # Figure composite
├── 01_overview.png                    # Vue d'ensemble
├── 02_coverage_by_kingdom.png         # Par royaume
├── 03_top_classes.png                 # Top 20 classes
├── 04_heatmap_coverage.png            # Carte de chaleur
├── 05_sequence_distribution.png       # Distribution
└── 06_genomes_by_class.png           # Génomes complets
```

## Limitations et considérations

### Limitations techniques

1. **Rate limiting**:
   - NCBI: ~3 requêtes/seconde (avec clé API: 10/sec)
   - BOLD: Pas de limite stricte, mais délai recommandé

2. **Exactitude taxonomique**:
   - Dépend de la nomenclature utilisée dans les bases
   - Certaines espèces peuvent avoir plusieurs noms
   - Mises à jour taxonomiques non synchronisées

3. **Qualité des données**:
   - NCBI: qualité variable (soumissions publiques)
   - BOLD: données plus standardisées mais moins diversifiées

### Limitations biologiques

1. **Biais taxonomiques**:
   - Vertébrés et plantes sur-représentés
   - Arthropodes: surtout insectes et araignées
   - Fungi et bactéries sous-représentés

2. **Biais géographiques**:
   - Données mondiales, pas spécifiques au Québec
   - Certaines espèces jamais séquencées
   - Priorités de séquençage (modèles, agriculture, santé)

3. **Types de données**:
   - BOLD: surtout COI (un seul gène)
   - NCBI: très variable selon les groupes
   - Génomes complets: minorité d'espèces

## Interprétation des résultats

### Catégories de couverture

- **Excellente** (>80%): Génome + multiples marqueurs
- **Bonne** (50-80%): Plusieurs séquences disponibles
- **Partielle** (10-50%): Quelques séquences, souvent un seul marqueur
- **Faible** (<10%): Très peu de données
- **Absente** (0%): Aucune donnée disponible

### Priorités de séquençage

Les espèces sans données génomiques peuvent être priorisées selon:
1. Importance écologique
2. Statut de conservation
3. Endémisme
4. Facilité d'échantillonnage
5. Questions de recherche spécifiques

## Prochaines étapes

1. **Analyse approfondie**:
   - Identifier les lacunes critiques
   - Comparer avec autres régions
   - Analyser les tendances temporelles

2. **Valorisation**:
   - Site web interactif
   - Rapports par groupe taxonomique
   - Recommandations pour nouveaux séquençages

3. **Mise à jour**:
   - Réexécution périodique (ex: annuelle)
   - Intégration d'autres bases (ENA, JGI, etc.)
   - Ajout de nouvelles espèces observées

## Scripts disponibles

1. **query_genomic_data.R**: Interrogation des APIs
2. **create_visualizations.R**: Génération des figures
3. **README.md**: Guide d'utilisation détaillé

## Références

- NCBI: https://www.ncbi.nlm.nih.gov/
- BOLD Systems: http://www.boldsystems.org/
- rentrez package: https://cran.r-project.org/package=rentrez
- BDQC: https://biodiversite-quebec.ca/
