# Protein Overlap Analysis Report

## Executive Summary

- **Total unique proteins**:  1047
- **Factors analyzed**:  6
- **Pairwise comparisons**:  15

## Factor Summary

- ** Relapse_Risk **:  475  proteins ( 45.37 % of total unique proteins)
- ** Sex **:  436  proteins ( 41.64 % of total unique proteins)
- ** Subject **:  303  proteins ( 28.94 % of total unique proteins)
- ** Days_Until_Relapse **:  189  proteins ( 18.05 % of total unique proteins)
- ** Age **:  81  proteins ( 7.74 % of total unique proteins)
- ** Treatment **:  11  proteins ( 1.05 % of total unique proteins)

## Top Overlaps (by count)

- ** Relapse_Risk  ∩  Sex **:  162  proteins ( 34.1 % of  Relapse_Risk ,  37.2 % of  Sex )
- ** Sex  ∩  Subject **:  79  proteins ( 18.1 % of  Sex ,  26.1 % of  Subject )
- ** Relapse_Risk  ∩  Subject **:  69  proteins ( 14.5 % of  Relapse_Risk ,  22.8 % of  Subject )
- ** Days_Until_Relapse  ∩  Relapse_Risk **:  58  proteins ( 30.7 % of  Days_Until_Relapse ,  12.2 % of  Relapse_Risk )
- ** Days_Until_Relapse  ∩  Sex **:  48  proteins ( 25.4 % of  Days_Until_Relapse ,  11 % of  Sex )
- ** Days_Until_Relapse  ∩  Subject **:  29  proteins ( 15.3 % of  Days_Until_Relapse ,  9.6 % of  Subject )
- ** Age  ∩  Relapse_Risk **:  24  proteins ( 29.6 % of  Age ,  5.1 % of  Relapse_Risk )
- ** Age  ∩  Sex **:  24  proteins ( 29.6 % of  Age ,  5.5 % of  Sex )
- ** Age  ∩  Subject **:  19  proteins ( 23.5 % of  Age ,  6.3 % of  Subject )
- ** Age  ∩  Days_Until_Relapse **:  6  proteins ( 7.4 % of  Age ,  3.2 % of  Days_Until_Relapse )

## Highest Similarity (by Jaccard Index)

- ** Relapse_Risk  ∩  Sex **: Jaccard =  0.216  ( 162  proteins)
- ** Sex  ∩  Subject **: Jaccard =  0.12  ( 79  proteins)
- ** Relapse_Risk  ∩  Subject **: Jaccard =  0.097  ( 69  proteins)
- ** Days_Until_Relapse  ∩  Relapse_Risk **: Jaccard =  0.096  ( 58  proteins)
- ** Days_Until_Relapse  ∩  Sex **: Jaccard =  0.083  ( 48  proteins)

## Specific Overlaps of Interest

- ** Sex and Age **:  24  proteins
- ** Sex and Days Until Relapse **:  48  proteins
- ** Sex and Relapse Risk **:  162  proteins
- ** Sex and Subject **:  79  proteins
- ** Sex and Treatment **:  3  proteins
- ** All Clinical **:  2  proteins
- ** Age and Relapse Risk **:  24  proteins
- ** Age and Days Until Relapse **:  6  proteins
- ** Relapse Risk and Days Until Relapse **:  58  proteins
- ** Core 3plus factors **:  67  proteins
- ** Core 4plus factors **:  6  proteins

## Core Proteins (Multi-Factor)

- ** 5 + factors**:  1  proteins
- ** 4 + factors**:  6  proteins
- ** 3 + factors**:  67  proteins

## Files Generated

### Visualizations:
- `factor_sizes.png`: Bar chart of proteins per factor
- `overlap_heatmap.png`: Heatmap of overlap counts
- `jaccard_heatmap.png`: Heatmap of Jaccard similarity indices
- `overlap_summary.png`: Scatter plot of overlap vs similarity

### Data Tables:
- `comprehensive_protein_table.csv`: All proteins with factor associations
- `pairwise_overlaps.csv`: Detailed pairwise overlap statistics
- `overlap_matrix.csv`: Overlap count matrix
- `factor_summary_stats.csv`: Summary statistics per factor
- `specific_overlaps_summary.csv`: Summary of specific overlap analyses
- `proteins_*.csv`: Individual protein lists for specific overlaps

## Key Biological Insights

1. **Sex effects** represent the largest single factor affecting protein expression
2. **Individual differences** (Subject) show substantial protein variation, indicating personalized signatures
3. **Clinical factors** (relapse risk, age) have meaningful overlaps suggesting shared pathways
4. **Core proteins** affected by multiple factors may represent central regulatory nodes
5. **Factor-specific proteins** indicate unique biological mechanisms

## Statistical Interpretation

- **Overlap counts** show absolute number of shared proteins
- **Jaccard indices** quantify proportional similarity (0-1 scale)
- **High Jaccard + high overlap** = strong biological relationship
- **High overlap + low Jaccard** = large factors with some shared biology
- **Low overlap** = factor-specific mechanisms

## Recommendations for Further Analysis

1. **Pathway analysis** on core proteins (3+ factors)
2. **Functional enrichment** of factor-specific proteins
3. **Network analysis** of highly overlapping factors
4. **Clinical validation** of sex-specific and relapse-associated proteins
5. **Mechanistic studies** of proteins unique to single factors
