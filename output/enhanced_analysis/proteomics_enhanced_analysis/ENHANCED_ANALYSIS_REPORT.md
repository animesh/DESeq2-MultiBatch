# Enhanced Proteomics Analysis Report

## Data Summary
- Total proteins analyzed: 3337
- Total samples: 92
- Batch 1 samples: 46
- Batch 2 samples: 46

## Batch Correction Results
- Mean absolute batch effect: 2.251 log2 units
- Correction method: DESeq2-inspired symmetric median centering
- Statistical analysis: Full DESeq2 pipeline with negative binomial modeling

## Comprehensive Annotation Analysis

### Categorical Variables:
- ** Batch **: 0 → 0 significant proteins
- ** Sex **: 240 → 179 significant proteins
- ** Treatment **: 65 → 36 significant proteins
- ** Location **: 0 → 0 significant proteins
- ** Relapse_Risk **: 225 → 133 significant proteins
- ** Subject **: 148 → 150 significant proteins

### Numerical Variables:
- ** Age **: 0 → 0 significant proteins
- ** Days_Until_Relapse **: 0 → 0 significant proteins
