## InferredGeneRegArchaicHominins

This code:

1. trains gene expression imputation models for non-GWARRs (GWARRs = genes without archaic regulatory regions) using only the individuals with no archaic hominin ancestry in the gene’s regulatory region 
2. applies these models, which are trained only on human-ancestry sequences, to individuals with archaic hominin ancestry in the regulatory region, and evaluates their performance

Archaic hominins include:

1. high-quality archaic genome from a 122,000-year-old Neanderthal individual found in the Altai Mountains
2. the 30x genome from a 72,000-year-old Denisovan from the Altai Mountains
3. the 30× coverage genome of an approximately 52,000-year-old Neanderthal from Croatia (the Vindija Neanderthal)

Dan Zhou & Eric Gamazon
