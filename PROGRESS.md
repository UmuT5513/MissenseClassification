# Project Progress

Last Updated: 2026-03-16

## Status Summary

- Current phase:
- Next milestone:
- Blocking issues:

## Milestones

- [ ] Confirm target genes list and panel membership
- [ ] Pin ClinVar filtering rules and label mapping
- [ ] Download and load ClinVar variant_summary.txt.gz
- [ ] Filter ClinVar by ReviewStatus, type, consequence, genes, labels
- [ ] Build minimal VCF (CHROM, POS, REF, ALT)
- [ ] Run offline VEP with dbNSFP and --af_gnomad
- [ ] Load and validate VEP output TSV
- [ ] Extract nucleotide context (±5 nt) from GRCh38 FASTA
- [ ] Extract amino-acid context (±5 AA) from idmapping_2026_03_16.fasta
- [ ] Compute biochemical features (delta scales, Grantham)
- [ ] Assemble feature table and handle missing values
- [ ] Merge all features and drop genomic address columns
- [ ] EDA: class balance and feature distributions
- [ ] Preprocessing: encoding, scaling, train/test split
- [ ] Train baseline model(s) and evaluate metrics
- [ ] Hyperparameter tuning and model selection
- [ ] Feature importance and interpretability (SHAP or built-in)
- [ ] Final report, model export, and notebook cleanup

## Work Log

- 2026-03-16:
  - 

## Decisions

- 

## Risks / Open Questions

- ClinVar strict ReviewStatus filters may yield a small dataset
  - Mitigation: track class counts after filtering; allow fallback to 2-star if counts are too low (only if approved).
- VEP offline setup and dbNSFP availability can delay progress
  - Mitigation: confirm cache + dbNSFP paths early; keep a small test VCF to validate VEP before full run.
- Protein position mapping errors for long proteins (e.g., BRCA2)
  - Mitigation: validate protein position parsing with spot checks; flag variants with unmapped positions.
- Missing gnomAD AF for some variants; need consistent imputation
  - Mitigation: impute Global MAF = 0 and add a boolean missing flag.
