# CYP2B6 Pharmacogenomics in the Cloud

Star-allele calling in the Luhya population using open 1000 Genomes data — implementation, results, and the ethics of publishing them.

**Computing Bioinformatics · Group 4A · Strathmore University**

**Repository:** https://github.com/Ilhaann/cyp2b6-cloud-pharmacogenomics

| Name | Admission No. |
| --- | --- |
| Benson Gitonga | 156740 |
| Elvis Wafuke | 154210 |
| Ilhan Hamud | 169004 |
| Joy Wanjohi | 159540 |
| Melissa Ndeti | 150851 |
| Muhammad Abdalla | 169962 |

---

## What this project does

We ran a reproducible pharmacogenomics pipeline entirely on cloud infrastructure to call *CYP2B6* star alleles (drug-metabolism variants) in the Luhya of Webuye (LWK) using open-access 1000 Genomes data, compared the frequencies against four other populations, translated the results into predicted clinical phenotypes, and examined the ethics of publishing population-level genetic findings.

*CYP2B6* clears drugs including efavirenz and nevirapine (HIV), cyclophosphamide, bupropion, and methadone, so variant frequency has direct dosing implications.

## Headline results

- **99 Luhya individuals** genotyped (198 haplotypes, 4,650 variants).
- **20.2% predicted poor metabolizers** vs 6.1% European and 2.9% Han Chinese.
- **60.6% carry reduced *CYP2B6* function** (poor + intermediate metabolizers).
- ***9** is the striking allele: 35.9% of Luhya haplotypes carry the 516G>T variant; the normal-function *1 reference allele accounts for only ~45%.
- ***18** is enriched in Africa, not Luhya-specific: absent from European and Han Chinese samples but lower in the Luhya (6.1%) than the Yoruba (11.6%).
- ***6** could not be separated from *9: one of its two defining variants (rs2279343, 785A>G) is missing from the reference panel in both GRCh37 and GRCh38.

## Pipeline

Four stages, each on a managed cloud service — no stage runs on local hardware:

1. **Data access** — stream open 1000 Genomes genotypes (no download of raw data to local disk).
2. **Region import** — `pypgx.import_variants` slices the *CYP2B6* region directly (no separate bcftools stage needed).
3. **Star-allele calling** — PyPGx calls diplotypes across 508 samples / five populations.
4. **Phenotype translation** — diplotypes mapped to CPIC clinical phenotypes.

The reference panel is already haplotype-phased, so the Beagle phasing step (and its Java runtime) drops out.

**Compute environment:** Google Colab (PaaS) free tier — managed Python runtime plus shell, with `pandas`, `scipy`, and `PyPGx`. Full run ≈ 6 minutes of compute (358 seconds of star-allele calling).

## Repository contents

```
.
├── README.md
├── pipeline.py            # reproduction of the documented pipeline
├── figures/               # result figures used in the presentation
│   ├── fig1_lwk_frequencies.png
│   ├── fig2_cross_population.png
│   ├── fig3_phenotypes.png
│   ├── fig4_runlog.png
│   └── fig5_pipeline.png
├── results/               # aggregate output tables (population-level only)
└── CYP2B6_Results_Cloud_Ethics.pptx   # presentation slides
```

## Reproducing the analysis

```bash
pip install pypgx pandas scipy
python pipeline.py
```

The pipeline is reproducible from a single script. Only aggregate, population-level data is published — no individual-level predictions — consistent with the project's ethics statement.

## Ethics summary

Genetic data is sensitive personal data under Kenya's Data Protection Act 2019. Because the 1000 Genomes data is open-access and non-identifying, the Act does not bind this project directly, but the project is designed to the same standard:

- Findings are framed as statistical distributions, never fixed properties of a community.
- Only aggregate allele frequencies and phenotype distributions are released; no per-person results.
- The findings must not be used for insurance underwriting, employment screening, or any individual clinical decision without qualified genetic counselling.
- Open access so Kenyan clinicians can reach the findings without a paywall; the Luhya community is acknowledged.

## Data source

1000 Genomes Project / International Genome Sample Resource (IGSR), open-access release. No patient-identifying information is used.
