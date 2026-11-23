"""
CYP2B6 pharmacogenomics pipeline — cloud reproduction script.

Calls CYP2B6 star alleles in the Luhya (LWK) and four comparator populations
from open 1000 Genomes genotypes, translates diplotypes into CPIC metabolizer
phenotypes, and writes aggregate (population-level only) result tables.

Environment: Google Colab (PaaS). Managed Python runtime + shell.
    pip install pypgx pandas scipy statsmodels

Data: 1000 Genomes Project / IGSR open-access release (GRCh38).
No individual-level output is written — only aggregate frequencies and
phenotype distributions, consistent with the project's ethics statement.
"""

from pathlib import Path
import pandas as pd
import pypgx
from statsmodels.stats.proportion import proportion_confint

GENE = "CYP2B6"
ASSEMBLY = "GRCh38"

# IGSR 1000 Genomes population codes used in this study (5 populations, 508 samples).
# Luhya plus two other African, one European, one East Asian population.
# Confirm these match the exact codes used in your run before publishing.
POPULATIONS = {
    "LWK": "Luhya (Webuye, Kenya)",
    "YRI": "Yoruba (Ibadan, Nigeria)",
    "ESN": "Esan (Nigeria)",
    "GBR": "British (England and Scotland)",
    "CHB": "Han Chinese (Beijing)",
}

OUT = Path("results")
OUT.mkdir(exist_ok=True)


def call_alleles(vcf_path: str):
    """Stage 1-3: import the CYP2B6 region, predict and call star alleles.

    pypgx.import_variants performs the region slice itself, so no separate
    bcftools stage is needed. The reference panel is already haplotype-phased,
    so the Beagle phasing step is skipped.
    """
    archive = pypgx.import_variants(GENE, vcf_path, assembly=ASSEMBLY)
    alleles = pypgx.predict_alleles(archive)
    genotypes = pypgx.call_genotypes(alleles=alleles)
    phenotypes = pypgx.call_phenotypes(genotypes)
    return genotypes, phenotypes


def wilson_ci(count, total):
    """95% Wilson confidence interval for a proportion."""
    low, high = proportion_confint(count, total, alpha=0.05, method="wilson")
    return round(low * 100, 1), round(high * 100, 1)


def allele_frequency_table(genotypes: pd.DataFrame, sample_pop: pd.Series):
    """Table: per-population star-allele frequencies with 95% Wilson CIs."""
    rows = []
    # genotypes carry two haplotype calls per sample; count alleles per population
    haplo = []
    for sample, row in genotypes.iterrows():
        pop = sample_pop.get(sample)
        for a in str(row["Genotype"]).split("/"):
            haplo.append((pop, a))
    hdf = pd.DataFrame(haplo, columns=["pop", "allele"])
    for pop, sub in hdf.groupby("pop"):
        n = len(sub)
        for allele, cnt in sub["allele"].value_counts().items():
            low, high = wilson_ci(cnt, n)
            rows.append({
                "population": pop,
                "star_allele": allele,
                "count": int(cnt),
                "haplotypes": n,
                "frequency_pct": round(cnt / n * 100, 1),
                "ci95_low": low,
                "ci95_high": high,
            })
    return pd.DataFrame(rows).sort_values(["population", "frequency_pct"],
                                          ascending=[True, False])


def phenotype_table(phenotypes: pd.DataFrame, sample_pop: pd.Series):
    """Table: per-population CPIC metabolizer-phenotype distribution."""
    df = phenotypes.copy()
    df["population"] = df.index.map(sample_pop)
    out = (df.groupby(["population", "Phenotype"]).size()
             .rename("count").reset_index())
    totals = out.groupby("population")["count"].transform("sum")
    out["percent"] = (out["count"] / totals * 100).round(1)
    return out


def main(vcf_path: str, panel_path: str):
    # panel_path: IGSR sample->population mapping (sample, pop, super_pop)
    panel = pd.read_csv(panel_path, sep="\t")
    sample_pop = panel.set_index("sample")["pop"]

    genotypes, phenotypes = call_alleles(vcf_path)

    freq = allele_frequency_table(genotypes, sample_pop)
    pheno = phenotype_table(phenotypes, sample_pop)

    freq.to_csv(OUT / "table1_allele_frequencies.tsv", sep="\t", index=False)
    pheno.to_csv(OUT / "table2_phenotype_distribution.tsv", sep="\t", index=False)

    # Limitations register — no-calls and indeterminate diplotypes are reported,
    # never dropped from the denominator.
    indeterminate = phenotypes[phenotypes["Phenotype"].str.contains(
        "Indeterminate", case=False, na=False)]
    indeterminate.to_csv(OUT / "table3_indeterminate_calls.tsv", sep="\t")

    print("Wrote aggregate tables to", OUT.resolve())
    print(freq.head(20).to_string(index=False))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True,
                    help="Phased 1000 Genomes VCF covering the CYP2B6 region")
    ap.add_argument("--panel", required=True,
                    help="IGSR integrated sample panel (sample/pop/super_pop TSV)")
    args = ap.parse_args()
    main(args.vcf, args.panel)
