# SNP Masking

Common germline SNPs under a primer or probe cause allele dropout and assay
failure. vflank masks them to `N` in the `Masked__` records so a downstream
designer avoids those positions. Both `small` and `fusion` share the same
options.

## Two backends — `--pop-source`

=== "Local VCFs (default)"

    ```bash
    vflank small run variants.maf -r GRCh37.fasta -g hg19 \
        --pop-source vcf --pop-vcf-dir gnomad_v2.1.1/
    ```

    Reproducible, offline, unlimited scale. Verify coverage first:

    ```bash
    vflank small list-vcf gnomad_v2.1.1/ -g hg19 --pop-data both
    ```

=== "gnomAD API (no download)"

    ```bash
    vflank small run variants.maf -r GRCh37.fasta -g hg19 --pop-source api
    ```

    Queries <https://gnomad.broadinstitute.org/api> per flank region — no data
    to download. Rate-limited to ~10 requests/min, so best for **small cohorts**;
    identical variants are cached to one request.

| | Local VCF | gnomAD API |
|---|---|---|
| Setup | download per-chromosome VCFs | none |
| Reproducible / offline | ✅ | ❌ (network; data updates) |
| Scale | unlimited | small cohorts (rate limit) |

## Which data — `--pop-data {genome,exome,both}`

Default **`genome`**. Flanks often fall in non-coding regions where only genomes
have data; `exome` adds power in coding regions; `both` masks the **union** (a
position is masked if it is a common SNP in *either* cohort).

A position is masked when a **single-base** gnomAD SNP there has
`max(AF, AF_grpmax) ≥ --af-threshold` (default `0.001` = 0.1%). Indels are not
masked. Requesting `exome`/`both` from the VCF backend without the exome files
fails fast — never a silent genome-only fallback.

## Genome builds

| `--genome-build` | Reference | gnomAD |
|---|---|---|
| `hg19` (default) | GRCh37 | v2.1.1 |
| `hg38` | GRCh38 | v4.1 |

gnomAD v4 is GRCh38-only — there is no v4 GRCh37.

## Downloading gnomAD VCFs (GRCh37 example)

```bash
mkdir -p gnomad_v2.1.1 && cd gnomad_v2.1.1
BASE=https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes
for CHR in {1..22} X; do
  wget -c "$BASE/gnomad.genomes.r2.1.1.sites.${CHR}.vcf.bgz"
  wget -c "$BASE/gnomad.genomes.r2.1.1.sites.${CHR}.vcf.bgz.tbi"
done
```

Exome files use `gnomad.exomes.*` (needed only for `--pop-data exome`/`both`).
Download only the chromosomes your variants are on to save space.
