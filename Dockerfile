# syntax=docker/dockerfile:1

# --- Build the wheel ---
FROM python:3.12-slim AS builder
WORKDIR /src
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir build && python -m build --wheel --outdir /dist

# --- Runtime ---
FROM python:3.12-slim
LABEL org.opencontainers.image.source="https://github.com/rhshah/vFlank" \
      org.opencontainers.image.description="Variant-aware flanking-sequence extraction and masking for ddPCR assay design" \
      org.opencontainers.image.licenses="Apache-2.0"

# pysam ships manylinux wheels with htslib bundled, so no system htslib needed.
COPY --from=builder /dist/*.whl /tmp/
RUN pip install --no-cache-dir /tmp/*.whl && rm -rf /tmp/*.whl

# Reference FASTAs / VCFs are mounted at runtime, e.g.:
#   docker run --rm -v "$PWD:/data" ghcr.io/rhshah/vflank \
#       small run /data/variants.maf -r /data/GRCh37.fasta -g hg19 -o /data/out.fasta
ENTRYPOINT ["vflank"]
CMD ["--help"]
