# Web tool

There is a hosted web front-end for `vflank` — an editable grid where you paste
one variant (or a small batch) and get back the masked target sequence, no
install required:

[Open the web app :material-open-in-new:](https://vflank-webapp.onrender.com/){ .md-button .md-button--primary }

It is a thin layer over this library: the browser app calls
`vflank.run_small` / `run_fusion` and renders the result. All the science —
flank extraction, junction building, SNP masking — lives here in `vflank`. The
web tool is scoped to the **no-download path** (UCSC reference + gnomAD APIs,
modes A/B): single variants or tiny batches (≤10), **no BAM, no PHI**. For BAM
consensus, large batches, or scripting, use the [CLI](reference/cli.md) or the
[Python API](reference/api.md).

!!! note "First load may be slow"
    The app is hosted on a free tier that sleeps after ~15 minutes of
    inactivity, so the first request after a while can take ~30 s to wake. Once
    it is up it is responsive.

## Small variant

Paste a variant (here, BRAF V600E on GRCh37), press **Run**, and `vflank`
pulls the reference from UCSC and common SNPs from gnomAD — nothing is
downloaded — and returns the flanking sequence with the variant and common
SNPs masked. Download the FASTA or a Primer3 input.

<div class="video-embed">
  <video controls preload="metadata" playsinline poster="">
    <source src="assets/small-variant-demo.mp4" type="video/mp4">
    Your browser does not support embedded video.
    <a href="assets/small-variant-demo.mp4">Download the demo (MP4)</a>.
  </video>
</div>

## Fusion

Switch to **Fusion** and enter the two breakpoints (here, EML4–ALK). `vflank`
builds the chimeric junction a probe spans, masks common SNPs in the flanks so
the probe avoids them, and returns the junction sequence to download.

<div class="video-embed">
  <video controls preload="metadata" playsinline poster="">
    <source src="assets/fusion-demo.mp4" type="video/mp4">
    Your browser does not support embedded video.
    <a href="assets/fusion-demo.mp4">Download the demo (MP4)</a>.
  </video>
</div>

---

The web app is a separate repository —
[rhshah/vFlank-webapp](https://github.com/rhshah/vFlank-webapp) — and pins a
released `vflank`. See its README to run it locally or deploy your own.
