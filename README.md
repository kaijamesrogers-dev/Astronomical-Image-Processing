# Astronomical Image Processing — Deep Galaxy Survey (B2)

This repository contains work for the **B2: Astronomical Image Processing** third-year lab. The goal is to process a deep optical **CCD** image of an extragalactic field, detect galaxies, perform simple photometry, convert fluxes to **apparent magnitudes**, and derive **galaxy number counts** (N(<m)) for comparison with published results. 

---

## Project aims

Using a provided reduced (science-ready) optical image:

* Read a **FITS** image and relevant header metadata (especially **MAGZPT** for AB magnitude calibration). 
* Detect astronomical sources (image segmentation) and build an object catalogue. 
* Measure source brightness via (initially) **aperture photometry** and convert to **AB magnitudes**. 
* Compute **cumulative number counts** and compare against the expected Euclidean trend:

\[
\log N(m)=0.6m+\text{constant}
\]

  then discuss deviations (evolution, cosmology, incompleteness, contamination). 
* Validate code rigorously using small cutouts and/or synthetic data before running on the full mosaic. 

---

## Dataset

* Deep optical imaging of an extragalactic field observed with the **Kitt Peak 4m telescope** using a **Sloan r-band** filter (central wavelength ~620 nm). 
* The field supports the **SWIRE** Spitzer survey; the provided mosaic has already been combined and processed (“reduced”), but still contains real-world artefacts (e.g., bright star bleeding/blooming and spatially varying noise). 

---

## Recommended tools

* FITS viewing: **SAOImage DS9** (recommended in the lab guide for inspecting scaling, artefacts, and background structure). 
* Python: `astropy` for FITS I/O + `numpy` for computation (recommended). 
* Alternatives: Matlab `fitsread`, or C/C++ via CFITSIO/CCfits if you prefer compiled workflows. 

---

## Typical workflow (what this repo implements)

1. **Inspect the image**

   * Use DS9 with `zscale` and log stretch to see faint sources and diagnose artefacts/background gradients. 

2. **Load FITS + header**

   * Read pixel data and extract **MAGZPT** (instrumental zero point for AB mags in this dataset). 

3. **Background + noise estimation**

   * Histogram pixel values: expect a Gaussian-ish sky/background peak with a positive tail from sources. Use this to motivate detection thresholds (e.g., “how many σ above background?”). 

4. **Source detection (segmentation)**

   * Thresholding / connected components / masking bright-star bleed regions, etc.

5. **Photometry**

   * Aperture flux (and sky subtraction) → magnitude conversion using the header calibration.

6. **Number counts + interpretation**

   * Build (N(<m)), compare slope to 0.6m expectation, and discuss:

     * incompleteness roll-over at faint end,
     * bright extended sources losing flux outside a fixed aperture,
     * foreground star contamination,
     * genuine cosmology/evolution effects. 

---

## Reproducibility checklist

* Record your **parameter choices** (threshold σ, aperture radius, masking rules, edge cuts) and justify them using image statistics and tests on subregions/synthetic images. 
* Keep intermediate outputs (masks, segmentation maps, catalogues) so results can be regenerated and audited.

---

## Assessment context (course logistics)

This experiment sits within the third-year lab cycle structure and contributes to assessed presentation/report components for the module. Ensure your written work is **independent** and properly referenced/cited. 

---

## Academic integrity

* Your report text must be written independently; do not share write-ups between partners. It’s fine to have identical figures/numerical results if you both did the same experiment, but your explanation/structure must be your own. 

---

## Acknowledgements

Experiment guide: *Computational Image Processing: A Deep Galaxy Survey* (Imperial College third-year lab, B2 Astronomical Image Processing).  
Third-year lab organisation/assessment guidance provided by the department. 
