#!/usr/bin/env python3
"""
simulation — Generate synthetic mosaic.fits + truth_catalogue.csv

Run this script first, then run src.py as normal.

Ground truth is printed to stdout when this script runs.
Every injected source is recorded in truth_catalogue.csv so you can
cross-match the src.py catalogue against the known input.
"""

import numpy as np
from astropy.io import fits

# =============================================================================
# Ground-truth parameters — change these to explore different scenarios
# =============================================================================

RNG_SEED    = 42              # fix for full reproducibility

# Image geometry
IMG_H       = 2000            # total image rows
IMG_W       = 1500            # total image cols
# Crop must match the slices in src.py (line 21-22)
CROP_R      = slice(1050, 1700)   # 650 rows
CROP_C      = slice(390,  1300)   # 910 cols
CROP_H      = CROP_R.stop - CROP_R.start   # 650
CROP_W      = CROP_C.stop - CROP_C.start   # 910

# Background — matched to the real mosaic.fits statistics (Lab Book, Lab 3)
BG_MU       = 3421.0          # mean background  (ADU)
BG_SIGMA    = 20.0            # background noise (ADU)

# PSF — circular Gaussian
PSF_SIGMA   = 1.5             # pixels  (FWHM ≈ 3.53 px)

# Photometric calibration (stored in the FITS header; read by src.py)
MAGZPT      = 26.0
MAGZRR      = 0.02

# Plate scale (must match src.py, line 331)
PIX_SCALE   = 0.258           # arcsec / pixel

# Number-count normalisation: log10[N(<m) / deg²] = 0.6·m + C_NC
# Chosen so that N(<18) ≈ 200 detectable sources in the cropped field.
C_NC        = -5.98

# Magnitude range of injected sources
MAG_MIN     = 10.0            # brightest
MAG_MAX     = 18.5            # fainter than detection limit (good stress test)

# Keep sources ≥ MARGIN pixels from crop edge so PSF tails stay inside
MARGIN      = 15

OUTPUT_FITS      = "sim_mosaic.fits"
OUTPUT_TRUTH_CSV = "truth_catalogue.csv"

# =============================================================================
# Derived quantities (printed as ground truth)
# =============================================================================

area_deg2  = CROP_H * CROP_W * PIX_SCALE**2 / 3600**2

# Total sources to inject (from the number-count law integrated to MAG_MAX)
N_total    = int(round(10**(0.6 * MAG_MAX + C_NC) * area_deg2))

# Detection thresholds that src.py will compute from the fit_gaussian() output
SEED_THR   = BG_MU + 3.0 * BG_SIGMA    # seed pixel must exceed this
RING_THR   = BG_MU + 2.5 * BG_SIGMA    # ring median must exceed this
MIN_RAD    = 2                          # min_radius in src.py

# Faint limit: ring at r=2 must be above RING_THR
#   peak_above_bg × <exp(-r²/2σ²)>_{ring r=2} > 2.5·σ_bg
#   <exp> over annulus 1<r≤2  ≈  ∫₁²exp(-r²/2σ²)·r dr / ∫₁²r dr
num2   = PSF_SIGMA**2 * (np.exp(-1/(2*PSF_SIGMA**2)) - np.exp(-4/(2*PSF_SIGMA**2)))
denom2 = 1.5    # ∫₁² r dr
ring2_factor = num2 / denom2
flux_min_r2 = 2.5 * BG_SIGMA / ring2_factor * (2*np.pi*PSF_SIGMA**2)
MAG_LIM    = MAGZPT - 2.5 * np.log10(flux_min_r2)

# =============================================================================
# Sample source magnitudes (inverse-CDF, exact 0.6 slope, no Poisson noise)
# =============================================================================

rng = np.random.default_rng(RNG_SEED)

def sample_mags(n, m_lo, m_hi):
    """Draw n magnitudes from dN/dm ∝ 10^(0.6·m)."""
    u  = rng.uniform(0, 1, n)
    lo = 10**(0.6 * m_lo)
    hi = 10**(0.6 * m_hi)
    return np.log10(u * (hi - lo) + lo) / 0.6

mags   = sample_mags(N_total, MAG_MIN, MAG_MAX)
fluxes = 10**((MAGZPT - mags) / 2.5)    # total counts per source

# Source positions: random inside the crop with MARGIN buffer
r_lo, r_hi = CROP_R.start + MARGIN, CROP_R.stop  - MARGIN
c_lo, c_hi = CROP_C.start + MARGIN, CROP_C.stop  - MARGIN

src_row = rng.integers(r_lo, r_hi, N_total)   # full-image row  → y in src.py
src_col = rng.integers(c_lo, c_hi, N_total)   # full-image col  → x in src.py

# =============================================================================
# Build image: background + PSF-convolved stars
# =============================================================================

image = rng.normal(BG_MU, BG_SIGMA, (IMG_H, IMG_W)).astype(np.float32)

PSF_R = int(np.ceil(5 * PSF_SIGMA))    # truncate PSF at 5σ (≈ 8 px)

for sr, sc, flux in zip(src_row, src_col, fluxes):
    r0, r1 = max(0, sr - PSF_R), min(IMG_H, sr + PSF_R + 1)
    c0, c1 = max(0, sc - PSF_R), min(IMG_W, sc + PSF_R + 1)
    rr, cc  = np.ogrid[r0:r1, c0:c1]
    psf     = np.exp(-((rr - sr)**2 + (cc - sc)**2) / (2 * PSF_SIGMA**2))
    psf    /= psf.sum()          # normalise → exactly 'flux' counts added
    image[r0:r1, c0:c1] += flux * psf

# =============================================================================
# Write FITS (drop-in replacement for mosaic.fits)
# =============================================================================

hdr = fits.Header()
hdr['MAGZPT']  = (MAGZPT,  'Photometric zero-point')
hdr['MAGZRR']  = (MAGZRR,  'Zero-point uncertainty (mag)')
hdr['BUNIT']   = ('ADU',   'Pixel units')
hdr['COMMENT'] = 'Simulated image - ground truth in truth_catalogue.csv'

fits.PrimaryHDU(image, header=hdr).writeto(OUTPUT_FITS, overwrite=True)
print(f"Written: {OUTPUT_FITS}  ({IMG_H} × {IMG_W} px, float32)")

# =============================================================================
# Write truth catalogue
# Coordinates are in the CROPPED frame (as src.py reports them)
# =============================================================================

x_crop = src_col - CROP_C.start   # col in cropped image  (src.py 'x')
y_crop = src_row - CROP_R.start   # row in cropped image  (src.py 'y')

truth = np.column_stack([x_crop, y_crop, mags, fluxes])
np.savetxt(
    OUTPUT_TRUTH_CSV, truth,
    header='x_crop  y_crop  mag_true  flux_true',
    fmt=['%d', '%d', '%.4f', '%.4f'],
    comments='# ',
)
print(f"Written: {OUTPUT_TRUTH_CSV}  ({N_total} sources)")

# =============================================================================
# Print full ground-truth summary
# =============================================================================

print()
print("=" * 65)
print("  SIMULATION GROUND TRUTH")
print("=" * 65)
print(f"  Random seed           : {RNG_SEED}")
print(f"  Image size            : {IMG_H} × {IMG_W} pixels")
print(f"  Crop applied in src.py: rows {CROP_R}, cols {CROP_C}")
print(f"                          → {CROP_H} × {CROP_W} px")
print(f"  Field area            : {area_deg2:.5f} deg²  "
      f"({CROP_H*CROP_W*PIX_SCALE**2:.0f} arcsec²)")
print()
print(f"  Background  μ         : {BG_MU} ADU")
print(f"  Background  σ         : {BG_SIGMA} ADU")
print(f"    fit_gaussian() range  [3300, 3550] spans ±6σ → perfect Gaussian fit")
print()
print(f"  PSF σ                 : {PSF_SIGMA} px  (FWHM = {2.355*PSF_SIGMA:.2f} px)")
print()
print(f"  MAGZPT                : {MAGZPT}")
print(f"  MAGZRR                : {MAGZRR}")
print(f"  Theoretical slope     : 0.6  (exact, by construction)")
print(f"  Normalisation C_NC    : {C_NC}")
print()
print(f"  Total sources injected: {N_total}  (m = {MAG_MIN} to {MAG_MAX})")
print()
print("  Cumulative number counts:")
print(f"  {'m_lim':>6}  {'N_injected':>11}  {'N_theory':>9}"
      f"  {'log10(N/deg²)_sim':>20}  {'log10(N/deg²)_theory':>21}")
print(f"  {'-'*6}  {'-'*11}  {'-'*9}  {'-'*20}  {'-'*21}")
for mlim in [14.0, 15.0, 16.0, 17.0, 17.5, 18.0, 18.5]:
    n_inj    = int(np.sum(mags <= mlim))
    n_theory = 10**(0.6 * mlim + C_NC) * area_deg2
    l10_sim  = np.log10(n_inj / area_deg2) if n_inj > 0 else float('nan')
    l10_th   = 0.6 * mlim + C_NC
    print(f"  {mlim:>6.1f}  {n_inj:>11d}  {n_theory:>9.1f}"
          f"  {l10_sim:>20.3f}  {l10_th:>21.3f}")
print()
print("  Detection thresholds (src.py will recover these from the image):")
print(f"    seed threshold  = μ + 3σ   = {SEED_THR:.1f} ADU  "
      f"→ flags candidate peaks")
print(f"    ring threshold  = μ + 2.5σ = {RING_THR:.1f} ADU  "
      f"→ stops ring expansion")
print(f"    min_radius      = {MIN_RAD} px    → rejects hot pixels")
print(f"    Estimated faint limit      ≈ m = {MAG_LIM:.2f}")
print()
print("  Expected peak value (centre pixel) for stars of various magnitudes:")
print(f"  {'mag':>5}  {'flux':>8}  {'peak_px':>9}  {'seed_ok':>8}  "
      f"{'ring_r2_ok':>11}  {'est_radius':>11}")
print(f"  {'-'*5}  {'-'*8}  {'-'*9}  {'-'*8}  {'-'*11}  {'-'*11}")
for m in [13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 18.5]:
    f       = 10**((MAGZPT - m) / 2.5)
    peak_ab = f / (2 * np.pi * PSF_SIGMA**2)
    peak_px = BG_MU + peak_ab
    seed_ok = peak_px > SEED_THR
    r2_ab   = peak_ab * ring2_factor
    r2_ok   = (BG_MU + r2_ab) > RING_THR
    # rough radius: find r where ring mean hits RING_THR
    #   peak_ab × <exp>_{ring r} = 2.5·σ_bg
    #   estimate using ring midpoint r ≈ r - 0.5
    est_rad = int(np.sqrt(-2*PSF_SIGMA**2 * np.log(
        2.5*BG_SIGMA / peak_ab / ring2_factor * 0.584)) + 1) if peak_ab > 0 else 0
    print(f"  {m:>5.1f}  {f:>8.0f}  {peak_px:>9.1f}  "
          f"{'Yes' if seed_ok else 'No':>8}  "
          f"{'Yes' if r2_ok else 'No':>11}  "
          f"{'~' + str(max(0, est_rad)) + ' px':>11}")
print()
print("  Note: magnitudes reported by src.py will be ~0.05–0.3 mag fainter")
print("  than mag_true due to aperture losses (flux outside the measured radius).")
print("  This is a systematic offset — the slope 0.6 is preserved.")
print("=" * 65)