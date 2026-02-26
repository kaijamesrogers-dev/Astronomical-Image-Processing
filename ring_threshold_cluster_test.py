"""Sweep ring-threshold sigma multipliers and quantify source merging in an ROI.

This script is standalone and does not modify src.py. It mirrors the source
detection logic used there, while exposing `ring_sigma` for repeated tests.
"""

import argparse
import csv
from pathlib import Path

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.ndimage import generate_binary_structure, label
from scipy.optimize import curve_fit


def parse_args():
    parser = argparse.ArgumentParser(
        description="Test whether lowering ring threshold merges dense ROI sources."
    )
    parser.add_argument("--fits", required=True, help="Path to FITS image")
    parser.add_argument("--x1", type=int, required=True, help="ROI x start (inclusive)")
    parser.add_argument("--x2", type=int, required=True, help="ROI x end (exclusive)")
    parser.add_argument("--y1", type=int, required=True, help="ROI y start (inclusive)")
    parser.add_argument("--y2", type=int, required=True, help="ROI y end (exclusive)")
    parser.add_argument(
        "--seed_sigma", type=float, required=True, help="Seed threshold multiple k_seed"
    )
    parser.add_argument(
        "--ring_sigmas",
        required=True,
        help="Comma-separated ring threshold multiples, e.g. 3.0,2.5,2.0,1.5,1.0",
    )
    parser.add_argument("--output_dir", required=True, help="Output directory")
    return parser.parse_args()


def parse_ring_sigmas(raw):
    vals = []
    for item in raw.split(","):
        s = item.strip()
        if not s:
            continue
        vals.append(float(s))
    if not vals:
        raise ValueError("No valid values parsed from --ring_sigmas")
    return vals


def validate_roi_bounds(shape, x1, x2, y1, y2):
    if len(shape) < 2:
        raise ValueError("FITS data is not a 2D image.")
    h, w = shape[:2]
    if not (0 <= x1 < x2 <= w and 0 <= y1 < y2 <= h):
        raise ValueError(
            f"Invalid ROI bounds for image shape (h={h}, w={w}): "
            f"x[{x1}:{x2}], y[{y1}:{y2}]"
        )


def fit_gaussian_like_src(image_data, xmin=3300, xmax=3550):
    """Background fit matching src.py behavior, without plotting."""
    masked_data = image_data[(image_data >= xmin) & (image_data <= xmax)]
    if masked_data.size == 0:
        raise ValueError(
            f"No pixels in background fit range [{xmin}, {xmax}] inside ROI."
        )

    bin_centres = np.arange(xmin, xmax + 1)
    bin_edges = np.arange(xmin - 0.5, xmax + 1.5, 1)
    counts, _ = np.histogram(masked_data, bins=bin_edges)
    if np.max(counts) <= 0:
        raise ValueError("Background histogram is empty or degenerate.")

    def gaussian(x, amp, mu, sigma):
        return amp * np.exp(-(x - mu) ** 2 / (2 * sigma**2))

    amp0 = counts.max()
    mu0 = bin_centres[np.argmax(counts)]
    sigma0 = float(np.std(masked_data))
    if sigma0 <= 0:
        sigma0 = 1.0

    popt, _ = curve_fit(gaussian, bin_centres, counts, p0=[amp0, mu0, sigma0], maxfev=10000)
    _, mu, sigma = popt
    return float(mu), float(abs(sigma))


def detect_sources_like_src(image_data, mu, sigma, seed_sigma=3.0, ring_sigma=2.0):
    """Mirror src.py detection logic, with ring_sigma exposed for sweeps."""
    seed_threshold = mu + seed_sigma * sigma
    ring_threshold = mu + ring_sigma * sigma
    min_radius = 1

    mask = np.zeros(image_data.shape, dtype=bool)
    sources = []

    height, width = image_data.shape
    while True:
        masked_data = np.copy(image_data).astype(float)
        masked_data[mask] = -np.inf

        max_value = np.max(masked_data)
        if max_value < seed_threshold:
            break

        max_index = np.argmax(masked_data)
        y, x = np.unravel_index(max_index, masked_data.shape)

        found_radius = 0
        for r in range(1, 101):
            ymin = max(0, y - r)
            ymax = min(height, y + r + 1)
            xmin = max(0, x - r)
            xmax = min(width, x + r + 1)

            cutout = masked_data[ymin:ymax, xmin:xmax]
            yy, xx = np.ogrid[ymin:ymax, xmin:xmax]
            dist = np.sqrt((xx - x) ** 2 + (yy - y) ** 2)

            in_ring = (dist > r - 1) & (dist <= r)
            if np.count_nonzero(in_ring) == 0:
                break

            ring_median = np.median(cutout[in_ring])
            if ring_median < ring_threshold:
                found_radius = r - 1
                break
        else:
            found_radius = r

        if found_radius <= min_radius:
            mask[y, x] = True
            continue

        radius = int(1.5 * found_radius)
        if radius <= 6:
            radius = 6

        sources.append((int(x), int(y), float(max_value), int(radius)))

        ymin2 = max(0, y - radius)
        ymax2 = min(height, y + radius + 1)
        xmin2 = max(0, x - radius)
        xmax2 = min(width, x + radius + 1)

        yy2, xx2 = np.ogrid[0 : (ymax2 - ymin2), 0 : (xmax2 - xmin2)]
        dist2 = np.sqrt((xx2 - (x - xmin2)) ** 2 + (yy2 - (y - ymin2)) ** 2)
        mask[ymin2:ymax2, xmin2:xmax2][dist2 <= radius] = True

    return sources


def sources_to_mask(shape, sources):
    mask = np.zeros(shape, dtype=bool)
    h, w = shape
    for x, y, _, r in sources:
        ymin = max(0, y - r)
        ymax = min(h, y + r + 1)
        xmin = max(0, x - r)
        xmax = min(w, x + r + 1)
        yy, xx = np.ogrid[ymin:ymax, xmin:xmax]
        dist = np.sqrt((xx - x) ** 2 + (yy - y) ** 2)
        mask[ymin:ymax, xmin:xmax][dist <= r] = True
    return mask


def mask_component_metrics(mask):
    total_masked = int(np.count_nonzero(mask))
    if total_masked == 0:
        return 0, 0.0

    structure = generate_binary_structure(2, 2)
    labels, n_components = label(mask, structure=structure)

    max_area = 0
    for comp_id in range(1, n_components + 1):
        area = int(np.count_nonzero(labels == comp_id))
        if area > max_area:
            max_area = area

    largest_fraction = float(max_area / total_masked)
    return int(n_components), largest_fraction


def save_overlay_png(roi, mask, ring_sigma, output_dir):
    p1, p99 = np.percentile(roi, [1, 99])
    if p99 <= p1:
        p1, p99 = float(np.min(roi)), float(np.max(roi) + 1.0)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.imshow(roi, cmap="gray", origin="lower", vmin=p1, vmax=p99)
    ax.imshow(mask.astype(float), cmap="autumn", origin="lower", alpha=0.35, vmin=0.0, vmax=1.0)
    ax.set_title(f"ROI + Mask Overlay (ring_sigma={ring_sigma:.3g})")
    ax.set_xlabel("x (ROI pixels)")
    ax.set_ylabel("y (ROI pixels)")
    fig.tight_layout()
    out_path = output_dir / f"overlay_ring_{str(ring_sigma).replace('.', 'p')}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def save_summary_csv(rows, output_dir):
    out_path = output_dir / "ring_threshold_summary.csv"
    fieldnames = [
        "ring_sigma",
        "N_sources_detected",
        "N_components",
        "largest_component_fraction",
        "masked_area_fraction",
        "total_masked_pixels",
        "roi_pixels",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def save_trend_plots(rows, output_dir):
    ring = np.array([r["ring_sigma"] for r in rows], dtype=float)
    nsrc = np.array([r["N_sources_detected"] for r in rows], dtype=float)
    ncomp = np.array([r["N_components"] for r in rows], dtype=float)
    largest = np.array([r["largest_component_fraction"] for r in rows], dtype=float)
    area_frac = np.array([r["masked_area_fraction"] for r in rows], dtype=float)

    order = np.argsort(ring)[::-1]
    ring, nsrc, ncomp, largest, area_frac = (
        ring[order],
        nsrc[order],
        ncomp[order],
        largest[order],
        area_frac[order],
    )

    fig1, ax1 = plt.subplots(figsize=(7, 5))
    ax1.plot(ring, nsrc, marker="o", color="tab:blue", label="N_sources_detected")
    ax1.set_xlabel("ring_sigma")
    ax1.set_ylabel("N_sources_detected", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax1.grid(True, alpha=0.3)
    ax2 = ax1.twinx()
    ax2.plot(ring, ncomp, marker="s", color="tab:red", label="N_components")
    ax2.set_ylabel("N_components", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    ax1.set_title("ring_sigma vs source/component counts")
    fig1.tight_layout()
    fig1.savefig(output_dir / "ring_vs_counts.png", dpi=150)
    plt.close(fig1)

    fig2, ax = plt.subplots(figsize=(7, 5))
    ax.plot(ring, largest, marker="o", label="largest_component_fraction")
    ax.plot(ring, area_frac, marker="s", label="masked_area_fraction")
    ax.set_xlabel("ring_sigma")
    ax.set_ylabel("Fraction")
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title("ring_sigma vs area fractions")
    fig2.tight_layout()
    fig2.savefig(output_dir / "ring_vs_fractions.png", dpi=150)
    plt.close(fig2)


def find_knee(rows):
    ordered = sorted(rows, key=lambda r: r["ring_sigma"], reverse=True)
    for i, cur in enumerate(ordered):
        k = cur["ring_sigma"]

        if cur["largest_component_fraction"] > 0.5:
            return (
                k,
                "largest_component_fraction > 0.5",
            )

        if i == 0:
            continue

        prev = ordered[i - 1]
        prev_nsrc = prev["N_sources_detected"]
        prev_ncomp = prev["N_components"]

        if prev_nsrc > 0:
            drop_sources = (prev_nsrc - cur["N_sources_detected"]) / prev_nsrc
            if drop_sources >= 0.2:
                return (
                    k,
                    f"N_sources_detected drop = {100*drop_sources:.1f}% vs next-higher ring_sigma",
                )

        if prev_ncomp > 0:
            drop_comp = (prev_ncomp - cur["N_components"]) / prev_ncomp
            if drop_comp >= 0.2:
                return (
                    k,
                    f"N_components drop = {100*drop_comp:.1f}% vs next-higher ring_sigma",
                )

    return None, None


def main():
    args = parse_args()
    ring_sigmas = parse_ring_sigmas(args.ring_sigmas)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with fits.open(args.fits) as hdul:
        data = np.asarray(hdul[0].data)

    validate_roi_bounds(data.shape, args.x1, args.x2, args.y1, args.y2)
    roi = data[args.y1 : args.y2, args.x1 : args.x2]
    roi_pixels = int(roi.size)

    try:
        mu, sigma = fit_gaussian_like_src(roi)
        bg_method = "gaussian_fit_3300_3550"
    except Exception:
        mu, sigma = float(np.mean(roi)), float(np.std(roi))
        if sigma <= 0:
            sigma = 1e-6
        bg_method = "mean_std_fallback"

    print(f"ROI shape: {roi.shape}, pixels={roi_pixels}")
    print(f"Background ({bg_method}): mu={mu:.3f}, sigma={sigma:.3f}")
    print(f"Seed sigma fixed at: {args.seed_sigma:.3f}")
    print(f"Ring sigma sweep: {ring_sigmas}")

    rows = []
    for ring_sigma in ring_sigmas:
        sources = detect_sources_like_src(
            roi,
            mu,
            sigma,
            seed_sigma=args.seed_sigma,
            ring_sigma=ring_sigma,
        )
        mask = sources_to_mask(roi.shape, sources)
        n_components, largest_fraction = mask_component_metrics(mask)
        total_masked = int(np.count_nonzero(mask))
        masked_area_fraction = float(total_masked / roi_pixels) if roi_pixels > 0 else 0.0

        row = {
            "ring_sigma": float(ring_sigma),
            "N_sources_detected": int(len(sources)),
            "N_components": int(n_components),
            "largest_component_fraction": float(largest_fraction),
            "masked_area_fraction": float(masked_area_fraction),
            "total_masked_pixels": int(total_masked),
            "roi_pixels": int(roi_pixels),
        }
        rows.append(row)
        save_overlay_png(roi, mask, ring_sigma, output_dir)

    save_summary_csv(rows, output_dir)
    save_trend_plots(rows, output_dir)

    knee_sigma, reason = find_knee(rows)
    if knee_sigma is None:
        print("Knee result: no merge knee detected by the requested criteria.")
    else:
        print(f"Knee result: merging starts at ring_sigma={knee_sigma:.3g} (criterion: {reason}).")

    print(f"Outputs written to: {output_dir}")


if __name__ == "__main__":
    main()
