"""
Cut individual crystals out of powder micrographs and emit them as sprites the
simulator can scatter.

    python make_sprites.py shot1.png shot2.png -o grains.js

Why sprites: procedurally drawn fragments always read as cartoons because real
milled crystals have conchoidal edges, cleavage steps and surface texture that
no small set of rules reproduces. Borrowing the real outlines fixes that.

How it works: the background in a powder micrograph is darker than the
particles, so an Otsu threshold separates them. Touching particles are split
with a distance-transform watershed. Each surviving particle is cut out with a
feathered alpha matte, trimmed, normalised for brightness and written as a
grayscale+alpha PNG — a quarter the bytes of RGBA, and powder is colourless
anyway.

Particles are rejected if they touch the image border (cut off), are too small
to carry detail, contain large dark pockets (background leaking through a hole
in the matte), or share most of their outline with a neighbour — the last one
matters because a watershed cut is a straight artificial edge, and a sprite with
one of those composites like a sliced-off crop rather than a crystal.

LICENSING: whatever goes in comes out embedded in the simulator and ships to
your students. Use micrographs you own or that carry a licence allowing it.
"""

from __future__ import annotations

import argparse
import base64
import io
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops
from skimage.morphology import closing, disk, remove_small_holes
from skimage.segmentation import watershed


def segment(gray: np.ndarray, min_area: int, split: bool) -> np.ndarray:
    """Label map of candidate particles."""
    thr = threshold_otsu(gray)
    fg = gray > thr
    fg = remove_small_holes(fg, max_size=max(32, min_area // 4))
    fg = closing(fg, disk(2))

    if not split:
        return label(fg)

    # particles touch, so split them where the distance transform says there
    # are separate cores
    dist = ndi.distance_transform_edt(fg)
    coords = peak_local_max(dist, labels=fg, min_distance=7,
                            footprint=np.ones((5, 5)), exclude_border=False)
    markers = np.zeros(dist.shape, dtype=int)
    for i, (r, c) in enumerate(coords, start=1):
        markers[r, c] = i
    if markers.max() == 0:
        return label(fg)
    return watershed(-dist, markers, mask=fg)


def feather(mask: np.ndarray, px: float) -> np.ndarray:
    """Soften the cut edge so sprites do not show a hard pixel staircase."""
    a = ndi.gaussian_filter(mask.astype(np.float32), px)
    a = np.clip((a - 0.35) / 0.4, 0, 1)
    return a


def extract(path: Path, *, min_area: int, max_particles: int,
            sprite_px: int, split: bool, crop_bottom: int,
            max_shared: float, max_dark: float) -> list:
    im = Image.open(path).convert("L")
    arr = np.array(im)
    if crop_bottom:                      # drop instrument banners
        arr = arr[:-crop_bottom]

    labels = segment(arr, min_area, split)
    h, w = arr.shape
    out = []

    for reg in regionprops(labels):
        if reg.area < min_area:
            continue
        r0, c0, r1, c1 = reg.bbox
        if r0 == 0 or c0 == 0 or r1 >= h or c1 >= w:
            continue                      # clipped by the frame
        if reg.solidity > 0.97:
            continue                      # suspiciously blobby
        if reg.solidity < 0.45:
            continue                      # a shredded mis-segmentation
        ar = reg.axis_major_length / max(reg.axis_minor_length, 1e-6)
        if ar > 7:
            continue                      # a smear, not a crystal

        mask = (labels[r0:r1, c0:c1] == reg.label)
        patch = arr[r0:r1, c0:c1].astype(np.float32)

        # A watershed cut is a straight artificial edge. Grow the mask and see
        # how much of the collar lands on a different particle: a lot means this
        # outline is mostly invention, so drop it.
        grown = ndi.binary_dilation(mask, disk(2))
        collar = grown & ~mask
        neigh = labels[r0:r1, c0:c1][collar]
        if collar.sum():
            shared = np.count_nonzero((neigh != 0) & (neigh != reg.label)) / collar.sum()
            if shared > max_shared:
                continue

        alpha = feather(mask, 0.8)
        if alpha.max() <= 0:
            continue

        # normalise brightness so sprites from different micrographs, taken at
        # different gains, sit together without one looking blown out
        inside = patch[mask]
        lo, hi = np.percentile(inside, (4, 98))
        if hi - lo < 12:
            continue
        norm = np.clip((patch - lo) / (hi - lo), 0, 1)

        # Dark pockets inside the outline are background showing through a hole
        # in the matte. They read as punched-out holes once composited.
        if np.count_nonzero(norm[mask] < 0.18) / max(mask.sum(), 1) > max_dark:
            continue

        norm = 0.18 + 0.82 * norm

        side = max(r1 - r0, c1 - c0)
        scale = min(1.0, sprite_px / side)
        tw, th = max(4, int((c1 - c0) * scale)), max(4, int((r1 - r0) * scale))

        la = np.dstack([(norm * 255).astype(np.uint8),
                        (alpha * 255).astype(np.uint8)])
        sprite = Image.fromarray(la, mode="LA").resize((tw, th), Image.LANCZOS)

        buf = io.BytesIO()
        sprite.save(buf, format="PNG", optimize=True)
        out.append({
            "uri": "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode(),
            "w": tw, "h": th,
            "area": int(reg.area),
            "aspect": round(float(ar), 2),
            "solidity": round(float(reg.solidity), 3),
            "bytes": buf.getbuffer().nbytes,
            "source": path.name,
        })

    # biggest first: the large ones carry the detail that sells the illusion
    out.sort(key=lambda s: -s["area"])
    return out[:max_particles]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("images", nargs="+", type=Path)
    ap.add_argument("-o", "--out", type=Path, default=Path("grains.js"))
    ap.add_argument("--min-area", type=int, default=400,
                    help="smallest particle to keep, in source pixels")
    ap.add_argument("--per-image", type=int, default=14,
                    help="how many particles to take from each micrograph")
    ap.add_argument("--sprite-px", type=int, default=104,
                    help="longest edge of a stored sprite")
    ap.add_argument("--no-split", action="store_true",
                    help="skip the watershed; use when particles are already separate")
    ap.add_argument("--crop-bottom", type=int, default=0,
                    help="rows to drop from the bottom, for instrument banners")
    ap.add_argument("--max-shared", type=float, default=0.18,
                    help="reject a particle if more than this fraction of its "
                         "outline is a cut against a touching neighbour")
    ap.add_argument("--max-dark", type=float, default=0.10,
                    help="reject a particle with more than this fraction of its "
                         "interior in deep shadow")
    a = ap.parse_args(argv)

    grains = []
    for p in a.images:
        if not p.exists():
            print(f"  {p}: not found", file=sys.stderr)
            continue
        got = extract(p, min_area=a.min_area, max_particles=a.per_image,
                      sprite_px=a.sprite_px, split=not a.no_split,
                      crop_bottom=a.crop_bottom, max_shared=a.max_shared,
                      max_dark=a.max_dark)
        print(f"  {p.name}: {len(got)} particles, "
              f"{sum(g['bytes'] for g in got) / 1024:.0f} kB")
        grains.extend(got)

    if not grains:
        print("No particles survived. Try lowering --min-area.", file=sys.stderr)
        return 1

    payload = [{"uri": g["uri"], "w": g["w"], "h": g["h"]} for g in grains]
    a.out.write_text(
        "/* Crystal sprites cut from powder micrographs by make_sprites.py.\n"
        f"   {len(payload)} particles. Regenerate with your own images rather\n"
        "   than editing by hand. */\n"
        "const GRAIN_SPRITES = " + json.dumps(payload, separators=(",", ":")) + ";\n",
        encoding="utf-8")
    total = sum(g["bytes"] for g in grains) / 1024
    print(f"\n{len(grains)} sprites -> {a.out}  ({a.out.stat().st_size / 1024:.0f} kB "
          f"of JS, {total:.0f} kB of PNG)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
