import numpy as np
import matplotlib.pyplot as plt

def annotate_plot(data, wcs, det_tbl, matches, out_png, subtitle=None):
    img = data.astype(float).copy()
    p1, p99 = np.nanpercentile(img, 1), np.nanpercentile(img, 99)
    img = (img - p1) / (p99 - p1 + 1e-9)
    img = np.clip(img, 0, 1)

    fig = plt.figure(figsize=(9, 7), dpi=150)
    ax = plt.subplot(projection=wcs)
    ax.imshow(img, origin='lower', cmap='gray')
    title = 'Detections and Catalog Matches'
    if subtitle: title += f"\n{subtitle}"
    ax.set_title(title); ax.set_xlabel('RA'); ax.set_ylabel('Dec')

    if len(det_tbl) > 0:
        ax.scatter(det_tbl['xcentroid'], det_tbl['ycentroid'], transform=ax.get_transform('pixel'),
                   s=20, marker='o', facecolors='none', edgecolors='C0', linewidths=1.2, label='Detections')

    if len(matches) > 0 and 'matched' in matches.colnames:
        matched = np.array(matches['matched'], dtype=bool)
        if matched.any():
            ax.scatter(det_tbl['xcentroid'][matched], det_tbl['ycentroid'][matched],
                       transform=ax.get_transform('pixel'),
                       s=30, marker='x', color='C1', linewidths=1.2, label='Matched')
    ax.legend(loc='upper right', frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, bbox_inches='tight')
    plt.close(fig)
