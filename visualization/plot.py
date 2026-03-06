import geopandas as gpd
import matplotlib.pyplot as plt
import rasterio
import rasterio.warp
import numpy as np
from matplotlib.colors import LogNorm, PowerNorm
from matplotlib.lines import Line2D
from pathlib import Path

# ==========================================================
# Scenario & CRS
# ==========================================================
scenario_name = "moderate_action"
PLOT_CRS = "EPSG:4326"

# ==========================================================
# Paths
# ==========================================================
load_zones_path = Path("visualization/load_zones/load_zones.shp")
demand_tif_path = Path(f"user_inputs/{scenario_name}/demand.tif")
sited_h2_path = Path(f"outputs/{scenario_name}/sited_h2_plants.gpkg")

# ==========================================================
# Technology Reference Capacities (MW)
# ==========================================================
ref_capacity = {
    "gas_smr": 150,
    "gas_smr_ccs": 150,
    "gas_atr_ccs": 205,
    "bio_smr": 150,
    "bio_smr_ccs": 150,
    "bio_atr_ccs": 205,
    "coal_gas": 205,
    "coal_gas_ccs": 205,
    "biomass_gas": 48,
}

# ==========================================================
# Load Vector Data
# ==========================================================
load_zones = gpd.read_file(load_zones_path).to_crs(PLOT_CRS)
plants = gpd.read_file(sited_h2_path).to_crs(PLOT_CRS)

# Use centroids for plotting
plants_pts = plants.copy()
plants_pts["geometry"] = plants.geometry.centroid

# Plot bounds with padding
minx, miny, maxx, maxy = load_zones.total_bounds
pad = 0.05  # 2% padding
dx = maxx - minx
dy = maxy - miny
minx -= pad * dx
maxx += pad * dx
miny -= pad * dy
maxy += pad * dy

# ==========================================================
# PLOT 1 — Hydrogen Demand (LOG SCALE, ZERO INCLUDED)
# ==========================================================
fig, ax = plt.subplots(figsize=(5, 12), dpi=150)

with rasterio.open(demand_tif_path) as src:
    # Reproject raster to EPSG:4326
    dst_transform, width, height = rasterio.warp.calculate_default_transform(
        src.crs, PLOT_CRS, src.width, src.height, *src.bounds
    )

    demand = np.zeros((height, width), dtype=np.float32)

    rasterio.warp.reproject(
        source=rasterio.band(src, 1),
        destination=demand,
        src_transform=src.transform,
        src_crs=src.crs,
        dst_transform=dst_transform,
        dst_crs=PLOT_CRS,
        resampling=rasterio.warp.Resampling.nearest,
    )

    # Replace nodata with zero
    demand = np.nan_to_num(demand, nan=0.0)

    # Compute raster extent
    extent = (
        dst_transform.c,
        dst_transform.c + dst_transform.a * width,
        dst_transform.f + dst_transform.e * height,
        dst_transform.f,
    )




from matplotlib.colors import SymLogNorm

# Small threshold for linear region around zero
linthresh = 1.0  # values smaller than this are mapped linearly

norm = SymLogNorm(linthresh=linthresh,
                  vmin=demand.min(),
                  vmax=1e6,
                  base=10)

img = ax.imshow(
    demand,
    extent=extent,
    cmap=plt.cm.viridis,
    norm=norm,
    interpolation='nearest',
    zorder=1
)

# Overlay load zones
load_zones.plot(
    ax=ax,
    facecolor="none",
    edgecolor="black",
    linewidth=0.5,
    zorder=2,
)

# Colorbar
cbar = fig.colorbar(img, ax=ax, fraction=0.035, pad=0.02)
cbar.set_label("Hydrogen Demand (kg/year)", fontsize=10)

# Axes formatting
ax.set_xlim(minx, maxx)
ax.set_ylim(miny, maxy)
ax.set_aspect('equal')
ax.set_xlabel("Longitude", fontsize=11)
ax.set_ylabel("Latitude", fontsize=11)
ax.set_title(f"Hydrogen Demand - {scenario_name}", fontsize=10, fontweight='bold')

plt.tight_layout()
plt.show()

# ==========================================================
# PLOT 2 — Hydrogen Production Plants (Centroids, Scaled by Capacity)
# ==========================================================
fig, ax = plt.subplots(figsize=(5, 12), dpi=150)

# Background
load_zones.plot(
    ax=ax,
    facecolor="#f0f0f0",
    edgecolor="#444444",
    linewidth=0.7,
    zorder=1,
)

ref_capacity = {
    "gas_smr": 150,
    "gas_smr_ccs": 150,
    "gas_atr_ccs": 205,
    "bio_smr": 150,
    "bio_smr_ccs": 150,
    "bio_atr_ccs": 205,
    "coal_gas": 205,
    "coal_gas_ccs": 205,
    "biomass_gas": 48,
}

# Color map for technologies
safe_colors = [
    "#E60000",  # red
    "#86B6F0",  # light blue
    "#3760EA",  # blue
    "#D3A97C",  # light bown
    "#957157",  # darker brown
    "#563610",  # dark brown
    "#313131",  # gray
    "#F0E442",  # yellow
    "#BABABA"   # different shade of gray
]

tech_colors = {
    tech: safe_colors[i % len(safe_colors)]
    for i, tech in enumerate(ref_capacity.keys())
}

# Scaling factor for marker size
size_scale = .2

legend_handles = []
capacity_handles = []

for tech in ref_capacity.keys():
    subset = plants_pts[plants_pts["prod_tech"] == tech]
    if subset.empty:
        continue

    # Marker sizes proportional to capacity
    sizes = subset["total_capacity_MW"].values * size_scale

    ax.scatter(
        subset.geometry.x,
        subset.geometry.y,
        s=sizes,
        c=[tech_colors[tech]],
        edgecolors="black",
        linewidths=.7,
        alpha=1,
        zorder=3,
    )

    # Tech legend handle
    legend_handles.append(
        Line2D(
            [0], [0],
            marker='o',
            color='w',
            markerfacecolor=tech_colors[tech],
            markeredgecolor='black',
            markersize=10,
            label=tech,
            linestyle='None'
        )
    )

# Capacity legend examples (based on marker area, not radius)
example_caps = [50, 150, 250]  # t/day
for c in example_caps:
    marker_size = np.sqrt(c * size_scale)  # Convert area to radius for markersize
    capacity_handles.append(
        Line2D(
            [0], [0],
            marker='o',
            color='w',
            markerfacecolor='gray',
            markeredgecolor='black',
            markersize=marker_size,
            label=f"{c} MW",
            linestyle='None'
        )
    )
# Tech legend
legend1 = ax.legend(
    handles=legend_handles,
    title="Production Technology",
    title_fontsize='small',      # Shrinks the title
    prop={'size': 8},            # Shrinks the labels
    loc='upper center',
    bbox_to_anchor=(0.5, -0.17),  # Adjusted closer to plot
    frameon=True,
    ncol=3,
    labelspacing=0.5,            # Vertical space between entries
    handletextpad=0.4            # Space between icon and text
)

# Capacity legend
legend2 = ax.legend(
    handles=capacity_handles,
    title="Plant Capacity",
    title_fontsize='small',
    prop={'size': 8},
    loc='upper center',
    bbox_to_anchor=(0.5, -0.5), # Moved up since tech legend is now smaller
    frameon=True,
    ncol=3,
    labelspacing=0.5
)

# Add the first legend back
ax.add_artist(legend1)

# Increase bottom margin so legends fit
plt.subplots_adjust(bottom=0.5)

# Axes and title
ax.set_xlim(minx, maxx)
ax.set_ylim(miny, maxy)
ax.set_aspect('equal')
ax.set_xlabel("Longitude", fontsize=10)
ax.set_ylabel("Latitude", fontsize=10)
ax.set_title(f"Hydrogen Plant Build-Out - {scenario_name}", fontsize=10, fontweight='bold')
plt.show()