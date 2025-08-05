import os
os.chdir('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
from pyhere import here
import pandas as pd
import numpy as np
import geopandas as gpd
import json
import os
import matplotlib.pyplot as plt
import tifffile
import session_info

sample_path = here('raw-data/sample_info/VistoSeg_samples_n24.csv')
tissue_positions_path = here('processed-data', '01_spaceranger', '{}', 'outs', 'spatial', 'tissue_positions.csv')
scale_factors_path = here('processed-data', '01_spaceranger', '{}', 'outs','spatial', 'scalefactors_json.json')
nuclear_path = here('processed-data', 'spaceranger_segment', '{}','outs', 'nucleus_segmentations.geojson')
out_path = here('processed-data', 'spaceranger_segment', '{}.csv')
plot_dir = here('plots', 'spaceranger_segment')

os.makedirs(os.path.join(plot_dir, 'segmentation_check'), exist_ok = True)
os.makedirs(os.path.join(plot_dir, 'alignment_check'), exist_ok = True)
os.makedirs(os.path.join(plot_dir, '5_nucleus_spot'), exist_ok = True)
os.makedirs(out_path.parent, exist_ok = True)

################################################################################
#   Functions
################################################################################

def overlay_plot(img, spot_gdf, nuc_gdf, start, box_len, plot_path):
    #   Subset image, spots, and nuclei to a small region
    img_small = np.flipud(
        img[start[1]:(start[1] + box_len), start[0]:(start[0] + box_len), :]
    )
    spot_gdf_small = spot_gdf.cx[
        start[0]:(start[0] + box_len), start[1]:(start[1] + box_len)
    ]
    nuc_gdf_small = nuc_gdf.cx[
        start[0]:(start[0] + box_len), start[1]:(start[1] + box_len)
    ]
    
    #   Overlay H&E, spots, and segmentations for a small region
    fig, ax = plt.subplots()
    spot_gdf_small['geometry'].plot(ax = ax, facecolor = 'none')
    nuc_gdf_small['geometry'].plot(
        ax = ax, facecolor = 'none', edgecolor = 'blue'
    )
    plt.imshow(
        img_small,
        extent = (start[0], start[0] + box_len, start[1], start[1] + box_len)
    )
    plt.gca().invert_yaxis()
    plt.savefig(plot_path)
    plt.close('all')

################################################################################
#   Main
################################################################################

sample_df = pd.read_csv(sample_path, header = None)

#   Get the sample ID and image path for this task
task_id = int(os.getenv('SLURM_ARRAY_TASK_ID'))
sample_id = sample_df.iloc[task_id, 0]
img_path=here(f'raw-data/images/{sample_id}.tif')

#   Fill in paths for this particular sample
tissue_positions_path = str(tissue_positions_path).format(sample_id)
nuclear_path = str(nuclear_path).format(sample_id)
scale_factors_path = str(scale_factors_path).format(sample_id)
out_path = str(out_path).format(sample_id)

#   Read in Visium spots as a GeoDataFrame consisting of Points
spot_df = pd.read_csv(tissue_positions_path)
spot_gdf = gpd.GeoDataFrame(
    spot_df,
    geometry = gpd.points_from_xy(
        spot_df['pxl_col_in_fullres'], spot_df['pxl_row_in_fullres']
    )
)

#   Make the Visium spots have the correct radius (give them nonzero area)
with open(scale_factors_path, 'r') as f:
    scale_factors = json.load(f)

spot_gdf['geometry'] = spot_gdf.buffer(
    scale_factors['spot_diameter_fullres'] / 2
)

#   For some reason, the default CRS is 'EPSG:4326', which seems inappropriate
#   for flat 2D data like a Visium capture area / image. Force the CRS to be 
#   None and read in the nuclear segmentations
nuc_gdf = gpd.read_file(nuclear_path).set_crs(None, allow_override = True)

#   Plot spots and nuclear segmentations spatially (full capture area) to
#   ensure coordinate systems are properly aligned
fig, ax = plt.subplots()
spot_gdf['geometry'].plot(ax = ax, facecolor = 'none')
nuc_gdf['geometry'].plot(ax = ax, facecolor = 'blue', edgecolor = 'blue')
plt.gca().invert_yaxis()
plt.savefig(os.path.join(plot_dir, 'alignment_check', f'{sample_id}.png'))
plt.close('all')

#   Define a subregion for plotting, centered at the center of all Visium spots
spot_bounds = spot_gdf.total_bounds
start = (
    int((spot_bounds[2] + spot_bounds[0] + 1000) / 2),
    int((spot_bounds[3] + spot_bounds[1] + 1000) / 2)
)

img = tifffile.imread(img_path)

overlay_plot(
    img = img,
    spot_gdf = spot_gdf,
    nuc_gdf = nuc_gdf,
    start = start,
    box_len = 1000,
    plot_path = os.path.join(plot_dir, 'segmentation_check', f'{sample_id}.png')
)

# Count nuclei whose centroids fall within each Visium spot
nuc_centroids = nuc_gdf.copy()
nuc_centroids['geometry'] = nuc_centroids.geometry.centroid

# enitre nuclei within
within_join = gpd.sjoin(nuc_gdf, spot_gdf, how='left', predicate='within')
within_counts = within_join.groupby('index_right').size().reindex(spot_gdf.index, fill_value=0)

# only nuclei centroid within
centroid_join = gpd.sjoin(nuc_centroids, spot_gdf, how='left', predicate='within')
centroid_counts = centroid_join.groupby('index_right').size().reindex(spot_gdf.index, fill_value=0)

# only on spot boundary
intersect_join = gpd.sjoin(nuc_gdf, spot_gdf, how='left', predicate='intersects')
intersect_counts = intersect_join.groupby('index_right').size().reindex(spot_gdf.index, fill_value=0)
# Assign
spot_gdf['num_nuclei_within'] = within_counts
spot_gdf['num_nuclei_centroid'] = centroid_counts
spot_gdf['num_nuclei_intersect'] = intersect_counts

#### proportion spot covered #####
# Step 0: Reset indexes so we can track source IDs
nuc_gdf = nuc_gdf.reset_index().rename(columns={"index": "nuc_index"})
spot_gdf = spot_gdf.reset_index().rename(columns={"index": "spot_index"})

# Step 1: Overlay to compute overlapping geometry
intersections = gpd.overlay(nuc_gdf, spot_gdf, how='intersection')

# Step 2: Area of nuclear segmentation within the spot
intersections['nuc_area_in_spot'] = intersections.geometry.area

# Step 3: Sum area per spot
area_per_spot = intersections.groupby('spot_index')['nuc_area_in_spot'].sum()

# Step 4: Assign to original spot_gdf
spot_gdf['nuclei_pixel_area'] = spot_gdf['spot_index'].map(area_per_spot).fillna(0)

# Step 5: Compute total spot area and coverage
spot_gdf['spot_area'] = spot_gdf.geometry.area
spot_gdf['nuclei_coverage'] = spot_gdf['nuclei_pixel_area'] / spot_gdf['spot_area']
#   To verify counting of nuclei per spot worked, examine the first spot with
#   5 reported nuclei
overlay_plot(
    img = img,
    spot_gdf = spot_gdf,
    nuc_gdf = nuc_gdf,
    start = (
        spot_gdf
            .loc[spot_gdf['num_nuclei_centroid'] == 5, :]
            .iloc[[0], :]
            .total_bounds[:2]
            .astype(int)
    ),
    box_len = int(np.ceil(scale_factors['spot_diameter_fullres'])),
    plot_path = os.path.join(plot_dir, '5_nucleus_spot', f'{sample_id}.png')
)

#   Clean up and export nuclei counts
spot_gdf['key'] = spot_gdf['barcode'] + '_' + sample_id
spot_gdf = spot_gdf.loc[
    : , ['key', 'num_nuclei_within', 'num_nuclei_intersect', 'num_nuclei_centroid', 'nuclei_coverage']
]
spot_gdf.to_csv(out_path, index = False)

session_info.show()
