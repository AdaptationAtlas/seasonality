#!/bin/bash

################################################################################
# GLASS NDVI HDF4 → GeoTIFF CONVERSION SCRIPT (NO SCALING)
#
# This script converts all .hdf files in the specified INPUT_DIR to Float32
# NDVI GeoTIFF files in OUTPUT_DIR using gdal_translate. Assumes the data is
# already scaled to NDVI float values (~ -0.2 to 1.0).
#
# 🛠 REQUIREMENTS:
# - GDAL with HDF4 support (including libgdal-hdf4 plugin)
# - Conda environment with gdal, hdf4, and libgdal-hdf4 installed
#
# ✅ SETUP INSTRUCTIONS:
# 1. Install Miniforge (lightweight Conda):
#    Apple Silicon:
#      curl -L -o Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
#    Intel Macs:
#      curl -L -o Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh
#    Then:
#      bash Miniforge3.sh
#      conda init zsh && exec zsh
#
# 2. Create HDF4-enabled GDAL environment:
#    conda create -n hdf_env -c conda-forge gdal hdf4 libgdal-hdf4 python
#    conda activate hdf_env
#
# 3. Confirm support:
#    gdalinfo --formats | grep -i hdf
#    # Should list: HDF4, HDF4Image
#
# 4. Run this script inside the environment:
#    chmod +x convert_glass_ndvi_direct.sh
#    ./convert_glass_ndvi_direct.sh
#
# 📦 OUTPUT:
# - GeoTIFFs preserving original NDVI float values, Float32 format
# - Fails logged to gdal_errors.log
################################################################################

# === USER CONFIGURATION ===
INPUT_DIR="/Volumes/clim_dat/climate_raw/glass_ndvi"
OUTPUT_DIR="/Volumes/clim_dat/climate_raw/glass_ndvi_tif"
LOG_FILE="gdal_errors.log"

# === SETUP ===
mkdir -p "$OUTPUT_DIR"
rm -f "$LOG_FILE"
touch "$LOG_FILE"

echo "🚀 Starting GLASS NDVI batch conversion at $(date)" | tee -a "$LOG_FILE"

count_total=0
count_success=0
count_failed=0

# === MAIN LOOP ===
for f in "$INPUT_DIR"/*.hdf; do
((count_total++))
base=$(basename "$f" .hdf)
out="$OUTPUT_DIR/${base}.tif"

echo "🔄 Processing: $base" | tee -a "$LOG_FILE"

if [ -f "$out" ]; then
echo "  ✅ Skipping (already exists)" | tee -a "$LOG_FILE"
continue
fi

# Convert using GDAL; preserve original NDVI float values
gdal_translate "$f" "$out" \
-a_nodata -9999 -ot Float32 2>>"$LOG_FILE"

if [ $? -eq 0 ]; then
echo "  ✅ Success" | tee -a "$LOG_FILE"
((count_success++))
else
  echo "  ❌ ERROR processing $base" | tee -a "$LOG_FILE"
((count_failed++))
fi
done

# === SUMMARY ===
echo "✅ Batch conversion complete at $(date)" | tee -a "$LOG_FILE"
echo "Total files: $count_total" | tee -a "$LOG_FILE"
echo "  Successful: $count_success" | tee -a "$LOG_FILE"
echo "  Failed:     $count_failed" | tee -a "$LOG_FILE"
