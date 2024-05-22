"""
Configuration / Settings for subbasins.py

Edit this file carefully before running `delineate()`.

See README for more information about the options and for instructions
on how to download the input data for delineating watersheds in areas
outside of the sample data provided for Iceland.

"""

# Directory containing the merged, basin-scale MERIT-Hydro flow direction rasters (.tif)
# Download from https://mghydro.com/watersheds/rasters
# For all paths, do not include a trailing slash.
MERIT_FDIR_DIR = "data/raster/flowdir_basins"
MERIT_FDIR_DIR = r"C:\Data\GIS\MERITHydro\flow_dir_basins"

# Directory containing the merged, basin-scale MERIT-Hydro flow accumulation rasters (.tif)
# Download from https://mghydro.com/watersheds/rasters
MERIT_ACCUM_DIR = "data/raster/accum_basins"
MERIT_ACCUM_DIR = r"C:\Data\GIS\MERITHydro\accum_basins"

# Set to True if you want the script to write status messages to the console
VERBOSE = True

# Set to True to make a bunch of plots of each watershed.
# (Just for debugging. Slows down the script a lot.)
PLOTS = False
NETWORK_DIAGRAMS = True

# Folder where you have stored the Merit-BASINS catchment shapefiles.
# These files need to be downloaded from: https://www.reachhydro.org/home/params/merit-basins
HIGHRES_CATCHMENTS_DIR = "data/shp/merit_catchments"
HIGHRES_CATCHMENTS_DIR = r"C:\Data\GIS\MERITBasins\catchments\src"

# Location of simplified unit catchment boundaries vector data (shapefiles)
# Download from: https://mghydro.org/watersheds/share/catchments_simplified.zip
LOWRES_CATCHMENTS_DIR = "data/shp/catchments_simplified"

# Folder where you have stored the MERIT-Basins River flowline shapefiles
# Download from: https://www.reachhydro.org/home/params/merit-basins
RIVERS_DIR = "C:/Data/GIS/MERITBasins/rivers"

# Folder where the script will write the output GeoJSON files or shapefiles
OUTPUT_DIR = "output"

# The file extension will determine the types of geodata files the script creates.
#   "gpkg" for GeoPackage (recommended)
#   "geojson" for GeoJSON files
#   "shp" for shapefile
# The list of possibilities depends one what is supported by GeoPandas.
# See: https://geopandas.org/en/stable/docs/user_guide/io.html#writing-spatial-data
OUTPUT_EXT = "geojson"

# Directory to store Python pickle files. It can be slow for Python to
# read shapefiles and create a GeoDataFrame. Once you have done this once, you
# can save time in future runs by storing the GeoDataFrame as a .pkl file.
# Enter a blank string, '' if you do NOT want the script to create .pkl files.
# Please note that these files can be large! (Up to around 1 GB for large basins.)
PICKLE_DIR = 'pkl'

# Threshold for number of upstream pixels that defines a stream
# These values worked will in my testing, but you might try changing if the
# outlet is not getting snapped to a river centerline properly
THRESHOLD_SINGLE = 500
THRESHOLD_MULTIPLE = 5000

SIMPLIFY = True
