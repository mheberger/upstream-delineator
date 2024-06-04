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
MERIT_FDIR_DIR = r"C:\Data\GIS\MERITHydro\flow_dir_basins"
# MERIT_FDIR_DIR = "data/raster/flowdir_basins"

# Directory containing the merged, basin-scale MERIT-Hydro flow accumulation rasters (.tif)
# Download from https://mghydro.com/watersheds/rasters
MERIT_ACCUM_DIR = r"C:\Data\GIS\MERITHydro\accum_basins"
# MERIT_ACCUM_DIR = "data/raster/accum_basins"

# Set to True if you want the script to write status messages to the console
VERBOSE = True

# Set to True to make a bunch of plots of each watershed.
# (Mostly for debugging. Slows down the script a lot.)
PLOTS = False

# Set to true to output a network diagram of the river network.
# This is a simplified view of the flow pathways.
# IMPORTANT: For this to work, you need to have GraphViz installed on your computer.
#  (Not just the graphviz Python library, which lets you access its functions.)
#  Download installers here: https://graphviz.org/download/
NETWORK_DIAGRAMS = False

# Folder where you have stored the Merit-BASINS unit catchment shapefiles.
# These files need to be downloaded from: https://www.reachhydro.org/home/params/merit-basins
CATCHMENTS_DIR = r"C:\Data\GIS\MERITBasins\catchments\src"
# CATCHMENTS_DIR = "data/shp/merit_catchments"


# Folder where you have stored the MERIT-Basins River flowline shapefiles
# Download from: https://www.reachhydro.org/home/params/merit-basins
RIVERS_DIR = "C:/Data/GIS/MERITBasins/rivers"

# Folder where the script will write the output geodata
OUTPUT_DIR = "output"

# The output file extension will determine the types of geodata files the script creates.
#   "gpkg" for GeoPackage (recommended)
#   "geojson" for GeoJSON files
#   "shp" for shapefile
# The list of possibilities depends one what is supported by GeoPandas.
# See: https://geopandas.org/en/stable/docs/user_guide/io.html#writing-spatial-data
OUTPUT_EXT = "gpkg"

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

# Simplify the output geodata? This will remove some vertices
# from the watershed boundary and river centerlines and produce smaller files.
SIMPLIFY = True

# If SIMPLIFY is True, set SIMPLIFY_TOLERANCE to a value in decimal degrees.
SIMPLIFY_TOLERANCE = 0.0008

# Do you wish to retun the whole watershed (a single polygon) for EACH individual the outlet point?
WATERSHEDS = False

# Output all of the river polylines (even small ones). Perhaps useful for display and mapping.
OUTPUT_ALL_RIVERS = False

# Watersheds created with Merit-Hydro data tend to have many "donut holes"
# ranging from one or two pixels to much larger. Setting FILL = True, will
# fill in these donut holes, and generally result in a better appearance
# and smaller output files.
FILL = True

# If FILL = True, you many choose to to fill only those donut holes that are smaller than
# a certain size. In other words, it will keep the big holes and fill in the small
# ones. This is roughly the number of pixels, on the 3 arcsecond grid.
# Set to 0 to fill ALL holes. Setting FILL_THRESHOLD = 100 will fill all the holes
# that are less than 100 pixels in size. (Little ones are usually minor topological
# errors in the input data, while larger holes *may* be more meaningful, reflecting
# surface drainage patterns.)
FILL_THRESHOLD = 100

# Consolidate the sub-basins to make a larger size? If set to true, the script will
# merge adjacent subbasins such that
# (a) subbasins do not exceed a max. threshold area and
# (b) the network topology is maintained (overall connectivity of the flow network)

CONSOLIDATE = True
MAX_AREA = 1500  # in km²

# MERGE tiny junction nodes? Only activated when CONSOLIDATE == True.
# TODO: NOT CURRENTLY IMPLEMENTED.
# These sometimes occur around confluences, where two tributaries
# join the mainstem close to one another. If you want to consider them joining at the same
# location, set the length below, in km. If you don't want to do this, set it to a big number like 99999
THRESHOLD_LENGTH = 2
