"""
Configuration / Settings for subbasins.py

Edit this file carefully before running `delineate()`.

See README for more information about the options and for instructions
on how to download the input data for delineating watersheds in areas
outside of the sample data provided for Iceland.

"""

# Set to True if you want the script to write status messages to the console
VERBOSE = True

# Directory (folder) containing the merged, basin-scale MERIT-Hydro flow direction rasters (.tif)
# Download from https://mghydro.com/watersheds/rasters
# For all paths, do not include a trailing slash.
#MERIT_FDIR_DIR = r"C:\Data\GIS\MERITHydro\flow_dir_basins"
MERIT_FDIR_DIR = "data/raster/flowdir_basins"

# Folder where you have stored the Merit-BASINS unit catchment shapefiles.
# These files need to be downloaded from: https://www.reachhydro.org/home/params/merit-basins
#CATCHMENTS_DIR = r"C:\Data\GIS\MERITBasins\catchments\src"
CATCHMENTS_DIR = "data/shp/merit_catchments"


# Folder where you have stored the MERIT-Basins River flowline shapefiles
# Download from: https://www.reachhydro.org/home/params/merit-basins
#RIVERS_DIR = "C:/Data/GIS/MERITBasins/rivers"
RIVERS_DIR = "data/shp/merit_rivers"

# Folder where the script will write the output geodata
OUTPUT_DIR = "output"

# The output file extension will determine the types of geodata files the script creates.
#   "gpkg" for GeoPackage (recommended)
#   "geojson" for GeoJSON files
#   "shp" for shapefile
# The list of possibilities depends on what is currently supported by GeoPandas.
# See: https://geopandas.org/en/stable/docs/user_guide/io.html#writing-spatial-data
OUTPUT_EXT = "gpkg"

# Set to True to make a bunch of plots of each watershed,  focused on the raster-based delineation.
# (Mostly for debugging. Slows down the script a lot.)
PLOTS = False

# Directory to put plots created by the script.
PLOTS_DIR = 'plots'

# Directory to store Python pickle files. It can be slow for Python to
# read shapefiles and create a GeoDataFrame. Once you have done this once, you
# can save time in future runs by storing the GeoDataFrame as a .pkl file.
# Enter a blank string, '' if you do NOT want the script to create .pkl files.
# Please note that these files can be large! (Up to around 1 GB for large basins.)
PICKLE_DIR = 'pkl'

# Threshold for number of upstream pixels that defines a stream
# These values worked will in my testing, but you might try changing if the
# outlet is not getting snapped to a river centerline properly
THRESHOLD_SINGLE = 500     # Where the outlet is in a unit catchment without upstream neighbors; finds smaller streams.
THRESHOLD_MULTIPLE = 5000  # Where the outlet is in a unit catchment with upstream neighbors; finds larger rivers.

# Simplify the output geodata? This will remove some vertices
# from the watershed boundary polygons and river reach centerlines and produce smaller files  
# However, it may also create topology problems -- misaligned edges, slivers, and dangles.
# But the appearance may be jagged when zoomed in. Better results may be obtained with GIS or mapshaper.
SIMPLIFY = False

# If SIMPLIFY is True, set SIMPLIFY_TOLERANCE to a value in decimal degrees.
# This is equivalent to the parameter epsilon in the Ramer–Douglas–Peucker algorithm
SIMPLIFY_TOLERANCE = 0.0008

# Do you wish to retun the whole watershed (a single polygon) for EACH individual outlet point?
# If True the script will create files with names beginning with watershed_ in the output directory
WATERSHEDS = True

# Output a separate geodata file with ALL of the river polylines (even small ones)? 
# This may be useful for display and mapping.
OUTPUT_ALL_RIVERS = True

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

# Consolidate the sub-basins to make them larger? If set to True, the script will
# selectively merge adjacent subbasins such that:
# (a) subbasins do not exceed a maximum area in MAX_AREA
# (b) the network topology is maintained (overall connectivity of the flow network)
# If you set MAX_AREA to a very high number, the network will be collapsed to the maximum
# extent possible while maintaining subbasins for your outlets and the necessary junctions.

CONSOLIDATE = False
MAX_AREA = 750  # in km²

# Output a network diagram of the river network?
# This is a simplified view of the flow pathways.
# IMPORTANT: For this to work, you need to have GraphViz installed on your computer.
#  (Not just the graphviz Python library, which lets you access its functions.)
#  Download installers here: https://graphviz.org/download/
NETWORK_DIAGRAMS = False

# See a list of available formats here: https://graphviz.org/docs/outputs/
DIAGRAM_FORMAT = 'pdf'

# Show the area of unit catchments on the network diagram (in the label, and size of bubble)?
# If you are plotting a large river basins, the diagram can be very large, so better to choose
SHOW_AREA = False

# Make the river network diagram vertical (top to bottom). If false, plot will be left to right.
VERTICAL_PLOT = True

# Do you want to export the river network graph?
SAVE_NETWORK = True

# What kind of file do you want for the network?
# 'pkl':  Python pickle file
# 'json': JSON file
# 'gml':  GML (Graph Modeling Language), a common graph file format.
# 'xml':  GraphML is an XML-based file format for graphs.
NETWORK_FILE_EXT = 'xml'
