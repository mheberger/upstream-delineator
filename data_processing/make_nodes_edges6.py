"""
Creates a set of points at the centroid of unit catchments,
and simple, linear links to connect them.
"""

import warnings
warnings.filterwarnings('ignore')
import psycopg2
import geopandas as gpd
from shapely.geometry import Point, LineString

# MEGABASIN
id = 77

# Connect to your PostgreSQL database
conn = psycopg2.connect(
    host="localhost",
    database="basins",
    user="postgres",
    password="dbpw"
)

# First, create the *nodes* as points at the catchments' centroids
sql_query = f"""
SELECT *
FROM merit_basins6_{id};
"""

gdf = gpd.read_postgis(sql_query, conn, geom_col='geom')

# Calculate the centroid
gdf['centroid'] = gdf['geom'].centroid

# For debugging, find rows where centroid is None (?)
null_centroids = gdf[gdf['centroid'].isnull()]

# Create a new geometry column with Point objects
gdf['geom'] = gdf['centroid'].apply(lambda x: Point(x))

# Drop the centroid helper column
gdf = gdf.drop(columns=['centroid'])

output_file = f"nodes6_{id}.gpkg"
gdf.to_file(output_file, driver="GPKG")

# Creat a set of dictionaries relating comid to point geometries, and to the to_node (next downstream reach)
gdf.set_index('comid3', inplace=True)

P = gdf['geom'].to_dict()
D = gdf['nextdown'].to_dict()

# We can simply reuse the existing GeoDataFrame.
# But we'll filter out the terminal nodes (which drain to ocean)
gdf = gdf[gdf['nextdown'] > 0]

# Iterate over the rows in the GeoDataFrame to construct the connectors
for comid in gdf.index:
    from_node = comid
    to_node = D[comid]
    try:
        from_point = P[from_node]
        to_point = P[to_node]
    except:
        print(f"Something went wrong! Check row with comid = {comid}")
    connector = LineString([from_point, to_point])
    gdf.loc[comid, 'geom'] = connector

output_file = f"edges6_{id}.gpkg"
gdf.to_file(output_file, driver="GPKG")