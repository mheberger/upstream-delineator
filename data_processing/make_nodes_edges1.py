"""
NODES & EDGES
Creates a set of points at the centroid of unit catchments,
and simple, linear links to connect them
"""

import warnings
warnings.filterwarnings('ignore')
import psycopg2
import geopandas as gpd
from shapely.geometry import Point, LineString


# Connect to your PostgreSQL database
conn = psycopg2.connect(
    host="localhost",
    database="basins",
    user="postgres",
    password="dbpw"
)

id = 77

# What this SQL query is doing: MERGES the unit catchments with sorder = 1
# when they drain to the same downstream
# For coastal catchments, where there is no upstream neighbor and which drain to the ocean
# (sorder = 1 AND nextdown = 0), we will exclude these.

sql_query = f"""
SELECT * 
FROM merit_basins_{id}
WHERE geom_simple IS NOT Null;
"""

gdf = gpd.read_postgis(sql_query, conn, geom_col='geom_simple')

# Turn the polygon into a point by calculating the centroid
gdf['centroid'] = gdf['geom_simple'].centroid

# Create a new geometry column with Point objects

gdf['geom'] = gdf['centroid'].apply(lambda x: Point(x))

# Drop the original geometry column and the helper column
gdf = gdf.drop(columns=['geom_simple'])
gdf = gdf.drop(columns=['centroid'])

# Update the GeoDataFrame's geometry to be the centroids
gdf = gdf.set_geometry('geom')

output_file = f"nodes{id}.gpkg"
gdf.to_file(output_file, driver="GPKG")

# Create a Dictionary relating comid: Point. We'll use this to construct edges.
gdf.set_index('comid', inplace=True)
P = gdf['geom'].to_dict()
D = gdf['nextdown'].to_dict()

# Now we are going to load the Rivers file, and replace the actual river geometry
# with a simple 2-point connector between the up_node and the down_node.
sql = f"""
SELECT comid,  sorder, nextdownid as nextdown, maxup as numup, uparea, up1, up2, up3, up4, geom_simple as geom
FROM merit_rivers_{id}
WHERE nextdownid > 0;
"""

rivers_gdf = gpd.read_postgis(sql, conn, geom_col='geom')
rivers_gdf['connector'] = None
rivers_gdf.set_index('comid', inplace=True)

# Iterate over the rows in the rivers GeoDataFrame to construct the connectors
for comid in rivers_gdf.index:
    from_node = comid
    to_node = D[comid]
    from_point = P[from_node]
    to_point = P[to_node]
    connector = LineString([from_point, to_point])
    rivers_gdf.loc[comid, 'geom'] = connector

output_file = f"edges{id}.gpkg"
rivers_gdf.to_file(output_file, driver="GPKG")