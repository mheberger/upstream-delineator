"""
NODES & EDGES #2
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
FROM merit_basins2_{id}
"""

gdf = gpd.read_postgis(sql_query, conn, geom_col='geom')

# Turn the polygon into a point by calculating the centroid
gdf['centroid'] = gdf['geom'].centroid

# Create a new geometry column with Point objects
gdf['geom'] = gdf['centroid'].apply(lambda x: Point(x))

# Drop the helper column
gdf = gdf.drop(columns=['centroid'])

output_file = f"nodes2_{id}.gpkg"
gdf.to_file(output_file, driver="GPKG")

# Create a Dictionary relating comid: Point. We'll use this to construct edges.
gdf.set_index('comid2', inplace=True)
P = gdf['geom'].to_dict()
D = gdf['nextdown'].to_dict()

# Remove any rows from the table that do not have a downstream
edges_gdf = gdf[gdf['nextdown'] !=0]

edges_gdf['connector'] = None

# Iterate over the rows in the rivers GeoDataFrame to construct the connectors
for comid in edges_gdf.index:
    from_node = comid
    to_node = D[comid]
    from_point = P[from_node]
    to_point = P[to_node]
    connector = LineString([from_point, to_point])
    edges_gdf.loc[comid, 'geom'] = connector

output_file = f"edges2_{id}.gpkg"
edges_gdf.to_file(output_file, driver="GPKG")