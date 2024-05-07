"""
Creates a set of points at the centroid of unit catchments,
and simple, linear links to connect them.
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

id = 27

# First, create the *nodes* as points at the catchments' centroids
sql_query = f"""
SELECT *
FROM merit_basins3_{id}
WHERE sorder > 0;
"""

gdf = gpd.read_postgis(sql_query, conn, geom_col='geom')

# Calculate the centroid
gdf['centroid'] = gdf['geom'].centroid

# Create a new geometry column with Point objects
gdf['geom'] = gdf['centroid'].apply(lambda x: Point(x))

# Drop the centroid helper column
gdf = gdf.drop(columns=['centroid'])

# Update the GeoDataFrame's geometry to be the centroids
gdf = gdf.set_geometry('geom')

output_file = f"nodes2_{id}.gpkg"
gdf.to_file(output_file, driver="GPKG")

# Get the comid: point as a dictionary
gdf.set_index('comid2', inplace=True)
P = gdf['geom'].to_dict()
D = gdf['nextdown'].to_dict()

# Now we are going to load the Rivers file, and replace the actual river geometry
# with a simple 2-point connector between the up_node and the down_node.
sql = f"""
SELECT comid, sorder, uparea, nextdownid as nextdown, geom_simple as geom
FROM merit_rivers2_{id}
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

output_file = f"edges2_{id}.gpkg"
rivers_gdf.to_file(output_file, driver="GPKG")