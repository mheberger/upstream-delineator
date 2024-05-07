"""
Final step in fixing the simplified geometries. I found out my previous method in QGIS was removing some small
features. So I redid the simplification using mapshaper. Now I am taking the geometry from the simplified
shapefiles and using it to update the database tables.

This demonstrates a really simple workflow for inserting geodata into a PostgreSQL database using GeoPandas. 

1. Read data using GeoPandas into a GeoDataFrame.
2. Use GeoPandas' `to_postgis()` method to write a table to the DB. 

So much simpler than all the other methods I had previously tried!!!

"""

import warnings
warnings.filterwarnings('ignore')
import psycopg2
import geopandas
from sqlalchemy import create_engine

basins = [11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44,
          45, 46, 47, 48, 51, 52, 53, 54, 55, 56, 57, 61, 62, 63, 64, 65, 66, 67, 71, 72, 73, 74, 75, 76, 77, 78,
          81, 82, 83, 84, 85, 86, 91]

basin = 77
print(basin)
# Step #1: Read the shapefile into a GeoDataFrame
shp = fr'C:\Data\GIS\MERITBasins\catchments\simplified\cat_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01.shp'
print("Reading shapefile")
gdf = geopandas.read_file(shp)
gdf.rename(columns={"COMID": "comid"}, inplace=True)

# Step #2: Write this file to my PostgreSQL database as a temporary table.
# Create a SQLAlchemy engine
host = "localhost"
database = "basins"
user = "postgres"
password = "dbpw"
port = 5432
engine = create_engine('postgresql://{}:{}@{}:{}/{}'.format(user, password, host, port, database))

print("Writing temp. table to Postgres")
gdf.to_postgis('a_temp', engine, if_exists='replace', index=False)


# Step #3: Run a query to (a) join the temp table to my catchments table and update the geometry.

conn = psycopg2.connect(
    host=host,
    database=database,
    user=user,
    password=password
)

cur = conn.cursor()

sql = """
ALTER TABLE a_temp
ALTER COLUMN comid TYPE INTEGER;
"""
cur.execute(sql)

sql = f"""
UPDATE merit_basins_{basin} as M                 
SET geom_simple = a_temp.geometry
FROM a_temp 
WHERE M.comid = a_temp.COMID; 
"""

print("Updating ")
cur.execute(sql)
conn.commit()

cur.close()
conn.close()
