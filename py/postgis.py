"""
Experiment on Nov. 23, 2023
Optionally, if the user has installed the PostgreSQL database system, plus the PostGIS
extension, the script can use that for much faster dissolve. 

I thought that this would be much faster than using GeoPandas to do the dissolve. 

But it turns out that my fast_dissolve() routine is actually faster (!).

So, abandoned. Was an interesting idea, but not worthwhile.
"""

import psycopg2
import os

from py.util import get_largest
os.environ['USE_PYGEOS'] = '0'
import geopandas
import matplotlib.pyplot as plt
from shapely import wkb
from shapely.geometry import Polygon

DBNAME = "basins"
DBUSER = "postgres"
DBPW = "dbpw"
HOST = "localhost"
PORT = "5432"


def postgis_installed() -> bool:
    """
    Checks whether the local machine has
    expects there to be global variables for the DBNAME, DBUSER, and DBPW
    :return: True or False
    """
    try:
        # Connect to PostgreSQL
        conn = psycopg2.connect(f"dbname={DBNAME} user={DBUSER} password={DBPW}")
        cursor = conn.cursor()

        # Check PostgreSQL version
        cursor.execute("SELECT version();")
        postgres_version = cursor.fetchone()
        print("PostgreSQL version:", postgres_version[0])

        # Check for PostGIS extension
        cursor.execute("SELECT postgis_full_version();")
        postgis_extension = cursor.fetchone()
        print(postgis_extension)
        if postgis_extension:
            print("PostGIS is installed.")
            bPostGIS = True
        else:
            print("PostGIS is not installed.")
            bPostGIS = False

        # Close connections
        cursor.close()
        conn.close()
        return bPostGIS

    except psycopg2.OperationalError as e:
        print("Error connecting to PostgreSQL:", e)
        return False


def dissolve_postgis(gdf: geopandas.GeoDataFrame) -> geopandas.GeoDataFrame:
    # Establish a connection to the PostgreSQL database
    # engine = create_engine(f"postgresql://{DBUSER}:{DBPW}@{HOST}:{PORT}/{DBNAME}")

    # Convert the GeoDataFrame's geometry to WKT format
    #gdf['geom_wkb'] = gdf['geometry'].apply(lambda x: x.wkb_hex)
    #gdf = gdf.assign(geom_wkb=gdf['geometry'].apply(lambda x: x.wkb_hex))

    #wkb_hex_list = gdf['geometry'].apply(lambda geom: geom.wkb_hex).tolist()
    geometries = gdf.geometry

    polygons = [geom.wkb_hex for geom in gdf['geometry'] if isinstance(geom, Polygon)]

    # Dissolve the geometries using a CTE in PostgreSQL
    conn = psycopg2.connect(
        dbname=DBNAME,
        user=DBUSER,
        password=DBPW,
        host=HOST,
        port=PORT
    )

    cursor = conn.cursor()


    # Prepare data for insertion into the query using psycopg2's execute_values
    query = "SELECT ST_Buffer(ST_Buffer(ST_Union(ARRAY[\n"
    for polygon in polygons:
        query += f"ST_SetSRID('{polygon}'::geometry, 4326),\n"
    query = query[:-2] + "]\n), 0.0004, 'join=mitre'), -0.0004, 'join=mitre') AS geom"

    # Execute the query using execute_values to insert the data
    cursor.execute(query)

    polygons_hex = cursor.fetchone()[0]

    shapely_polygons = wkb.loads(bytes.fromhex(polygons_hex))

    # Sometimes the result of this is a MultiPolygon, with a little dangling one-pixel size polygon
    # This is pretty easy to fix.
    # Here we should determine if it is a MutliPolygon, and if so, discard all but the largest piece
    poly = get_largest(shapely_polygons)
    polygon_str = poly.wkb.hex()

    sql = f"""
    SELECT ST_MakePolygon(geom) as poly FROM (
        SELECT ST_ExteriorRing(geom) as geom FROM (
            SELECT ST_SetSRID('{polygon_str}'::geometry, 4326) as geom
        ) as ring
    ) as poly;
    """
    cursor.execute(sql)
    dissolved_geometry = cursor.fetchone()[0]

    shapely_geometry = wkb.loads(bytes.fromhex(dissolved_geometry))

    dissolved_gdf = geopandas.GeoDataFrame({'geometry': [shapely_geometry]})
    dissolved_gdf.set_geometry('geometry', inplace=True)
    return dissolved_gdf


if __name__ == "__main__":
    # is_installed = postgis_installed()
    # print(is_installed)

    # Test the fast PostGIS-enabled dissolve function
    # First, load some sample data.
    fname = r"C:\Data\GIS\MERITBasins\catchments\simplified\cat_pfaf_27_MERIT_Hydro_v07_Basins_v01.shp"

    gdf = geopandas.read_file(fname)

    # For debugging, let's just do a small selection of unit catchments.
    #comid_list = [27001951,27001973,27001950,27001890,27001969,27001946,27001966,27001943,27001885,27001942,27001941,27001954,27001881,27001938,27001880,27001957,27001878,27001935,27001875,27001932,27001931,27001960,27001930,27001705,27001961,27001929,27001702,27001928,27001927,27001925,27001924,27001963,27001861,27001923,27001858,27001922,27001691,27001958,27001856,27001921,27001920,27001780,27001684,27001683,27001846,27001842,27001676,27001898,27001897,27001674,27001896,27001839,27001895,27001671,27001837,27001836,27001670,27001835,27001669,27001833,27001831,27001664,27001662,27001822,27001658,27001657,27001820,27001655,27001915,27001742,27001653,27001818,27001812,27001576,27001810,27001651,27001630,27001808,27001807,27001805,27001804,27001646,27001803,27001801,27001800,27001643,27001797,27001796,27001795,27001638,27001637,27001750,27001604,27001747,27001745,27001744,27001899,27001559,27001558,27001553,27001546,27001545,27001544,27001543,27001542,27001541,27001540,27001539,27001538]
    #comid_list = [27001650, 27001811, 27001852, 27001863, 27001867, 27001870]
    #selection = gdf[gdf['COMID'].isin(comid_list)]

    selection = gdf
    #selection.plot()
    #plt.show()

    dissolved = dissolve_postgis(selection)
    print(len(dissolved))
    dissolved.plot()
    plt.show()
