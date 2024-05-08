"""
Inserts the lat, lng coordinates for the *basin outlet* into the MERIT-Basins river reaches tables.
This is simply the end point of the river reach.

Note that the rivers are ALWAYS backwards, and so their endpoint with PostGIS function `ST_StartPoint()`

I simply write the information back to the table.
"""

import psycopg2

# Connect to your PostgreSQL database
conn = psycopg2.connect(
    host="localhost",
    database="basins",
    user="postgres",
    password="dbpw"
)

# 27
basins = [11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 28, 29, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44,
          45, 46, 47, 48, 51, 52, 53, 54, 55, 56, 57, 61, 62, 63, 64, 65, 66, 67, 71, 72, 73, 74, 75, 76, 77, 78,
          81, 82, 83, 84, 85, 86, 91]

cur = conn.cursor()

for basin in basins:
    print(basin)
    sql = f"""
    ALTER TABLE merit_rivers_{basin} ADD COLUMN end_lat double precision;
    ALTER TABLE merit_rivers_{basin} ADD COLUMN end_lng double precision;

    UPDATE merit_rivers_{basin} 
    SET 
    end_lat = ST_Y(ST_StartPoint(geom_detail)),
    end_lng = ST_X(ST_StartPoint(geom_detail));
    """
    cur.execute(sql)

conn.commit()
cur.close()
conn.close()
