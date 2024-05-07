import psycopg2
from beep import beep

# Connect to your PostgreSQL database
conn = psycopg2.connect(
    host="localhost",
    database="basins",
    user="postgres",
    password="dbpw"
)

#27
basins = [11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 28, 29, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44,
          45, 46, 47, 48, 51, 52, 53, 54, 55, 56, 57, 61, 62, 63, 64, 65, 66, 67, 71, 72, 73, 74, 75, 76, 77, 78,
          81, 82, 83, 84, 85, 86, 91]

cur = conn.cursor()

for basin in basins:
    print(basin)
    sql = f"""
    DROP TABLE IF EXISTS merit_outlets_{basin};
    
    CREATE TABLE merit_outlets_{basin} (
        comid integer,  
        geom geometry(Point, 4326)  
    );
    
    INSERT INTO merit_outlets_{basin} (comid, geom)
    SELECT comid, ST_StartPoint(geom_detail) as geom
    FROM merit_rivers_{basin};
    
    CREATE INDEX idx_outlets_{basin} ON merit_outlets_{basin} USING GIST (geom);
    """
    cur.execute(sql)

conn.commit()
cur.close()
conn.close()



beep()
