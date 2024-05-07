"""Removing the column `basin` from my unit catchments in PostgreSQL database
"""
import psycopg2

# Connect to your PostgreSQL database
conn = psycopg2.connect(
    host="localhost",
    database="basins",
    user="postgres",
    password="dbpw"
)

#11
basins = [12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44,
          45, 46, 47, 48, 51, 52, 53, 54, 55, 56, 57, 61, 62, 63, 64, 65, 66, 67, 71, 72, 73, 74, 75, 76, 77, 78,
          81, 82, 83, 84, 85, 86, 91]

cur = conn.cursor()

for basin in basins:
    sql = f"""
    ALTER TABLE merit_basins_{basin}
        ADD COLUMN numup smallint,
        ADD COLUMN up1 integer,
        ADD COLUMN up2 integer,
        ADD COLUMN up3 integer,
        ADD COLUMN up4 integer,
        ADD COLUMN nextdown integer,
        ADD COLUMN sorder smallint;
    """

    cur.execute(sql)
    conn.commit()

    sql = f"""
    UPDATE merit_basins_{basin} as mb                 
    SET numup = mr.maxup,
        up1 = mr.up1,
        up2 = mr.up2,
        up3 = mr.up3,
        up4 = mr.up4,
        nextdown = mr.nextdownid,
        sorder = mr.sorder
    FROM merit_rivers_{basin} as mr
    WHERE mb.comid = mr.comid; 
    """

    cur.execute(sql)
    conn.commit()

cur.close()
conn.close()
