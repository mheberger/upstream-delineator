"""
Remove the columns end_lat and end_lng from merit_basins tables.
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
    ALTER TABLE merit_basins_{basin}
        DROP COLUMN end_lat,
        DROP COLUMN end_lng
    """

    cur.execute(sql)

conn.commit()
cur.close()
conn.close()
