"""
Step #1

Merge unit catchments with Strahler stream order = 1 with their immediate downstream neighbors.

Practically, the way to do this is to assign each unit catchment the field `mergeid` then do a dissolve.

Updated so that if the unit catchment's area is in the 98th percentile or up, we do NOT merge it.
This helps keep the distribution of catchment areas slightly more homogeneous.

"""

import psycopg2

# Connect to your PostgreSQL database
conn = psycopg2.connect(
    host="localhost",
    database="basins",
    user="postgres",
    password="dbpw"
)
cur = conn.cursor()

id = 11

# New, let's get the 99th percentile of the areas so we can exclude those!
sql = f"""
SELECT PERCENTILE_CONT(0.99) WITHIN GROUP (ORDER BY unitarea) AS percentile_99
FROM merit_basins_{id};
"""
cur.execute(sql)
p99 = cur.fetchone()[0]


# What this SQL query should do:
# Finds all the MERGES the unit catchments with sorder = 1
# when they drain to the same downstream
# For coastal catchments, where there is no upstream neighbor and which drain to the ocean
# (sorder = 1 AND nextdown = 0), we will exclude these.

sql = f"DROP TABLE IF EXISTS merit_basins2_{id};"
cur.execute(sql)
conn.commit()

sql_query = f"""
CREATE TABLE merit_basins2_{id} AS
SELECT 
    CASE 
        WHEN (sorder = 1 
          AND unitarea < {p99}) THEN nextdown
        ELSE comid
    END AS comid2,
    ST_Union(geom_simple) as geom, 
    COUNT(comid) as num_merged,
    MAX(sorder) as sorder, 
    SUM(unitarea) as unitarea
FROM merit_basins_{id} 
WHERE nextdown > 0 OR sorder > 1
GROUP BY comid2;
"""
cur.execute(sql_query)
conn.commit()

# I could not figure out a reliable way to keep the nextdown field after the dissolve
# so we'll do a table join and extract the info from the original table.
sql = f"""
ALTER TABLE merit_basins2_{id}
ADD COLUMN nextdown INTEGER;
"""
cur.execute(sql)
conn.commit()

sql = f"""
UPDATE merit_basins2_{id} as m2
SET nextdown = mb.nextdown
FROM merit_basins_{id} as mb
WHERE m2.comid2 = mb.comid; 
"""
cur.execute(sql)
conn.commit()

# Inconvenient and confusing to have field named comid2
sql = f"ALTER TABLE merit_basins2_{id} RENAME COLUMN comid2 TO comid;"
cur.execute(sql)

# Add the spatial index!
sql = f"""
    CREATE INDEX merit_basins2_{id}_idx
    ON merit_basins2_{id}
    USING GIST(geom);
"""

# Execute the SQL query
cur.execute(sql)

RIVERS = True
if RIVERS:
    # Rivers
    sql_query = f"""
    DROP TABLE IF EXISTS merit_rivers2_{id};
    CREATE TABLE merit_rivers2_{id} AS
    SELECT r.comid, r.lengthkm, r.uparea, r.geom_simple
    FROM merit_rivers_{id} r
    INNER JOIN merit_basins2_{id} b ON r.comid = b.comid;
    """
    cur = conn.cursor()
    cur.execute(sql_query)

# Close the database connection
conn.commit()
cur.close()
conn.close()

