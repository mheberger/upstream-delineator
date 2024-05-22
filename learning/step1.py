"""
Step #1

Merge unit catchments with Strahler stream order = 1 with their immediate downstream neighbors.

Practically, the way to do this is to assign each unit catchment the field `mergeid` then do a dissolve.

Updated so that if the unit catchment's area is in the 98th percentile or up, we do NOT merge it.
This helps keep the distribution of catchment areas slightly more homogeneous.

"""
from update_network_data import update_network
from util.db import cursor, db_close


def step1(basin: int, max_area_threshold: int):
    """
    First step in reducing the number of subcatchments, while preserving the total area and the
    connectivity of the river network.

    Finds all the MERGES the unit catchments with sorder = 1 which drain to same downstream neighbor
    For coastal catchments, where there is no upstream neighbor, and they drain to the ocean
    (sorder = 1 AND nextdown = 0), we will exclude these.

    Creates the temporary tables `basins1` and `river1`
    """

    sql = f"DROP TABLE IF EXISTS basins1;"
    cursor.execute(sql)

    sql_query = f"""
    CREATE TABLE basins1 AS
    SELECT 
        CASE 
            WHEN (sorder = 1 
              AND unitarea < {max_area_threshold}) THEN nextdown
            ELSE comid
        END AS comid2,
        ST_Union(geom_simple) as geom, 
        COUNT(comid) as num_merged,
        MAX(sorder) as sorder, 
        SUM(unitarea) as unitarea
    FROM merit_basins_{basin} 
    WHERE nextdown > 0 OR sorder > 1
    GROUP BY comid2;
    """
    cursor.execute(sql_query)

    # I could not figure out a reliable way to keep the nextdown field after the dissolve
    # so we'll do a simple table join and extract the info from the original table.
    sql = f"""
    ALTER TABLE basins1
    ADD COLUMN nextdown INTEGER;
    """
    cursor.execute(sql)

    sql = f"""
    UPDATE basins1 as m2
    SET nextdown = mb.nextdown
    FROM merit_basins_{basin} as mb
    WHERE m2.comid2 = mb.comid; 
    """
    cursor.execute(sql)

    # I created the new field `comid2` to avoid conflict with `comid`,
    # but it's inconvenient and confusing, so rename it back to `comid`.
    sql = f"ALTER TABLE basins1 RENAME COLUMN comid2 TO comid;"
    cursor.execute(sql)

    # Rivers
    sql_query = f"""
    DROP TABLE IF EXISTS rivers1;
    CREATE TABLE rivers1 AS
    SELECT r.comid, r.lengthkm, r.uparea, r.geom_simple
    FROM merit_rivers_{basin} r
    INNER JOIN basins1 b ON r.comid = b.comid;
    """
    cursor.execute(sql_query)

    # Close the database connection
    db_close()


if __name__ == "__main__":
    step1(11, 150)
    # Every time I run one of the steps, I have to update the river network data (shreve, strahler, numup)
    update_network("basins1")
