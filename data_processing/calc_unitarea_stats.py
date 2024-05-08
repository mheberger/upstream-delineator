"""

"""

import psycopg2
from util.db import connection, db_close
import pandas
from util.beep import beep

# Connect to your PostgreSQL database

basins = [11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44,
          45, 46, 47, 48, 51, 52, 53, 54, 55, 56, 57, 61, 62, 63, 64, 65, 66, 67, 71, 72, 73, 74, 75, 76, 77, 78,
          81, 82, 83, 84, 85, 86, 91]

results_df = pandas.DataFrame(columns=['basin', 'p01', 'p02', 'p25','p50','p75','p98','max_area_threshold'])

for basin in basins:
    print(basin)
    sql = f"""
    SELECT 
        {basin} as basin,
        PERCENTILE_CONT(0.01) WITHIN GROUP (ORDER BY unitarea) AS p01,
        PERCENTILE_CONT(0.02) WITHIN GROUP (ORDER BY unitarea) AS p02,
        PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY unitarea) AS p25,
        PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY unitarea) AS p50,
        PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY unitarea) AS p75,
        PERCENTILE_CONT(0.98) WITHIN GROUP (ORDER BY unitarea) AS p98,
        PERCENTILE_CONT(0.99) WITHIN GROUP (ORDER BY unitarea) AS max_area_threshold
    FROM merit_basins_{basin};
    """
    # Execute the query and fetch the result
    result = pandas.read_sql_query(sql, connection)

    # Append the result as a new row to the DataFrame
    if basin == 11:
        results_df = result
    else:
        results_df = pandas.concat([results_df, result], ignore_index=True)

results_df.to_csv('area_stats.csv', index=False)

db_close()

beep()
