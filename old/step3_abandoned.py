"""
Step #3 We just did some merging. Now we need to 'invert' the network.
For every node, we have the `nextdown` which is the to_node.
However, for each node, we also need to know the number of upstream nodes, `numup`
and have a list of the upstream nodes. Since there are a maximum of 4 in our dataset,
I found it convenient to put these in fields up1, up2, up3, and up4.

This is important later on for the delineation algorithms.

TODO: I think this could be made more efficient by using the EDGES from the previous step.
  Just create a data structure that has `from_node` and `to_node`
  And then we can do a graph 'reverse' - I think there is a built-in function for this.

TODO: update the Strahler and Shreve stream orders. Because we can basically just repeat step 2 (right?).

"""

import warnings
warnings.filterwarnings('ignore')
import psycopg2
import geopandas as gpd
import pandas as pd
from sqlalchemy import create_engine
from graph_tools import calculate_strahler_stream_order, calculate_shreve_stream_order

id = 27

host = "localhost"
database = "basins"
user = "postgres"
password = "dbpw"
port = 5432

# Connect to your PostgreSQL database
conn = psycopg2.connect(
    host=host,
    database=database,
    user=user,
    password=password
)

# Create a SQLAlchemy engine
engine = create_engine('postgresql://{}:{}@{}:{}/{}'.format(user, password, host, port, database))

# Define your SQL query
sql_query = f"""
SELECT *
FROM merit_basins2_{id};
"""

# Read the result of the SQL query into a GeoDataFrame
gdf = gpd.read_postgis(sql_query, conn, geom_col='geom')
gdf.set_index('comid2', inplace=True)

# Create a dictionary relating comid: nextdown
# D for down
D = gdf['nextdown'].to_dict()

# Now we will build an inverted network dictionary.
# R for reversed
R = {}
for node, ds in D.items():
    if ds == 0:
        continue
    if ds not in R:
        R[ds] = [node]
    else:
        mylist = R[ds]
        mylist.append(node)
        R[ds] = mylist

# Now, update the table with this new information
C = {}
U1 = {}
U2 = {}
U3 = {}
U4 = {}

for k, mylist in R.items():
    C[k] = len(mylist)
    if len(mylist) > 0:
        U1[k] = mylist[0]
    if len(mylist) > 1:
        U2[k] = mylist[1]
    if len(mylist) > 2:
        U3[k] = mylist[2]
    if len(mylist) > 3:
        U4[k] = mylist[3]
    if len(mylist) > 4:
        print("Not expecting more than 4 upstream nodes. Investigate.")
        raise Exception

gdf = gdf.assign(numup=pd.Series(C))
gdf['numup'] = gdf['numup'].fillna(0)
gdf = gdf.assign(up1=pd.Series(U1))
gdf = gdf.assign(up2=pd.Series(U2))
gdf.reset_index(inplace=True)

gdf.to_postgis('merit_basins3_27', engine, if_exists='replace')
