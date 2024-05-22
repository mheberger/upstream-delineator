"""
Step #4: Prune the network

We are going to merge adjacent nodes whenever there is only
one from_node and one to_node. This makes the network more sparse and topologically correct.
It also has the desired effect of making the subcatchments larger.

Later learned this is called network pruning, more specifically has been called
a "linear branch collapser" pruner. Citaton: https://doi.org/10.1186/s12859-015-0486-3
(paper has nothing to do with hydrology or physical sciences!)


TODO: This one needs a lot of work! This is a kind of tricky problem!

The current mering algorithm is not aggressive enough.
    What we need to do is find and remove the very small branch nodes,
    but we have to do it intelligently.
  is not preserved. Need to go back and figure this out!!! Obviously extremely important!
  I think there should be separate logic for branch nodes. Can search for branch nodes where
  numup > 1 and unitarea < p10 and automatically merge them with their downstream neighbor

"""

import warnings
warnings.filterwarnings('ignore')
import psycopg2
import geopandas as gpd
from sqlalchemy import create_engine

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

# Let's figure out the size of a really small basin (<5%-ile) and a rather large basin (>80%-ile)
# so that we can do another round of merging, intelligently.
cur = conn.cursor()
sql = f"""
SELECT PERCENTILE_CONT(0.10) WITHIN GROUP (ORDER BY unitarea) AS percentile_05
FROM merit_basins3_{id};
"""
cur.execute(sql)
p05 = cur.fetchone()[0]

sql = f"""
SELECT PERCENTILE_CONT(0.90) WITHIN GROUP (ORDER BY unitarea) AS percentile_80
FROM merit_basins_{id};
"""
cur.execute(sql)
p80 = cur.fetchone()[0]


# Create a SQLAlchemy engine
engine = create_engine('postgresql://{}:{}@{}:{}/{}'.format(user, password, host, port, database))

# Read in data from Step #3
sql_query = f"""
SELECT *
FROM merit_basins3_{id};
"""

# Read the result of the SQL query into a GeoDataFrame
gdf = gpd.read_postgis(sql_query, conn, geom_col='geom')

# Add a new column where we will put the comid of the branch node with which a node shall be merged.
gdf['mergeid'] = gdf['comid2']

gdf.set_index('comid2', inplace=True)

numup = gdf['numup'].to_dict()
order = gdf['sorder'].to_dict()
to_node = gdf['nextdown'].to_dict()
U1 = gdf['up1'].to_dict()
D = gdf['nextdown'].to_dict()
A = gdf['unitarea'].to_dict()

# Get the first set of internal or "stem" nodes:  ◙---◙---o
#stem_nodes = gdf[ (gdf['numup'] < 2) & (gdf['nextdown'] !=0) ].index.tolist()
#stem_nodes.sort()
#print(len(gdf))
#print(len(stem_nodes))

# The full set of target nodes to dissolve will include the stem nodes
#target_nodes = stem_nodes.copy()
#print(len(target_nodes))

# For each of the target nodes, find the downstream node that is a branch, or outlet.

def get_mergeid(node):
    # Here is the logic to decide whether each unit catchment (node) should get
    cumulative_area = 0
    node_area = A[node]
    cumulative_area = cumulative_area + node_area

    while True:
        # find the next downstream node
        ds_node = D[node]
        # Find out whether the downstream node is a branch, i.e. has more than one piece flowing into it.
        n_up = numup[ds_node]

        if n_up > 1:
            # If the number of upstream nodes > 1, that means it's a branch, and we should NOT merge.
            return node

        # Consider whether the cumulative area will be to big.
        dsnode_area = A[node]
        cumulative_area = cumulative_area + dsnode_area

        # Some logic based on the area of the nodes.
        # If the node is large, don't continue, UNLESS the other is very small, merge.
        if (cumulative_area > p80) and not (node_area < p05 or dsnode_area < p05):
            return node

        if D[ds_node] == 0:
            # If this is zero, it means that we have hit a terminal (ocean) node, so we can stop
            return ds_node

        node = ds_node


# Iterate over each of the nodes that is a candidate for merging with its downstream neighbor.
# Only do the merge if it or its neighbor are really small, OR if merging them would not create
# an abnormally large catchment!

nodes = gdf[gdf['nextdown'] > 0].index.tolist()

for node in nodes:
    if node == 27001184:
        print(1)
    merge_node = get_mergeid(node)
    gdf.loc[node, 'mergeid'] = merge_node
    gdf.loc[node, 'nextdown'] = D[merge_node]


# TODO: This is mostly working the way I want it, but I am losing the correct nextdown value
#  Better will be to load a temp table into Postgres DB, then use PostGIS commands to do the dissolve. 

dissolved_gdf = gdf.dissolve(by='mergeid', aggfunc='sum')
print(len(dissolved_gdf))

dissolved_gdf.reset_index(inplace=True)
dissolved_gdf.rename(columns={'mergeid': 'comid2'}, inplace=True)
dissolved_gdf.to_postgis('merit_basins4_27', engine, if_exists='replace')