"""
Testing the different dissolve methods to find out which is fastest

1. Regular geopandas dissolve
2. My custom fast_dissolve() function
3. PostGIS-enabled dissolve
"""

import geopandas
import time

import matplotlib.pyplot as plt

from postgis import dissolve_postgis
from fast_dissolve import dissolve_geopandas


def timer(func):
    """
    Decorator function to time another Python function
    Usage:
       @timer
       def myfunction()

    timer(myfunction)
    """
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = round(end_time - start_time, 3)
        print(f"Execution time of '{func.__name__}': {execution_time} seconds")
        return execution_time, result
    return wrapper


@timer
def dissolver(gdf: geopandas.GeoDataFrame, opt: int) -> geopandas.GeoDataFrame:
    if opt == 1:
        return geopandas.geoseries.GeoSeries([geom for geom in gdf.unary_union.geoms])
    elif opt == 2:
        return dissolve_geopandas(gdf)
    else:
        return dissolve_postgis(gdf)

print("Reading shapefile")
fname = r"C:\Data\GIS\MERITBasins\catchments\simplified\cat_pfaf_63_MERIT_Hydro_v07_Basins_v01.shp"
gdf = geopandas.read_file(fname)

methods = ["GeoPandas", "fast_dissolve()", "PostGIS"]
for i in range(1, 3):
    print(f"Dissolving with {methods[i]}")
    t, dissolved_gdf = dissolver(gdf, i)
    dissolved_gdf.plot()
    plt.title(methods[i])
    fname = f"diss{i}.jpg"
    plt.savefig(fname)
