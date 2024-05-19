import geopandas as gpd
from shapely.geometry import Point, MultiPolygon

# Create GeoDataFrame with Points
points_data = {
    'geometry': [Point(1, 1), Point(2, 2), Point(3, 3)],
    'data': [1, 2, 3]
}
gdf_points = gpd.GeoDataFrame(points_data, crs="EPSG:4326")

# Create GeoDataFrame with MultiPolygons
multipolygons_data = {
    'geometry': [MultiPolygon([Point(0, 0).buffer(2).envelope, Point(5, 5).buffer(2).envelope])],
    'data': ['A']
}
gdf_multipolygons = gpd.GeoDataFrame(multipolygons_data, crs="EPSG:4326")

# Create GeoDataFrame with Points
points_data = {
    'geometry': [Point(1, 1), Point(2, 2), Point(3, 3)],
    'data': [1, 2, 3]
}
gdf_points = gpd.GeoDataFrame(points_data, crs="EPSG:4326")

# Create GeoDataFrame with MultiPolygons
multipolygons_data = {
    'geometry': [MultiPolygon([Point(0, 0).buffer(2).envelope, Point(5, 5).buffer(2).envelope])],
    'data': ['A']
}
gdf_multipolygons = gpd.GeoDataFrame(multipolygons_data, crs="EPSG:4326")

# Perform the overlay analysis
result = gpd.overlay(gdf_points, gdf_multipolygons, how='intersection')

# Display the result
print(result)
