import geopandas as gpd
import pandas as pd
import ast
import sys

def mph_to_ms(speed):
    # converts mph to m/s
    return (speed * 1609.344) / (60 * 60)

def calculate_time(speed, distance):
    # speed is in mph, distance is in meters
    speed_in_ms = mph_to_ms(speed)
    return distance / speed_in_ms

'''
The 'Intersections' and 'Roads' attributes in the Roads and Intersections Geodataframes, respectively, are represented as a nested list.
If they were a one-dimensional list, GeoPandas would think that they represent geometry, which they don't.
However, when reading these attributes from the dataframe, Python interprets them as strings.
This is an example: '[[180, 240, 360]]'
This function converts this string representation of a nested list into a 1-dimensional Python list.
The example output would be [180, 240, 360], where each value is an integer.

Parameter:
    bad_string - a string representation of a nested list

Return Value:
    new_list[0] - a Python list representation of the input
'''
def convert_string_to_list(bad_string):
    new_list = ast.literal_eval(bad_string)
    return new_list[0]

# roads = gpd.read_file(r"Plymouth_Roads.geojson")

# intersections = gpd.read_file(r"Plymouth_Intersections.geojson")

# print(roads.head())
# print(intersections.head())

# Find the roads that intersect with a given intersection
def map_function(intersection, road_geodataframe):
    indices = [[]]
    for road in road_geodataframe.itertuples():
        if intersection['geometry'].intersects(road.geometry):
            indices[0].append(road.Index)
    return indices

# Find the intersections that intersection with a given road
def map_function2(road, intersection_geodataframe):
    indices = [[]]
    for intersection in intersection_geodataframe.itertuples():
        if road.geometry.intersects(intersection.geometry):
            indices[0].append(intersection.Index)
    return indices

# Find the intersections that are one road segment away from the current intersection
def map_function3(intersection, road_gdf):
    connecting_intersections = [[]]
    roads = intersection['Roads']
    roads = convert_string_to_list(roads)
    for road in roads:
        new_intersections = road_gdf.loc[road]['Intersections']
        new_intersections = convert_string_to_list(new_intersections)
        for int in new_intersections:
            if (int != intersection.name) and (int not in connecting_intersections[0]):
                connecting_intersections[0].append((int, float(road_gdf.loc[road]['TimeToTravel'])))
    return connecting_intersections
