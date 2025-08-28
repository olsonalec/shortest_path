def mph_to_ms(speed):
    # converts mph to m/s
    return (speed * 1609.344) / (60 * 60)

def calculate_time(speed, distance):
    # speed is in mph, distance is in meters
    speed_in_ms = mph_to_ms(speed)
    return distance / speed_in_ms

# Find the roads that intersect with a given intersection
def map_function(intersection, road_geodataframe):
    # For some reason, geopandas doesn't like single-dimensional lists. I think it's because the library interprets them as a list of coordinates.
    # To get around this, the list of indices has to be a two-dimensional list, where the first element is the actual list.
    indices = [[]]
    for road in road_geodataframe.itertuples():
        if intersection['geometry'].intersects(road.geometry):
            indices[0].append(road.Index)
    return indices

# Find the intersections that intersection with a given road
def map_function2(road, intersection_geodataframe):
    # For some reason, geopandas doesn't like single-dimensional lists. I think it's because the library interprets them as a list of coordinates.
    # To get around this, the list of indices has to be a two-dimensional list, where the first element is the actual list.
    indices = [[]]
    for intersection in intersection_geodataframe.itertuples():
        if road.geometry.intersects(intersection.geometry):
            indices[0].append(intersection.Index)
    return indices

# Find the intersections that are one road segment away from the current intersection
def map_function3(intersection, road_gdf):
    connecting_intersections = [[{}]]
    roads = intersection['Roads'][0]
    # roads = convert_string_to_list(roads)
    for road in roads:
        new_intersections = road_gdf.loc[road]['Intersections'][0]
        # new_intersections = convert_string_to_list(new_intersections)
        for int in new_intersections:
            if (int != intersection.name) and (int not in connecting_intersections[0]):
                connecting_intersections[0][0][int] = float(road_gdf.loc[road]['TimeToTravel'])
    return connecting_intersections
