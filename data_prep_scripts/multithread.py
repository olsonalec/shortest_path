import geopandas as gpd
import pandas as pd
import data_prep
import time
from pandarallel import pandarallel

'''
Function to convert speed in mph to m/s.

Parameters:
    speed (int) - speed limit of the road in mph

Return value:
    (float) the speed limit in m/s
'''
def mph_to_ms(speed):
    return (speed * 1609.344) / (60 * 60)


'''
Function to calculate the approximate time it takes to travel a road segment.

Parameters:
    speed (int) - speed limit of the road in mph
    distance (float) - distance of the road segment in meters

Return Value:
    (float) time in seconds taken to travel the road segment
'''
def calculate_time(speed, distance):
    speed_in_ms = mph_to_ms(speed)
    return distance / speed_in_ms


if __name__ == '__main__':

    roads = gpd.read_file('../data/Hennepin_Roads.geojson')
    intersections = gpd.read_file('../data/Hennepin_Intersections.geojson')

    roads_speedlimit = []

    # Add speedlimits to roads that don't currently have one
    i = 0
    for road in roads.itertuples():
        speedlimit = road.ROUTESPEED
        if speedlimit == 0:
            roads.loc[i, 'ROUTESPEED'] = 30
        i += 1

    # Calculate approximate time to travel each road segment
    roads_time = [] # will store the amount of time it takes to travel each road segment
    for road in roads.itertuples():
        road_speed = road.ROUTESPEED
        road_distance = road.Shape_Length
        roads_time.append(calculate_time(road_speed, road_distance))

    roads_time_df = gpd.GeoDataFrame(pd.Series(roads_time))
    roads['TimeToTravel'] = roads_time_df

    start_time = time.time()

    pandarallel.initialize(progress_bar=True)

    # Find the roads that intersect with a given intersection
    intersections['Roads'] = intersections.parallel_apply(data_prep.map_function, axis=1, args=(roads,))

    # Find the intersections that intersection with a given road
    roads['Intersections'] = roads.parallel_apply(data_prep.map_function2, axis=1, args=(intersections,))

    # Find the intersections that are one road segment away from a given intersection, and calculate the time it takes to reach each one
    intersections['NeighboringIntersections'] = intersections.parallel_apply(data_prep.map_function3, axis=1, args=(roads,))

    # Calculating the speedlimit and time to travel could also be parallelized

    end_time = time.time()
    print(f'Total time: {end_time - start_time}')

    intersections.to_file('../data/Hennepin_Intersections_Prepped2.geojson')
    roads.to_file('../data/Hennepin_Roads_Prepped2.geojson')
