import geopandas as gpd
import geojson
import json
import shapely


roads_gpd = gpd.read_file('new_data/Seven_Counties_Roads_v2.geojson')


for i in range(roads_gpd.shape[0]): # in roads_gpd.itertuples():
    # print(type(tuple.geometry))
    if not isinstance(roads_gpd.loc[i, 'geometry'], shapely.geometry.linestring.LineString):
        print(f'beginning: {(roads_gpd.loc[i, 'geometry'])}')
        coordinates = shapely.get_coordinates(roads_gpd.loc[i, 'geometry'])
        new_coordinates = {}
        # print(coordinates)
        # print()
        for j in range(coordinates.shape[0]):
            if coordinates[j][0] not in new_coordinates:
                new_coordinates[coordinates[j][0]] = [coordinates[j][1]]
            elif coordinates[j][1] not in new_coordinates[coordinates[j][0]]:
                new_coordinates[coordinates[j][0]].append(coordinates[j][1])

        new_list = []
        for key, value in new_coordinates.items():
            for j in range(len(value)):
                new_list.append([round(float(key), 4), round(float(value[j]), 4)])
        # print(new_list)

        new_coords = shapely.LineString(new_list)
        # print(new_coords)

        # print(type(new_coords))

        roads_gpd.loc[i, 'geometry'] = new_coords

        print(f'end: {(roads_gpd.loc[i])}')


# roads_gpd.to_file('new_data/Seven_Counties_Roads_v2.geojson')



    # if type(tuple.geometry) != 'shapely.geometry.linestring.LineString':
        # print(tuple.geometry)