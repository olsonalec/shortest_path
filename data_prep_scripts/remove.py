import geopandas as gpd

x = 0

if x == 1:
    roads = gpd.read_file('../data/Hennepin_Roads_Prepped.geojson')

    print(roads.info())

    roads = roads.drop('CTU_NAME_L', axis=1)
    roads = roads.drop('CTU_NAME_R', axis=1)
    roads = roads.drop('ROUTE_ID', axis=1)
    roads = roads.drop('OBJECTID', axis=1)
    roads = roads.drop('FID', axis=1)
    roads = roads.drop('StartX', axis=1)
    roads = roads.drop('EndX', axis=1)
    roads = roads.drop('StartY', axis=1)
    roads = roads.drop('EndY', axis=1)

    print(roads.info())

    roads.to_file('data/Hennepin_Roads_Prepped3.geojson')
else:
    intersections = gpd.read_file('../data/Hennepin_Intersections_Prepped.geojson')

    print(intersections.info())

    intersections = intersections.drop('OBJECTID', axis=1)
    intersections = intersections.drop('POINT_X', axis=1)
    intersections = intersections.drop('POINT_Y', axis=1)

    print(intersections.info())

    intersections.to_file('data/Hennepin_Intersections_Prepped3.geojson')