import json
import geojson


with open('new_data/Hennepin_Roads.geojson', 'r') as fp:
    contents = json.loads(fp.read())

    intersections = []

    # extract starting and ending coordinates of each road segment
    for i in range(len(contents['features'])):
        # these two variables are lists of length 2; index 0 is the x coordinate, and index 1 is the y coordinate
        intersection_1 = contents['features'][i]['geometry']['coordinates'][0]
        intersection_2 = contents['features'][i]['geometry']['coordinates'][1]
        
        # make sure there are no duplicates because multiple road segments will have the same starting or ending coordaintes
        if intersection_1 not in intersections:
            intersections.append(intersection_1)
        if intersection_2 not in intersections:
            intersections.append(intersection_2)
    

# create new geojson file of intersections:

# convert each intersection, which is a list of two coordinates, into a Point object
intersection_Points = []
for i in range(len(intersections)):
    intersection_Points.append(geojson.Point(intersections[i]))

# write the intersections to a geojson file
name = 'Hennepin_Intersections'
with open(f'new_data/{name}.geojson', 'w') as fp:
    # header
    fp.write(
        '{\n' \
        '\"type\": \"FeatureCollection\",\n' \
        '\"name\": \"' + name + '\",\n' \
        '\"crs\": { \"type\": \"name\", \"properties\": { \"name\": \"urn:ogc:def:crs:EPSG::26915\" } },\n' \
        '\"features\": [\n' 
    )

    # content
    length = len(intersection_Points)
    for i in range(length):
        fp.write(f'{{ \"type\": \"Feature\", \"properties\": {intersection_Points[i]}}}')
        if i != length - 1:
            fp.write(',\n')

    # footer
    fp.write('\n]' \
             '\n}')
