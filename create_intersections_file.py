import json
import geojson
import time
import sys

start = time.time()
with open('new_data/Seven_Counties_Roads_v2.geojson', 'r') as fp:
    contents = json.loads(fp.read())

    intersections = {}      # keys are x coordinates; values are lists of y coordinates

    # extract starting and ending coordinates of each road segment
    for i in range(len(contents['features'])):
        if contents['features'][i]['geometry']['type'] == 'MultiLineString':
            print(contents['features'][i])
        else:
        
            # these two variables are lists of length 2; index 0 is the x coordinate, and index 1 is the y coordinate
            intersection_1 = contents['features'][i]['geometry']['coordinates'][0]      # starting coordinates of road segment
            intersection_2 = contents['features'][i]['geometry']['coordinates'][1]      # ending coordinates of road segment

            # make sure there are no duplicates because multiple road segments will have the same starting or ending coordaintes
            if intersection_1[0] not in intersections:      # check if x coordinate is in the dictionary
                intersections[intersection_1[0]] = [intersection_1[1]]
            else:
                if intersection_1[1] not in intersections[intersection_1[0]]:   # if the x coordinate is already there, check to see if the y coordinate is there as well
                    intersections[intersection_1[0]].append(intersection_1[1])
            if intersection_2[0] not in intersections:
                intersections[intersection_2[0]] = [intersection_2[1]]
            else:
                if intersection_2[1] not in intersections[intersection_2[0]]:       # if the x coordinate is already there, check to see if the y coordinate is there as well
                    intersections[intersection_2[0]].append(intersection_2[1])
    
# create new geojson file of intersections:

# convert each intersection, which is a list of two coordinates, into a Point object
intersection_Points = []
for key, value in intersections.items():
    for i in range(len(value)):
        intersection_Points.append(geojson.Point([key, value[i]]))
# for i in range(len(intersections)):
    # intersection_Points.append(geojson.Point(intersections[i]))

# write the intersections to a geojson file
name = 'Seven_Counties_Intersections_v2'
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

end = time.time()
print(end - start)