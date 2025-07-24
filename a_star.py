import geopandas as gpd
import matplotlib.pyplot as plt
import ast
import time
import math


roads_gpd = gpd.read_file('Ramsey_Roads_Prepped.geojson')
intersections_gpd = gpd.read_file('Ramsey_Intersections_Prepped.geojson')


source = 2     # an intersection
dest = 2000       # an intersection

class Road:
    def __init__(self, index, cost, speed, start_x, start_y, end_x, end_y):
        self.index = index      # the index into the roads_gpd where this road is stored
        self.travel_time = cost        # the cost to travel this road segment
        self.speed = speed      # speed limit on this road
        self.start_x = start_x  # x coordinate of the starting point
        self.start_y = start_y  # x coordinate of the ending point
        self.end_x = end_x      # x coordinate of the ending point
        self.end_y = end_y      # y coordinate of the ending point

class Vertex:
    def __init__(self, index, x, y):
        self.x = x      # x coordinate of the intersection
        self.y = y      # y coordinate of the intersection
        self.cost = float("inf")        # the total cost of the path to reach this vertex
        self.prev = None    # the previous Vertex in the path
        self.index = index      # the index into the intersections_gpd where this vertex is stored
        self.connections = {}   # a dict of Vertex objects that can be reached by this Vertex - i.e. intersections that are one road segment away; keys are Vertex objects and values are Road objects

'''
A function to intialize the Vertex objects.

Parameters:
    gdf (Geodataframe) - a GeoPandas Geodataframe containing intersection points
    source (int) - the starting vertex; this is an index into the Geodataframe

Return Value:
    vertex_list (list) - a list of Vertex objects
'''
def initialize_vertices(gdf, source):
    vertex_list = []        # a list to store a Vertex object associated with each intersection
    for intersection in gdf.itertuples():
        index = intersection.Index
        new_vertex = Vertex(index, intersection.geometry.x, intersection.geometry.y)
        vertex_list.append(new_vertex)

    # Change the starting vertex's cost to 0
    vertex_list[source].cost = 0

    return vertex_list

'''
A function to initialize the connection attribute for each Vertex object.
The connection attribute stores the index of other intersections that are one road segment away from a given intersection.
Parameters:
    int_gdf - a GeoPandas Geodataframe containing intersection points
    road_gdf - a GeoPandas Geodataframe containing roads
    vertices - a list of Vertex objects
'''
def initialize_ints(int_gdf, road_gdf, vertices):
    for vertex in vertices:
        roads = int_gdf.loc[vertex.index]['Roads']
        roads = convert_string_to_list(roads)

        # find the intersections that are one road segment away from the current intersection
        connections = {}
        for road in roads:
            new_intersections = road_gdf.loc[road]['Intersections']
            new_intersections = convert_string_to_list(new_intersections)
            for intersection in new_intersections:
                if (intersection != vertex.index) and (intersection not in connections):
                    linestring_coordinates = list(road_gdf.loc[road].geometry.coords)
                    start_x = linestring_coordinates[0][0]
                    start_y = linestring_coordinates[0][1]
                    end_x = linestring_coordinates[-1][0]
                    end_y = linestring_coordinates[-1][1]
                    connections[vertices[intersection]] = Road(road, road_gdf.loc[road]['TimeToTravel'], road_gdf.loc[road]['ROUTESPEED'], start_x, start_y, end_x, end_y)
                    # connections.append(vertices[intersection])
            vertex.connections = connections
    
'''
A function to find the index of the minimum element in a list of unsorted objects.

Parameter:
    list_of_vertices (list) - a list of Vertex objects

Return Value:
    min_vertex (Vertex) - the Vertex object that has the lowest travel cost in the list
'''
def find_min(list_of_vertices):
    min_value = float("inf")
    min_vertex = None
    n = len(list_of_vertices)
    for i in range(n):
        if list_of_vertices[i].cost < min_value:
            min_value = list_of_vertices[i].cost
            min_vertex = list_of_vertices[i]
    return min_vertex

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

'''
Converts seconds into minutes and seconds.

Parameter:
    seconds (float) - the value of seconds

Return Values:
    min (int) - an integer representing the number of minutes
    sec (int) - an integer representing the number of seconds
'''
def convert_sec_to_min(seconds):
    return math.ceil(seconds / 60)
    # min = int(seconds // 60)
    # sec = int(seconds % 60)
    # return min, sec

def add_buffer_time(minutes):
    return math.ceil(minutes * 1.1)  

def calculate_distance(start_x, start_y, end_x, end_y):
    return ((((start_x - end_x) **2) + (start_y - end_y) ** 2) ** 0.5)

'''
A function to find the shortest path between two intersections.

Parameters:
    graph - a list containing all the unvisited vertices in the graph; this is the list returned by the initialize() function
    dest - the index into the intersections Geodataframe of the destination vertex

Return Value:
    visited_intersections (list) - a list of all the Vertex objects that were visited by the algorithm
                                    This list contains Vertex objects that are on the shortest path as well as Vertex objects that are not in the shortest path
                                    The prev attribute of each Vertex object will be used to determine the shortest path. After this function returns, simply find
                                        the ending vertex, and use the prev attribute to work your way backwards until you reach the starting vertex. This process
                                        returns the shortest path.
'''
def a_star(graph, dest, dest_x, dest_y):
    visited_vertices = []       # a list containing all the indices of all the vertices that have been visited so far
    visited_intersections = []  # a list containing all the Vertex objects that have been visited so far
    chosen_roads = []           # a list containing indicies for all the roads for the route
    while dest not in visited_vertices:
        min_vertex = find_min(graph)        # find the Vertex in the graph with the lowest travel cost

        # update distances to the neighbor nodes
        neighbors = min_vertex.connections
        for neighbor, road in neighbors.items():
            g = min_vertex.cost + road.travel_time

            # calculate h
            # find out if the starting or ending point of the road is where the intersection is at
            if neighbor.x == road.start_x and (neighbor.y == road.start_y):
                # calculate distance between destination vertex and the road segment's ending coordinates
                h = calculate_distance(road.end_x, dest_x, road.end_y, dest_y) * road.speed
            else:
                # use the roads starting coordinates
                h = calculate_distance(road.start_x, dest_x, road.start_y, dest_y) * road.speed

            new_cost = g + h
            # new_cost = min_vertex.cost + road.travel_time
            if new_cost < neighbor.cost:
                neighbor.cost = g
                neighbor.prev = min_vertex
        
        # add the current node to the set of visited nodes
        visited_vertices.append(min_vertex.index)
        visited_intersections.append(min_vertex)
        chosen_roads.append(min_vertex.connections)

        # remove the current vertex from the set of unvisited nodes
        graph.remove(min_vertex)
    
    return visited_intersections, visited_vertices
    

list_of_vertices = initialize_vertices(intersections_gpd, source)
initialize_ints(intersections_gpd, roads_gpd, list_of_vertices)

start_time = time.time()
dest_x = intersections_gpd.loc[dest, 'geometry'].x
dest_y = intersections_gpd.loc[dest, 'geometry'].y
visited_intersections, visited_vertices = a_star(list_of_vertices, dest, dest_x, dest_y)
end_time = time.time()

# Find the shortest path
chosen_vertices = []        # a list of Vertex objects that compose the shortest path
chosen_vertices_idxs = [dest]   # a list of indices into the intersections_gdb that compose the shortest path; this list is just used for graphing the vertices
dest_idx = -1
for i in range(len(visited_intersections)):
    if visited_intersections[i].index == dest:
        dest_idx = i
        break

prev_intersection = visited_intersections[dest_idx].prev
total_travel_time = visited_intersections[dest_idx].cost
minutes = convert_sec_to_min(total_travel_time)

# Start with the destination Vertex and walk backwards until you've reached the source vertex
while prev_intersection != None:
    chosen_vertices.append(prev_intersection)
    chosen_vertices_idxs.append(prev_intersection.index)
    prev_intersection = prev_intersection.prev

print(f'The time it will take to travel this route is approximately {add_buffer_time(minutes)} minutes.')
print(f'Time taken to run the A* Algorithm: {end_time - start_time} seconds\n')

print(f'Total number of intersections: {intersections_gpd.shape[0]}')
print(f'Total number of intersections visited by the A* algorithm: {len(visited_vertices)}\n')

fig, ax = plt.subplots()

roads_gpd.plot(ax=ax)

# Plot the intersections that were chosen to create the fastest route
intersections_gpd.loc[chosen_vertices_idxs, 'geometry'].plot(ax=ax, color='r')

# Plot all the intersections that were visited by the A* algorithm
# intersections_gpd.loc[visited_vertices, 'geometry'].plot()

plt.show()