import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import ast
import time
import math
import sys

start_time = time.time()
roads_gdf = gpd.read_file('data/Plymouth_Roads_Prepped.geojson')
intersections_gdf = gpd.read_file('data/Plymouth_Intersections_Prepped2.geojson')
end_time = time.time()
print(f'Time to read the geojson files: {end_time - start_time}')


source = 0     # an intersection
dest = 50       # an intersection

class Vertex:
    def __init__(self, index):
        self.cost = float("inf")        # the total cost of the path to reach this vertex
        self.prev = None    # the previous Vertex in the path
        self.index = index      # the index into the intersections_gdf where this vertex is stored

    
'''
A function to find the index of the minimum element in a list of unsorted objects.

Parameter:
    list_of_vertices (list) - a list of Vertex objects

Return Value:
    min_vertex (Vertex) - the Vertex object that has the lowest travel cost in the list
'''
def find_min(list_of_vertices, visited_vertices):
    start = time.time()
    min_value = float("inf")
    min_vertex = None
    n = len(list_of_vertices)
    for i in range(n):
        if (list_of_vertices[i].cost < min_value) and (list_of_vertices[i].index not in visited_vertices):
            min_value = list_of_vertices[i].cost
            min_vertex = list_of_vertices[i]
    stop = time.time()
    search_time = stop - start
    return min_vertex, search_time

'''
The 'Intersections' and 'Roads' attributes in the Roads and Intersections Geodataframes, respectively, are represented as a nested list.
If they were a one-dimensional list, GeoPandas would think that they represent geometry, which they don't.
However, when reading these attributes from the dataframe, Python interprets them as strings.
This is an example: '[[180, 240, 360]]'
This function converts this string representation of a nested list into a 1-dimensional Python list.
The example output would be [180, 240, 360], where each element is an integer.

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

def add_buffer_time(minutes):
    return math.ceil(minutes * 1.1)

'''
A function to find the shortest path between two intersections.

Parameters:
    graph - a list of Vertex objects that contains all the unvisited vertices in the graph
    dest - the index into the intersections Geodataframe of the destination vertex

Return Value:
    visited_intersections (list) - a list of all the Vertex objects that were visited by the algorithm
                                    This list contains Vertex objects that are on the shortest path as well as Vertex objects that are not in the shortest path
                                    The prev attribute of each Vertex object will be used to determine the shortest path. After this function returns, simply find
                                        the ending vertex, and use the prev attribute to work your way backwards until you reach the starting vertex. This process
                                        returns the shortest path.
'''
def Dijkstra(graph, array, dest):
    visited_vertices = []       # a list containing all the indices of all the vertices that have been visited so far
    visited_intersections = []  # a list containing all the Vertex objects that have been visited so far
    global total_search_time
    total_search_time = 0
    additional_time = 0
    while dest not in visited_vertices:
        min_vertex, search_time = find_min(array, visited_vertices)        # find the Vertex in the graph with the lowest travel cost
        total_search_time += search_time
        start2 = time.time()
        # update distances to the neighbor nodes
        min_vertex_cost = min_vertex.cost
        neighbors = graph[min_vertex.index]     # neighbors is a list of tuples: [(index, weight), (index, weight), ...]
        for neighbor in neighbors:
            new_cost = min_vertex_cost + neighbor[1]
            if new_cost < array[neighbor[0]].cost:
                array[neighbor[0]].cost = new_cost
                array[neighbor[0]].prev = min_vertex
        
        # add the current node to the set of visited nodes
        visited_vertices.append(min_vertex.index)
        visited_intersections.append(min_vertex)
        # chosen_roads.append(min_vertex.connections)
        end2 = time.time()
        additional_time += (end2 - start2)
        # remove the current vertex from the set of unvisited nodes
        # array.remove(min_vertex)
    
    return visited_intersections, visited_vertices, additional_time
    

'''
Function to initialize the graph.

Parameters:
    graph (dictionary) - A representation of the network of intersections; each key is an index into the intersections GeoDataFrame.
                            Each value is a list of tuples. The first element of each tuple is the index of a neighboring intersection.
                            The second element of each tuple is the cost to travel from the key to the neighboring intersection.
    int_gdf (GeoDataFrame) - A GeoDataFrame that stores information for each intersection.

Return value:
    graph (dictionary) - The initialized graph
'''
def initialize_graph(graph, int_gdf):
    length = int_gdf.shape[0]       # Find the total number of intersections
    for i in range(length):
        neighbors = convert_string_to_list(intersections_gdf.loc[i]['NeighboringIntersections'])
        graph[i] = neighbors

    return graph


'''
Function to intialize the array that stores intersections and the current cost to get there from the starting node.

Parameters:
    int_gdf (GeoDataFrame) - A GeoDataFrame that stores information for each intersection.
    start (int) - The index into the GeoDataFrame where the starting intersection is located.

Return value:
    vertex_array (list) - The initialized list. Each element is a Vertex object that represents one intersection.
                            Each Vertex has a starting cost of infinity, except for the starting Vertex, which has a cost of zero.
                            Each Vertex object also stores the index of the corresponding intersection in the GeoDataFrame.
                            Initially, this index is equivalent to the Vertex's index in the list.
'''
def initialize_array(int_gdf, start):
    n = int_gdf.shape[0]            # Find the total number of intersections
    vertex_array = list(np.zeros(n))
    for i in range(n):
        vertex_array[i] = Vertex(i)     # this initializes each Vertex object to have a cost of infinity
    vertex_array[start].cost = 0        # Set the cost of the source Vertex to be zero

    return vertex_array




start_time = time.time()
graph = initialize_graph({}, intersections_gdf)
vertex_array = initialize_array(intersections_gdf, source)
visited_intersections, visited_vertices, add_time = Dijkstra(graph, vertex_array, dest)
end_time = time.time()
print(f'Additional time: {add_time}')
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
print(f'Time taken to run Dijkstra\'s Algorithm: {end_time - start_time} seconds\n')
print(f'Time spent finding minimum value in list: {total_search_time}')


fig, ax = plt.subplots()

roads_gdf.plot(ax=ax)

# Plot the intersections that were chosen to create the fastest route
intersections_gdf.loc[chosen_vertices_idxs, 'geometry'].plot(ax=ax, color='r')

# Plot all the intersections that were visited by Dijkstra's algorithm
# intersections_gdf.loc[visited_vertices, 'geometry'].plot()

print(f'Total number of intersections: {intersections_gdf.shape[0]}')
print(f'Number of intersections visited by Dijkstra\'s algorithm: {len(visited_vertices)}\n')

plt.show()