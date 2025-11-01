import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import ast
import time
import math


start_time = time.time()
roads_gpd = gpd.read_file('data/Ramsey_Roads_Prepped.geojson')
intersections_gpd = gpd.read_file('data/Ramsey_Intersections_Prepped.geojson')
end_time = time.time()
print(f'Time to read the geojson files: {end_time - start_time}')


source = 2     # an intersection
dest = 8000       # an intersection

class Road:
    def __init__(self, index, cost):
        self.index = index      # the index into the roads_gpd where this road is stored
        self.travel_time = cost        # the cost to travel this road segment

class Vertex:
    def __init__(self, index):
        self.cost = float("inf")        # the total cost of the path to reach this vertex
        self.prev = None    # the previous Vertex in the path
        self.index = index      # the index into the intersections_gpd where this vertex is stored
        self.visited = False    # a boolean representing if Dijkstra's algorithm has visited this node yet
        self.neighbors = [] # a list of pointers to Vertex objects that can be reached by this Vertex - i.e. intersections that are one road segment away

# class Node:
#     def __init__(self):
#         self.visited = False
#         self.prev = None
#         self.cost = float("inf")

    
'''
A function to find the index of the minimum element in a list of unsorted objects.

Parameter:
    list_of_vertices (list) - a list of Vertex objects

Return Value:
    idx (int) - the index into the array where the Vertex object with the lowest travel cost is stored
'''
def find_min(list_of_vertices):
    min_value = float("inf")
    idx = float("inf")
    n = len(list_of_vertices)
    for i in range(n):
        if (list_of_vertices[i].cost < min_value) and not (list_of_vertices[i].visited):
            min_value = list_of_vertices[i].cost
            idx = i
    return idx

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
    # min = int(seconds // 60)
    # sec = int(seconds % 60)
    # return min, sec

def add_buffer_time(minutes):
    return math.ceil(minutes * 1.1)


# Updates in late August:
# initialize graph
# The graph has this structure: {index: {neighboring_index: cost_to_reach, neighboring_index: cost_to_reach},
#                                index: {neighboring_index: cost_to_reach}, etc.}
# Each index is an index into the intersections GeoDataframe
# The indices in the graph are intersections, and the weights are the costs to reach neighboring intersections.
graph = {}
length = intersections_gpd.shape[0]
for i in range(length):
    neighboring_intersections = convert_string_to_list(intersections_gpd.loc[i]['NeighboringIntersections'])
    graph[i] = neighboring_intersections
# print(graph)

# Create a array. Each element is a Vertex object.
# Each node also has pointers to its neighboring nodes.
my_array = list(np.zeros(length))

# Create a Vertex object for each intersection. Right now, the neighbor attribute of each Vertex object is empty
for i in range(length):
    my_array[i] = Vertex(i)

# Set the cost for the starting intersection to be zero
my_array[source].cost = 0

# Update the neighbors attribute for each Vertex in the heap
for i in range(length):
    neighboring_intersections = convert_string_to_list(intersections_gpd.loc[i]['NeighboringIntersections'])
    # neighboring_intersections is a list of tuples: [(intersection, cost), (intersection, cost), (intersection, cost)]. We need to extract the intersections from this list.
    try:
        neighboring_intersections = list(neighboring_intersections[0].keys())
    except:     # for some reason, there are some intersections that don't have any neighbors
        neighboring_intersections = []
    for intersection in neighboring_intersections:
        my_array[i].neighbors.append(my_array[intersection])

def dijkstra_update(graph, dest, vertex_array):
    visited_vertices = []       # a list containing all the indices of all the vertices that have been visited so far
    visited_intersections = []  # a list containing all the Vertex objects that have been visited so far

    while dest not in visited_vertices:
        # extract the vertex with lowest cost from the heap
        min_vertex_idx = find_min(vertex_array)

        # add the vertex's index to the list of visited vertices
        # visited_vertices.append(min_vertex_idx)

        # update distances to the neighboring vertices
        for neighbor in vertex_array[min_vertex_idx].neighbors:       # each neighbor is a Vertex object
            if neighbor.visited == False:
                # new_cost = current cost + cost to reach the neighboring Vertex
                new_cost = vertex_array[min_vertex_idx].cost + graph[min_vertex_idx][0][neighbor.index]
                if new_cost < neighbor.cost:        # check to see if the new route to reach this Vertex is faster than the previous fastest route
                    neighbor.cost = new_cost
                    neighbor.prev = vertex_array[min_vertex_idx]
        
        # add the vertex's index to the list of visited vertices
        visited_vertices.append(min_vertex_idx)
        visited_intersections.append(vertex_array[min_vertex_idx])

        vertex_array[min_vertex_idx].visited = True

        # remove the current vertex from the heap
        # vertex_array.remove(min_vertex_idx)
    
    return visited_intersections, visited_vertices



start_time = time.time()
visited_intersections, visited_vertices = dijkstra_update(graph, dest, my_array)
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
print(f'Time taken to run Dijkstra\'s Algorithm: {end_time - start_time} seconds\n')
# print(f'Time spent finding minimum value in list: {total_search_time}')

fig, ax = plt.subplots()

roads_gpd.plot(ax=ax)

# Plot the intersections that were chosen to create the fastest route
intersections_gpd.loc[chosen_vertices_idxs, 'geometry'].plot(ax=ax, color='r')

# Plot all the intersections that were visited by Dijkstra's algorithm
# intersections_gpd.loc[visited_vertices, 'geometry'].plot()

print(f'Total number of intersections: {intersections_gpd.shape[0]}')
print(f'Number of intersections visited by Dijkstra\'s algorithm: {len(visited_vertices)}\n')

plt.show()