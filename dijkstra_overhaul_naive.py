import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import ast
import time
import math
import sys

start_time = time.time()
roads_gpd = gpd.read_file('data/Ramsey_Roads_Prepped.geojson')
intersections_gpd = gpd.read_file('data/Ramsey_Intersections_Prepped.geojson')
end_time = time.time()
print(f'Time to read the geojson files: {end_time - start_time}')


source = 2     # an intersection
dest = 13001       # an intersection

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
        new_vertex.connections2 = convert_string_to_list(intersection.NeighboringIntersections)
        new_vertex = Vertex(index)
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
                    connections[vertices[intersection]] = Road(road, road_gdf.loc[road]['TimeToTravel'])
                    # connections.append(vertices[intersection])
            vertex.connections = connections
    
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
def Dijkstra(graph, dest):
    visited_vertices = []       # a list containing all the indices of all the vertices that have been visited so far
    visited_intersections = []  # a list containing all the Vertex objects that have been visited so far
    chosen_roads = []           # a list containing indicies for all the roads for the route
    global total_search_time
    total_search_time = 0
    while dest not in visited_vertices:
        min_vertex, search_time = find_min(graph)        # find the Vertex in the graph with the lowest travel cost
        total_search_time += search_time
        # update distances to the neighbor nodes
        neighbors = min_vertex.connections
        for neighbor, road in neighbors.items():
            new_cost = min_vertex.cost + road.travel_time
            if new_cost < neighbor.cost:        # check to see if the new route to reach this node is faster than the previous fastest route
                neighbor.cost = new_cost
                neighbor.prev = min_vertex
        
        # add the current node to the set of visited nodes
        visited_vertices.append(min_vertex.index)
        visited_intersections.append(min_vertex)
        chosen_roads.append(min_vertex.connections)

        # remove the current vertex from the set of unvisited nodes
        graph.remove(min_vertex)
    
    return visited_intersections, visited_vertices
    
# start_time = time.time()
# list_of_vertices = initialize_vertices(intersections_gpd, source)
# initialize_ints(intersections_gpd, roads_gpd, list_of_vertices)
# end_time = time.time()
# print(f'Time to initialize vertices: {end_time - start_time}\n')


'''
start = time.time()
adjacency_list = list(np.zeros(intersections_gpd.shape[0]))
# find the intersections that are one road segment away from the current intersection
i = 0
for intersection in intersections_gpd.itertuples():
    roads = intersection.Roads
    roads = convert_string_to_list(roads)
    neighboring_intersections = {}  # keys are indices; values are the time it takes to travel there (i.e. the weight of the road segment in between the two vertices)
    for road in roads:
        new_intersections = roads_gpd.loc[road]['Intersections']
        new_intersections = convert_string_to_list(new_intersections)
        for int in new_intersections:
            if (int != i) and (int not in neighboring_intersections):
                neighboring_intersections[int] = roads_gpd.loc[road]['TimeToTravel']
    adjacency_list[i] = neighboring_intersections
    i += 1

end = time.time()
print(adjacency_list)
print(f'Time to build the adjacency list: {end - start}')
'''


# Updates in late August:
# initialize graph
# The graph has this structure: {index: {neighboring_index: cost_to_reach, neighboring_index: cost_to_reach},
#                                index: {neighboring_index: cost_to_reach}, etc.}
# The indices in the graph are intersections, and the weights are the costs to reach neighboring intersections.
graph = {}
length = intersections_gpd.shape[0]
for i in range(length):
    neighboring_intersections = convert_string_to_list(intersections_gpd.loc[i]['NeighboringIntersections'])
    graph[i] = neighboring_intersections
#    print(neighboring_intersections)
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


sys.exit()






###################################################################################################
###################################################################################################


def dijkstra_new(graph, dest, adj_list):
    visited_vertices = []       # a list containing all the indices of all the vertices that have been visited so far
    visited_intersections = []  # a list containing all the Vertex objects that have been visited so far
    global total_search_time
    total_search_time = 0
    while dest not in visited_vertices:
        min_vertex, search_time = find_min(graph)        # find the Vertex in the graph with the lowest travel cost
        total_search_time += search_time
        # update distances to the neighbor nodes
        # neighbors = min_vertex.connections
        neighbors = adj_list[min_vertex.index]      # find the neighboring intersections

        for neighbor in neighbors:
            # each neighbor is an index into the graph
            new_cost = min_vertex.cost + 0


        for neighbor, road in neighbors.items():
            new_cost = min_vertex.cost + road.travel_time
            if new_cost < neighbor.cost:        # check to see if the new route to reach this node is faster than the previous fastest route
                neighbor.cost = new_cost
                neighbor.prev = min_vertex
        
        # add the current node to the set of visited nodes
        visited_vertices.append(min_vertex.index)
        visited_intersections.append(min_vertex)

        # remove the current vertex from the set of unvisited nodes
        graph.remove(min_vertex)
    
    return visited_intersections, visited_vertices

sys.exit()

start_time = time.time()
visited_intersections, visited_vertices = Dijkstra(list_of_vertices, dest)
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
print(f'Time spent finding minimum value in list: {total_search_time}')


fig, ax = plt.subplots()

roads_gpd.plot(ax=ax)

# Plot the intersections that were chosen to create the fastest route
intersections_gpd.loc[chosen_vertices_idxs, 'geometry'].plot(ax=ax, color='r')

# Plot all the intersections that were visited by Dijkstra's algorithm
# intersections_gpd.loc[visited_vertices, 'geometry'].plot()

print(f'Total number of intersections: {intersections_gpd.shape[0]}')
print(f'Number of intersections visited by Dijkstra\'s algorithm: {len(visited_vertices)}\n')

plt.show()