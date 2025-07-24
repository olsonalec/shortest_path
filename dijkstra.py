import geopandas as gpd
import matplotlib.pyplot as plt
import ast
import time
import math
import sys

start_time = time.time()
roads_gpd = gpd.read_file('data/Plymouth_Roads_Prepped.geojson')
intersections_gpd = gpd.read_file('data/Plymouth_Intersections_Prepped.geojson')
end_time = time.time()
print(f'Time to read the geojson files: {end_time - start_time}')


source = 0     # an intersection
dest = 2000       # an intersection

class Road:
    def __init__(self, index, cost):
        self.index = index      # the index into the roads_gpd where this road is stored
        self.travel_time = cost        # the cost to travel this road segment

class Vertex:
    def __init__(self, index):
        self.cost = float("inf")        # the total cost of the path to reach this vertex
        self.prev = None    # the previous Vertex in the path
        self.index = index      # the index into the intersections_gpd where this vertex is stored
        self.connections = {}   # a dict of Vertex objects that can be reached by this Vertex - i.e. intersections that are one road segment away; keys are Vertex objects and values are Road objects


'''
Binary Min-Heap Functions (implemented with a list)
'''

class MinHeap:
    def __init__(self, heap, size):
        self.heap = heap        # a list representing a mininum priority binary heap; each element is a Vertex object
        self.size = size        # the number of elements in the heap

    # i is an index into the list
    def parent(self, i):
        return (i - 1) // 2

    def left_child(self, i):
        return i * 2 + 1

    def right_child(self, i):
        return i * 2 + 2

    def swap(self, i, j):
        temp = self.heap[i]
        self.heap[i] = self.heap[j]
        self.heap[j] = temp

    # move a node up the tree until it is in place
    def heapify(self, i):
        parent = self.parent(i)
        if (parent >= 0) and (self.heap[i].cost < self.heap[parent].cost):
            self.swap(i, parent)

    # move a node down the tree until it is in place
    def min_heapify(self, i):
        l = self.left_child(i)
        r = self.right_child(i)
        smallest = -1
        if (l < self.size) and (self.heap[l].cost < self.heap[i].cost):
            smallest = l
        else:
            smallest = i
        if (r < self.size) and (self.heap[r].cost < self.heap[smallest].cost):
            smallest = r
        if smallest != i:
            print(self.heap[i].cost, self.heap[smallest].cost)
            self.swap(i, smallest)
            print(self.heap[i].cost, self.heap[smallest].cost)
            self.min_heapify(smallest)

    def extract_min(self):
        min = self.heap[0]
        self.heap[0] = self.heap[self.size]
        self.min_heapify(0)
        return min

    def decrease_key(self, index, key):
        if key > self.heap[index].cost:
            pass        # in this case, the new key is larger than the old key
        else:
            self.heap[index].cost = key
            # move heap[i] up the heap until it is in place
            while (i > 0) and (self.heap[self.parent(i)].cost > self.heap[i].cost):
                self.swap(self.heap, index, self.parent(i))
                i = self.parent(i)


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
    min_vertex (Vertex) - the Vertex object that has the lowest travel cost in the list
'''
def find_min(list_of_vertices):
    start = time.time()
    min_value = float("inf")
    min_vertex = None
    n = len(list_of_vertices)
    for i in range(n):
        if list_of_vertices[i].cost < min_value:
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
            # print(neighbor.cost)
            new_cost = min_vertex.cost + road.travel_time
            if new_cost < neighbor.cost:        # check to see if the new route to reach this node is faster than the previous fastest route
                neighbor.cost = new_cost
                neighbor.prev = min_vertex
            # print(neighbor.cost)
            # print()
        
        # add the current node to the set of visited nodes
        visited_vertices.append(min_vertex.index)
        visited_intersections.append(min_vertex)
        chosen_roads.append(min_vertex.connections)

        # remove the current vertex from the set of unvisited nodes
        graph.remove(min_vertex)
    
    return visited_intersections, visited_vertices
    
start_time = time.time()
list_of_vertices = initialize_vertices(intersections_gpd, source)
initialize_ints(intersections_gpd, roads_gpd, list_of_vertices)
with open('output.txt', 'w') as fp:
    for vertex in list_of_vertices:
        fp.write(f'{str(vertex)}')
        fp.write('\n\n')
    for vertex in list_of_vertices:
        fp.write(f'{str(vertex.connections)}\n')
    # fp.write([str(x) for x in list_of_vertices])
    # fp.write('\n\n\n')
    # fp.write([x.connections for x in list_of_vertices])

end_time = time.time()
print(f'Time to initialize vertices: {end_time - start_time}\n')


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