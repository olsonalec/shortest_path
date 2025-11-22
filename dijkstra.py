import geopandas as gpd
import matplotlib.pyplot as plt
import ast
import time
import math


start_idx = 20000
end_idx = 0

print('Loading data...')
roads_gpd = gpd.read_file('data/Hennepin_Roads_Prepped.geojson')
intersections_gpd = gpd.read_file('data/Hennepin_Intersections_Prepped.geojson')


##########################################################################
# Updates November 2025

class Vertex:
    def __init__(self, index, neighbors, coordinates):
        self.index = index
        self.neighbors = neighbors
        self.coordinates = coordinates

class Node:
    def __init__(self, neighbors, coordinates):
        self.visited = False
        self.cost = float("inf")
        self.prev = None
        self.neighbors = neighbors      # {<vertex_of_neighboring_intersection>: <cost_to_reach_that_intersection>}
        self.coordinates = coordinates  # [<x_coordinate>, <y_coordinate>]

# Each coordinate is in UTM, so we can calculate Euclidean distance directly without having to convert to a different coordiante system.
# (Ex. if the coordinates were in degrees, we'd have to convert them first.)
# This is the heuristic function for A*
def calculate_distance(start_x, start_y, end_x, end_y):
    return ((((start_x - end_x) **2) + (start_y - end_y) ** 2) ** 0.5)

def convert_meters_to_miles(meters):
    return meters / 1609.344

def heuristic(start_x, start_y, end_x, end_y, max_speed):
    distance = convert_meters_to_miles(calculate_distance(start_x, start_y, end_x, end_y))
    return distance / max_speed

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


##### NEED TO SET UP vertex AND graph PROPERLY BEFORE THIS FUNCTION WILL WORK
# vertex is a Vertex object
def expand_vertex(vertex, graph):
    expanded_states = []
    # vertex.neighbors is a dictionary
    neighbors = list(vertex.neighbors.keys())
    for neighbor in neighbors:
        expanded_states.append(Vertex(neighbor, graph[neighbor].neighbors, graph[neighbor].coordinates))
     
    return expanded_states

    # extract neighbor intersection indices from the neighbors attribute, which is a dictionary
    vertex_neighbors = list(vertex.neighbors.keys())
    for i in range(len(vertex_neighbors)):
        neighbor = vertex_neighbors[i]
        expanded_states.append(Vertex(neighbor, graph[neighbor], [int_gdf.iloc[neighbor].POINT_X, int_gdf.iloc[neighbor].POINT_Y]))

    return expanded_states

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


# The following functions implement standard heap operations.

########## ELEMENTS IN THE HEAP SHOULD BE Vertex OBJECTS, NOT Node OBJECTS ##########

# get the index of the parent node given a starting index i
def parent(i):
    return (i - 1) // 2

# get the index of the left child node given a starting index i
def left_child(i):
    return i * 2 + 1

# get the index of the right child node given a starting index i
def right_child(i):
    return i * 2 + 2

# push an item onto the heap
# runs in O(log2(n)) time because Python's append() operator runs in O(1) amortized time and heapify() runs in O(log2(n)) time
#   item is a State object
#   heap is the minimum binary heap, represented as a list
#   lookup is a dictionary that is used to map State objects to their costs
#       the keys are strings representing room numbers, and the values are Node objects
def push(item, heap, lookup):
    heap.append(item)
    heap = heapify(-1, heap, lookup)
    
    return heap

# swap two items in the heap
# runs in O(1) time
#   i and j are indices
#   heap is the minimum binary heap, represented as a list
def swap(i, j, heap):
    temp = heap[i]
    heap[i] = heap[j]
    heap[j] = temp

    return heap

# move a node up the heap until it is in place
# runs in O(log2(n)) time because the heap is a binary tree
#   i is an index
#   heap is the minimum binary heap, implemented as a list
#   lookup is a dictionary that is used to map State objects to their costs
#       the keys are strings representing room numbers, and the values are Node objects
def heapify(i, heap, lookup):
    parent_idx = parent(i)

    # lookup[heap[x].index] is a Node that corresponds to a State
    # lookup[heap[x].index].cost is the f value for that Node
    if (parent_idx >= 0) and (lookup[heap[i].index].cost < lookup[heap[parent_idx].index].cost):
        heap = swap(i, parent_idx, heap)
        heap = heapify(parent_idx, heap, lookup)

    return heap

# move a node down the heap until it is in place
# runs in O(log2(n)) time because the heap is a binary tree
#   i is an index
#   heap is the minimum binary heap, implemented as a list
#   lookup is a dictionary that is used to map State objects to their costs
#       the keys are strings representing room numbers, and the values are Node objects
def min_heapify(i, heap, lookup):
    l = left_child(i)
    r = right_child(i)
    smallest = -1

    # lookup[heap[x].index] is a Node that corresponds to a State
    # lookup[heap[x].index].cost is the f value for that Node

    if (l < len(heap)) and (lookup[heap[l].index].cost < lookup[heap[i].index].cost):
        smallest = l
    else:
        smallest = i

    if (r < len(heap)) and (lookup[heap[r].index].cost < lookup[heap[smallest].index].cost):
        smallest = r

    if smallest != i:
        heap = swap(i, smallest, heap)
        heap = min_heapify(smallest, heap, lookup)
    
    return heap

# remove the minimum element from the heap (it will be at the front because this is a minimum binary heap)
# runs in O(log2(n)) time because Python's pop() function runs in O(1) time when removing the last element of a list, and min_heapify() runs in O(log2(n)) time
#   heap is the minimum binary heap, implemented as a list
#   lookup is a dictionary that is used to map State objects to their costs
#       the keys are strings representing room numbers, and the values are Node objects
def extract_min(heap, lookup):
    if len(heap) == 0:
        return None
    
    min = heap[0]

    if len(heap) > 1:
        '''
        If there are multiple items in the heap:
            - copy the last element to the front
            - remove the element at the end
            - call min_heapify to move this element to the correct place in the heap
        '''
        heap[0] = heap[-1]
        heap.pop(-1)
        heap = min_heapify(0, heap, lookup)

    return min, heap


def a_star(node_dict, starting_intersection, goal_intersection):
    # This dictionary contains Node objects, each of which keeps track of a corresponding State object.
    # The keys are strings representing room numbers; values are Node objects
    # Node objects keep track of their State's cost, f value, whether or not A* has visited the State, and, if it has been visited, the previous State in the shortest path
    # node_dict = {}
    
    # initialize the maze by creating a Node object for each state (room) in the maze
    # for room in maze.keys():
        # node_dict[room] = Node()
    # retrieve the room number of the starting state
    start_position = starting_intersection.index

    # retrieve the coordinates of the starting and goal intersections
    start_x = starting_intersection.coordinates[0]
    start_y = starting_intersection.coordinates[1]
    goal_x = goal_intersection.coordinates[0]
    goal_y = goal_intersection.coordinates[1]

    # initialize the Node object associated with the starting state
    node_dict[start_position].visited = True        # by default, the starting state has already been visited
    node_dict[start_position].cost = 0              # by default, the starting state has cost 0
    
    # This minimum binary heap will contain states as they are generated. Initially, only the starting state is in the heap.
    heap = [starting_intersection]                         # each element in the heap is a State object
    
    # Keep track of how many intersections have been explored. This information will be used to evaluate the runtime of this algorithm.
    intersections_visited = 1

    while len(heap) > 0:
        # Choose the state with lowest estimated cost.
        state, heap = extract_min(heap, node_dict)      # extract_min() chooses the state with lowest f cost [f(x) = g(x) + h(x)]

        # Terminate A* if the algorithm has reached the goal state.
        if state.index == goal_intersection.index:
            print(f'Total number of intersections visited: {intersections_visited}.')
            return node_dict

        # Generate the states that can be reached from the given state.
        expanded_states = expand_vertex(state, node_dict)         # expanded_states is a list of State objects

        # Each neighbor is a State object.
        for neighbor in expanded_states:
            intersections_visited += 1
            # retrieve the room number of the neighbor
            index = neighbor.index

            # retrieve coordinates
            # neighbor_x = neighbor.coordinates[0]
            # neighbor_y = neighbor.coordinates[1]
            
            # get cost of traveling from the current node to this neighbor node
            neighbor_cost = node_dict[state.index].neighbors[index]

            # h_neighbor = heuristic(neighbor_x, neighbor_y, goal_x, goal_y, 70)           # the estimated cost to reach the goal state from the neighboring state
            g_neighbor = node_dict[state.index].cost + neighbor_cost                 # The cost to reach the neighboring state from the starting state. 
            #                                                                 # node_dict[state.index].cost is the cost to reach the previous node. Therefore, the cost to reach the neighboring node is node_dict[state.index].cost +
            # f_neighbor = g_neighbor + h_neighbor

            # Check to make sure that this state hasn't already been visited.
            # If it hasn't been visited, then check to see if the new cost to reach the state is less than the current best-known cost to reach it.
            if (node_dict[index].visited == False) and (node_dict[index].cost > g_neighbor):
                # update the true cost and f(x) to reach this node
                # remember, State objects are immutable, so we must update the corresponding Node object, which is stored in the node_dict dictionary
                node_dict[index].cost = g_neighbor

                # add the neighboring state to the heap and update its prev attribute to keep track of the previous state in the path
                heap = push(neighbor, heap, node_dict)
                node_dict[index].prev = state.index       # state.index is the room number of the state that Dijkstra is currently expanding

        # after expanding all the states of this node, mark the current state as visited
        node_dict[state.index].visited = True
        intersections_visited += 1
    # if we reach this point, A* has failed to find a path from the starting state to the goal state
    return None


# each key is an index into the intersections GeoDataFrame
# each value is a Node object
# IGNORE THIS COMMENT each value is a dictionary in the form {<vertex_of_neighboring_intersection>: <cost_to_reach_that_intersection>}
# IGNORE THIS COMMENT each value is a list of strings that are themselves indices into the intersections GeoDataFrame (i.e. intersections that are one road segment away from the key intersection)
print('Constructing road network...')
road_network = {}
for intersection in intersections_gpd.itertuples():
    index = intersection.Index
    neighbors = convert_string_to_list(intersection.NeighboringIntersections)   # neighbors is now a list with one element, which is a dictionary
        # the dictionary is structured {<vertex_of_neighboring_intersection>: <cost_to_reach_that_intersection>}
    neighbors = neighbors[0]    # get the dictionary out of the list
    # road_network[index] = neighbors
    # neighbor_indices = list(neighbors.keys())
    # road_network[index] = neighbor_indices
    road_network[index] = Node(neighbors, [intersections_gpd.iloc[index].geometry.x, intersections_gpd.iloc[index].geometry.y])

# initialize the start and stop intersections
start_intersection = Vertex(start_idx, road_network[start_idx].neighbors, [intersections_gpd.iloc[start_idx].geometry.x, intersections_gpd.iloc[start_idx].geometry.y])
destination_intersection = Vertex(end_idx, road_network[end_idx].neighbors, [intersections_gpd.iloc[end_idx].geometry.x, intersections_gpd.iloc[end_idx].geometry.y])

print('Running A*...')
start_time = time.time()
nodes = a_star(road_network, start_intersection, destination_intersection)
end_time = time.time()

if nodes == None:
    print('Failed to find a path from the starting state to the goal state.')
else:
    # using the prev attribute of each Node, construct the shortest path from the goal state to the beginning state, then reverse it to find the true path
    path = []
    path.append(destination_intersection.index)
    prev = nodes[destination_intersection.index].prev

    # find time to travel this path
    total_travel_time = convert_sec_to_min(nodes[destination_intersection.index].cost)

    while prev != None:
        path.append(prev)
        prev = nodes[prev].prev

    print('Path: { ', end='')
    for i in range(len(path) - 1, -1, -1):
        print(f'{path[i]}', end=' ')
    print('}')

print(f'Number of intersections selected: {len(path)}.')
print(f'The time it will take to travel this route is approximately {add_buffer_time(total_travel_time)} minutes.')
print(f'Time taken by A* to find the shortest path: {end_time - start_time} seconds.')


fig, ax = plt.subplots()

roads_gpd.plot(ax=ax)

# Plot the intersections that were chosen to create the fastest route
intersections_gpd.loc[path, 'geometry'].plot(ax=ax, color='r')

# Plot all the intersections that were visited by the A* algorithm
# intersections_gpd.loc[visited_vertices, 'geometry'].plot()

plt.show()