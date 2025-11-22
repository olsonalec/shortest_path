'''
SCROLL DOWN PAST ALL THE CLASS AND FUNCTION DEFINITIONS IN ORDER TO REACH THE SECTION WHERE YOU CAN MANIPULATE THE STARTING AND ENDING ROOM NUMBERS.
'''

import csv

'''
Each State object represents a position in the maze where the A* algorithm is.
State objects, once generated, never change.
State objects can be generated multiple times. For example, if rooms X and Y both are adjacent to Z, then Z could be generated twice, once when X is visited by A* and once when Y is visited.
Each unique State object has a corresponding Node object, described below. The Node objects are mutable.
'''
class State:
    '''
    Parameters:
        - position is a two-character string representing a room in the maze
        - neighbors is a list of strings
            each string is a room that neighbors the position room and is accessible (i.e. there is not a wall between the two rooms)
    '''
    def __init__(self, position, neighbors):
        self.position = position
        self.neighbors = neighbors
        
'''
Each State might be generated multiple times. Therefore, each unique State object has a corresponding Node object that keeps track of the following properties:
    visited (Boolean) - represents if the A* algorithm has expanded this State
    cost (int) - the cost to reach this state from the starting state; this is equal to g(x)
    f (float) - f(x) = g(x) + h(x); the cost to reach this state from the starting state (g) + the estimated cost to reach the goal state (h)
    prev (string) - the room number of the previous State in the path (used to construct the shortest path once A* is done running)
'''
class Node:
    def __init__(self):
        self.visited = False
        self.cost = float("inf")
        self.f = float("inf")
        self.prev = None

'''
Description of Function:
    This function calculates the Euclidean distance between the input state and the goal state.
    The Manhatten distance between any two nodes in the maze is always less than or equal to the Euclidean distance between those nodes.
    Because we navigate the maze using Manhatten distance, this function never overestimates the number of states on a path between the input state and the goal state.
    Therefore, for every node n, with an intermediate node n' which is somewhere on the path between n and the goal node, h(n) <= (true cost between n and n') + h(n')
    This means that this heuristic function is consistent, which means that it is also admissible.  
    As shown in class, this means that the A* algorithm will always generate an optimal solution.
Parameter:
    current_coordinates - a string of length two representing the starting room, ex. '01'
    goal_coordinates - a string of length two representing the goal room, ex. '23'
Return Value:
    distance - a float that estimates the cost to travel to the goal state from the current state
'''
def heuristic(current_coordinates, goal_coordinates):
    start_x, start_y = int(current_coordinates[0]), int(current_coordinates[1])
    end_x, end_y = int(goal_coordinates[0]), int(goal_coordinates[1])

    # calculate the Euclidean distance between these two rooms
    distance = ((start_x - end_x) ** 2 + (start_y - end_y) ** 2) ** 0.5

    return distance

'''
Description of Function:
    For each accessible neighboring room of the input state (i.e. there is not a wall between the two rooms), this function creates a new State object.
    Each of these new State objects are appended to the expanded_states list, which is then returned
Parameters:
    s - a State object
    maze - a dictionary where each key is a room and each value is a list of rooms that are accessible by the key room
Return Value:
    expanded_states is a list of State objects that are reachable from the input state
'''
def expand_state(s, maze):
    expanded_states = []
    for i in range(len(s.neighbors)):
        neighbor = s.neighbors[i]               # neighbor is a string representing a room number in the maze
        expanded_states.append(State(neighbor, maze[neighbor]))         # maze[neighbor] is a list of room numbers, representing the rooms that are adjacent to neighbor
    
    return expanded_states

'''
A* will use a minimum binary heap to manage states.
This heap will be implemented as a list.
Each item in the heap is a State object.
The following functions implement standard heap operations.
'''
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

    # lookup[heap[x].position] is a Node that corresponds to a State
    # lookup[heap[x].position].f is the f value for that Node
    if (parent_idx >= 0) and (lookup[heap[i].position].f < lookup[heap[parent_idx].position].f):
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

    # lookup[heap[x].position] is a Node that corresponds to a State
    # lookup[heap[x].position].f is the f value for that Node

    if (l < len(heap)) and (lookup[heap[l].position].f < lookup[heap[i].position].f):
        smallest = l
    else:
        smallest = i

    if (r < len(heap)) and (lookup[heap[r].position].f < lookup[heap[smallest].position].f):
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

'''
The main A* function that finds the shortest route between two rooms in the maze.
Parameters:
    - maze: a dictionary that represents the maze; each key is a room number, and each value is a list of rooms that are accessible from the key room
    - starting_state: a State object representing the starting room
    - goal_state: a State object representing the goal room
Return Value:
    - node_dict: A dictionary that contains information about each state in the maze.
                    The keys are strings that represent room numbers.
                    The values are Node objects that will be used to construct the shortest path between the starting and ending states once A* returns.
'''
def a_star(maze, starting_state, goal_state):
    # This dictionary contains Node objects, each of which keeps track of a corresponding State object.
    # The keys are strings representing room numbers; values are Node objects
    # Node objects keep track of their State's cost, f value, whether or not A* has visited the State, and, if it has been visited, the previous State in the shortest path
    node_dict = {}
    
    # initialize the maze by creating a Node object for each state (room) in the maze
    for room in maze.keys():
        node_dict[room] = Node()

    # retrieve the room number of the starting state
    start_position = starting_state.position

    # initialize the Node object associated with the starting state
    node_dict[start_position].visited = True        # by default, the starting state has already been visited
    node_dict[start_position].cost = 0              # by default, the starting state has cost 0
    node_dict[start_position].f = heuristic(start_position, goal_state.position)    # calculate the estimated cost from the starting state to the goal state
    
    # This minimum binary heap will contain states as they are generated. Initially, only the starting state is in the heap.
    heap = [starting_state]                         # each element in the heap is a State object
    
    while len(heap) > 0:
        # Choose the state with lowest estimated cost.
        state, heap = extract_min(heap, node_dict)      # extract_min() chooses the state with lowest f cost [f(x) = g(x) + h(x)]

        # Terminate A* if the algorithm has reached the goal state.
        if state.position == goal_state.position:
            return node_dict
        
        # Generate the states that can be reached from the given state.
        expanded_states = expand_state(state, maze)         # expanded_states is a list of State objects
        
        # Each neighbor is a State object.
        for neighbor in expanded_states:
            # retrieve the room number of the neighbor
            position = neighbor.position
            
            h_neighbor = heuristic(position, goal_state.position)           # the estimated cost to reach the goal state from the neighboring state
            g_neighbor = node_dict[state.position].cost + 1                 # The cost to reach the neighboring state from the starting state. The cost to travel between two adjacent rooms (states) is always 1.
                                                                            # node_dict[state.position].cost is the cost to reach the previous node. Therefore, the cost to reach the neighboring node is node_dict[state.position].cost + 1
            f_neighbor = g_neighbor + h_neighbor

            # Check to make sure that this state hasn't already been visited.
            # If it hasn't been visited, then check to see if the new cost to reach the state is less than the current best-known cost to reach it.
            if (node_dict[position].visited == False) and (node_dict[position].f > f_neighbor):
                # update the true cost and f(x) to reach this node
                # remember, State objects are immutable, so we must update the corresponding Node object, which is stored in the node_dict dictionary
                node_dict[position].cost = g_neighbor
                node_dict[position].f = f_neighbor

                # add the neighboring state to the heap and update its prev attribute to keep track of the previous state in the path
                heap = push(neighbor, heap, node_dict)
                node_dict[position].prev = state.position       # state.position is the room number of the state that A* is currently expanding

        # after expanding all the states of this node, mark the current state as visited
        node_dict[state.position].visited = True

    # if we reach this point, A* has failed to find a path from the starting state to the goal state
    return None

'''
maze is a dictionary to store the maze
each key is a two-character string representing the room number
each value is a list of two-character strings representing the rooms that are reachable from the key
'''
maze = {}

# read the maze file
with open('maze.csv', 'r') as fp:
    csvreader = csv.reader(fp, delimiter=',', quotechar='\"', skipinitialspace=True)
    
    # skip the header
    next(csvreader)

    # create the maze data structure
    for row in csvreader:
        room = row[0]
        accessible_neighbors = row[1].split(',')
        maze[room] = accessible_neighbors


########## EDIT THE STARTING AND ENDING ROOMS HERE ##########

# initialize the start and goal states
start_state = State('03', maze['03'])
goal_state = State('30', maze['30'])

################# DO NOT EDIT ANYWHERE ELSE #################


nodes = a_star(maze, start_state, goal_state)

if nodes == None:
    print('Failed to find a path from the starting state to the goal state.')
else:
    # using the prev attribute of each Node, construct the shortest path from the goal state to the beginning state, then reverse it to find the true path
    path = []
    path.append(goal_state.position)
    prev = nodes[goal_state.position].prev

    while prev != None:
        path.append(prev)
        prev = nodes[prev].prev

    print('Path: { ', end='')
    for i in range(len(path) - 1, -1, -1):
        print(f'{path[i]}', end=' ')
    print('}')