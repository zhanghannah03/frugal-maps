#!/usr/bin/env python3

from cgitb import small
import typing
from util import read_osm_data, great_circle_distance, to_local_kml_url

# NO ADDITIONAL IMPORTS!


ALLOWED_HIGHWAY_TYPES = {
    'motorway', 'trunk', 'primary', 'secondary', 'tertiary', 'unclassified',
    'residential', 'living_street', 'motorway_link', 'trunk_link',
    'primary_link', 'secondary_link', 'tertiary_link',
}


DEFAULT_SPEED_LIMIT_MPH = {
    'motorway': 60,
    'trunk': 45,
    'primary': 35,
    'secondary': 30,
    'residential': 25,
    'tertiary': 25,
    'unclassified': 25,
    'living_street': 10,
    'motorway_link': 30,
    'trunk_link': 30,
    'primary_link': 30,
    'secondary_link': 30,
    'tertiary_link': 25,
}


def build_internal_representation(nodes_filename, ways_filename):
    """
    Create any internal representation you you want for the specified map, by
    reading the data from the given filenames (using read_osm_data)
    """
    nodes_data = read_osm_data(nodes_filename)
    ways_data = read_osm_data(ways_filename)
    # Prune the invalid ways and invalid nodes.
    valid_nodes, valid_ways = pruning(ways_data)
    # Create the internal representation.
    nodes_dict = {}
    for node in nodes_data:
        if node['id'] in valid_nodes:
            nodes_dict[node['id']] = {
                'lat': node['lat'],
                'lon': node['lon'],
                'connect': set(),
                'ways': [], # List of tuples in the form of (way id, node id, node id)
            }
    # Determine which nodes are connected to which nodes.
    for way_id in valid_ways.keys():
        for index in range(len(valid_ways[way_id]['nodes'])):
            node = valid_ways[way_id]['nodes'][index]
            # If the way is oneway, then add the next node along the way.
            if valid_ways[way_id]['oneway']:
                if index != len(valid_ways[way_id]['nodes']) - 1:
                    nodes_dict[node]['connect'].add(valid_ways[way_id]['nodes'][index + 1])
                    nodes_dict[node]['ways'].append((way_id, valid_ways[way_id]['nodes'][index + 1]))
            # If the way is two-way, then add the previous and next node along the way.
            else:
                if index == 0:
                    if len(valid_ways[way_id]['nodes']) > 1:
                        nodes_dict[node]['connect'].add(valid_ways[way_id]['nodes'][1])
                        nodes_dict[node]['ways'].append((way_id, valid_ways[way_id]['nodes'][1]))
                elif index == len(valid_ways[way_id]['nodes']) - 1:
                    if len(valid_ways[way_id]['nodes']) > 1:
                        nodes_dict[node]['connect'].add(valid_ways[way_id]['nodes'][index - 1])
                        nodes_dict[node]['ways'].append((way_id, valid_ways[way_id]['nodes'][index - 1]))
                else:
                    nodes_dict[node]['connect'].add(valid_ways[way_id]['nodes'][index - 1])
                    nodes_dict[node]['connect'].add(valid_ways[way_id]['nodes'][index + 1])
                    nodes_dict[node]['ways'].append((way_id, valid_ways[way_id]['nodes'][index - 1], valid_ways[way_id]['nodes'][index + 1]))
    return [nodes_dict, valid_ways] 


def pruning(ways_data):
    """
    Given a dataset of the ways in the form of a generator, return a tuple of
    the following objects:
    1. A set of the ids of the valid nodes
    2. A set of the ids of the valid ways
    """
    valid_nodes = set()
    valid_ways = {}
    for way in ways_data:
        # A way is considered to be valid if and only if it has a highway tag 
        # and its tag is an allowed highway type.
        if 'highway' in way['tags'].keys():
            if way['tags']['highway'] in ALLOWED_HIGHWAY_TYPES:
                for node in way['nodes']:
                    valid_nodes.add(node)
                valid_ways[way['id']] = {
                    'nodes': way['nodes'],
                    'oneway': False,
                }
                if 'oneway' in way['tags'].keys():
                    if way['tags']['oneway'] == 'yes':
                        valid_ways[way['id']]['oneway'] = True
                if 'maxspeed_mph' in way['tags'].keys():
                    valid_ways[way['id']]['maxspeed'] = way['tags']['maxspeed_mph']
                else:
                    valid_ways[way['id']]['maxspeed'] = DEFAULT_SPEED_LIMIT_MPH[way['tags']['highway']]
    return [valid_nodes, valid_ways]


def find_short_path_nodes(map_rep, node1, node2):
    """
    Return the shortest path between the two nodes

    Parameters:
        map_rep: the result of calling build_internal_representation
        node1: node representing the start location
        node2: node representing the end location

    Returns:
        a list of node IDs representing the shortest path (in terms of
        distance) from node1 to node2
    """
    nodes_dict = map_rep[0]
    agenda = Queue()
    expanded = set()
    # If the start and end locations are the same, return a list of the start node.
    if node1 == node2:
        return [node1]
    # For all valid nodes, initialize its distance from source node as infinity.
    for node in nodes_dict.keys():
        nodes_dict[node]['distance'] = float('inf')
    # Initialize the source node's distance as 0.
    nodes_dict[node1]['distance'] = 0
    agenda.add(node1, 0)
    # Initialize the parent pointers dictionary.
    parent = {}
    parent[node1] = 'start'
    # Implement Dijkstra's algorithm.
    while agenda.length() != 0:
        current = agenda.pop_min()
        if current == node2:
            break 
        if current in expanded:
            continue 
        expanded.add(current)
        current_lat = nodes_dict[current]['lat']
        current_lon = nodes_dict[current]['lon']
        current_loc = (current_lat, current_lon)
        for neighbor in nodes_dict[current]['connect']:
            # Add each neighbor and its distance value to the queue.
            agenda.add(neighbor, nodes_dict[neighbor]['distance'])
            neighbor_lat = nodes_dict[neighbor]['lat']
            neighbor_lon = nodes_dict[neighbor]['lon']
            neighbor_loc = (neighbor_lat, neighbor_lon)
            # Edge relaxation technique
            new_distance = nodes_dict[current]['distance'] + great_circle_distance(current_loc, neighbor_loc)
            if new_distance < nodes_dict[neighbor]['distance']:
                agenda.change_distance(neighbor, new_distance)
                nodes_dict[neighbor]['distance'] = new_distance
                parent[neighbor] = current
    # If the terminal node is unreachable from the source node, then return None.
    if nodes_dict[node2]['distance'] == float('inf'):
        return None 
    result = [node2]
    state = node2
    # Return the list of ids along the shortest (distance) path.
    while parent[state] != 'start':
        result.append(parent[state])
        state = parent[state]
    result.reverse()
    return result
    

class Queue():
    def __init__(self):
        self.H = {}
    def add(self, node_id, distance):
        self.H[node_id] = distance
    # Deletes the key-value pair with the minimum distance and returns the id.
    def pop_min(self):
        min_distance = float('inf')
        min_id = None 
        for id in self.H.keys():
            if self.H[id] < min_distance:
                min_distance = self.H[id] 
                min_id = id 
        del self.H[min_id]
        return min_id 
    def change_distance(self, node_id, new_distance):
        self.H[node_id] = new_distance 
    def length(self):
        return len(self.H.keys())


def find_short_path(map_rep, loc1, loc2):
    """
    Return the shortest path between the two locations

    Parameters:
        map_rep: the result of calling build_internal_representation
        loc1: tuple of 2 floats: (latitude, longitude), representing the start
              location
        loc2: tuple of 2 floats: (latitude, longitude), representing the end
              location

    Returns:
        a list of (latitude, longitude) tuples representing the shortest path
        (in terms of distance) from loc1 to loc2.
    """
    result = []
    nodes_dict = map_rep[0]
    node1 = closest_node(nodes_dict, loc1)
    node2 = closest_node(nodes_dict, loc2)
    path = find_short_path_nodes(map_rep, node1, node2)
    if path is None:
        return None 
    for node in path:
        result.append((nodes_dict[node]['lat'], nodes_dict[node]['lon']))
    return result


def closest_node(nodes_dict, loc):
    """
    Return the closest node (in terms of distance) to the given location.
    """
    smallest_dist = float('inf')
    closest_node = None 
    # For each node, it calculates the distance to the given location and compares them.
    for node in nodes_dict:
        node_loc = (nodes_dict[node]['lat'], nodes_dict[node]['lon'])
        distance = great_circle_distance(loc, node_loc)
        if distance < smallest_dist:
            smallest_dist = distance
            closest_node = node 
    return closest_node
        

def find_fast_path(map_rep, loc1, loc2):
    """
    Return the shortest path between the two locations, in terms of expected
    time (taking into account speed limits).

    Parameters:
        map_rep: the result of calling build_internal_representation
        loc1: tuple of 2 floats: (latitude, longitude), representing the start
              location
        loc2: tuple of 2 floats: (latitude, longitude), representing the end
              location

    Returns:
        a list of (latitude, longitude) tuples representing the shortest path
        (in terms of time) from loc1 to loc2.
    """
    nodes_dict = map_rep[0]
    valid_ways = map_rep[1]
    agenda = Queue()
    expanded = set()
    # Find the closest nodes to the start and end locations.
    node1 = closest_node(nodes_dict, loc1)
    node1_loc = (nodes_dict[node1]['lat'], nodes_dict[node1]['lon'])
    node2 = closest_node(nodes_dict, loc2)
    node2_loc = (nodes_dict[node2]['lat'], nodes_dict[node2]['lon'])
    # If the start and end locations are the same, return a list of the start and end nodes.
    if node1 == node2:
        return [node1_loc]
    # For all valid nodes, initialize its time from source node as infinity.
    for node in nodes_dict.keys():
        nodes_dict[node]['time'] = float('inf')
    # Initialize the source node's time as 0.
    nodes_dict[node1]['time'] = 0
    agenda.add(node1, 0)
    # Initialize the parent pointers dictionary.
    parent = {}
    parent[node1_loc] = 'start'
    # Implement Dijkstra's algorithm.
    while agenda.length() != 0:
        current = agenda.pop_min()
        if current == node2:
            break 
        if current in expanded:
            continue 
        expanded.add(current)
        current_lat = nodes_dict[current]['lat']
        current_lon = nodes_dict[current]['lon']
        current_loc = (current_lat, current_lon)
        for neighbor in nodes_dict[current]['connect']:
            # Add each neighbor and its distance value to the queue.
            agenda.add(neighbor, nodes_dict[neighbor]['time'])
            neighbor_lat = nodes_dict[neighbor]['lat']
            neighbor_lon = nodes_dict[neighbor]['lon']
            neighbor_loc = (neighbor_lat, neighbor_lon)
            # Determine the way id.
            max_speed = -1 * float('inf')
            way_id = None 
            for value in nodes_dict[current]['ways']:
                if neighbor == value[1]:
                    if valid_ways[value[0]]['maxspeed'] > max_speed:
                        max_speed = valid_ways[value[0]]['maxspeed']
                        way_id = value[0]
                if len(value) == 3:
                    if neighbor == value[2]:
                        if valid_ways[value[0]]['maxspeed'] > max_speed:
                            max_speed = valid_ways[value[0]]['maxspeed']
                            way_id = value[0]
            # Edge relaxation technique
            speed = valid_ways[way_id]['maxspeed']
            new_time = nodes_dict[current]['time'] + great_circle_distance(current_loc, neighbor_loc) / speed
            if new_time < nodes_dict[neighbor]['time']:
                agenda.change_distance(neighbor, new_time)
                nodes_dict[neighbor]['time'] = new_time
                parent[neighbor_loc] = current_loc 
    # If the terminal node is unreachable from the source node, then return None.
    if nodes_dict[node2]['time'] == float('inf'):
        return None 
    result = [node2_loc]
    state = node2_loc 
    # Return the list of ids along the shortest (distance) path.
    while parent[state] != 'start':
        result.append(parent[state])
        state = parent[state]
    result.reverse()
    print(len(result))
    return result 


if __name__ == '__main__':
    pass
