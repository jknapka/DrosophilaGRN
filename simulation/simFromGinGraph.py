# Build a simulated transcriptome expression dataset from a
# GINsim simulation graph.
from __future__ import print_function
import re
import sys
import random

def valuesFromName(name):
    ''' Split a name like s0123 into a list of integer values
    like [0,1,2,3] '''
    # One value per character.
    while name[0] not in "0123456789":
        name = name[1:]
    values = [0] * len(name)
    for ii in range(len(values)):
        values[ii] = int(name[ii])
    return values

class Node:

    def __init__(self,name,levels,genes):
        self.name = name
        self.vals = levels
        self.genes = genes
        self.edges = []

    def addTransition(self,destNode,kind):
        self.edges.append((destNode,kind))

NODE_ORDER_RE = re.compile('nodeorder="(([^ ]+) ?)+"')
NODE_RE = re.compile('<node id="([^"]*)">')
EDGE_RE = re.compile('<edge .* from="([^"]+)" to="([^"]+)">')

def parseGraph(lines):
    ''' Parse a GINsim XML file and derive the network structure. '''
    node_order = []
    nodes = {}
    for line in lines:
        line = line.strip()
        m = NODE_ORDER_RE.search(line)
        if m is not None:
            print('Got NODE_ORDER: %s'%line)
            node_order = m.group(1).split()
            print('  node_order: %s'%(node_order,))
            continue
        m = NODE_RE.search(line)
        if m is not None:
            print('Got NODE: %s'%line)
            node = Node(m.group(1),valuesFromName(m.group(1)),node_order)
            nodes[node.name] = node
            print('  node: %s [%s]'%(node.name,node.vals))
            continue
        m = EDGE_RE.search(line)
        if m is not None:
            print('Got EDGE: %s'%line)
            nFrom = nodes[m.group(1)]
            nTo = nodes[m.group(2)]
            nFrom.edges.append(nTo)
            print(' edge from: %s to: %s'%(nFrom.name,nTo.name))
    return nodes

def pathFromNode(node,nPoints):
    path = [node]
    def chooseEdge(node):
        return node.edges[random.randint(0,len(node.edges)-1)]
    while len(path) < nPoints and node.edges:
        node = chooseEdge(node)
        path.append(node)
    return path

def simFromNode(node,nPoints):
    path = [node]

if __name__ == '__main__':
    lines = open(sys.argv[1])
    nodes = parseGraph(lines)
    for node in nodes.values():
        print('node %s -> %s'%(node.name,[n.name for n in node.edges]))
    if len(sys.argv) > 1:
        pathLen = int(sys.argv[2])
        path = pathFromNode(nodes.values()[random.randint(0,len(nodes.values())-1)],pathLen)
        print('path: %s'%([node.name for node in path],))
