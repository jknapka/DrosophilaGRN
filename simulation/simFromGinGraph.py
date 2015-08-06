# Build a simulated transcriptome expression dataset from a
# GINsim simulation graph.
from __future__ import print_function
import re
import sys
import random
import numpy.random as nrnd

def houseNoise(y,y1,theta,Q):
    y,y1,theta,Q = map(float,(y,y1,theta,Q))
    if y != y1:
        return ( (1.0 - ((abs( y1-y ) )/( sum( [abs(d-y) for d in range(int(Q))] ) ) ) ) * (theta / (Q-1.0)) ) * (1.0 - theta) + theta/Q
    else:
        return ((theta/(Q-1.0)) +1.0 -theta) * (1.0 - theta) + (theta/Q)

def houseNoiseDistribution(y,theta,Q):
    ''' Build a probability distribution for the possible
    values of y1 between 0 and Q-1, given the value y,
    the noise level = theta, and the number of levels of
    y = Q '''
    dist = [houseNoise(y,y1,theta,Q) for y1 in range(Q)]
    #error = sum(dist) - 1.0
    #error = error / float(Q)
    #dist = map(lambda x: x-error,dist)
    return dist

def applyHouseNoise(y,theta,Q):
    ''' Build the probability distribution for y based on 
    the house noise model and the given value of Y, and
    return the new y with noise applied. '''
    dist = houseNoiseDistribution(y,theta,Q)
    accum = 0.0
    choice = random.random()
    for y1 in range(Q):
        accum += dist[y1]
        if accum >= choice:
            return y1
    # We should never get here, but if we do,
    # we should return Q-1.
    print("UNEXPECTED: choice value %s > sum of discrete distribution %s."%(choice,sum(dist)))
    return Q-1

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
    ''' Represents a graph node. '''
    def __init__(self,name,levels,genes):
        self.name = name
        self.levels = levels
        self.genes = genes
        self.edges = []

    def addTransition(self,destNode,kind):
        self.edges.append((destNode,kind))

NODE_ORDER_RE = re.compile('nodeorder="([^"]*)"')
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
            node = Node(m.group(1),valuesFromName(m.group(1)),node_order[:])
            nodes[node.name] = node
            print('  node: %s [%s]'%(node.name,node.levels))
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

def floatIdentity(x):
    return float(x)

def constantFn(x):
    return lambda y: x

def simFromNode(node,pathLen,nSamples,activationFns={'DEFAULT':floatIdentity},noiseFns={'DEFAULT': constantFn(0.5)}):
    ''' Create a list of samples of gene activity along
    a state path in the network starting at the given node.
    Each sampling operation will generate an expression
    level for each gene at a given node. '''
    path = pathFromNode(node,pathLen)
    return simForPath(path,nSamples)

def simForPath(path,nSamples,activationFns={'DEFAULT':floatIdentity},noiseFns={'DEFAULT': constantFn(0.5)},randomize=False):

    sampFn = doSample
    if randomize:
        sampFn = randomSample

    # Figure out how many samples to take at each point.
    # MOSTLY DONE: parameterize the relationship between activation
    # levels and expression levels. This might differ for each
    # gene so a map is required.
    samplesPerPoint = nSamples / len(path)
    extraSamples = nSamples % len(path)

    sampleList = []

    print("Samples per point: %s; extra samples %s"%(samplesPerPoint,extraSamples))

    for nd in path:
        for samp in range(samplesPerPoint):
            sampleList = sampleList + sampFn(nd,samp,activationFns,noiseFns)
        if extraSamples > 0:
            sampleList = sampleList + sampFn(nd,samp+1,activationFns,noiseFns)
            extraSamples -= 1

    return sampleList

def getFunctionForGene(gname,fmap):
    fn = fmap.get(gname,None)
    if fn is None:
        fn = fmap.get('DEFAULT',None)
        if fn is None:
            fn = constantFn(-42)
    return fn

def doSample(node,replicate,activationFns,noiseFns):
    ''' Sample all genes at the given node. '''
    result = []
    print("Sampling genes: %s"%(node.genes,))
    for gidx in range(len(node.genes)):
        gene = node.genes[gidx]

        # Compose a sampler function from the activation and noise
        # functions for the gene, or the default ones if no gene-
        # specific ones are supplied.
        sampler = lambda x: getFunctionForGene(gene,activationFns)(x) + getFunctionForGene(gene,noiseFns)(x)

        level = node.levels[gidx]
        sample = (gene,'%s_%s'%(node.name,replicate),max(sampler(level),0.0))
        print("Adding sample: %s"%(sample,))
        result.append(sample)

    return result

def randomSample(node,replicate,activationFns,noiseFns):
    result = []
    print("Randomly sampling genes: %s"%(node.genes,))
    for gidx in range(len(node.genes)):
        gene = 'RG%s'%random.randint(0,10000000) 

        # Compose a sampler function from the activation and noise
        # functions for the gene, or the default ones if no gene-
        # specific ones are supplied.
        sampler = lambda x: getFunctionForGene(gene,activationFns)(x) + getFunctionForGene(gene,noiseFns)(x)

        level = node.levels[random.randint(0,len(node.genes)-1)]
        sample = (gene,'%s_%s'%(node.name,replicate),max(sampler(level),0.0))
        print("Adding sample: %s"%(sample,))
        result.append(sample)

    return result

def linearActivation(factor):
    def activation(level):
        return level * factor
    return activation

def squareActivation(factor):
    def activation(level):
        return (level*level*factor)
    return activation

def expActivation(base,factor):
    def activation(level):
        return base**level * factor
    return activation

def normalNoise(scale):
    def noise(x):
        return nrnd.normal(scale=scale)
    return noise

def scaledNormalNoise(factor,scale):
    def noise(x):
        if x <= 0.0:
            x = factor
        return nrnd.normal(scale=x*scale/factor)
    return noise

def printSamples(samples):
    genes = {}
    conditions = {}
    for (gene,condition,value) in samples:
        genes[gene] = True
        conditions[condition] = True
    genes = genes.keys()
    genes.sort()
    conditions = conditions.keys()
    conditions.sort()

    condList = [None] * (1+len(conditions))
    matrix = [None] * (1+len(genes))
    for ii in range(len(matrix)):
        matrix[ii] = condList[:]
        if ii>0:
            matrix[ii][0] = genes[ii-1]

    matrix[0][0] = 'Gene\Condition'
    matrix[0][1:] = conditions[:]

    for (gene,condition,value) in samples:
        geneIdx = [c[0] for c in matrix].index(gene)
        condIdx = matrix[0].index(condition)
        matrix[geneIdx][condIdx] = value

    for line in matrix:
        print(','.join(map(str,line)))

if __name__ == '__main__':
    if sys.argv[1] == 'noiseTest':
        dist = houseNoiseDistribution(int(sys.argv[2]),float(sys.argv[4]),int(sys.argv[3]))
        print("House noise distribution: %s sums to %s"%(dist,sum(dist)))
        sys.exit(0)
    lines = open(sys.argv[1])
    nodes = parseGraph(lines)
    for node in nodes.values():
        print('node %s -> %s'%(node.name,[n.name for n in node.edges]))
    if len(sys.argv) > 2:
        pathLen = int(sys.argv[2])
        path = pathFromNode(nodes.values()[random.randint(0,len(nodes.values())-1)],pathLen)
        print('path: %s'%([node.name for node in path],))

        activationFns = {'DEFAULT':squareActivation(3.0)}
        noiseFns = {'DEFAULT':scaledNormalNoise(3.0,10.0) }

        samples = simForPath(path,8,activationFns,noiseFns)
        print('%s samples: %s'%(len(samples),samples,))

        if len(sys.argv)>3:
            nRandom = int(sys.argv[3])
            samples = samples + simForPath(path,8,activationFns,noiseFns,randomize=True)

        printSamples(samples)


