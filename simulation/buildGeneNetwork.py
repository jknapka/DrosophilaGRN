''' Build the gene network implied by the input list
of binary gene interactions (source,target). '''
import sys

def buildNetwork(inf):
    interactionMatrix = {}
    for line in inf.readlines():
        src,dst = line.strip().split(',')
        targets = interactionMatrix.setdefault(src,[])
        opTargets = interactionMatrix.setdefault(dst,[])
        if src not in opTargets:
            targets.append(dst)
    return interactionMatrix

def allPathsFrom(gene,netw):
    return map(lambda x: [gene] + x,map(lambda g: allPathsFrom(g,netw),netw[gene]))

if __name__ == '__main__':
    netw = buildNetwork(sys.stdin)
    import random
    srcGenes = netw.keys()
    nextList = srcGenes
    pathLen = 0
    while pathLen < 20 and len(nextList) > 0:
        pathLen += 1
        nextGene = nextList[random.randint(0,len(nextList)-1)]
        print nextGene
        nextList = netw[nextGene]
    print('%s'%(allPathsFrom(srcGenes[random.randint(0,len(srcGenes)-1)],netw),))


