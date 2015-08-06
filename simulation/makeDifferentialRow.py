import sys
import random
import optparse

from simFromGinGraph import *

colNames = ['FBgn',
    'CG number',
    'Gene ID',
    'EC_R1_1',
    'EC_R1_2',
    'EC_R2_1',
    'EC_R2_2',
    'EC_R3_1',
    'EC_R3_2',
    'EC_R4_1',
    'EC_R4_2',
    'EC_R5_1',
    'EC_R5_2',
    'EB_R1_1',
    'EB_R1_2',
    'EB_R2_1',
    'EB_R2_2',
    'EB_R3_1',
    'EB_R3_2',
    'EB_R4_1',
    'EB_R4_2',
    'EB_R5_1',
    'EB_R5_2',
    'EE_R1_1',
    'EE_R1_2',
     'EE_R2_1',
     'EE_R2_2',
     'EE_R3_1',
     'EE_R3_2',
     'EE_R4_1',
     'EE_R4_2',
     'EE_R5_1',
     'EE_R5_2',
     'ISC_R1_1',
     'ISC_R1_2',
     'ISC_R2_1',
     'ISC_R2_2',
     'ISC_R3_1',
     'ISC_R3_2',
     'ISC_R4_1',
     'ISC_R4_2',
     'ISC_R5_1',
     'ISC_R5_2',
]

def setColNames(replicates):
    global colNames
    colNames = ['FBgn','CG number','Gene ID']
    for cellType in ['EC','EB','EE','ISC']:
        for region in ['R1','R2','R3','R4','R5']:
            for rr in range(1,replicates+1):
                colNames.append('%s_%s_%s'%(cellType,region,rr))

def loadGeneStats(fname):
    # "fbid","gene","level","min","mean","max","sd"
    result = {}
    with open(fname) as inf:
        lines = inf.readlines()
        for line in lines:
            line = line.strip()
            toks = line.split(',')
            if toks[0] == '"fbid"':
                continue
            fbid = toks[0][1:-1]
            gnm = toks[1][1:-1]
            if toks[6] == 'NA':
                toks[6] = '0.001'
            ldict = {'fbid':fbid,'gene':gnm,'level':int(toks[2]),'min':float(toks[3]),'mean':float(toks[4]),'max':float(toks[5]),'sd':float(toks[6])}
            levelMap = result.setdefault(fbid,{})
            levelMap[ldict['level']] = ldict
    return result

def getColIndexes(cat,reverse=False):
    def checker(x):
        if x in ['FBgn','CG number','Gene ID']:
            return False
        if reverse:
            return cat not in x
        return cat in x
    cols = [nm for nm in colNames if checker(nm)]
    return map(lambda cn: colNames.index(cn),cols)

geneStats = None

def randomExpValue(fbid,level):
    gDict = geneStats[fbid]
    lDict = gDict[level]
    expVal = lDict['min'] + random.random()*(lDict['max']-lDict['min'])
    #print "Choosing value for gene %s level %s min %s max %s --> %s"%(fbid,level,lDict['min'],lDict['max'],expVal)
    return expVal

interactions = None

def loadInteractions(fname):
    result = []
    with open(fname) as inf:
        lines = inf.readlines()
        for line in lines:
            inter = line.strip().split(',')
            if len(inter) == 2:
                result.append(inter)
    return result

def chooseRandomInteraction():
    global interactions
    if interactions is None:
        interactions = loadInteractions('data/interactions-present-in-data.txt')
    fb1,fb2 = random.choice(interactions)
    return fb1,fb2

def mapLevels(g1dat,g2dat,lvls):
    #print 'Mapping levels %s'%(lvls,)
    result = lvls[:]
    maxG1Lvl = float(max([x for (x,y) in lvls]))
    maxG2Lvl = float(max([y for (x,y) in lvls]))
    for ii in range(len(lvls)):
        mapEntry = []
        for jj in range(2):
            lvl = lvls[ii][jj]
            mxLvl = [maxG1Lvl,maxG2Lvl][jj]
            gdat = (g1dat,g2dat)[jj]
            levels = gdat.keys()
            levels.sort()
            mxTarget = float(levels[-1])
            howFar = (float(lvl)-1.0)/(mxLvl-1.0)
            target = int(1.0 + (mxTarget-1.0)*howFar)
            mapEntry.append(target)
        result[ii] = tuple(mapEntry)
    #print 'result: %s'%(result,)
    return result
        
def buildDiffRows(g1,g2,ebLvls,ebNz,ecLvls,ecNz,eeLvls,eeNz,iscLvls,iscNz):
    global geneStats
    if geneStats is None:
        geneStats = loadGeneStats('data/geneStats.csv')

    g1Dat = geneStats[g1][1]
    g2Dat = geneStats[g2][1]

    def levelMapper(lvls):
        return mapLevels(geneStats[g1],geneStats[g2],lvls)
    ebLvls,ecLvls,eeLvls,iscLvls = map(levelMapper,[ebLvls,ecLvls,eeLvls,iscLvls])

    def concatWithNone(x,y):
        if x is None: x=[]
        if y is None: y=[]
        return x+y

    allLevels = reduce(concatWithNone,(ebLvls,ecLvls,eeLvls,iscLvls))
    Q = max([y for (x,y) in allLevels])
    xQ = max([x for (x,y) in allLevels])
    categories = [('EB',ebLvls,ebNz),('EC',ecLvls,ecNz),('EE',eeLvls,eeNz),('ISC',iscLvls,iscNz)]

    gResult1 = [0.0] * len(colNames)
    gResult2 = [0.0] * len(colNames)

    gResult1[0] = g1
    gResult2[0] = g2
    gResult1[1] = 'CG%s'%g1[4:]
    gResult2[1] = 'CG%s'%g2[4:]
    gResult1[2] = g1Dat['gene']
    gResult2[2] = g2Dat['gene']

    for (cat,valMap,noise) in categories:
        catcols = getColIndexes(cat)

        for col in catcols:
            if valMap is None:
                valMap = [(random.randint(1,xQ),random.randint(1,Q))]
            level1,level2 = random.choice(valMap)
            level1 = applyHouseNoise(level1-1,noise,xQ)+1
            level2 = applyHouseNoise(level2-1,noise,Q)+1
            #print "Chose levels %s=%s, %s=%s"%(g1,level1,g2,level2)
            gResult1[col] = randomExpValue(g1,level1)
            gResult2[col] = randomExpValue(g2,level2)

    return gResult1,gResult2

def toMap(vals,dflt):
    if vals is None and dflt is None:
        return None
    result = []
    if type(vals) == str:
        vals = map(int,vals.split(','))
    if type(dflt) == str:
        dflt = map(int,dflt.split(','))
    if vals is None:
        vals = dflt
    for ii in range(0,len(vals),2):
        result.append((vals[ii],vals[ii+1]))
    return result

def parseOptions(args):
    psr = optparse.OptionParser(conflict_handler='resolve')
    psr.add_option('--f1',dest='fbid1',help='First FBid',action='store',type='string',default=None)
    psr.add_option('--f2',dest='fbid2',help='First FBid',action='store',type='string',default=None)
    psr.add_option('--g1',dest='gene1',help='First gene',action='store',type='string',default=None)
    psr.add_option('--g2',dest='gene2',help='Second gene',action='store',type='string',default=None)
    psr.add_option('--EB',dest="eb_levels",help="Enteroblast levels g1,g2,g1,g2...",action="store",type='string',default=None)
    psr.add_option('--EC',dest='ec_levels',help='Entrocyte levels g1,g2,g1,g2...',action='store',type='string',default=None)
    psr.add_option('--EE',dest='ee_levels',help='Entroendocrine levels g1,g2,g1,g2...',action='store',type='string',default=None)
    psr.add_option('--ISC',dest='isc_levels',help='Intestinal step cell levels g1,g2,g1,g2...',action='store',type='string',default=None)
    psr.add_option('-o','--other',dest='olevels',help='Default levels g1,g2,g1,g2...',action='store',type='string',default=None)
    psr.add_option('--neb',help='House noise level for EB cells 0 .. 1.0',dest='houseNoiseEB',action='store',type='float',default=None)
    psr.add_option('--nec',help='House noise level for EC cells 0 .. 1.0',dest='houseNoiseEC',action='store',type='float',default=None)
    psr.add_option('--nee',help='House noise level for EE cells 0 .. 1.0',dest='houseNoiseEE',action='store',type='float',default=None)
    psr.add_option('--nisc',help='House noise level ISCs 0 .. 1.0',dest='houseNoiseISC',action='store',type='float',default=None)
    psr.add_option('-n','--n-other',help='House noise level default 0 .. 1.0',dest='houseNoise',action='store',type='float',default=0.0)
    psr.add_option('-r','--replicates',help='Replicates in each cell type X region',dest='replicates',default=2,type='int')
    return psr.parse_args()

if __name__ == '__main__':
    opts,args = parseOptions(sys.argv)
    loadGeneStats('data/geneStats.csv')
    if opts.fbid1 == '_':
        opts.fbid1,opts.fbid2 = chooseRandomInteraction()
    if (opts.fbid1 is None):
        opts.fbid1 = geneToFbid(opts.gene1)
    if (opts.fbid2 is None):
        opts.fbid2 = geneToFbid(opts.gene2)
    fbid1 = opts.fbid1
    fbid2 = opts.fbid2

    if opts.houseNoiseEC is None: opts.houseNoiseEC = opts.houseNoise
    if opts.houseNoiseEE is None: opts.houseNoiseEE = opts.houseNoise
    if opts.houseNoiseEB is None: opts.houseNoiseEB = opts.houseNoise
    if opts.houseNoiseISC is None: opts.houseNoiseISC = opts.houseNoise

    opts.eb_levels = toMap(opts.eb_levels,opts.olevels)
    opts.ec_levels = toMap(opts.ec_levels,opts.olevels)
    opts.ee_levels = toMap(opts.ee_levels,opts.olevels)
    opts.isc_levels = toMap(opts.isc_levels,opts.olevels)

    setColNames(opts.replicates)
    r1,r2 = buildDiffRows(fbid1,fbid2,opts.eb_levels,opts.houseNoiseEB,opts.ec_levels,opts.houseNoiseEC,\
            opts.ee_levels,opts.houseNoiseEE,opts.isc_levels,opts.houseNoiseISC)

    print ','.join(colNames)
    print ','.join(map(str,r1))
    print ','.join(map(str,r2))

