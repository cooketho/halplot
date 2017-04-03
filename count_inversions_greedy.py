import newick
import sys
from collections import namedtuple
from itertools import permutations
from itertools import combinations
from itertools import combinations_with_replacement
from itertools import product
from scipy.misc import comb
import pickle
from numpy import random

random.seed(0)
HalEdge = namedtuple('HalEdge', 'parentGenome parentSeq parentStart parentEnd childGenome childSeq childStart childEnd rev')

class SeqPair:
    def __init__(self, parentGenome, childGenome, parentSeq, childSeq, bp=0, inversions=(0, 0), crosses={}):
        self.parentGenome = parentGenome
        self.parentSeq = parentSeq
        self.childGenome = childGenome
        self.childSeq = childSeq
        self.bp = bp
        self.inversions = inversions
        self.crosses = crosses
    def as_tuple(self, partner=None):
        return (self.parentGenome, self.parentSeq, self.childGenome, self.childSeq)
    def as_seq_tuple(self):
        return (self.parentSeq, self.childSeq)

script, tree, edges = sys.argv

def get_config(node, config, seqPairs):
    # Recursively calculate best configuration for nodes in the tree.
    inversionScore = 0
    if node.descendants:
        # Sum of sizes of edges mapped to each sequence in genome, keyed by sequence name.
        edgesPerSeq = {}
        # Lists of SeqPair objects, keyed by parent sequence name.
        pairsByParentSeq = {}
        for child in node.descendants:
            # Call self recursively to get child configurations, which must be fixed before proceeding.
            inversionScore2, config2 = get_config(child, config, seqPairs)
            config = config2
            inversionScore += inversionScore2
            # Get all SeqPair objects shared between parent and child genomes.
            try:
                pairsByChild = seqPairs[(node.name, child.name)].values()
            except:
                pairsByChild = []
            for p in pairsByChild:
                try:
                    edgesPerSeq[p.parentSeq] += p.bp
                    pairsByParentSeq[p.parentSeq].append(p)
                except:
                    edgesPerSeq[p.parentSeq] = p.bp
                    pairsByParentSeq[p.parentSeq] = [p]
        # Order in which sequences are added is random, and weighted by the total bp of edges mapped to them.
        totalSize = sum(edgesPerSeq.values())
        seqSamplingProb = [float(edgesPerSeq[s]) / totalSize for s in edgesPerSeq.keys()]
        # Generate a random order of the parental sequences.
        sequences = list(random.choice(edgesPerSeq.keys(), size=len(edgesPerSeq), replace=False, p=seqSamplingProb))
        # Numeric index of each parent seq in parent genome.
        seqIndex = {}
        # List of pairs to iterate over when calculating inversions.
        pairs = []
        # Initialize loop variable to store sequence orders as sequences are added.
        seqOrder = []
        seqFlip = []
        # Add sequences to index one at a time until it is full.
        inversionsLast = 0
        while (len(sequences) > 0):
            # Sample the next parent sequence.
            seq = sequences.pop(0)
            # Get all SeqPair objects mapping to this sequence.
            pairs += pairsByParentSeq[seq]
            # Inversion counts, keyed by tuple of (index, flip) tuples for each seq.
            inversions = {}
            for i in range(len(seqOrder) + 1):
                for flip in (0, 1):
                    seqOrderNext = seqOrder[:]
                    seqOrderNext.insert(i, seq)
                    seqFlipNext = seqFlip[:]
                    seqFlipNext.insert(i, flip)
                    # Update the configuration.
                    configNext = zip(range(len(seqOrderNext)), seqFlipNext)
                    config[node.name] = dict(zip(seqOrderNext, configNext))
                    # Count the number of inversions for this configuration.
                    inversions[(tuple(seqOrderNext), tuple(seqFlipNext))] = count_inversions(pairs, config)
            #for i in inversions:
            #    print i, inversions[i]
            #print
            # Configuration is chosen randomly with probabilities weighted by reciprocal of number of inversions.
            #inversionReciprocalSum = sum([1 / float(i - inversionsLast + 1) for i in inversions.values()])
            #orderSamplingProb = [1 / (float(i - inversionsLast + 1) * inversionReciprocalSum) for i in inversions.values()]
            #i = random.choice(range(len(inversions)), size=1, p=orderSamplingProb)
            #inversionsLast = inversions.values()[i]
            #seqOrder = list(inversions.keys()[i][0])
            #seqFlip = list(inversions.keys()[i][1])

            bestConfig = min(inversions, key=inversions.get)
            seqOrder = list(bestConfig[0])
            seqFlip = list(bestConfig[1])

            # Update config with the best sequence order found.
            configNext = zip(range(len(seqOrder)), seqFlip)
            config[node.name] = dict(zip(seqOrder, configNext))
        inversionScore += count_inversions(pairs, config)
    else:
        pass
    return (inversionScore, config)

def count_inversions(pairs, index):
    inversions = 0
    for pair1, pair2 in combinations_with_replacement(pairs, 2):
        if pair1.childGenome != pair2.childGenome:
            continue
        # Get parent and child seq indices in this configuration.
        pInd1 = index[pair1.parentGenome][pair1.parentSeq][0]
        pInd2 = index[pair2.parentGenome][pair2.parentSeq][0]
        cInd1 = index[pair1.childGenome][pair1.childSeq][0]
        cInd2 = index[pair2.childGenome][pair2.childSeq][0]
        # Get parent and child seq orientations in this configuration.
        pFlip = index[pair1.parentGenome][pair1.parentSeq][1]
        cFlip = index[pair1.childGenome][pair1.childSeq][1]
        if pInd1 < pInd2:
            if cInd1 > cInd2:
                inversions += pair1.bp * pair2.bp
            elif cInd1 == cInd2:
                if cFlip:
                    inversions += pair1.crosses[(pair2.parentSeq, pair2.childSeq)][1]
                else:
                    inversions += pair1.crosses[(pair2.parentSeq, pair2.childSeq)][0]
        elif pInd1 > pInd2:
            if cInd1 < cInd2:
                inversions += pair1.bp * pair2.bp
            elif cInd1 == cInd2:
                if cFlip:
                    inversions += pair1.crosses[(pair2.parentSeq, pair2.childSeq)][0]
                else:
                    inversions += pair1.crosses[(pair2.parentSeq, pair2.childSeq)][1]
        elif pInd1 == pInd2:
            if cInd1 < cInd2:
                if pFlip:
                    inversions += pair1.crosses[(pair2.parentSeq, pair2.childSeq)][1]
                else:
                    inversions += pair1.crosses[(pair2.parentSeq, pair2.childSeq)][0]
            elif cInd1 > cInd2:
                if pFlip:
                    inversions += pair1.crosses[(pair2.parentSeq, pair2.childSeq)][0]
                else:
                    inversions += pair1.crosses[(pair2.parentSeq, pair2.childSeq)][1]
            elif cInd1 == cInd2:
                if pFlip == cFlip:
                    inversions += pair1.inversions[0]
                else:
                    inversions += pair1.inversions[1]
    return inversions

# Dict of dict of list of edges, keyed by (parentGenome, childGenome), (parentSeq, childSeq).
edgeLists = {}
# Sequence names, keyed by genome name.
sequences = {}
edges = open(edges, 'r')
for line in edges:
    line = line.strip()
    edge = HalEdge(*line.split())
    parentGenome, parentSeq, parentStart, parentEnd, childGenome, childSeq, childStart, childEnd, rev = line.split()
    parentStart = int(parentStart)
    parentEnd = int(parentEnd)
    childStart = int(childStart)
    childEnd = int(childEnd)
    rev = bool(rev)
    edge = HalEdge(parentGenome, parentSeq, parentStart, parentEnd, childGenome, childSeq, childStart, childEnd, rev)
    try:
        edgeLists[(parentGenome, childGenome)][(parentSeq, childSeq)].append(edge)
    except:
        try:
            edgeLists[(parentGenome, childGenome)][(parentSeq, childSeq)] = [edge]
        except:
            edgeLists[(parentGenome, childGenome)] = {(parentSeq, childSeq):[edge]}
edges.close()

# Make dict of lists of SeqPair objects, keyed by (parentGenome, childGenome), (parentSeq, childSeq).
seqPairs = {}
# Loop over all parent-child genome pairs.
for genomePair in edgeLists:
    seqPairs[genomePair] = {}
    for seqPair in edgeLists[genomePair]:
        bp = 0
        inv0 = 0
        inv1 = 0
        for i in range(len(edgeLists[genomePair][seqPair])):
            edge1 = edgeLists[genomePair][seqPair][i]
            bp1 = abs(edge1.parentEnd - edge1.parentStart)
            bp += bp1
            # Count inversions due to inverted edges.
            if edge1.rev:
                inv0 += int(comb(bp1, 2))
            else:
                inv1 += int(comb(bp1, 2))
            # Count inversions due to crossed edges.
            for j in range(i + 1, len(edgeLists[genomePair][seqPair])):
                edge2 = edgeLists[genomePair][seqPair][j]
                bp2 = abs(edge2.parentEnd - edge2.parentStart)
                if edge1.parentStart < edge2.parentStart:
                    if edge1.childStart > edge2.childStart:
                        inv0 += bp1 * bp2
                    else:
                        inv1 += bp1 * bp2
                elif edge1.parentStart > edge2.parentStart:
                    if edge1.childStart < edge2.childStart:
                        inv0 += bp1 * bp2
                    else:
                        inv1 += bp1 * bp2
        seqPairs[genomePair][seqPair] = SeqPair(parentGenome=genomePair[0],
                                                childGenome=genomePair[1],
                                                parentSeq=seqPair[0],
                                                childSeq=seqPair[1],
                                                bp=bp, inversions=(inv0, inv1),
                                                crosses={})
    # Loop over all parent-child pairs of sequences.
    for pair1, pair2 in combinations(seqPairs[genomePair].values(), 2):
        # Inversions when pair2 seq is to the right of pair1 seq.
        cross1 = 0
        # Inversions when pair2 seq is to the left of pair1 seq.
        cross2 = 0
        # Make lists of edges to iterate over.
        edges1 = edgeLists[genomePair][(pair1.parentSeq, pair1.childSeq)]
        edges2 = edgeLists[genomePair][(pair2.parentSeq, pair2.childSeq)]
        # If parent seq or child seq is shared.
        if pair1.parentSeq == pair2.parentSeq or pair1.childSeq == pair2.childSeq:
            # Count crosses.
            for edge1 in edges1:
                for edge2 in edges2:
                    bp1 = abs(edge1.parentEnd - edge1.parentStart)
                    bp2 = abs(edge2.parentEnd - edge2.parentStart)
                    # Parent seq is the same.
                    if pair1.parentSeq == pair2.parentSeq:
                        if edge2.parentStart < edge1.parentStart:
                            cross1 += bp1 * bp2
                        else:
                            cross2 += bp1 * bp2
                    # Child seq is the same.
                    elif pair1.childSeq == pair2.childSeq:
                        if edge2.childStart < edge1.childStart:
                            cross1 += bp1 * bp2
                        else:
                            cross2 += bp1 * bp2
        pair1.crosses[(pair2.parentSeq, pair2.childSeq)] = (cross1, cross2)
        pair2.crosses[(pair1.parentSeq, pair1.childSeq)] = (cross2, cross1)

#pickle.dump(seqPairs, open('seqpairs.p', 'wb'))

#seqPairs = pickle.load(open('seqpairs.p', 'rb'))

tree = newick.read(tree)
config = {'medaka':{'NC_019878.1':(0, 1)},
          'anole':{'NW_003338722.1':(0, 1)},
          'budgie':{'NW_004848279.1':(0, 0)},
          'human':{'NT_008705.17':(0, 0)}, 
          'opossum':{'NW_001582021.1':(0, 1)}}

bestScore = None
bestConfig = {}
for i in range(50):
    score, conf = get_config(tree[0], config, seqPairs)
    if not bestScore or score < bestScore:
        bestScore = score
        bestConfig = conf

#print bestScore
for genome in bestConfig:
    for seq in bestConfig[genome]:
        print "%s\t%s\t%i\t%i" % (genome, seq, bestConfig[genome][seq][0], bestConfig[genome][seq][1])

