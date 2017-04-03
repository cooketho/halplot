#include <iostream>
#include <cstdlib>
#include <hal.h>
#include <fstream>
#include "./intervaltree/IntervalTree.h"

using namespace std;
using namespace hal;

class HalEdge{
    // Holds information about an edge connecting two segments in a HAL graph.
public:
    string parentGenome;
    string parentSeq;
    int parentStart;
    int parentEnd;
    string childGenome;
    string childSeq;
    int childStart;
    int childEnd;
    bool reversed;
    set<string> features;
    HalEdge(const BottomSegment* parent, const TopSegment* child)
        : parentGenome(parent->getGenome()->getName())
        , parentSeq(parent->getSequence()->getName())
        , parentStart(parent->getStartPosition())
        , parentEnd(parent->getEndPosition())
        , childGenome(child->getGenome()->getName())
        , childSeq(child->getSequence()->getName())
        , childStart(child->getStartPosition())
        , childEnd(child->getEndPosition())
        , reversed(child->getParentReversed())
    {}
    void print(ostream& os, string feature) const
    {
        os << parentGenome << "\t";
        os << parentSeq << "\t";
        os << parentStart << "\t";
        os << parentEnd + 1 << "\t";
        os << childGenome << "\t";
        os << childSeq << "\t";
        os << childStart << "\t";
        os << childEnd + 1 << "\t";
        os << reversed << "\t";
        os << feature << endl;
    }
};
class HalPath{
    // Holds information about a chain of edges connecting two segments in a HAL graph.
public:
    vector<HalEdge> edges;
    string parent;
    string child;
    vector<string> features;
};

static void printAllEdges(ostream& os, AlignmentConstPtr alignment);
static void printSubgraphByNode(ostream& os, AlignmentConstPtr alignment, string rootGenome,
                                map<string, map<string, IntervalTree<string> > > intervalTrees,
                                float overlapA, float overlapB);
static void printSubgraphByLeaves(ostream& os, AlignmentConstPtr alignment,
                                  map<string, map<string, IntervalTree<string> > > intervalTrees);
static void getIntervals(ostream&os, string bedFile, 
                         map<string, map<string, IntervalTree<string> > >& intervalTrees);
static vector<HalPath> traverseDownRecursive(ostream& os, AlignmentConstPtr alignment,
                                          BottomSegmentIteratorPtr parent, HalPath parentPath,
                                          float overlapA, float overlapB);
static void traverseUpRecursive(ostream& os, AlignmentConstPtr alignment,
                                TopSegmentIteratorPtr child,
                                map<string, set<HalEdge> >& edges);

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->setDescription("Print HAL graph edges that connect user-specified intervals in two or more child genomes through the user-specified root genome");
  optionsParser->addArgument("halFile", "path to hal file to analyze");
  optionsParser->addOption("bed", "BED file of intervals "
                           "in bed file", "\"\"");
  optionsParser->addOption("root", "Name of root genome "
                           , "\"\"");

  string path;
  string bedFile;
  string rootGenome;
  try
  {
    optionsParser->parseOptions(argc, argv);
    path = optionsParser->getArgument<string>("halFile");
    bedFile = optionsParser->getOption<string>("bed");
    rootGenome = optionsParser->getOption<string>("root");

    size_t optCount = 0;
    if (bedFile != "\"\"") ++optCount;
    if (optCount > 1)
    {
        throw hal_exception("hmm");
    }
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
  }
  try
  {
    AlignmentConstPtr alignment = openHalAlignmentReadOnly(path, optionsParser);
    if (bedFile !=  "\"\"" && rootGenome != "\"\"")
        {
            map<string, map<string, IntervalTree<string> > > featureTrees;
            getIntervals(cout, bedFile, featureTrees);
            printSubgraphByNode(cout, alignment, rootGenome, featureTrees, 1.0, 0.0);
            //printSubgraphByLeaves(cout, alignment, featureTrees);
        }
    else
        {
            printAllEdges(cout, alignment);
        }
  }
  catch(hal_exception& e)
  {
    cerr << "hal exception caught: " << e.what() << endl;
    return 1;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    return 1;
  }
  
  return 0;
}

void printAllEdges(ostream& os, AlignmentConstPtr alignment)
{
    // Print all edges in the HAL graph encoded by alignment.
    set<const Genome*> genomes;
    const Genome* root = alignment->openGenome(alignment->getRootName());
    getGenomesInSubTree(root, genomes);
    for (set<const Genome*>::iterator i = genomes.begin(); i != genomes.end();
         ++i)
        {
            TopSegmentIteratorPtr top;
            TopSegmentIteratorConstPtr topEnd;
            BottomSegmentIteratorPtr bot;
            const Genome* parent = (*i)->getParent();
            top = (*i)->getTopSegmentIterator();
            topEnd = (*i)->getTopSegmentEndIterator();
            if (parent == NULL)
                {
                    continue;
                }
            bot = parent->getBottomSegmentIterator();
            while(top != topEnd)
                {
                    if (top->hasParent())
                        {
                            bot->toParent(top);
                            HalEdge edge = HalEdge(bot->getBottomSegment(), top->getTopSegment());
                            edge.print(os, "");
                        }
                    top->toRight();
                }
        }
}

void printSubgraphByNode(ostream& os, AlignmentConstPtr alignment, string rootGenome,
                         map<string, map<string, IntervalTree<string> > > intervalTrees,
                         float overlapA, float overlapB)
{
    // Print all edges connected to a given node.
    // Node is root right now, but plan to generalize soon.
    const Genome* root = alignment->openGenome(rootGenome);
    vector<string> leafNames = alignment->getLeafNamesBelow(alignment->getRootName());
    // Only output paths that span this set of leaves through root.
    set<string> leafSet(leafNames.begin(), leafNames.end());
    BottomSegmentIteratorPtr bot = root->getBottomSegmentIterator();
    BottomSegmentIteratorConstPtr botEnd = root->getBottomSegmentEndIterator();
    while (bot != botEnd)
        {
            vector<HalPath> paths;
            vector<HalPath> outputPaths;
            HalPath path;
            path.parent = root->getName();
            paths = traverseDownRecursive(os, alignment, bot, path, overlapA, overlapB);
            bot->toRight();
            // Set of child genomes containing a target feature that is linked to root via edges in a HalPath.
            set<string> leaves;
            for (vector<HalPath>::iterator i = paths.begin(); i != paths.end(); ++i)
                {
                    // Get terminal edge in the HalPath.
                    HalEdge lastEdge = *((*i).edges.end() - 1);
                    // Get interval tree of target features in the sequence at the terminus of HalPath.
                    IntervalTree<string> targetFeatures;
                    if (intervalTrees.find(lastEdge.childGenome) != intervalTrees.end())
                        {
                            if (intervalTrees[lastEdge.childGenome].find(lastEdge.childSeq) != intervalTrees[lastEdge.childGenome].end())
                                {
                                    targetFeatures = intervalTrees[lastEdge.childGenome][lastEdge.childSeq];
                                }
                        }
                    // Get vector of target features that overlap with the top segment of the terminal edge.
                    vector<Interval<string> > features;
                    targetFeatures.findOverlapping(lastEdge.childStart, lastEdge.childEnd, features);
                    // Assign vector of feature names to HalPath.
                    for (vector<Interval<string> >::iterator feature = features.begin(); feature != features.end(); ++feature)
                        {
                            (*i).features.push_back((*feature).value);
                        }
                    // Don't consider HalPaths that don't overlap with a target feature.
                    if (features.size() > 0)
                        {
                            leaves.insert(lastEdge.childGenome);
                            outputPaths.push_back(*i);
                        }
                }
            // Only print HalPaths that connect target features in more than one genome.
            if (leaves.size() > 1)
                {
                    for (vector<HalPath>::iterator i = outputPaths.begin(); i != outputPaths.end(); ++i)
                        {
                            for (vector<HalEdge>::iterator j = (*i).edges.begin(); j != (*i).edges.end(); ++j)
                                {
                                    for (vector<string>::iterator f = (*i).features.begin(); f != (*i).features.end(); ++f)
                                        {
                                            (*j).print(os, *f);
                                        }
                                }
                        }
                }
        }
}

void printSubgraphByLeaves(ostream& os, AlignmentConstPtr alignment,
                           map<string, map<string, IntervalTree<string> > > intervalTrees)
{
    // Print all edges connected to any one of a set of leaves.
    // Function is unfinished--doesn't do anything yet!
    set<const Genome*> genomes;
    const Genome* root = alignment->openGenome(alignment->getRootName());
    // Get leaf genomes to iterate over.
    getGenomesInSubTree(root, genomes);
    for (set<const Genome*>::iterator i = genomes.begin(); i != genomes.end(); ++i)
        {
            os << (*i)->getName() << endl;
            TopSegmentIteratorPtr top = (*i)->getTopSegmentIterator();
            TopSegmentIteratorConstPtr topEnd = (*i)->getTopSegmentEndIterator();
            while (top != topEnd)
                {
                    top->toRight();
                }
        }
}

void getIntervals(ostream&os, string bedFile,
                  map<string, map<string, IntervalTree<string> > >& intervalTrees)
{
    ifstream fin(bedFile.c_str());
    string genomeName;
    string seq;
    int start;
    int end;
    string feature;
    map<string, map<string, vector<Interval<string> > > > intervalVectors;
    while (fin >> genomeName >> seq >> start >> end >> feature)
        {
            end = end - 1;
            intervalVectors[genomeName][seq].push_back(Interval<string>(start, end, feature));
        }
    map<string, map<string, vector<Interval<string> > > >::iterator i = intervalVectors.begin();
    while (i != intervalVectors.end())
        {
            map<string, vector<Interval<string> > >::iterator j = i->second.begin();
            while (j != i->second.end())
                {
                    intervalTrees[i->first][j->first] = IntervalTree<string>(j->second);
                    ++j;
                }
            ++i;
        }
}


vector<HalPath> traverseDownRecursive(ostream& os, AlignmentConstPtr alignment,
                                      BottomSegmentIteratorPtr parent,
                                      HalPath parentPath,
                                      float overlapA, float overlapB)
{
    vector<HalPath> paths;
    // Get names of the child genomes.
    vector<string> children = alignment->getChildNames(parent->getBottomSegment()->getGenome()->getName());
    for (vector<string>::iterator i = children.begin(); i != children.end(); ++i)
        {
            const Genome* childGenome = alignment->openGenome(*i);
            // Make iterator for top segments in child genome.
            TopSegmentIteratorPtr top = childGenome->getTopSegmentIterator();
            TopSegmentIteratorConstPtr topEnd = childGenome->getTopSegmentEndIterator();
            // If parent has a child segment in this genome, move top iterator to it.
            if (parent->hasChildG(childGenome))
                {
                    top->toChildG(parent, childGenome);
                }
            else
                {
                    // Otherwise set top iterator to end.
                    top = topEnd;
                }
            while (top != topEnd)
                {
                    HalEdge edge = HalEdge(parent->getBottomSegment(), top->getTopSegment());
                    // Make iterator for bottom segments in child genome.
                    BottomSegmentIteratorPtr nextParent;
                    BottomSegmentIteratorConstPtr nextParentEnd;
                    nextParentEnd = childGenome->getBottomSegmentEndIterator();
                    // If top segment has a parse edge to a bottom segment, move bottom iterator to it.
                    if (top->hasParseDown())
                        {
                            nextParent = childGenome->getBottomSegmentIterator(top->getBottomParseIndex());
                        }
                    else
                        {
                            // Otherwise set iterator to end.
                            nextParent = nextParentEnd;
                            // Recursion stops here for top segment, so add terminal path to the vector to return.
                            HalPath terminalPath = parentPath;
                            terminalPath.child = childGenome->getName();
                            terminalPath.edges.push_back(edge);
                            paths.push_back(terminalPath);
                        }
                    while (nextParent != nextParentEnd)
                        {
                            int overlap = nextParent->getBottomSegment()->getEndPosition() - top->getTopSegment()->getStartPosition();
                            int lenB = nextParent->getBottomSegment()->getLength();
                            // Require some minimum overlap with the next bottom segment to continue the recursion.
                            if (float(overlap) / float(lenB) >= overlapB)
                                {
                                    HalPath nextParentPath = parentPath;
                                    nextParentPath.child = childGenome->getName();
                                    nextParentPath.edges.push_back(edge);
                                    // Search recursively for more edges connected to child bottom segment.
                                    vector<HalPath> childPaths = traverseDownRecursive(os, alignment, nextParent, nextParentPath, overlapA, overlapB);
                                    // Gather recursion output.
                                    paths.insert(paths.end(), childPaths.begin(), childPaths.end());
                                }
                            // Move the bottom segment one to the right.
                            nextParent->toRight();
                            if (nextParent->getTopParseIndex() != top->getTopSegment()->getArrayIndex())
                                {
                                    // If the start coordinate of the bottom segment isn't contained in the top,
                                    // then no more overlaps exist.
                                    break;
                                }
                        }
                    if (top->getTopSegment()->hasNextParalogy())
                        {
                            // Move top to next paralogous segment in child genome.
                            top->toNextParalogy();
                            if (top->getTopSegment()->isCanonicalParalog())
                                {
                                    // Only visit the `canonical` paralog once (otherwise we go in a circle).
                                    break;
                                }
                        }
                    else
                        {
                            // If no paralogs exist, move on.
                            break;
                        }
                }
        }
    return paths;
}

void traverseUpRecursive(ostream& os, AlignmentConstPtr alignment,
                         TopSegmentIteratorPtr child,
                         map<string, set<HalEdge> >& edges)
{
    // Function is unfinished--doesn't do anything yet!
}

