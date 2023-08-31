/**
 * generateKVertexCriticalGraphs.c
 * 
 * Author: Jorik Jooken (jorik.jooken@kuleuven.be)
 * 
 */

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "bitset.h"
#include "planarity.h" //Includes nauty.h
#include "read_graph/readGraph6.h"

#define USAGE \
"\nUsage:./genK2Hypohamiltonian [-b] [-p] [-g#] [-d#] [-D#] [-X#] [-e] [-h] n [res/mod]\n\n"

#define HELPTEXT \
"Generate all pairwise non-isomorphic K2-hypohamiltonian graphs of order `n`.\n\
Graphs are sent to stdout in graph6 format. For more information on the format,\n\
see http://users.cecs.anu.edu.au/~bdm/data/formats.txt.\n\
\n\
The `res/mod` argument, should always appear after the specified order `n`.\n\
Otherwise, the order in which the arguments appear does not matter. Be careful\n\
not to put an argument immediately after one with an option. E.g. -g#b will not\n\
recognise the -b argument.\n\
\n\
Mandatory arguments to long options are mandatory for short options too.\n\
    -b, --bipartite             only generate bipartite K2-hypohamiltonian\n\
                                 graphs\n\
    -p, --planar                only generate planar K2-hypohamiltonian graphs\n\
    -g, --girth=GIRTH           only generate K2-hypohamiltonian graphs of\n\
                                 girth at least GIRTH\n\
    -d, --minimum-degree=DEG    only generate K2-hypohamiltonian graphs with\n\
                                 minimum degree at least DEG\n\
    -D, --maximum-degree=DEG    only generate K2-hypohamiltonian graphs with\n\
                                 maximum degree at most DEG\n\
    -X, --splitlevel=LVL        set the splitlevel to LVL; used for\n\
                                 parallellization; a higher level means more\n\
                                 uniformity in running times across parts at \n\
                                 the cost of longer runtime\n\
    -e, --edges-count           print table of how many intermediate graphs had\n\
                                 a certain size\n\
    -h, --help                  print help message\n\
    res/mod                     split the generation in mod (not necessarily\n\
                                 equally big) parts. Here part res will be \n\
                                 executed. Splitting can cause a lot of extra\n\
                                 overhead and duplicate graphs. Make sure to\n\
                                 filter isomorphic graphs afterwards!\n"

//  Used in parallellisation
#define SPLITLEVEL 1

//  Macro's for nauty representation
#define FOREACH(element,nautySet)\
 for(int element = nextelement((nautySet),MAXM,-1); (element) >= 0;\
 (element) = nextelement((nautySet),MAXM,(element)))
#define REMOVEONEEDGE(g,i,j,MAXM)\
 DELELEMENT(GRAPHROW(g,i,MAXM),j); DELELEMENT(GRAPHROW(g,j,MAXM),i)

// Graph structure containing nauty and bitset representations.
struct graph {
    graph nautyGraph[MAXN];
    int numberOfVertices;
    int numberOfEdges;
    bitset* adjacencyList;
    bitset* forbiddenEdges;
    bitset* verticesOfDeg;
}; 

//  Struct for passing options.
struct options {
    int minimalGirth;
    int minimumDegree;
    int maximumDegree;
    int pathLength;
    int nbColors;
    bool planarFlag;
    bool bipartiteFlag;
    int modulo;
    int remainder;
    int splitLevel;
    int splitCounter;
};

//  Struct keeping counters
struct counters {
    long long unsigned int *nOfNonIsoGraphsWithEdges;
    long long unsigned int nOfGraphsFound;
    long long unsigned int nOfTimesCheckedGirth;
    long long unsigned int nOfTimesHadForbiddenGirth;
    long long unsigned int nOfTimesCheckedHamiltonicity;
    long long unsigned int nOfTimesWasHamiltonian;
    long long unsigned int nOfTimesCheckedIsomorphism;
    long long unsigned int nOfTimesWasIsomorphic;
    long long unsigned int nOfTimesCheckedTypeA;
    long long unsigned int nOfTimesContainedTypeA;
    long long unsigned int nOfTimesCheckedTypeB;
    long long unsigned int nOfTimesContainedTypeB;
    long long unsigned int nOfTimesCheckedTypeC;
    long long unsigned int nOfTimesContainedTypeC;
    long long unsigned int nOfTimesCheckedDegreeObstruction;
    long long unsigned int nOfTimesContainedDegreeObstruction;
    long long unsigned int nOfTimesCheckedTriangleObstruction;
    long long unsigned int nOfTimesContainedTriangleObstruction;
    long long unsigned int nOfTimesCheckedSquareObstruction;
    long long unsigned int nOfTimesContainedSquareObstruction;
    long long unsigned int nOfTimesCheckedArrowObstruction;
    long long unsigned int nOfTimesContainedArrowObstruction;
    long long unsigned int nOfTimesCheckedStarObstruction;
    long long unsigned int nOfTimesContainedStarObstruction;
    long long unsigned int nOfTimesNoObstructionApplied;
    long long unsigned int nOfTimesObstructionFound;
    long long unsigned int nOfTimesTypeAObstructionChosen;
    long long unsigned int nOfTimesTypeBObstructionChosen;
    long long unsigned int nOfTimesTypeCObstructionChosen;
    long long unsigned int nOfTimesDegreeObstructionChosen;
    long long unsigned int nOfTimesTriangleObstructionChosen;
    long long unsigned int nOfTimes4CycleObstructionChosen;
    long long unsigned int nOfTimesGeneral4CycleObstructionChosen;
    long long unsigned int nOfTimesGeneralTriangleObstructionChosen;
    long long unsigned int nOfTimesCheckedMaximumDegree;
    long long unsigned int nOfTimesMaximumDegreeExceeded;
    long long unsigned int nOfTimesCheckedPlanarity;
    long long unsigned int nOfTimesWasNonPlanar;

    long long unsigned int *nOfTimesPoorVertex;

    long long unsigned int *nOfTimesSimilarVertices;
    long long unsigned int *nOfTimesSimilarVerticesWithHull;
    long long unsigned int *nOfTimesSimilarEdges;
    long long unsigned int *nOfTimesSimilarEdgesWithHull;
    long long unsigned int *nOfTimesSimilarTriangle;
    long long unsigned int *nOfTimesSimilarP3;
    long long unsigned int *nOfTimesSimilarDiamond;
    long long unsigned int *nOfTimesSimilarP4;
    long long unsigned int *nOfTimesSimilarC4;
    long long unsigned int *nOfTimesSimilarK13;
    long long unsigned int *nOfTimesSimilarK4;
    long long unsigned int *nOfTimesSimilarComplP1PlusP3;
    long long unsigned int *nOfTimesCutVertex;
    long long unsigned int *nOfTimesNoSimilar;

    long long unsigned int *nOfTimesIsomorphismChecked;
    long long unsigned int *nOfNonTerminatingGraphs;
    long long unsigned int *nOfTimesNotC5Colorable;
    long long unsigned int *nOfTimesC5VertexCritical;
};


struct similarGraphsList
{
    int** firstList;
    int** secondList;
    long long size;
    long long capacity;
};

#define MAXFORDER 10

bitset adjacencyListF[MAXFORDER];

char * graphString = NULL;
size_t size;
int nVerticesF;

int arrayToCheckForInducedF[MAXFORDER];
bitset usedInInducedF;

#define INITIALCAPACITY 1
#define MAXSUBGRAPHORDER 7

// state related to C5-colorability
#define MAXVERTICES 64
#define MAXEDGES MAXVERTICES*MAXVERTICES

#define MAXPATHLENGTH 14
#define MAXNUMBERPATHS 20

int pathAsArrayQueue[MAXVERTICES][MAXNUMBERPATHS][MAXPATHLENGTH];
bitset pathAsBitsetQueue[MAXVERTICES][MAXNUMBERPATHS];
int pathQueueStart[MAXVERTICES];
int pathQueueEnd[MAXVERTICES];

int currentPathAsArray[2*MAXPATHLENGTH+5];
bitset currentPathAsBitset;
int currentPathStartIdx;
int currentPathEndIdx;

static bitset oldAdjacencyList[MAXVERTICES];

static bitset oldAvailableColorsAtIteration[MAXVERTICES][MAXVERTICES]; // might be wrong if n is close to MAXVERTICES
static bitset addedToQueueBitsetInIteration[MAXVERTICES]; // might be wrong if n is close to MAXVERTICES

static bitset availableColors[MAXVERTICES];
static int colorOfVertex[MAXVERTICES];
//static const int C5=5; // implement program for C_k colorability
#define C5 5
static bitset adjacentColorBitset[100];

static int vertexQueue[MAXVERTICES];
static int qStart;
static int qFinish;
static bitset canAddToQueueBitset;
static bitset notColoredYetBitset;
static int nbVerticesColored;

// similar graphs
int graph1[MAXSUBGRAPHORDER];
int graph2[MAXSUBGRAPHORDER];
int subgraphAdjacencyMatrix[MAXSUBGRAPHORDER][MAXSUBGRAPHORDER];

// hull
bitset hull[MAXVERTICES];

int max_palette_size[MAXVERTICES];
int colorUsed[100];

// cutvertex
bool articulation_vertex[MAXVERTICES];
int dfs_low[MAXVERTICES];
int dfs_num[MAXVERTICES];
int dfsNumberCounter;
#define DFS_WHITE -1
int dfs_parent[MAXVERTICES];
int dfsRoot;
int rootChildren;

/************************************************************************************
 * 
 * 
 *                      Macro's for dealing with graphs
 * 
 * 
 ************************************************************************************/

//  Initializer for empty graph.
#define emptyGraph(g) EMPTYGRAPH((g)->nautyGraph, (g)->numberOfVertices, MAXM);\
 (g)->verticesOfDeg[0] = compl(EMPTY,(g)->numberOfVertices);\
 (g)->numberOfEdges = 0;\
 for(int i = 0; i < (g)->numberOfVertices; i++) {(g)->adjacencyList[i] = EMPTY;\
 (g)->forbiddenEdges[i] = EMPTY;} 

//  Add one edge. Assume that every edge has already deg2 or more.
#define addEdge(g,i,j) {ADDONEEDGE((g)->nautyGraph, (i), (j), MAXM);\
 add((g)->adjacencyList[i], j); add((g)->adjacencyList[j],i);\
 (g)->numberOfEdges++;\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[i]) - 1], i);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[i])], i);\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[j]) - 1], j);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[j])], j);\
}

//  Remove one edge.
#define removeEdge(g,i,j) {REMOVEONEEDGE((g)->nautyGraph, (i), (j), MAXM);\
 removeElement((g)->adjacencyList[i], j); removeElement((g)->adjacencyList[j],i);\
 (g)->numberOfEdges--;\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[i]) + 1], i);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[i])], i);\
 removeElement((g)->verticesOfDeg[size((g)->adjacencyList[j]) + 1], j);\
 add((g)->verticesOfDeg[size((g)->adjacencyList[j])], j);\
}

//  Add edge to g->forbiddenEdges.
#define forbidEdge(g,i,j) {add((g)->forbiddenEdges[i],j);\
 add((g)->forbiddenEdges[j],i);}

//  Remove edge from g->forbiddenEdges.
#define permitEdge(g,i,j) {removeElement((g)->forbiddenEdges[i],j);\
 removeElement((g)->forbiddenEdges[j],i);}

#define areNeighbours(g,i,j) contains((g)->adjacencyList[(i)], j)


/************************************************************************************
 * 
 * 
 *                      Definitions for Nauty's splay tree
 * 
 * 
 ************************************************************************************/
typedef struct SPLAYNODE {
    graph* canonForm;
    struct SPLAYNODE *left, *right, *parent;
} SPLAYNODE;

SPLAYNODE *splayTreeArray[MAXN*(MAXN-1)/2] = {NULL};

#define SCAN_ARGS

#define ACTION(p) 

#define INSERT_ARGS , graph gCan[], int numberOfVertices,\
 bool *isPresent

int compareSplayNodeToGraph(SPLAYNODE* p, graph gCan[], int numberOfVertices);

#define COMPARE(p) compareSplayNodeToGraph(p, gCan, numberOfVertices);

#define PRESENT(p) {(*isPresent) = true;}

#define NOT_PRESENT(p) {p->canonForm = gCan; (*isPresent) = false;}

#define LOOKUP_ARGS , graph gCan[], int numberOfVertices

#include "splay.c"


/************************************************************************************
 * 
 * 
 *                          Functions for printing data 
 * 
 * 
 ************************************************************************************/

void printGraph(struct graph *g) {
    for(int i = 0; i < g->numberOfVertices; i++) {
        fprintf(stderr, "%d:", i);
        FOREACH(neighbour, GRAPHROW(g->nautyGraph, i, MAXM)) {
            fprintf(stderr, " %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr,"\n");   
}

void printNautyGraph(graph g[], int numberOfVertices) {
    for(int i = 0; i < numberOfVertices; i++) {
        fprintf(stderr, "%d:", i);
        FOREACH(neighbour, GRAPHROW(g, i, MAXM)) {
            fprintf(stderr, " %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr,"\n");   
}

void printAdjacencyList(struct graph *g) {
    for(int i = 0; i < g->numberOfVertices; i++) {
        fprintf(stderr, "%d:", i);
        forEach(neighbour, g->adjacencyList[i]) {
            fprintf(stderr, " %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

void printForbiddenEdges(struct graph *g) {
    for(int i = 0; i< g->numberOfVertices; i++) {
        forEach(nbr, g->forbiddenEdges[i]) {
            fprintf(stderr, " (%d,%d)\n", i, nbr);
        }
    }
}

//  Print gCan to stdout in graph6 format. 
void writeToG6(graph gCan[], int numberOfVertices) {
    char graphString[8 + numberOfVertices*(numberOfVertices - 1)/2];
    int pointer = 0;

    //  Save number of vertices in the first one, four or 8 bytes.
    if(numberOfVertices <= 62) {
        graphString[pointer++] = (char) numberOfVertices + 63;
    }
    else if(numberOfVertices <= 258047) {
        graphString[pointer++] = 63 + 63;
        for(int i = 2; i >= 0; i--) {
            graphString[pointer++] = (char) (numberOfVertices >> i*6) + 63;
        }
    }
    else if(numberOfVertices <= 68719476735) {
        graphString[pointer++] = 63 + 63;
        graphString[pointer++] = 63 + 63;
        for(int i = 5; i >= 0; i--) {
            graphString[pointer++] = (char) (numberOfVertices >> i*6) + 63;
        }
    }
    else {
        fprintf(stderr, "Error: number of vertices too large.\n");
        exit(1);
    }

    // Group upper triangle of adjacency matrix in groups of 6. See B. McKay's 
    // graph6 format.
    int counter = 0;
    char charToPrint = 0;
    for(int i = 1; i < numberOfVertices; i++) {
        for(int j = 0; j < i; j++) {
            charToPrint = charToPrint << 1;
            if(ISELEMENT(GRAPHROW(gCan, i, MAXM), j)) {
                charToPrint |= 1;
            }
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
                charToPrint = 0;
                counter = 0;
            }
        }
    }

    //  Pad final character with 0's.
    if(counter != 0) {
        while(counter < 6) {
            charToPrint = charToPrint << 1;
            if(++counter == 6) {
                graphString[pointer++] = charToPrint + 63;
            }
        }
    }

    //  End with newline and end of string character.
    graphString[pointer++] = '\n';
    graphString[pointer++] = '\0';
    printf("%s", graphString);
}

/************************************************************************************/

void initializeSimilarGraphsList(struct similarGraphsList *sgl)
{
    sgl->size=0;
    sgl->capacity=INITIALCAPACITY;
    sgl->firstList=malloc(INITIALCAPACITY*sizeof(int *));
    sgl->secondList=malloc(INITIALCAPACITY*sizeof(int *));
}

void addSimilarGraphs(struct similarGraphsList *sgl)
{
    if(sgl->size==sgl->capacity)
    {
        sgl->capacity*=2;
        sgl->firstList=realloc(sgl->firstList,sizeof(int *)*sgl->capacity);
        sgl->secondList=realloc(sgl->secondList,sizeof(int *)*sgl->capacity);
    }
    int *newGraph1=malloc(sizeof(int)*MAXSUBGRAPHORDER);
    memcpy(newGraph1,graph1,sizeof(int)*MAXSUBGRAPHORDER);
    sgl->firstList[sgl->size]=newGraph1;

    int *newGraph2=malloc(sizeof(int)*MAXSUBGRAPHORDER);
    memcpy(newGraph2,graph2,sizeof(int)*MAXSUBGRAPHORDER);
    sgl->secondList[sgl->size]=newGraph2;

    sgl->size++;
}

void freeSimilarGraphsList(struct similarGraphsList *sgl)
{
    for(int i=0; i<sgl->size; i++)
    {
        free(sgl->firstList[i]);
        free(sgl->secondList[i]);
    }
    free(sgl->firstList);
    free(sgl->secondList);
}

void printSimilarGraphsList(struct similarGraphsList *sgl)
{
    for(int i=0; i<sgl->size; i++)
    {
        for(int j=0; j<MAXSUBGRAPHORDER; j++)
        {
            fprintf(stderr,"%d ",sgl->firstList[i][j]);
        }
        fprintf(stderr,"\n");
        fprintf(stderr,"Is similar with:\n");
        for(int j=0; j<MAXSUBGRAPHORDER; j++)
        {
            fprintf(stderr,"%d ",sgl->secondList[i][j]);
        }
        fprintf(stderr,"\n");
    }
}

//  Uses Nauty to get canonical form.
void createCanonicalForm(graph g[], graph gCan[], int numberOfVertices) {

    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    statsblk stats;

    densenauty(g, lab, ptn, orbits, &options, &stats, MAXM,
     numberOfVertices, gCan);
}

//  Splaynode contains the canonical form of a graph checked earlier.
int compareSplayNodeToGraph(SPLAYNODE* p, graph gCan[], int numberOfVertices) {
    return memcmp(p->canonForm, gCan, numberOfVertices * sizeof(graph));
}

//  Check recursively whether a path can be extended to a cycle of length
//  smaller than minimalGirth.
bool canBeForbiddenCycle(struct graph *g, bitset nbrsOfStart, int minimalGirth,
 int secondToLast, int last, int pathLength) {
    if(!isEmpty(intersection(nbrsOfStart, g->adjacencyList[last]))
     && pathLength + 1 < minimalGirth) {
        return true;
    }
    if(pathLength + 1 >= minimalGirth) {
        return false;
    }
    forEach(neighbourOfLast,
     difference(g->adjacencyList[last], singleton(secondToLast))) {
        if(canBeForbiddenCycle(g, nbrsOfStart, minimalGirth, last,
         neighbourOfLast, pathLength + 1)) {
            return true;
        }
    }
    return false;
}

int isC5Colorable(struct graph* g, int n, int omittedVertex, int nbColors, int iterationCounter, int highestColorUsed)
{
    if(notColoredYetBitset==0) return 1;
    // search for vertex with fewest number of colors available
    int min_nb_colors_available=nbColors+5;
    int argmin=-1;
    forEach(u,notColoredYetBitset)
    {
        int amount=size(availableColors[u]);
        if(amount<min_nb_colors_available)
        {
            min_nb_colors_available=amount;
            argmin=u;
        }
    }
    // try to color that vertex
    forEach(color,availableColors[argmin])
    {
        if(color>=highestColorUsed+2) return 0;
        //fprintf(stderr,"Trying to color vertex %d with color %d\n",argmin,color);
        bitset oldAvailableColorsArgmin=availableColors[argmin];
        availableColors[argmin]=(1LL<<color);
        notColoredYetBitset=(notColoredYetBitset^(1LL<<argmin));
        // update neighbor availabilites        
        forEach(neigh,(g->adjacencyList[argmin]&notColoredYetBitset))
        {
            oldAvailableColorsAtIteration[iterationCounter][neigh]=availableColors[neigh];
            availableColors[neigh]=(availableColors[neigh]&adjacentColorBitset[color]);
        }
        int newHighestColorUsed=highestColorUsed;
        if(color>highestColorUsed) newHighestColorUsed=color;
        if(isC5Colorable(g,n,omittedVertex,nbColors,iterationCounter+1,newHighestColorUsed)) return 1;
        // restore neighbor availabilities
        forEach(neigh,(g->adjacencyList[argmin]&notColoredYetBitset))
        {
            availableColors[neigh]=oldAvailableColorsAtIteration[iterationCounter][neigh];
        }
        availableColors[argmin]=oldAvailableColorsArgmin;
        notColoredYetBitset=(notColoredYetBitset^(1LL<<argmin));
    }
    return 0;
}

int isC5ColorableHelper(struct graph* g, int n, int omittedVertex, int nbColors)
{
    if(n<=2) return 1;
    for(int i=0;i<nbColors;i++)
    {
        adjacentColorBitset[i]=compl((1LL<<i),nbColors);
    }
    for(int i=0; i<n; i++) availableColors[i]=(1LL<<nbColors)-1;

    int startVertex=0; // change me later; e.g. choose vertex with smallest degree, largest degree or some other heuristic
    if(omittedVertex==0) startVertex=1;

    int bestDegree=n+5;
    int argBest=-1;
    for(int i=0; i<n; i++)
    {
        if(i==omittedVertex) continue;
        int deg=size(g->adjacencyList[i]);
        if(deg<bestDegree)
        {
            bestDegree=deg;
            argBest=i;
        }
    }
    if(argBest!=-1) startVertex=argBest;

    availableColors[startVertex]=(1LL<<0);
    notColoredYetBitset=((1LL<<n)-1)^(1LL<<startVertex);
    if(0<=omittedVertex && omittedVertex<n) notColoredYetBitset=(notColoredYetBitset^(1LL<<omittedVertex));
    forEach(neigh,(g->adjacencyList[startVertex]&notColoredYetBitset))
    {
        availableColors[neigh]=(availableColors[neigh]&adjacentColorBitset[0]);
    }
    int ret=isC5Colorable(g,n,omittedVertex,nbColors,1,0);
    return ret;
}

/*void generateAllC5Colorings(struct graph* g, int n, bitset omittedVertices, int nbColors)
{
    if(notColoredYetBitset==0)
    {
        for(int i=0; i<n; i++)
        {
            if((1LL<<i)&omittedVertices) continue;
            for(int j=0; j<n; j++)
            {
                if(i==j) continue;
                if((1LL<<j)&omittedVertices) continue;
                int c1=colorOfVertex[i];
                int c2=colorOfVertex[j];
                if(!(c1!=c2))
                {
                    hull[i]=(hull[i]&compl((1LL<<j),n));  
                }
            }
        }
        return;
    }
    bitset allExceptOmitted=compl(omittedVertices,n);

    if(qStart>qFinish)
    {
        bitset oldAddedToQueueBitsetInIteration=addedToQueueBitsetInIteration[nbVerticesColored];
        bitset oldCanAddToQueueBitset=canAddToQueueBitset;
        forEach(neigh,canAddToQueueBitset)
        {
            addedToQueueBitsetInIteration[nbVerticesColored]=(addedToQueueBitsetInIteration[nbVerticesColored]|(1LL<<neigh));
            canAddToQueueBitset=(canAddToQueueBitset^(1LL<<neigh));
            vertexQueue[qFinish+1]=neigh;
            qFinish++;
            break;
        }
        generateAllC5Colorings(g,n,omittedVertices,nbColors);
        addedToQueueBitsetInIteration[nbVerticesColored]=oldAddedToQueueBitsetInIteration;
        canAddToQueueBitset=oldCanAddToQueueBitset;
        qFinish--;
        return;
    }
    else
    {
        int currentVertex=vertexQueue[qStart];
    
        notColoredYetBitset=(notColoredYetBitset^(1LL<<currentVertex));
        qStart++;
        nbVerticesColored++;
        forEach(validColor, availableColors[currentVertex])
        {
            // color the vertex
            colorOfVertex[currentVertex]=validColor;
            int recurse=1;
            int oldQFinish=qFinish;
            // iterate over all uncolored neighbors
            int break_neigh=-1;
            forEach(neigh,((g->adjacencyList[currentVertex]&notColoredYetBitset)&allExceptOmitted))
            {
                oldAvailableColorsAtIteration[nbVerticesColored][neigh]=availableColors[neigh];
                availableColors[neigh]=(availableColors[neigh]&adjacentColorBitset[validColor]);
                if(isEmpty(availableColors[neigh])) // backtrack
                {
                    recurse=0;
                    break_neigh=neigh;
                    break;
                }
                if(canAddToQueueBitset&(1LL<<neigh))
                {
                    addedToQueueBitsetInIteration[nbVerticesColored]=(addedToQueueBitsetInIteration[nbVerticesColored]|(1LL<<neigh));
                    canAddToQueueBitset=(canAddToQueueBitset^(1LL<<neigh));
                    vertexQueue[qFinish+1]=neigh;
                    qFinish++;
                }
            }
            if(recurse)
            {
                generateAllC5Colorings(g,n,omittedVertices,nbColors);
            }

            // restore state
            forEach(neigh,((g->adjacencyList[currentVertex]&notColoredYetBitset)&allExceptOmitted))
            {
                availableColors[neigh]=oldAvailableColorsAtIteration[nbVerticesColored][neigh];
                if(neigh==break_neigh) break;
            }
            forEach(neigh,addedToQueueBitsetInIteration[nbVerticesColored])
            {
                canAddToQueueBitset=(canAddToQueueBitset^(1LL<<neigh));
            }
            qFinish=oldQFinish;
            addedToQueueBitsetInIteration[nbVerticesColored]=0;   
        }
        notColoredYetBitset=(notColoredYetBitset^(1LL<<currentVertex));
        qStart--;
        nbVerticesColored--;
        return;
    }
}

void generateAllC5ColoringsHelper(struct graph* g, int n, bitset omittedVertices, int nbColors)
{
    for(int i=0;i<nbColors;i++)
    {
        adjacentColorBitset[i]=compl((1LL<<i),nbColors);
    }
    for(int i=0; i<n; i++) availableColors[i]=(1LL<<nbColors)-1;
    for(int i=0; i<n+1; i++)
    {
        addedToQueueBitsetInIteration[i]=0LL;
    }
    qStart=0;
    qFinish=0;
    int startVertex=0; // change me later; e.g. choose vertex with smallest degree, largest degree or some other heuristic
    while(true)
    {
        if(omittedVertices&(1LL<<startVertex))
        {
            startVertex++;
        }
        else break;
    }
    vertexQueue[0]=startVertex;
    availableColors[startVertex]=(1LL<<0);
    nbVerticesColored=0;
    notColoredYetBitset=(1LL<<n)-1;
    notColoredYetBitset=(notColoredYetBitset^omittedVertices);
    canAddToQueueBitset=notColoredYetBitset;
    canAddToQueueBitset=(canAddToQueueBitset^(1LL<<startVertex));
    return generateAllC5Colorings(g,n,omittedVertices,nbColors);
}*/

void generateAllC5Colorings(struct graph* g, int n, bitset omittedVertices, int nbColors, int iterationCounter, int calcPoorVertex, int highestColorUsed)
{
    if(notColoredYetBitset==0)
    {
        if(calcPoorVertex==1)
        {
            int omitted=-1;
            forEach(el,omittedVertices)
            {
                omitted=el;
                break;
            }
            for(int i=0; i<nbColors; i++)
            {
                colorUsed[i]=0;
            }
            forEach(neigh,g->adjacencyList[omitted])
            {
                colorUsed[colorOfVertex[neigh]]+=1;
            }
            int nonZero=0;
            for(int i=0; i<nbColors; i++)
            {
                if(colorUsed[i]!=0) nonZero+=1;
            }
            if(nonZero>max_palette_size[omitted])
            {
                max_palette_size[omitted]=nonZero;
            }
        }
        else
        {
            for(int i=0; i<n; i++)
            {
                if((1LL<<i)&omittedVertices) continue;
                for(int j=0; j<n; j++)
                {
                    if(i==j) continue;
                    if((1LL<<j)&omittedVertices) continue;
                    int c1=colorOfVertex[i];
                    int c2=colorOfVertex[j];
                    if(!(c1!=c2))
                    {
                        hull[i]=(hull[i]&compl((1LL<<j),n));  
                    }
                }
            }
        }
        return;
    }

    // search for vertex with fewest number of colors available
    /*int min_nb_colors_available=nbColors+5;
    int argmin=-1;
    forEach(u,notColoredYetBitset)
    {
        int amount=size(availableColors[u]);
        if(amount<min_nb_colors_available)
        {
            min_nb_colors_available=amount;
            argmin=u;
        }
    }*/

    /*
    // search for vertex with most number of colors available
    int min_nb_colors_available=-1;
    int argmin=-1;
    forEach(u,notColoredYetBitset)
    {
        int amount=size(availableColors[u]);
        if(amount>min_nb_colors_available)
        {
            min_nb_colors_available=amount;
            argmin=u;
        }
    }*/

    /*// search for vertex with most colored neighbors
    int most_colored_neighbors=-1;
    int argmin=-1; // I call it argmin so I do not have to rename this variable below, but it is actually an argmax!
    forEach(u,notColoredYetBitset)
    {
        int amount=size(g->adjacencyList[u]&(compl(notColoredYetBitset,n)));
        if(amount>most_colored_neighbors)
        {
            most_colored_neighbors=amount;
            argmin=u;
        }
    }*/

    // search for vertex with most colored neighbors
    // use the number of available colours as a tiebreak
    int most_colored_neighbors=-1;
    int min_nb_colors_available=nbColors+5;
    int argmin=-1; // I call it argmin so I do not have to rename this variable below, but it is actually an argmax!
    forEach(u,notColoredYetBitset)
    {
        int amount=size(g->adjacencyList[u]&(compl(notColoredYetBitset,n)));
        int amount2=size(availableColors[u]);
        if(amount>most_colored_neighbors || (amount==most_colored_neighbors && amount2<min_nb_colors_available))
        {
            most_colored_neighbors=amount;
            min_nb_colors_available=amount2;
            argmin=u;
        }
    }
    
    // try to color that vertex
    forEach(color,availableColors[argmin])
    {
        if(color>=highestColorUsed+2) return;
        bitset oldAvailableColorsArgmin=availableColors[argmin];
        availableColors[argmin]=(1LL<<color);
        colorOfVertex[argmin]=color;
        notColoredYetBitset=(notColoredYetBitset^(1LL<<argmin));
        // update neighbor availabilites        
        forEach(neigh,(g->adjacencyList[argmin]&notColoredYetBitset))
        {
            oldAvailableColorsAtIteration[iterationCounter][neigh]=availableColors[neigh];
            availableColors[neigh]=(availableColors[neigh]&adjacentColorBitset[color]);
        }
        int newHighestColorUsed=highestColorUsed;
        if(color>highestColorUsed) newHighestColorUsed=color;
        generateAllC5Colorings(g,n,omittedVertices,nbColors,iterationCounter+1,calcPoorVertex,newHighestColorUsed);
        if(calcPoorVertex==1)
        {
            int omitted=-1;
            forEach(el,omittedVertices)
            {
                omitted=el;
                break;
            }
            int degOmitted=size(g->adjacencyList[omitted]);
            // max_palette_size cannot increase further
            if(max_palette_size[omitted]==nbColors || max_palette_size[omitted]==degOmitted)
            {
                return;
            }
        }
        // restore neighbor availabilities
        forEach(neigh,(g->adjacencyList[argmin]&notColoredYetBitset))
        {
            availableColors[neigh]=oldAvailableColorsAtIteration[iterationCounter][neigh];
        }
        availableColors[argmin]=oldAvailableColorsArgmin;
        colorOfVertex[argmin]=-1;
        notColoredYetBitset=(notColoredYetBitset^(1LL<<argmin));
    }
}

void generateAllC5ColoringsHelper(struct graph* g, int n, bitset omittedVertices, int nbColors, int calcPoorVertex)
{
    /*if(n==10 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0)
    && ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))>0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))>0) && ((g->adjacencyList[9]&(1LL<<6))==0) && ((g->adjacencyList[9]&(1LL<<7))>0) && ((g->adjacencyList[9]&(1LL<<8))==0))
        {
            fprintf(stderr,"entered here!\n");
            fprintf(stderr,"omittedVertices: %ld\n",omittedVertices);
            fprintf(stderr,"calcPoorVertex: %d\n",calcPoorVertex);
        }*/
    if(calcPoorVertex==1)
    {
        for(int i=0; i<n; i++)
        {
            max_palette_size[i]=-1;
        }
    }
    for(int i=0;i<nbColors;i++)
    {
        adjacentColorBitset[i]=compl((1LL<<i),nbColors);
    }
    for(int i=0; i<n; i++) availableColors[i]=(1LL<<nbColors)-1;
    int startVertex=0; // change me later; e.g. choose vertex with smallest degree, largest degree or some other heuristic
    while(true)
    {
        if(omittedVertices&(1LL<<startVertex))
        {
            startVertex++;
        }
        else break;
    }
    availableColors[startVertex]=(1LL<<0);
    colorOfVertex[startVertex]=0;
    notColoredYetBitset=(1LL<<n)-1;
    notColoredYetBitset=(notColoredYetBitset^omittedVertices);
    notColoredYetBitset=(notColoredYetBitset^(1LL<<startVertex));
    forEach(neigh,(g->adjacencyList[startVertex]&notColoredYetBitset))
    {
        availableColors[neigh]=(availableColors[neigh]&adjacentColorBitset[0]);
    }
    return generateAllC5Colorings(g,n,omittedVertices,nbColors,2,calcPoorVertex,0);
}

int can_add_vertices_right(int amount, struct graph *g, int n, bitset *allowedOtherSide, int currVertex)
{
    if(amount==0)
    {
        // add current path to queue

        // queue full, remove first element
        if(pathQueueEnd[n]!=-1 && ((pathQueueEnd[n]+1)%MAXNUMBERPATHS==pathQueueStart[n]))
        {
            pathQueueStart[n]=(pathQueueStart[n]+1)%MAXNUMBERPATHS;
        }
        pathQueueEnd[n]=(pathQueueEnd[n]+1)%MAXNUMBERPATHS;
        pathAsBitsetQueue[n][pathQueueEnd[n]]=currentPathAsBitset;
        for(int i=currentPathStartIdx; i<=currentPathEndIdx; i++)
        {
            pathAsArrayQueue[n][pathQueueEnd[n]][i-currentPathStartIdx]=currentPathAsArray[i];
        }
        return 1;
    }
    bitset allowedNeighbours=intersection(*allowedOtherSide,g->adjacencyList[currVertex]);
    forEach(neighbour,allowedNeighbours)
    {
        bitset oldOtherSide=(*allowedOtherSide);
        (*allowedOtherSide)=intersection(*allowedOtherSide,compl(g->adjacencyList[currVertex],n-1));
        currentPathAsBitset=(currentPathAsBitset^(1LL<<neighbour));
        currentPathEndIdx++;
        currentPathAsArray[currentPathEndIdx]=neighbour;
        if(can_add_vertices_right(amount-1, g, n, allowedOtherSide, neighbour)) return 1;
        (*allowedOtherSide)=oldOtherSide;
        currentPathAsBitset=(currentPathAsBitset^(1LL<<neighbour));
        currentPathEndIdx--;
    }
    return 0;
}

int can_add_vertices_left(int amount, struct graph *g, int n, bitset *allowedOneSide, bitset *allowedOtherSide, int currVertex)
{
    // option 1: add vertices to the other side
    if(currVertex != n-1)
    {
        (*allowedOtherSide)=intersection(*allowedOtherSide,compl(g->adjacencyList[currVertex],n-1));
    }
    if(can_add_vertices_right(amount,g,n,allowedOtherSide,n-1)) return 1;
    
	// option 2: add vertices to this side
	int neighbour = -1;
	bitset allowedNeighbours=intersection(*allowedOneSide,g->adjacencyList[currVertex]);
	forEach(neighbour,allowedNeighbours)
	{
	    bitset oldOneSide=(*allowedOneSide);
	    bitset oldOtherSide=(*allowedOtherSide);
	    (*allowedOneSide)=intersection(*allowedOneSide,compl(g->adjacencyList[currVertex],n-1));
	    if(currVertex!=n-1)
	    {
		(*allowedOtherSide)=intersection(*allowedOtherSide,compl(g->adjacencyList[currVertex],n-1));
	    }
        currentPathAsBitset=(currentPathAsBitset^(1LL<<neighbour));
        currentPathStartIdx--;
        currentPathAsArray[currentPathStartIdx]=neighbour;
	    if(can_add_vertices_left(amount-1, g, n, allowedOneSide, allowedOtherSide, neighbour)) return 1;
	    (*allowedOneSide)=oldOneSide;
	    (*allowedOtherSide)=oldOtherSide;
        currentPathAsBitset=(currentPathAsBitset^(1LL<<neighbour));
        currentPathStartIdx++;
	}
    return 0;
}

int has_induced_path_of_length(int len, struct graph *g, int currVertex, struct options *options)
{
    if(len>currVertex+1) return 0;
    // first check if any of the saved paths are still induced paths or not
    if(pathQueueEnd[currVertex+1]!=-1) // queue is non-empty
    {
        bool loop_already_entered=false;
        for(int i=pathQueueStart[currVertex+1]; !loop_already_entered || (((i-1+MAXNUMBERPATHS)%MAXNUMBERPATHS)!= pathQueueEnd[currVertex+1]); i=(i+1)%MAXNUMBERPATHS)
        {
            loop_already_entered=true;
            bool is_induced_path=true;
            for(int j=0; j<options->pathLength && is_induced_path; j++)
            {
                int vertex=pathAsArrayQueue[currVertex+1][i][j];
                bitset shouldBeNeighbours=EMPTY;
                if(j-1>=0) shouldBeNeighbours=(shouldBeNeighbours|(1LL<<pathAsArrayQueue[currVertex+1][i][j-1]));
                if(j+1<options->pathLength) shouldBeNeighbours=(shouldBeNeighbours|(1LL<<pathAsArrayQueue[currVertex+1][i][j+1]));
                if(((g->adjacencyList[vertex])&(pathAsBitsetQueue[currVertex+1][i]))!=shouldBeNeighbours)
                {
                    is_induced_path=false;
                }
            }
            if(is_induced_path)
            {
                return 1;
            }
        }
    }    

    bitset allowedOneSide=compl(EMPTY,currVertex);
    bitset allowedOtherSide=compl(EMPTY,currVertex);

    currentPathAsBitset=singleton(currVertex);
    currentPathStartIdx=MAXPATHLENGTH+1;
    currentPathEndIdx=MAXPATHLENGTH+1;
    currentPathAsArray[currentPathStartIdx]=currVertex;

    return can_add_vertices_left(len-1,g,currVertex+1,&allowedOneSide,&allowedOtherSide,currVertex);
}

int recursive_check_has_induced_F(struct graph *g, int currVertex, int alreadyPlacedIdx, int currIdx)
{
    if(currIdx==alreadyPlacedIdx) return recursive_check_has_induced_F(g,currVertex,alreadyPlacedIdx,currIdx+1);
    if(currIdx==nVerticesF) return 1;
    forEach(u,compl(usedInInducedF,currVertex+1))
    {
        if(size(g->adjacencyList[u])<size(adjacencyListF[currIdx])) continue;
        int preservesAllAdjacencies=1;
        for(int i=0; i<currIdx; i++)
        {
            int v=arrayToCheckForInducedF[i];

            int u_prime=currIdx;
            int v_prime=i;

            if(g->adjacencyList[u]&(1LL<<v))
            {
                if(adjacencyListF[u_prime]&(1LL<<v_prime))
                {
                }
                else
                {
                    preservesAllAdjacencies=0;
                    break;
                }
            }
            else
            {
                if(adjacencyListF[u_prime]&(1LL<<v_prime))
                {                    
                    preservesAllAdjacencies=0;
                    break;
                }
            }
        }
        // also compare with element at alreadyPlacedIdx
        int v=arrayToCheckForInducedF[alreadyPlacedIdx];

        int u_prime=currIdx;
        int v_prime=alreadyPlacedIdx;

        if(g->adjacencyList[u]&(1LL<<v))
        {
            if(adjacencyListF[u_prime]&(1LL<<v_prime))
            {
            }
            else
            {
                preservesAllAdjacencies=0;
            }
        }
        else
        {
            if(adjacencyListF[u_prime]&(1LL<<v_prime))
            {                    
                preservesAllAdjacencies=0;
            }
        }

        if(preservesAllAdjacencies==1)
        {
            arrayToCheckForInducedF[currIdx]=u;
            usedInInducedF=(usedInInducedF^(1LL<<u));
            if(recursive_check_has_induced_F(g,currVertex,alreadyPlacedIdx,currIdx+1)) 
            {
                arrayToCheckForInducedF[currIdx]=-1;
                usedInInducedF=(usedInInducedF^(1LL<<u));
                return 1;
            }
            usedInInducedF=(usedInInducedF^(1LL<<u));
            arrayToCheckForInducedF[currIdx]=-1;
        }
    }
    return 0;
}

int has_induced_F(struct graph *g, int currVertex)
{
    if(nVerticesF>currVertex+1) return 0;
    for(int i=0; i<nVerticesF; i++)
    {
        if(size(g->adjacencyList[currVertex])<size(adjacencyListF[i])) continue;
        arrayToCheckForInducedF[i]=currVertex;
        usedInInducedF=(1LL<<currVertex);
        if(recursive_check_has_induced_F(g,currVertex,i,0)) 
        {
            arrayToCheckForInducedF[i]=-1;
            return 1;
        }
        arrayToCheckForInducedF[i]=-1;
    }
    return 0;    
}

// complexity O(|V|+|E|)
// stops as soon as some cutvertex was found (if all cutvertices are needed, function needs to be rewritten!)
bool articulationPointAndBridge(struct graph* g, int u)
{
    dfs_low[u] = dfs_num[u] = dfsNumberCounter++; //dfs_low[u] <= dfs_num[u]
    forEach(v,g->adjacencyList[u])
    {
       if(dfs_num[v] == DFS_WHITE) //tree edge
       {
            dfs_parent[v] = u;
            if(u==dfsRoot) rootChildren++;
            bool alreadyFound=articulationPointAndBridge(g, v);
            if(alreadyFound) return true;
            if(dfs_low[v] >= dfs_num[u]) //articulation point
            {
                articulation_vertex[u] = true;
                return true;
            }
            /*if(dfs_low[v->first] > dfs_num[u]) // articulation bridge
            {
               printf("Edge %d %d is a bridge\n", u, v->first);
            }*/
            if(dfs_low[v]<dfs_low[u])
            {
                dfs_low[u] = dfs_low[v]; 
            }              
       }
       else if(v != dfs_parent[u])
       {
            if(dfs_num[v]<dfs_low[u])
            {
                dfs_low[u]=dfs_num[v];
            }
       }
    }    
}

void addVertex(struct graph* g,  struct similarGraphsList* sgl, int currVertex, int numberOfVertices, int *ctr, struct options *options, struct counters *counters);

void compute_hull(struct graph* g, int currVertex, bitset omittedVertices, int nbColors)
{
    for(int i=0; i<currVertex; i++)
    {
        // edge between every pair of nodes and no self-loops
        hull[i]=(((1LL<<currVertex)-1)^(1LL<<i));
        // no edges incident with omitted vertices
        forEach(omittedVertex,omittedVertices)
        {
            if(i==omittedVertex)
            {
                hull[i]=0;
                break;
            }
            hull[i]=(hull[i]^(1LL<<omittedVertex));
        }
    }
    generateAllC5ColoringsHelper(g,currVertex,omittedVertices,nbColors,0);
}


void addEdgesInAllPossibleWays(struct graph* g,  struct similarGraphsList* sgl, int currVertex, int numberOfVertices, int *ctr, struct options *options, int otherEndPointAfter, bitset *validEndPoints, int numberSimilarVertices, struct counters *counters, int poor_vertex, int max_palette_size_G_minus_u)
{
    /*if(currVertex+1==6)
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
    }*/
    /*if(currVertex+1==7 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0))
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
    }*/
    /*if(currVertex+1==7 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0))
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
        fprintf(stderr,"numSimVert: %d\n",numberSimilarVertices);
        fprintf(stderr,"poor_vertex: %d\n",poor_vertex);
    }*/
    /*if(currVertex+1==8 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0))
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
    }*/

    /*if(currVertex+1==9 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0))
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
    }*/

    /*if(currVertex+1==10 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0))
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
    }*/

    /*if(currVertex+1==11 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0)
    && ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))>0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))>0) && ((g->adjacencyList[9]&(1LL<<6))==0) && ((g->adjacencyList[9]&(1LL<<7))>0) && ((g->adjacencyList[9]&(1LL<<8))==0))
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
    }*/

    /*if(currVertex+1==10 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0)
    && ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))>0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))>0) && ((g->adjacencyList[9]&(1LL<<6))==0) && ((g->adjacencyList[9]&(1LL<<7))>0) && ((g->adjacencyList[9]&(1LL<<8))==0))
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
        fprintf(stderr,"numSimVert: %d\n",numberSimilarVertices);
        fprintf(stderr,"poor_vertex: %d\n",poor_vertex);
    }*/

    /*if(currVertex+1==10 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0)
    && ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))==0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))==0) && ((g->adjacencyList[9]&(1LL<<6))>0) && ((g->adjacencyList[9]&(1LL<<7))==0) && ((g->adjacencyList[9]&(1LL<<8))>0))
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
        fprintf(stderr,"numSimVert: %d\n",numberSimilarVertices);
        fprintf(stderr,"poor_vertex: %d\n",poor_vertex);
    }*/

    /*if(currVertex+1==11 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0)
    && ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))==0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))==0) && ((g->adjacencyList[9]&(1LL<<6))>0) && ((g->adjacencyList[9]&(1LL<<7))>0) && ((g->adjacencyList[9]&(1LL<<8))>0))
        {
            printNautyGraph(g->nautyGraph,currVertex+1);
            fprintf(stderr,"numSimVert: %d\n",numberSimilarVertices);
            fprintf(stderr,"poor_vertex: %d\n",poor_vertex);
        }*/

    /*if(currVertex+1==10 && ((g->adjacencyList[1]&(1LL<<0))>0) && 
        ((g->adjacencyList[2]&(1LL<<0))==0) && ((g->adjacencyList[2]&(1LL<<1))>0) && 
        ((g->adjacencyList[3]&(1LL<<0))==0) && ((g->adjacencyList[3]&(1LL<<1))==0)
      && ((g->adjacencyList[3]&(1LL<<2))>0) && 
        ((g->adjacencyList[4]&(1LL<<0))>0) && ((g->adjacencyList[4]&(1LL<<1))==0) && ((g->adjacencyList[4]&(1LL<<2))==0) && ((g->adjacencyList[4]&(1LL<<3))>0) && 
        ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0) 
    && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<4))==0) && ((g->adjacencyList[6]&(1LL<<5))>0) 
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0))
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
        fprintf(stderr,"numSimVert: %d\n",numberSimilarVertices);
        fprintf(stderr,"poor_vertex: %d\n",poor_vertex);
    }*/

    // TODO: try changing the position of this check
    bool has_path=false;
    if(has_induced_path_of_length(options->pathLength,g,currVertex,options))
    {
        has_path=true; // graph should be Pt-free
    }
    bool has_F=false;
    if(!has_path && has_induced_F(g,currVertex))
    {
        has_F=true; // graph should be F-free
    }


    int nbColors=options->nbColors;
    bool graphColorable=true;
    bool allInducedSubgraphsColorable=true;
    // option 1 : graph is finished
    /*if(!isC5ColorableHelper(g,currVertex+1,123,nbColors)) // graph cannot be colored
    {
        graphColorable=false;
        counters->nOfTimesNotC5Colorable[currVertex]++;
        for(int omit=0; omit<currVertex && allInducedSubgraphsColorable; omit++)
        {
            if(!isC5ColorableHelper(g,currVertex+1,omit,nbColors)) allInducedSubgraphsColorable=false;
        }
    }*/ 

    if(!has_path && !has_F && !isC5ColorableHelper(g,currVertex+1,123,nbColors)) // graph cannot be colored
    {
        graphColorable=false;
        counters->nOfTimesNotC5Colorable[currVertex]++;
        for(int omit=0; omit<currVertex && allInducedSubgraphsColorable; omit++)
        {
            if(!isC5ColorableHelper(g,currVertex+1,omit,nbColors)) allInducedSubgraphsColorable=false;
        }
    }

    bool isolated=false;
    if(isEmpty(g->adjacencyList[currVertex]))
    {
        isolated=true; // graph should be connected
    }

    /*if(!isPresent && !isolated && !graphColorable && allInducedSubgraphsColorable && !has_F && !has_path)
    {
        counters->nOfTimesC5VertexCritical[currVertex]++;
        //writeToG6(g->nautyGraph,currVertex+1);
        //printNautyGraph(g->nautyGraph,currVertex+1);
        (*ctr)++;
    }*/

    //if(currVertex+1==10) isPresent=false; // remove me later!

    bool structureBroken1=false;
    if(poor_vertex<0)
    {
        structureBroken1=false;
        for(int i=0; i<numberSimilarVertices && !structureBroken1; i++)
        {
            if((g->adjacencyList[currVertex]&(1LL<<sgl->firstList[sgl->size-1][i]))>0 && (g->adjacencyList[currVertex]&(1LL<<sgl->secondList[sgl->size-1][i]))==0)
            {
                structureBroken1=true;
            }
        }

        if(numberSimilarVertices<=1) structureBroken1=true;
        // check that structure is really broken, there can't be a hidden edge between the newly added vertex and the non-edge similar vertex
        // TODO: here we compute the hull, which is quite expensive. Perhaps we can place this check at another point.
        ///if(structureBroken1 && !isolated && graphColorable && ! has_F && !has_path && !isPresent)
        if(structureBroken1 && !isolated && graphColorable && ! has_F && !has_path)
        {
            if(numberSimilarVertices>=1)
            {
                bitset omittedVertices=0;
                for(int i=0; i<numberSimilarVertices; i++)
                {
                    omittedVertices=(omittedVertices|(1LL<<(sgl->firstList[sgl->size-1][i])));
                }
                compute_hull(g,currVertex+1,omittedVertices,nbColors);
                int someOk=0;
                for(int i=0; i<numberSimilarVertices && (someOk==0); i++)
                {
                    if(((g->adjacencyList[sgl->firstList[sgl->size-1][i]])&(1LL<<currVertex))>0 && (((hull[sgl->secondList[sgl->size-1][i]])&(1LL<<currVertex))==0)) someOk=1;
                }
                if(someOk==0) structureBroken1=false;
            }
        }
    }

    bool structureBroken2=true;
    //if(poor_vertex>=0 && !isolated && graphColorable && ! has_F && !has_path && !isPresent)
    if(poor_vertex>=0 && !isolated && graphColorable && ! has_F && !has_path)
    {
        bitset omittedVertices=(1LL<<poor_vertex);
        generateAllC5ColoringsHelper(g,currVertex+1,omittedVertices,nbColors,1);
        int max_palette_size_H_minus_u=max_palette_size[poor_vertex];
        /*if(currVertex+1==10 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0)
    && ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))>0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))>0) && ((g->adjacencyList[9]&(1LL<<6))==0) && ((g->adjacencyList[9]&(1LL<<7))>0) && ((g->adjacencyList[9]&(1LL<<8))==0))
        {
            fprintf(stderr,"max palette size H minus u: %d\n",max_palette_size_H_minus_u);
            fprintf(stderr,"max palette size G minus u: %d\n",max_palette_size_G_minus_u);
        }*/
        /*if(currVertex+1==10 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0)
    && ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))==0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))==0) && ((g->adjacencyList[9]&(1LL<<6))>0) && ((g->adjacencyList[9]&(1LL<<7))==0) && ((g->adjacencyList[9]&(1LL<<8))>0))
        {
            fprintf(stderr,"max palette size H minus u: %d\n",max_palette_size_H_minus_u);
            fprintf(stderr,"max palette size G minus u: %d\n",max_palette_size_G_minus_u);
        }*/

        /*if(currVertex+1==10 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0)
    && ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))>0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))>0) && ((g->adjacencyList[9]&(1LL<<6))==0) && ((g->adjacencyList[9]&(1LL<<7))>0) && ((g->adjacencyList[9]&(1LL<<8))==0))
        {
            fprintf(stderr,"poor_vertex: %d\n",poor_vertex);
            fprintf(stderr,"max palette size H minus u: %d\n",max_palette_size_H_minus_u);
            fprintf(stderr,"max palette size G minus u: %d\n",max_palette_size_G_minus_u);
        }*/

        // palette size did not increase
        if(max_palette_size_H_minus_u<=max_palette_size_G_minus_u) structureBroken2=false;
    }

    
    //fprintf(stderr,"isPresent: %d\n", isPresent);

    /*if(currVertex+1==10 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0)
    && ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))>0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))>0) && ((g->adjacencyList[9]&(1LL<<6))==0) && ((g->adjacencyList[9]&(1LL<<7))>0) && ((g->adjacencyList[9]&(1LL<<8))==0))
        {
            printNautyGraph(g->nautyGraph,currVertex+1);
            fprintf(stderr,"isPresent %d\n",isPresent);
            fprintf(stderr,"isolated %d\n",isolated);
            fprintf(stderr,"graphColorable %d\n",graphColorable);
            fprintf(stderr,"structureBroken2 %d\n",structureBroken2);
        }*/

    /*if(currVertex+1==10 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0)
    && ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))==0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))==0) && ((g->adjacencyList[9]&(1LL<<6))>0) && ((g->adjacencyList[9]&(1LL<<7))==0) && ((g->adjacencyList[9]&(1LL<<8))>0))
        {
            printNautyGraph(g->nautyGraph,currVertex+1);
            fprintf(stderr,"isPresent %d\n",isPresent);
            fprintf(stderr,"isolated %d\n",isolated);
            fprintf(stderr,"graphColorable %d\n",graphColorable);
            fprintf(stderr,"structureBroken2 %d\n",structureBroken2);
            fprintf(stderr,"has_path %d\n",has_path);
            fprintf(stderr,"has_F %d\n",has_F);
        }*/

        /*if(currVertex+1==10 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0)
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0)
    && ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))==0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))==0) && ((g->adjacencyList[9]&(1LL<<6))>0) && ((g->adjacencyList[9]&(1LL<<7))>0) && ((g->adjacencyList[9]&(1LL<<8))>0))
        {
            printNautyGraph(g->nautyGraph,currVertex+1);
            fprintf(stderr,"isPresent %d\n",isPresent);
            fprintf(stderr,"isolated %d\n",isolated);
            fprintf(stderr,"graphColorable %d\n",graphColorable);
            fprintf(stderr,"structureBroken2 %d\n",structureBroken2);
            fprintf(stderr,"has_path %d\n",has_path);
            fprintf(stderr,"has_F %d\n",has_F);
        }*/

    /*if(isTheOne && entered)
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
        fprintf(stderr,"isPresent %d\n",isPresent);
        fprintf(stderr,"isolated %d\n",isolated);
        fprintf(stderr,"graphColorable %d\n",graphColorable);
        fprintf(stderr,"structureBroken1 %d\n",structureBroken1);
        fprintf(stderr,"structureBroken2 %d\n",structureBroken2);
        fprintf(stderr,"has_path %d\n",has_path);
        fprintf(stderr,"has_F %d\n",has_F);
        fprintf(stderr,"numSimVert: %d\n",numberSimilarVertices);
        fprintf(stderr,"poor_vertex: %d\n",poor_vertex);
    }*/

    /*if(currVertex+1==10 && ((g->adjacencyList[1]&(1LL<<0))>0) && 
        ((g->adjacencyList[2]&(1LL<<0))==0) && ((g->adjacencyList[2]&(1LL<<1))>0) && 
        ((g->adjacencyList[3]&(1LL<<0))==0) && ((g->adjacencyList[3]&(1LL<<1))==0)
      && ((g->adjacencyList[3]&(1LL<<2))>0) && 
        ((g->adjacencyList[4]&(1LL<<0))>0) && ((g->adjacencyList[4]&(1LL<<1))==0) && ((g->adjacencyList[4]&(1LL<<2))==0) && ((g->adjacencyList[4]&(1LL<<3))>0) && 
        ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0) 
    && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<4))==0) && ((g->adjacencyList[6]&(1LL<<5))>0) 
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0) &&
    ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))>0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))>0) && ((g->adjacencyList[9]&(1LL<<6))==0) && ((g->adjacencyList[9]&(1LL<<7))>0) && ((g->adjacencyList[9]&(1LL<<8))==0))
    {

        printNautyGraph(g->nautyGraph,currVertex+1);
        fprintf(stderr,"isPresent %d\n",isPresent);
        fprintf(stderr,"isolated %d\n",isolated);
        fprintf(stderr,"graphColorable %d\n",graphColorable);
        fprintf(stderr,"structureBroken1 %d\n",structureBroken1);
        fprintf(stderr,"structureBroken2 %d\n",structureBroken2);
        fprintf(stderr,"has_path %d\n",has_path);
        fprintf(stderr,"has_F %d\n",has_F);
        fprintf(stderr,"numSimVert: %d\n",numberSimilarVertices);
        fprintf(stderr,"poor_vertex: %d\n",poor_vertex);
    }*/

    /*if(currVertex+1==10 && ((g->adjacencyList[1]&(1LL<<0))>0) && 
        ((g->adjacencyList[2]&(1LL<<0))==0) && ((g->adjacencyList[2]&(1LL<<1))>0) && 
        ((g->adjacencyList[3]&(1LL<<0))==0) && ((g->adjacencyList[3]&(1LL<<1))==0)
      && ((g->adjacencyList[3]&(1LL<<2))>0) && 
        ((g->adjacencyList[4]&(1LL<<0))>0) && ((g->adjacencyList[4]&(1LL<<1))==0) && ((g->adjacencyList[4]&(1LL<<2))==0) && ((g->adjacencyList[4]&(1LL<<3))>0) && 
        ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0) 
    && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<4))==0) && ((g->adjacencyList[6]&(1LL<<5))>0) 
    && ((g->adjacencyList[7]&(1LL<<0))==0) && ((g->adjacencyList[7]&(1LL<<1))>0) && ((g->adjacencyList[7]&(1LL<<2))>0) && ((g->adjacencyList[7]&(1LL<<3))>0) && ((g->adjacencyList[7]&(1LL<<4))>0) && ((g->adjacencyList[7]&(1LL<<5))>0) && ((g->adjacencyList[7]&(1LL<<6))==0) && ((g->adjacencyList[8]&(1LL<<0))>0) && ((g->adjacencyList[8]&(1LL<<1))>0) && ((g->adjacencyList[8]&(1LL<<2))>0) && ((g->adjacencyList[8]&(1LL<<3))>0) && ((g->adjacencyList[8]&(1LL<<4))==0) && ((g->adjacencyList[8]&(1LL<<5))>0) && ((g->adjacencyList[8]&(1LL<<6))>0) && ((g->adjacencyList[8]&(1LL<<7))==0) &&
    ((g->adjacencyList[9]&(1LL<<0))==0) && ((g->adjacencyList[9]&(1LL<<1))==0) && ((g->adjacencyList[9]&(1LL<<2))>0) && ((g->adjacencyList[9]&(1LL<<3))>0) && ((g->adjacencyList[9]&(1LL<<4))>0) && ((g->adjacencyList[9]&(1LL<<5))==0) && ((g->adjacencyList[9]&(1LL<<6))>0) && ((g->adjacencyList[9]&(1LL<<7))>0) && ((g->adjacencyList[9]&(1LL<<8))>0))
    {
        printNautyGraph(g->nautyGraph,currVertex+1);
        fprintf(stderr,"isPresent %d\n",isPresent);
        fprintf(stderr,"isolated %d\n",isolated);
        fprintf(stderr,"graphColorable %d\n",graphColorable);
        fprintf(stderr,"structureBroken1 %d\n",structureBroken1);
        fprintf(stderr,"structureBroken2 %d\n",structureBroken2);
        fprintf(stderr,"has_path %d\n",has_path);
        fprintf(stderr,"has_F %d\n",has_F);
        fprintf(stderr,"numSimVert: %d\n",numberSimilarVertices);
        fprintf(stderr,"poor_vertex: %d\n",poor_vertex);
    }*/
    //if(!isPresent &&(!isolated || currVertex==0) && graphColorable && ((structureBroken1 && poor_vertex<0) || (structureBroken2 && poor_vertex>=0) || poor_vertex==-2) && !has_F && !has_path) {
    if((!isolated || currVertex==0) && ((structureBroken1 && poor_vertex<0) || (structureBroken2 && poor_vertex>=0) || poor_vertex==-2) && !has_F && !has_path) {

        bool isPresent=false;
        graph* gCan = malloc(sizeof(graph)*(currVertex+1));
        if(gCan == NULL) {
            fprintf(stderr, "Error: out of memory\n");
            exit(1);
        }
        counters->nOfTimesIsomorphismChecked[currVertex]++;
        createCanonicalForm(g->nautyGraph, gCan, currVertex+1);
        splay_insert(&splayTreeArray[g->numberOfEdges], gCan, currVertex+1, &isPresent);
        if(isPresent)
        {
            free(gCan);
        }
        if(!isPresent && !graphColorable && allInducedSubgraphsColorable)
        {
            counters->nOfTimesC5VertexCritical[currVertex]++;
            writeToG6(g->nautyGraph,currVertex+1);
            //printNautyGraph(g->nautyGraph,currVertex+1);
            (*ctr)++;
        }

        if(graphColorable && !isPresent)
        {
            counters->nOfNonTerminatingGraphs[currVertex]++;
            addVertex(g,sgl,currVertex+1,numberOfVertices,ctr,options,counters);
        }
    }

    // adding more edges will not change the situation
    //if(!graphColorable && !allInducedSubgraphsColorable) return; // consider removing this check and testing first if the graph is H-free or not

    // option 2: graph is not finished
    forEachAfterIndex(nextEndPoint,*validEndPoints,otherEndPointAfter)
    {
        //fprintf(stderr,"%d %d\n",nextEndPoint,otherEndPointAfter);
        // use current edge
        addEdge(g,currVertex,nextEndPoint);
        bitset oldValidEndPoints=(*validEndPoints);
        //(*validEndPoints)=intersection(*validEndPoints,compl(g->adjacencyList[nextEndPoint],currVertex)); // graph should be triangle-free
        (*validEndPoints)=intersection(*validEndPoints,compl(0LL,currVertex));
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,nextEndPoint,validEndPoints,numberSimilarVertices,counters,poor_vertex,max_palette_size_G_minus_u);
        (*validEndPoints)=oldValidEndPoints;
        removeEdge(g,currVertex,nextEndPoint);
    }
}

bool is_subset(bitset b1, bitset b2)
{
    return (b1&b2)==b1;
}

bool canBeExtended(struct graph* g, struct similarGraphsList* sgl, int currVertex, int subgraphOrder, int subgraphCurrentVertex, bitset verticesUsed, bool useHull, int nbColors)
{
    if(subgraphCurrentVertex==subgraphOrder)
    {
        bitset V1AsBitset=EMPTY;
        //bitset V2AsBitset=EMPTY;
        for(int i=0; i<subgraphOrder; i++)
        {
            V1AsBitset=(V1AsBitset|(1LL<<graph1[i]));
            //V2AsBitset=(V2AsBitset|(1LL<<graph2[i]));
        }
        if(useHull)
        {
            compute_hull(g,currVertex,V1AsBitset,nbColors);
        }
        else
        {
            for(int i=0; i<currVertex; i++)
            {
                hull[i]=0;
            }
        }
        // the neighbours of a vertex v in V1 should be a subset of the neighbours of the corresponding vertex from V2 in the graph G[V-V1]
        for(int i=0; i<subgraphOrder; i++)
        {
            if(!is_subset(g->adjacencyList[graph1[i]]&compl(V1AsBitset,g->numberOfVertices), ((g->adjacencyList[graph2[i]]&compl(V1AsBitset,g->numberOfVertices))|hull[graph2[i]]))) return false;
        }  
        // a vertex v in V2 cannot be connected to a non-neighbour of the corresponding vertex from V1 in the graph G[V1]
        for(int i=0; i<subgraphOrder; i++)
        {
            if(!is_subset(g->adjacencyList[graph2[i]]&V1AsBitset, ((g->adjacencyList[graph2[i]]&V1AsBitset)|hull[graph2[i]]))) return false;
        }
        addSimilarGraphs(sgl);
        return true;
    }
    forEach(i,compl(verticesUsed,currVertex))
    {
        bool keepsEdgesAndNonEdgesIntact=true;
        for(int prev=0; prev<subgraphCurrentVertex && keepsEdgesAndNonEdgesIntact; prev++)
        {
            bool connected1=false;
            if((g->adjacencyList[i]&(1LL<<graph1[prev]))>0) connected1=true;
            bool connected2=false;
            if(subgraphAdjacencyMatrix[prev][subgraphCurrentVertex]==1) connected2=true;
            if(connected1!=connected2) keepsEdgesAndNonEdgesIntact=false;
        }
        if(!keepsEdgesAndNonEdgesIntact) continue;
        forEach(j,compl((verticesUsed|(1LL<<i)),currVertex))
        {
            bool edgesAreSubset=true;
            for(int prev=0; prev<subgraphCurrentVertex && edgesAreSubset; prev++)
            {
                bool connected2=false;
                if(subgraphAdjacencyMatrix[prev][subgraphCurrentVertex]==1) connected2=true;
                if(!connected2) continue;
                bool connected1=false;
                if((g->adjacencyList[j]&(1LL<<graph2[prev]))>0) connected1=true;
                if(!connected1) edgesAreSubset=false;
            }
            if(edgesAreSubset)
            {
                graph1[subgraphCurrentVertex]=i;
                graph2[subgraphCurrentVertex]=j;
                bool can=canBeExtended(g,sgl,currVertex,subgraphOrder,subgraphCurrentVertex+1,(verticesUsed|(1LL<<i)|(1LL<<j)),useHull,nbColors);
                if(can) return true;
            }
        }
    }
    return false;
}

// find vertex sets graph1 and graph2 such that g[graph1] is isomorphic with the graph defined by subgraphAdjacencyMatrix and
// there is a subgraph of g[graph2] which is isomorphic with g[graph1] and
// graph1 and graph2 are disjoint and
// the neighbours of a vertex v in graph1 should be a subset of the neighbours of the corresponding vertex from graph2 in the graph G[V-graph1]
// and
// a vertex v in graph2 cannot be connected to a non-neighbour of the corresponding vertex from graph1 in the graph G[graph1]
bool findSimilarGraphs(struct graph* g, struct similarGraphsList* sgl, int currVertex, int subgraphOrder, bool useHull, int nbColors)
{
    for(int i=0; i<currVertex; i++)
    {
        for(int j=0; j<currVertex; j++)
        {
            if(i==j) continue;
            bitset verticesUsed=((1LL<<i)|(1LL<<j));
            graph1[0]=i;
            graph2[0]=j;
            bool can=canBeExtended(g,sgl,currVertex,subgraphOrder,1,verticesUsed,useHull,nbColors);
            if(can) return true;
        }
    }
    return false;
}

int maxReached=-1;

void addVertex(struct graph* g, struct similarGraphsList* sgl, int currVertex, int numberOfVertices, int *ctr, struct options *options, struct counters *counters)
{
    if(currVertex==0)
    {
        addVertex(g,sgl,currVertex+1,numberOfVertices,ctr,options,counters);
        return;
    }
    int nbColors=options->nbColors;
    if(currVertex>maxReached) maxReached=currVertex;
    if(currVertex==numberOfVertices) return;
    bitset validEndPoints=compl(EMPTY,currVertex);


    // look for similar vertices; complexity might be improved by doing proper bookkeeping
    // order in which vertices are traversed by the loops might be important
    int sm=-1;
    int la=-1;
    int bestDegree=100;
    for(int smallestVertex=currVertex-1; smallestVertex>=0; smallestVertex--)
    {
        for(int largestVertex=currVertex-1; largestVertex>=0; largestVertex--)
        {
            if(smallestVertex==largestVertex) continue;
            // neighbours of a are subset of neighbours of b
            if(is_subset(g->adjacencyList[smallestVertex], g->adjacencyList[largestVertex]) && size(g->adjacencyList[sm])<bestDegree)
            {
                sm=smallestVertex;
                la=largestVertex;
                bestDegree=size(g->adjacencyList[sm]);
            }
        }
    }
    if(sm!=-1)
    {

        /*if(currVertex==7 && ((g->adjacencyList[5]&(1LL<<0))>0) && ((g->adjacencyList[5]&(1LL<<1))>0) && ((g->adjacencyList[5]&(1LL<<2))>0) && 
      ((g->adjacencyList[5]&(1LL<<3))==0) && ((g->adjacencyList[5]&(1LL<<4))==0)
      && ((g->adjacencyList[6]&(1LL<<0))>0) && ((g->adjacencyList[6]&(1LL<<1))>0) && ((g->adjacencyList[6]&(1LL<<2))>0) && ((g->adjacencyList[6]&(1LL<<3))>0) && ((g->adjacencyList[6]&(1LL<<5))>0) && ((g->adjacencyList[6]&(1LL<<4))==0))
        {
            fprintf(stderr,"Sim vertex pair %d %d!\n",sm,la);
        }*/
        graph1[0]=sm;
        graph2[0]=la;
        int oldSglSize=sgl->size;
        addSimilarGraphs(sgl);
        counters->nOfTimesSimilarVertices[currVertex]++;
        addEdge(g,currVertex,sm);
        validEndPoints=intersection(validEndPoints,compl((1LL<<sm),currVertex));
        //validEndPoints=intersection(validEndPoints,compl(g->adjacencyList[sm],currVertex));
        validEndPoints=intersection(validEndPoints,compl((1LL<<la),currVertex));
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,1,counters,-1,-1);
        removeEdge(g,currVertex,sm);
        sgl->size=oldSglSize;
        return;
    }

    /*
    // look for similar vertices; complexity might be improved by doing proper bookkeeping
    // order in which vertices are traversed by the loops might be important
    int sm=-1;
    int la=-1;
    for(int smallestVertex=currVertex-1; smallestVertex>=0 && sm==-1; smallestVertex--)
    {
        for(int largestVertex=currVertex-1; largestVertex>=0 && sm==-1; largestVertex--)
        {
            if(smallestVertex==largestVertex) continue;
            // neighbours of a are subset of neighbours of b
            if(is_subset(g->adjacencyList[smallestVertex], g->adjacencyList[largestVertex]))
            {
                sm=smallestVertex;
                la=largestVertex;
            }
        }
    }
    if(sm!=-1)
    {
        graph1[0]=sm;
        graph2[0]=la;
        int oldSglSize=sgl->size;
        addSimilarGraphs(sgl);
        counters->nOfTimesSimilarVertices[currVertex]++;
        addEdge(g,currVertex,sm);
        validEndPoints=intersection(validEndPoints,compl((1LL<<sm),currVertex));
        //validEndPoints=intersection(validEndPoints,compl(g->adjacencyList[sm],currVertex));
        validEndPoints=intersection(validEndPoints,compl((1LL<<la),currVertex));
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,1,counters,-1,-1);
        removeEdge(g,currVertex,sm);
        sgl->size=oldSglSize;
        return;
    }*/

    // the higher this number, the more powerful the program, but also it could become much slower sometimes
    // of course, making this number higher could also turn a non-terminating program into a terminating one, so there is a trade-off.
    int hullConstant=3;
    if(currVertex<=hullConstant) // only compute when the current graph is small enough
    {
        // look for similar vertices where the hull is also taken into account
        sm=-1;
        la=-1;
        for(int smallestVertex=currVertex-1; smallestVertex>=0 && sm==-1; smallestVertex--)
        {
            // compute hull
            bitset omittedVertices=(1LL<<smallestVertex);
            compute_hull(g, currVertex, omittedVertices,nbColors);
            for(int largestVertex=currVertex-1; largestVertex>=0 && sm==-1; largestVertex--)
            {
                if(smallestVertex==largestVertex) continue;
                // neighbours of a are subset of neighbours of b
                if(is_subset(g->adjacencyList[smallestVertex], (g->adjacencyList[largestVertex]|hull[largestVertex])))
                //if(((g->adjacencyList[smallestVertex])&(g->adjacencyList[largestVertex]))==g->adjacencyList[smallestVertex])
                {
                    sm=smallestVertex;
                    la=largestVertex;
                }
            }
        }
        if(sm!=-1)
        {
            graph1[0]=sm;
            graph2[0]=la;
            int oldSglSize=sgl->size;
            addSimilarGraphs(sgl);
            counters->nOfTimesSimilarVerticesWithHull[currVertex]++;
            addEdge(g,currVertex,sm);
            validEndPoints=intersection(validEndPoints,compl((1LL<<sm),currVertex));
            //validEndPoints=intersection(validEndPoints,compl(g->adjacencyList[sm],currVertex));
            validEndPoints=intersection(validEndPoints,compl((1LL<<la),currVertex));
            addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,1,counters,-1,-1);
            removeEdge(g,currVertex,sm);
            sgl->size=oldSglSize;
            return;
        }
    }

    // this part about poor vertex might contain a mistake somewhere: find it later!
    // in fact: I am not sure if it contains a mistake. If we start from a single vertex, the output seems correct, so maybe it is another thing that's wrong!

    // look for poor vertex
    int poor_vertex=-1;
    int max_palette_size_G_minus_u=-1;
    for(int maybe_poor=0; maybe_poor<currVertex; maybe_poor++)
    {
        bitset omittedVertices=(1LL<<maybe_poor);
        generateAllC5ColoringsHelper(g,currVertex,omittedVertices,nbColors,1);
        max_palette_size_G_minus_u=max_palette_size[maybe_poor];
        if(max_palette_size_G_minus_u<nbColors)
        {
            poor_vertex=maybe_poor;
            break;
        }
    }
    if(poor_vertex != -1)
    {
        counters->nOfTimesPoorVertex[currVertex]++;
        addEdge(g,currVertex,poor_vertex);
        validEndPoints=intersection(validEndPoints,compl((1LL<<poor_vertex),currVertex));
  addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,-1,counters,poor_vertex,max_palette_size_G_minus_u);
        removeEdge(g,currVertex,poor_vertex);
        return;
    }

    // look for similar edges
    int subgraphOrder=2;
    for(int i=0; i<subgraphOrder; i++)
    {
        for(int j=0; j<subgraphOrder; j++)
        {
            if(i==j) subgraphAdjacencyMatrix[i][j]=0;
            else subgraphAdjacencyMatrix[i][j]=1;
        }
    }
    int oldSglSize=sgl->size;
    bool similarEdgesFound=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,false,nbColors);

    if(similarEdgesFound) // similar edges were found
    {
        counters->nOfTimesSimilarEdges[currVertex]++;
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
        sgl->size=oldSglSize;
        return;
    }

    // look for similar edges where the hull is also taken into account
    if(currVertex<=hullConstant) // only compute when the current graph is small enough
    {
        oldSglSize=sgl->size;
        similarEdgesFound=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,true,nbColors);
        if(similarEdgesFound) // similar edges were found
        {
            counters->nOfTimesSimilarEdgesWithHull[currVertex]++;
            addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
            sgl->size=oldSglSize;
            return;
        }
    }

    // look for similar triangle
    subgraphOrder=3;
    for(int i=0; i<subgraphOrder; i++)
    {
        for(int j=0; j<subgraphOrder; j++)
        {
            if(i!=j) subgraphAdjacencyMatrix[i][j]=1;
            else subgraphAdjacencyMatrix[i][j]=0;
        }
    }
    oldSglSize=sgl->size;
    bool similarTriangleFound=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,false,nbColors);
    if(similarTriangleFound) // similar triangle was found
    {
        //fprintf(stderr,"Entered!\n");
        counters->nOfTimesSimilarTriangle[currVertex]++;
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
        sgl->size=oldSglSize;
        return;
    }
    // look for similar triangle where the hull is also taken into account
    /*oldSglSize=sgl->size;
    similarTriangleFound=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,true,nbColors);
    if(similarTriangleFound) // similar triangle was found
    {
        //fprintf(stderr,"Entered!\n");
        counters->nOfTimesSimilarTriangle[currVertex]++;
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
        sgl->size=oldSglSize;
        return;
    }*/
    
    // look for similar P3
    /*subgraphOrder=3;
    for(int i=0; i<subgraphOrder; i++)
    {
        for(int j=0; j<subgraphOrder; j++)
        {
            if(i-j==1 || j-i==1) subgraphAdjacencyMatrix[i][j]=1;
            else subgraphAdjacencyMatrix[i][j]=0;
        }
    }
    oldSglSize=sgl->size;
    bool similarP3Found=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,false,nbColors);
    if(similarP3Found) // similar P3 was found
    {
        //fprintf(stderr,"Entered!\n");
        counters->nOfTimesSimilarP3[currVertex]++;
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
        sgl->size=oldSglSize;
        return;
    }*/
    
    // look for similar diamond
    subgraphOrder=4;
    for(int i=0; i<subgraphOrder; i++)
    {
        for(int j=0; j<subgraphOrder; j++)
        {
            if(i!=j && !(i==0 && j==1)) subgraphAdjacencyMatrix[i][j]=1;
            else subgraphAdjacencyMatrix[i][j]=0;
        }
    }
    oldSglSize=sgl->size;
    bool similarDiamondFound=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,false,nbColors);
    if(similarDiamondFound) // similar diamond was found
    {
        counters->nOfTimesSimilarDiamond[currVertex]++;
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
        sgl->size=oldSglSize;
        return;
    }
    // look for similar diamond where the hull is also taken into account
    /*oldSglSize=sgl->size;
    similarDiamondFound=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,true,nbColors);
    if(similarDiamondFound) // similar diamond was found
    {
        counters->nOfTimesSimilarDiamond[currVertex]++;
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
        sgl->size=oldSglSize;
        return;
    }*/

    
    // look for similar P4
    /*subgraphOrder=4;
    for(int i=0; i<subgraphOrder; i++)
    {
        for(int j=0; j<subgraphOrder; j++)
        {
            if(i-j==1 || j-i==1) subgraphAdjacencyMatrix[i][j]=1;
            else subgraphAdjacencyMatrix[i][j]=0;
        }
    }
    oldSglSize=sgl->size;
    bool similarP4Found=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,false,nbColors);
    if(similarP4Found) // similar P4 was found
    {
        counters->nOfTimesSimilarP4[currVertex]++;
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
        sgl->size=oldSglSize;
        return;
    }

    // look for similar C4
    subgraphOrder=4;
    for(int i=0; i<subgraphOrder; i++)
    {
        for(int j=0; j<subgraphOrder; j++)
        {
            if(i-j==1 || j-i==1) subgraphAdjacencyMatrix[i][j]=1;
            else subgraphAdjacencyMatrix[i][j]=0;
        }
    }
    subgraphAdjacencyMatrix[0][3]=1;
    subgraphAdjacencyMatrix[3][0]=1;
    oldSglSize=sgl->size;
    bool similarC4Found=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,false,nbColors);
    if(similarC4Found) // similar C4 was found
    {
        counters->nOfTimesSimilarC4[currVertex]++;
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
        sgl->size=oldSglSize;
        return;
    }

    // look for similar K13
    subgraphOrder=4;
    for(int i=0; i<subgraphOrder; i++)
    {
        for(int j=0; j<subgraphOrder; j++)
        {
            subgraphAdjacencyMatrix[i][j]=0;
        }
    }
    subgraphAdjacencyMatrix[0][1]=1;
    subgraphAdjacencyMatrix[1][0]=1;
    subgraphAdjacencyMatrix[0][2]=1;
    subgraphAdjacencyMatrix[2][0]=1;
    subgraphAdjacencyMatrix[0][3]=1;
    subgraphAdjacencyMatrix[3][0]=1;
    oldSglSize=sgl->size;
    bool similarK13Found=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,false,nbColors);
    if(similarK13Found) // similar K13 was found
    {
        counters->nOfTimesSimilarK13[currVertex]++;
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
        sgl->size=oldSglSize;
        return;
    }

    // look for similar K4
    subgraphOrder=4;
    for(int i=0; i<subgraphOrder; i++)
    {
        for(int j=0; j<subgraphOrder; j++)
        {
            if(i!=j) subgraphAdjacencyMatrix[i][j]=1;
            else subgraphAdjacencyMatrix[i][j]=0;
        }
    }
    oldSglSize=sgl->size;
    bool similarK4Found=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,false,nbColors);
    if(similarK4Found) // similar K4 was found
    {
        //fprintf(stderr,"Entered!\n");
        counters->nOfTimesSimilarK4[currVertex]++;
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
        sgl->size=oldSglSize;
        return;
    }

    // look for similar ComplP1PlusP3
    subgraphOrder=4;
    for(int i=0; i<subgraphOrder; i++)
    {
        for(int j=0; j<subgraphOrder; j++)
        {
            if(i==0 || j==0)
            {
                if(i+j==1) subgraphAdjacencyMatrix[i][j]=1;
                else subgraphAdjacencyMatrix[i][j]=0;
                continue;
            }
            if(i!=j) subgraphAdjacencyMatrix[i][j]=1;
            else subgraphAdjacencyMatrix[i][j]=0;
        }
    }
    oldSglSize=sgl->size;
    bool similarComplP1PlusP3Found=findSimilarGraphs(g,sgl,currVertex,subgraphOrder,false,nbColors);
    if(similarComplP1PlusP3Found) // similar ComplP1PlusP3 was found
    {
        //fprintf(stderr,"Entered!\n");
        counters->nOfTimesSimilarComplP1PlusP3[currVertex]++;
        addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,subgraphOrder,counters,-1,-1);
        sgl->size=oldSglSize;
        return;
    }*/

    // look for cut vertex
    // attention, do not uncomment! This is not implemented correctly at the moment
    // the newly added vertex u should be adjacent to some vertex in V(G)\x, but it can also still be adjacent to x.
    // the algorithm was never able to apply the lemma anyways, so I will not implement this logic for now.
    /*
    dfsNumberCounter = 0;
    for(int i=0; i<currVertex; i++)
    {
        dfs_num[i]=DFS_WHITE;
        articulation_vertex[i]=false;
    }
    for(int i=0; i<currVertex; i++)
    {
        if(dfs_num[i] == DFS_WHITE)
        {
            dfsRoot = i; rootChildren=0;
            bool found=articulationPointAndBridge(g,i);
            articulation_vertex[dfsRoot] = (rootChildren > 1); //special case
            if(articulation_vertex[dfsRoot]) found=true;
            if(found) break;          
        }
    }
    for(int i=0; i<currVertex; i++)
    {
        if(articulation_vertex[i])
        {
            counters->nOfTimesCutVertex[currVertex]++;
            //validEndPoints=intersection(validEndPoints,compl((1LL<<i),currVertex));
            addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,0,counters,-1,-1);
            return;
        }
    }*/

    // no similarity lemma could be used
	counters->nOfTimesNoSimilar[currVertex]++;
	addEdgesInAllPossibleWays(g,sgl,currVertex,numberOfVertices,ctr,options,-1,&validEndPoints,0,counters,-1,-1);
}

//  Generate K2hypohamiltonian graphs of order numberOfVertices.
void generateK2HypohamiltonianGraphs(int numberOfVertices, 
    struct counters *counters, struct options *options) {
    struct graph g = {.numberOfVertices = numberOfVertices};
    g.adjacencyList = malloc(sizeof(bitset)*numberOfVertices);
    g.forbiddenEdges = malloc(sizeof(bitset)*numberOfVertices);
    g.verticesOfDeg = malloc(sizeof(bitset)*numberOfVertices);
    for(int i = 0; i < g.numberOfVertices; i++) {
        g.verticesOfDeg[i] = EMPTY;
    }
    fprintf(stderr,"Number of vertices: %d\n",numberOfVertices);
    int ctr=0;
    struct similarGraphsList sgl;
    initializeSimilarGraphsList(&sgl);

    // comment later!
    //doIt2(&g);

    /*addEdge(&g,0,1);
    addEdge(&g,1,2);
    addEdge(&g,2,3);
    addEdge(&g,3,4);*/

    //addEdge(&g,4,0);

    /*addEdge(&g,4,5);
    addEdge(&g,5,6);
    addEdge(&g,6,0);*/

    //(struct graph* g,  struct similarGraphsList* sgl, int currVertex, int numberOfVertices, int *ctr, struct options *options, int otherEndPointAfter, bitset *validEndPoints, int numberSimilarVertices, struct counters *counters, int poor_vertex, int max_palette_size_G_minus_u)

    // make clique on nbColors+1 vertices
    if(numberOfVertices>=options->nbColors+1)
    {
        for(int i=0; i<MAXVERTICES; i++)
        {
            pathQueueStart[i]=0;
            pathQueueEnd[i]=-1;
        }
        emptyGraph(&g);
        for(int i=0; i<options->nbColors+1; i++)
        {
            for(int j=i+1; j<options->nbColors+1; j++)
            {
                addEdge(&g,i,j);
            }
        }
        bitset validEndPoints=0LL;
        initializeSimilarGraphsList(&sgl);
        //addVertex(&g,&sgl,options->nbColors+1,numberOfVertices,&ctr,options,counters);
        addEdgesInAllPossibleWays(&g,&sgl,options->nbColors,numberOfVertices,&ctr,options,numberOfVertices,&validEndPoints,0,counters,-2,-2);
    }

    // make odd hole on 2s+1 vertices
    for(int s=2; s<=(options->pathLength-1)/2 && 2*s+1<=numberOfVertices; s++)
    //for(int s=2; s<=2 && 2*s+1<=numberOfVertices; s++)
    {
        initializeSimilarGraphsList(&sgl);
        for(int i=0; i<MAXVERTICES; i++)
        {
            pathQueueStart[i]=0;
            pathQueueEnd[i]=-1;
        }
        emptyGraph(&g);
        for(int i=0; i<2*s+1; i++)
        {
            addEdge(&g,i,(i+1)%(2*s+1));
        }
        //addVertex(&g,&sgl,2*s+1,numberOfVertices,&ctr,options,counters);
        bitset validEndPoints=0LL;
        addEdgesInAllPossibleWays(&g,&sgl,2*s,numberOfVertices,&ctr,options,numberOfVertices,&validEndPoints,0,counters,-2,-2);
    }

    
    // make odd antihole on 2s+1 vertices
    for(int s=3; s<=options->nbColors && 2*s+1<=numberOfVertices; s++)
    //for(int s=options->nbColors; s<=options->nbColors && 2*s+1<=numberOfVertices; s++)
    {
        initializeSimilarGraphsList(&sgl);
        for(int i=0; i<MAXVERTICES; i++)
        {
            pathQueueStart[i]=0;
            pathQueueEnd[i]=-1;
        }
        emptyGraph(&g);
        for(int i=0; i<2*s+1; i++)
        {
            for(int j=i+2; j<2*s+1; j++)
            {
                if(((i+1)%(2*s+1))==j) continue;
                if(((j+1)%(2*s+1))==i) continue;
                addEdge(&g,i,j);
            }
        }
        //addVertex(&g,&sgl,2*s+1,numberOfVertices,&ctr,options,counters);
        bitset validEndPoints=0LL;
        addEdgesInAllPossibleWays(&g,&sgl,2*s,numberOfVertices,&ctr,options,numberOfVertices,&validEndPoints,0,counters,-2,-2);
    }

    // make empty graph
    /*
    for(int i=0; i<MAXVERTICES; i++)
    {
        pathQueueStart[i]=0;
        pathQueueEnd[i]=-1;
    }
    emptyGraph(&g);
    //addVertex(&g,&sgl,2*s+1,numberOfVertices,&ctr,options,counters);
    bitset validEndPoints=0LL;
    addEdgesInAllPossibleWays(&g,&sgl,0,numberOfVertices,&ctr,options,numberOfVertices,&validEndPoints,0,counters,-2,-2);*/


    fprintf(stderr,"ctr: %d\n",ctr);
    fprintf(stderr,"maxReached: %d\n",maxReached);
    freeSimilarGraphsList(&sgl);
    free(g.adjacencyList);
    free(g.forbiddenEdges);
    free(g.verticesOfDeg);
}

int main(int argc, char ** argv) {
    /*Initialize the optional command line flags.*/

    int opt;
    bool haveModResPair = false;
    bool printCountsPerEdgeNumber = false;
    struct options options = 
        {.minimalGirth = -1,
         .minimumDegree = -1,
         .maximumDegree = -1,
         .pathLength = -1,
         .nbColors = -1,
         .planarFlag = false,
         .bipartiteFlag = false,
         .splitLevel = SPLITLEVEL,
         .splitCounter = 0,
         .remainder = 0,
         .modulo = 1};
    char *graphClass = "";
    while(1) {
        int option_index = 0;
        static struct option long_options[] =
        {
            {"bipartite", no_argument, NULL, 'b'},
            {"minimum-degree", required_argument, NULL, 'd'},
            {"maximum-degree", required_argument, NULL, 'D'},
            {"path-length", required_argument, NULL, 'l'},
            {"number-colors", required_argument, NULL, 'k'},
            {"edges-count", no_argument, NULL, 'e'},
            {"girth", required_argument, NULL, 'g'},
            {"help", no_argument, NULL, 'h'},
            {"planar", no_argument, NULL, 'p'},
            {"splitlevel", required_argument, NULL, 'X'}
        };
        opt = getopt_long(argc, argv, "bcd:D:l:k:eg:hpX:", long_options, &option_index);
        if(opt == -1) break;
        char *ptr;
        switch(opt) {
            case 'b':
                options.bipartiteFlag = true;
                graphClass = "bipartite ";
                break;
            case 'd':
                options.minimumDegree = strtol(optarg, &ptr, 10);
                break;
            case 'D':
                options.maximumDegree = strtol(optarg, &ptr, 10);
                break;
            case 'l':
                options.pathLength = strtol(optarg, &ptr, 10);
                break;
            case 'k':
                options.nbColors = strtol(optarg, &ptr, 10);
                break;
            case 'e':
                printCountsPerEdgeNumber = true;
                break;
            case 'h':
                fprintf(stderr, "%s", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'g': 
                options.minimalGirth = strtol(optarg, &ptr, 10);
                break;
            case 'p':
                options.planarFlag = true;
                graphClass = "planar ";
                break;
            case 'X':
                options.splitLevel = strtol(optarg, &ptr, 10);
                break;
        }
    }
    //  Check if there is a non-option argument.
    if (optind >= argc) {
        fprintf(stderr, "Error: add number of vertices.");
        fprintf(stderr, "%s", USAGE);
        return 1;
    }
    //  First non-option argument should be number of vertices.
    char* endptr;
    int numberOfVertices = strtol(argv[optind++], &endptr, 10);

    if(numberOfVertices <= 2 || numberOfVertices > MAXN) {
        fprintf(stderr, "Error: n needs to be a number between 3 and %d.",
         MAXN);
        fprintf(stderr, "%s",USAGE);
        return 1;
    }

    //  Check for other non-option arguments.
    while (optind < argc) {
        bool pairIsInvalid = false;
        options.remainder = strtol(argv[optind], &endptr, 10);
        if( !endptr || *endptr != '/' || *(endptr+1) == '\0') {
            pairIsInvalid = true;
        }
        options.modulo = strtol(endptr+1, &endptr, 10);
        if( !endptr || *endptr != '\0') {
            pairIsInvalid = true;
        }
        if(options.modulo <= options.remainder) {
            pairIsInvalid = true;
        }
        if(haveModResPair) {
            fprintf(stderr,
             "Error: You can only add one res/mod pair as an argument.\n");
            fprintf(stderr, "%s\n", USAGE);
            fprintf(stderr,
             "Use ./genK2Hypohamiltonian --help for more detailed instructions.\n");
            return 1;
        }
        haveModResPair = true;

        if(pairIsInvalid) {
            fprintf(stderr,
                 "Error: Invalid res/mod pair: '%s'.\n", argv[optind]);
            fprintf(stderr, "%s\n", USAGE);
            fprintf(stderr,
             "Use ./genK2Hypohamiltonian --help for more detailed instructions.\n");
            return 1;
        }
        fprintf(stderr, "Class=%d/%d. Splitlevel = %d.\n",
         options.remainder, options.modulo, options.splitLevel);
        optind++;
    }

    /*Do some correctness checks.*/

    //  Cubic K2-hypohamiltonian graphs have girth at least 5.
    if(options.maximumDegree == 3) {
        if(options.minimalGirth < 5) {
            options.minimalGirth = 5;
        }
    }

    //  Bipartite K2-hypohamiltonian graphs can only have even girth.
    nauty_check(WORDSIZE, SETWORDSNEEDED(numberOfVertices),
     numberOfVertices, NAUTYVERSIONID);

    /*Initialize all counters to 0.*/
    struct counters counters = {};
    counters.nOfNonIsoGraphsWithEdges = calloc(sizeof(long long unsigned int),
     numberOfVertices*(numberOfVertices - 1)/2+1);
    counters.nOfTimesPoorVertex = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarVertices = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarVerticesWithHull = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarEdges = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarEdgesWithHull = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarTriangle = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarP3 = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarDiamond = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarP4 = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarC4 = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarK13 = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarK4 = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesSimilarComplP1PlusP3 = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesCutVertex = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesNoSimilar = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);

    counters.nOfTimesIsomorphismChecked = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfNonTerminatingGraphs = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesNotC5Colorable = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    counters.nOfTimesC5VertexCritical = calloc(sizeof(long long unsigned int),
     numberOfVertices+1);
    /*Main routine*/
    clock_t start = clock();
    // read the graph F (the generator generates (P_t,F)-free obstructions
    // The input should contain a single graph F. (I.e. the algorithm does currently not support (P_t,F_1,F_2)-free for example.
    while(getline(&graphString, &size, stdin) != -1) {
		nVerticesF = getNumberOfVertices(graphString);
		loadGraph(graphString, nVerticesF, adjacencyListF);
    
        fprintf(stderr,"Number of vertices in F: %d\n",nVerticesF);
        for(int i=0; i<nVerticesF; i++)
        {
            for(int j=0; j<nVerticesF; j++)
            {
                if((1LL<<i)&adjacencyListF[j]) fprintf(stderr,"1");
                else fprintf(stderr,"0");
            }
            fprintf(stderr,"\n");
        }
        generateK2HypohamiltonianGraphs(numberOfVertices, &counters, &options);
    }
    clock_t end = clock();

    /*Output results*/
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    if(options.minimalGirth != -1) {
        fprintf(stderr,
         "\rGenerated %lld %sgraphs of order %d with girth at least %d in %f seconds.\n",
         counters.nOfGraphsFound, graphClass, numberOfVertices, options.minimalGirth, time_spent);
    }
    else {
        fprintf(stderr,
         "\rGenerated %lld %sgraphs of order %d in %f seconds.\n", 
         counters.nOfGraphsFound, graphClass, numberOfVertices, time_spent);
    }
    if(options.minimumDegree != -1) {
        fprintf(stderr, "Their minimum degree is %d.\n", options.minimumDegree);
    }
    if(options.maximumDegree != -1) {
        fprintf(stderr, "Their maximum degree is %d.\n", options.maximumDegree);
    }
    if(options.planarFlag) {
        fprintf(stderr, "Times graph was non-planar: %llu (%.2f%%)\n",
         counters.nOfTimesWasNonPlanar, 
         (double) 100*counters.nOfTimesWasNonPlanar/counters.nOfTimesCheckedPlanarity);
    }
    fprintf(stderr, "---\n");
    fprintf(stderr, "Number of non-isomorphic graphs checked: %llu\n",
     counters.nOfTimesCheckedIsomorphism - counters.nOfTimesWasIsomorphic);
    fprintf(stderr, "Times contained type A: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n",
     counters.nOfTimesContainedTypeA, 
     (double) 100*counters.nOfTimesContainedTypeA/counters.nOfTimesCheckedTypeA,
     counters.nOfTimesTypeAObstructionChosen, 
     (double) 100*counters.nOfTimesTypeAObstructionChosen/counters.nOfTimesObstructionFound);
    // fprintf(stderr, "Times contained type B: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n",
    //  counters.nOfTimesContainedTypeB, 
    //  (double) 100*counters.nOfTimesContainedTypeB/counters.nOfTimesCheckedTypeB,
    //  counters.nOfTimesTypeBObstructionChosen, 
    //  (double) 100*counters.nOfTimesTypeBObstructionChosen/counters.nOfTimesObstructionFound);
    fprintf(stderr, "Times contained type C: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n",
     counters.nOfTimesContainedTypeC, 
     (double) 100*counters.nOfTimesContainedTypeC/counters.nOfTimesCheckedTypeC,
     counters.nOfTimesTypeCObstructionChosen, 
     (double) 100*counters.nOfTimesTypeCObstructionChosen/counters.nOfTimesObstructionFound);
    fprintf(stderr, "Times contained bad degree vertex: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n",
     counters.nOfTimesContainedDegreeObstruction, 
     (double) 100*counters.nOfTimesContainedDegreeObstruction/counters.nOfTimesCheckedDegreeObstruction,
     counters.nOfTimesDegreeObstructionChosen, 
     (double) 100*counters.nOfTimesDegreeObstructionChosen/counters.nOfTimesObstructionFound);
    if(options.minimalGirth < 3) {
        fprintf(stderr,
         "Times contained general triangle obstruction: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n", 
         counters.nOfTimesContainedStarObstruction, 
         (double) 100*counters.nOfTimesContainedStarObstruction/counters.nOfTimesCheckedStarObstruction,
         counters.nOfTimesGeneralTriangleObstructionChosen, 
         (double) 100*counters.nOfTimesGeneralTriangleObstructionChosen/counters.nOfTimesObstructionFound);   
        fprintf(stderr,
         "Times contained deg 3 triangle vertex: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n", 
         counters.nOfTimesContainedTriangleObstruction, 
         (double) 100*counters.nOfTimesContainedTriangleObstruction/counters.nOfTimesCheckedTriangleObstruction,
         counters.nOfTimesTriangleObstructionChosen, 
         (double) 100*counters.nOfTimesTriangleObstructionChosen/counters.nOfTimesObstructionFound);
    }
    if(options.minimalGirth < 4) {
        fprintf(stderr,
         "Times contained deg 3 4-cycle vertex: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n", 
         counters.nOfTimesContainedSquareObstruction, 
         (double) 100*counters.nOfTimesContainedSquareObstruction/counters.nOfTimesCheckedSquareObstruction,
         counters.nOfTimes4CycleObstructionChosen, 
         (double) 100*counters.nOfTimes4CycleObstructionChosen/counters.nOfTimesObstructionFound);
        fprintf(stderr,
         "Times contained general 4-cycle obstruction: %llu (%.2f%%)\n\tChosen: %llu (%.2f%%)\n", 
         counters.nOfTimesContainedArrowObstruction, 
         (double) 100*counters.nOfTimesContainedArrowObstruction/counters.nOfTimesCheckedArrowObstruction,
         counters.nOfTimesGeneral4CycleObstructionChosen, 
         (double) 100*counters.nOfTimesGeneral4CycleObstructionChosen/counters.nOfTimesObstructionFound); 
    }   
    fprintf(stderr, "Times no obstruction applied: %llu (%.2f%% of non-isomorphic graphs)\n",
     counters.nOfTimesNoObstructionApplied, 
     (double) 100*counters.nOfTimesNoObstructionApplied / 
     (counters.nOfTimesCheckedIsomorphism - counters.nOfTimesWasIsomorphic));
    fprintf(stderr, "---\n");
    fprintf(stderr, "Times adding successor edge made graph hamiltonian: %llu (%.2f%%)\n",
     counters.nOfTimesWasHamiltonian, 
     (double) 100*counters.nOfTimesWasHamiltonian/counters.nOfTimesCheckedHamiltonicity);
    if(options.minimalGirth >= 4) {
        fprintf(stderr, "Times adding successor edge gave forbidden girth: %llu (%.2f%%)\n",
         counters.nOfTimesHadForbiddenGirth, 
         (double) 100*counters.nOfTimesHadForbiddenGirth/counters.nOfTimesCheckedGirth);
    }
    if(options.maximumDegree != -1) {
        fprintf(stderr, "Times adding successor edge gave vertex forbidden degree: %llu (%.2f%%)\n",
         counters.nOfTimesMaximumDegreeExceeded, 
         (double) 100*counters.nOfTimesMaximumDegreeExceeded/counters.nOfTimesCheckedMaximumDegree);
    }
    fprintf(stderr, "Size of graphs in splay trees: %.3f gigabytes\n", 
     (double) (sizeof(graph)*numberOfVertices *
     (counters.nOfTimesCheckedIsomorphism - counters.nOfTimesWasIsomorphic))/1000000000);

    if(printCountsPerEdgeNumber) {
        for(int i = numberOfVertices+3; i <= numberOfVertices*(numberOfVertices - 1)/2; i++) {
            fprintf(stderr, "%8llu non-isomorphic graphs with %d edges generated.\n",
             counters.nOfNonIsoGraphsWithEdges[i], i);
            if(counters.nOfNonIsoGraphsWithEdges[i]==0) {
                fprintf(stderr,"No graphs with more edges generated.\n");
                break;
            }
        }
    }

    fprintf(stderr,"Start of statistics\n");
    fprintf(stderr,"Results for C5-vertex-critical P%d-free graphs on %d vertices\n",options.pathLength,numberOfVertices);
    fprintf(stderr,"The maximum recursion depth reached by the algorithm was %d\n",maxReached);
    fprintf(stderr,"\n");

    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar vertices lemma was applied %llu times on level %d\n",counters.nOfTimesSimilarVertices[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar vertices lemma with hull was applied %llu times on level %d\n",counters.nOfTimesSimilarVerticesWithHull[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Poor vertex lemma was applied %llu times on level %d\n",counters.nOfTimesPoorVertex[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar edges lemma was applied %llu times on level %d\n",counters.nOfTimesSimilarEdges[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar edges lemma with hull was applied %llu times on level %d\n",counters.nOfTimesSimilarEdgesWithHull[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar Triangle lemma was applied %llu times on level %d\n",counters.nOfTimesSimilarTriangle[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar P3 lemma was applied %llu times on level %d\n",counters.nOfTimesSimilarP3[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar Diamond lemma was applied %llu times on level %d\n",counters.nOfTimesSimilarDiamond[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar P4 lemma was applied %llu times on level %d\n",counters.nOfTimesSimilarP4[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar C4 lemma was applied %llu times on level %d\n",counters.nOfTimesSimilarC4[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar K13 lemma was applied %llu times on level %d\n",counters.nOfTimesSimilarK13[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar K4 lemma was applied %llu times on level %d\n",counters.nOfTimesSimilarK4[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Similar ComplP1PlusP3 lemma was applied %llu times on level %d\n",counters.nOfTimesSimilarComplP1PlusP3[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Cutvertex lemma was applied %llu times on level %d\n",counters.nOfTimesCutVertex[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"No similarity lemma was applied %llu times on level %d\n",counters.nOfTimesNoSimilar[i],i);
    }
    fprintf(stderr,"\n");

    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"Isomorphism was checked %llu times on level %d\n",counters.nOfTimesIsomorphismChecked[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"There were %llu non-terminating graphs on level %d\n",counters.nOfNonTerminatingGraphs[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"There were %llu not-C5-colorable graphs on level %d\n",counters.nOfTimesNotC5Colorable[i],i);
    }
    fprintf(stderr,"\n");
    for(int i=0; i<numberOfVertices; i++)
    {
        fprintf(stderr,"There were %llu C5-vertex-critical graphs on level %d\n",counters.nOfTimesC5VertexCritical[i],i);
    }
    fprintf(stderr,"\n");
    return 0;
}
