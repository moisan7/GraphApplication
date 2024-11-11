#include "pch.h"
#include "Graph.h"
#include <set>
#include <stack>
#include <limits>
#include <list>
#include <cmath>

// Declaration of global variables
CVertex* destinationVertex;                      // Target destination vertex for the TSP
CTrack optimalSolution(NULL);                    // Stores the optimal solution path
CVisits visitList(NULL);                         // List of vertices to visit
double shortestPathLength = std::numeric_limits<double>::max(); // Shortest path length found
double currentPathLength = 0.0;                  // Current path length
std::list<CVertex*> orderedVertices;             // Ordered list of vertices
std::list<CVertex*> unorderedVertices;           // Unordered list of vertices

// Struct to represent a node in the path
struct PathNode {
    CEdge* edge;                                 // Edge to the current vertex
    PathNode* previousNode;                      // Previous node in the path
};

// Forward declarations of functions
void recursiveBacktracking(PathNode* previousNode, CVertex* currentVertex);
bool allVerticesVisited();

// Function to check if all vertices have been visited
bool allVerticesVisited() {
    for (CVertex* vertex : visitList.m_Vertices) {
        bool isVertexVisited = false;
        for (CEdge* edge : vertex->m_Edges) {
            if (edge->m_Processed) {
                isVertexVisited = true;
                break;
            }
        }
        if (!isVertexVisited) return false;
    }
    return true;
}

// Backtracking solution for the Traveling Salesman Problem (TSP)
CTrack SalesmanTrackBacktracking(CGraph& graph, CVisits& visits) {
    visitList = visits;
    currentPathLength = 0.0;
    shortestPathLength = std::numeric_limits<double>::max();

    // Initialize source and destination vertices
    CVertex* startVertex = visitList.m_Vertices.front();
    destinationVertex = visitList.m_Vertices.back();
    visitList.m_Vertices.pop_back();  // Remove destination vertex from visit list

    // Mark all edges as unprocessed
    for (CEdge& edge : graph.m_Edges) {
        edge.m_Processed = false;
    }

    // Begin backtracking
    recursiveBacktracking(nullptr, startVertex);
    return optimalSolution;
}

// Recursive backtracking function to explore possible paths
void recursiveBacktracking(PathNode* previousNode, CVertex* currentVertex) {
    if (currentVertex == destinationVertex && allVerticesVisited()) {
        if (currentPathLength < shortestPathLength) {
            optimalSolution.Clear();
            while (previousNode) {
                optimalSolution.m_Edges.push_front(previousNode->edge);
                previousNode = previousNode->previousNode;
            }
            shortestPathLength = currentPathLength;
        }
        return;
    }

    PathNode pathNode;
    pathNode.previousNode = previousNode;

    for (CEdge* edge : currentVertex->m_Edges) {
        if (!edge->m_Processed && currentPathLength + edge->m_Length < shortestPathLength) {
            edge->m_Processed = true;
            pathNode.edge = edge;
            currentPathLength += edge->m_Length;

            recursiveBacktracking(&pathNode, edge->m_pDestination);

            // Backtrack to explore other paths
            currentPathLength -= edge->m_Length;
            edge->m_Processed = false;
        }
    }
}

// ==============================================
// Greedy + Backtracking TSP Solution
// ==============================================

// Struct to hold distance information for greedy search
struct DistanceMatrixElement {
    double distance;                              // Distance between vertices
    std::list<CEdge*> pathEdges;                  // Edges forming the path
};

// Greedy search variables
DistanceMatrixElement distanceMatrix[20][20];    // Distance matrix
int pathOrder[20];                               // Current path order
int bestPathOrder[20];                           // Best path order found
double bestGreedyPathLength = std::numeric_limits<double>::max(); // Best path length found
double currentGreedyPathLength = 0.0;            // Current path length

// Recursive function to find the shortest path using a greedy approach
void greedyBacktrackingRec(int depth) {
    if (currentGreedyPathLength >= bestGreedyPathLength) return;

    if (depth == orderedVertices.size() - 1) {
        currentGreedyPathLength += fabs(distanceMatrix[pathOrder[depth - 1]][pathOrder[depth]].distance);
        if (currentGreedyPathLength < bestGreedyPathLength) {
            bestGreedyPathLength = currentGreedyPathLength;
            std::copy(std::begin(pathOrder), std::end(pathOrder), std::begin(bestPathOrder));
        }
        currentGreedyPathLength -= fabs(distanceMatrix[pathOrder[depth - 1]][pathOrder[depth]].distance);
    } else {
        for (int j = 1; j < orderedVertices.size() - 1; j++) {
            bool isInPath = false;
            for (int k = 1; k < depth && !isInPath; k++) {
                isInPath = (j == pathOrder[k]);
            }
            if (!isInPath) {
                pathOrder[depth] = j;
                currentGreedyPathLength += fabs(distanceMatrix[pathOrder[depth - 1]][pathOrder[depth]].distance);
                greedyBacktrackingRec(depth + 1);
                currentGreedyPathLength -= fabs(distanceMatrix[pathOrder[depth - 1]][pathOrder[depth]].distance);
            }
        }
    }
}

// Greedy approach to initialize and find the approximate TSP solution
CTrack SalesmanTrackBacktrackingGreedy(CGraph& graph, CVisits& visits) {
    orderedVertices = visits.m_Vertices;
    unorderedVertices = orderedVertices;

    // Initialize distance matrix using Dijkstra's algorithm for all vertices
    int i = 0;
    for (CVertex* startVertex : orderedVertices) {
        Dijkstra(graph, startVertex);
        int j = 0;
        for (CVertex* endVertex : unorderedVertices) {
            distanceMatrix[i][j].distance = endVertex->m_DijkstraDistance;
            while (endVertex->m_pDijkstraPrevious) {
                distanceMatrix[i][j].pathEdges.push_front(endVertex->m_pDijkstraPrevious);
                endVertex = endVertex->m_pDijkstraPrevious->m_pOrigin;
            }
            j++;
        }
        i++;
    }

    bestGreedyPathLength = std::numeric_limits<double>::max();
    currentGreedyPathLength = 0.0;
    pathOrder[0] = 0;   
    pathOrder[orderedVertices.size() - 1] = orderedVertices.size() - 1;
    greedyBacktrackingRec(1);

    CTrack result(&graph);
    for (int m = 0; m < orderedVertices.size() - 1; m++) {
        while (!distanceMatrix[bestPathOrder[m]][bestPathOrder[m + 1]].pathEdges.empty()) {
            result.m_Edges.push_back(distanceMatrix[bestPathOrder[m]][bestPathOrder[m + 1]].pathEdges.front());
            distanceMatrix[bestPathOrder[m]][bestPathOrder[m + 1]].pathEdges.pop_front();
        }
    }
    return result;
}
