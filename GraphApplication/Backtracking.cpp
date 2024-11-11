#include "pch.h"
#include "Graph.h"
#include <set>
#include <stack>
#include <limits>
#include <list>
#include <cmath>

// Global variables for storing the state of the TSP problem
CVertex* destinationVertex;                      // Destination vertex for the TSP path
CTrack optimalSolution(NULL);                    // To store the optimal TSP solution path
CVisits visitList(NULL);                         // List of vertices to visit
double shortestPathLength = std::numeric_limits<double>::max(); // Shortest path length found
double currentPathLength = 0.0;                  // Length of the current path being explored
std::list<CVertex*> orderedVertices;             // Ordered list of vertices for TSP solution
std::list<CVertex*> unorderedVertices;           // Unordered list of vertices for TSP exploration

// Struct to represent a node in the backtracking search path
struct PathNode {
    CEdge* edge;                                 // Edge leading to the current vertex
    PathNode* previousNode;                      // Previous node in the path (for backtracking)
};

// Function declarations for the backtracking solution
void recursiveBacktracking(PathNode* previousNode, CVertex* currentVertex);
bool allVerticesVisited();

// Function to check if all vertices have been visited by the current path
bool allVerticesVisited() {
    // Iterate over each vertex in the visit list
    for (CVertex* vertex : visitList.m_Vertices) {
        bool isVertexVisited = false;

        // Check if there is any processed edge for this vertex
        for (CEdge* edge : vertex->m_Edges) {
            if (edge->m_Processed) {
                isVertexVisited = true;
                break;
            }
        }

        // If any vertex has not been visited, return false
        if (!isVertexVisited) return false;
    }

    return true; // All vertices have been visited
}

// Main backtracking function to solve the Traveling Salesman Problem (TSP)
CTrack SalesmanTrackBacktracking(CGraph& graph, CVisits& visits) {
    visitList = visits;
    currentPathLength = 0.0;
    shortestPathLength = std::numeric_limits<double>::max();

    // Initialize the start and destination vertices for the TSP
    CVertex* startVertex = visitList.m_Vertices.front();
    destinationVertex = visitList.m_Vertices.back();
    visitList.m_Vertices.pop_back();  // Remove destination vertex from visit list

    // Mark all edges as unprocessed
    for (CEdge& edge : graph.m_Edges) {
        edge.m_Processed = false;
    }

    // Start the recursive backtracking search from the start vertex
    recursiveBacktracking(nullptr, startVertex);
    return optimalSolution;  // Return the optimal solution found
}

// Recursive backtracking function to explore paths
void recursiveBacktracking(PathNode* previousNode, CVertex* currentVertex) {
    // Check if the current vertex is the destination and all vertices have been visited
    if (currentVertex == destinationVertex && allVerticesVisited()) {
        // If the current path length is shorter than the previously found path, update the solution
        if (currentPathLength < shortestPathLength) {
            optimalSolution.Clear();

            // Reconstruct the optimal path by following the previous nodes
            while (previousNode) {
                optimalSolution.m_Edges.push_front(previousNode->edge);
                previousNode = previousNode->previousNode;
            }
            shortestPathLength = currentPathLength; // Update the shortest path length
        }
        return;
    }

    PathNode pathNode;
    pathNode.previousNode = previousNode;

    // Explore all edges connected to the current vertex
    for (CEdge* edge : currentVertex->m_Edges) {
        // If the edge has not been processed and the resulting path is promising, explore it
        if (!edge->m_Processed && currentPathLength + edge->m_Length < shortestPathLength) {
            edge->m_Processed = true;  // Mark the edge as processed
            pathNode.edge = edge;      // Update the current edge in the path
            currentPathLength += edge->m_Length; // Add the edge length to the current path length

            // Recurse to the next vertex
            recursiveBacktracking(&pathNode, edge->m_pDestination);

            // Backtrack: undo the changes to explore other paths
            currentPathLength -= edge->m_Length;
            edge->m_Processed = false;  // Mark the edge as unprocessed
        }
    }
}

// ==============================================
// Greedy + Backtracking TSP Solution
// ==============================================

// Struct to store the distance and edges information for the greedy approach
struct DistanceMatrixElement {
    double distance;                              // Distance between vertices
    std::list<CEdge*> pathEdges;                  // Edges forming the path between vertices
};

// Greedy search variables for solving TSP
DistanceMatrixElement distanceMatrix[20][20];    // Matrix to store distances between vertices
int pathOrder[20];                               // Current order of vertices in the path
int bestPathOrder[20];                           // Best order of vertices found so far
double bestGreedyPathLength = std::numeric_limits<double>::max(); // Length of the best greedy path found
double currentGreedyPathLength = 0.0;            // Length of the current greedy path being explored

// Recursive function to explore paths using a greedy approach
void greedyBacktrackingRec(int depth) {
    // Stop recursion if the current path length is already greater than the best found so far
    if (currentGreedyPathLength >= bestGreedyPathLength) return;

    // If all vertices have been considered, check the path length
    if (depth == orderedVertices.size() - 1) {
        currentGreedyPathLength += fabs(distanceMatrix[pathOrder[depth - 1]][pathOrder[depth]].distance);
        if (currentGreedyPathLength < bestGreedyPathLength) {
            bestGreedyPathLength = currentGreedyPathLength;
            // Store the best path order found
            std::copy(std::begin(pathOrder), std::end(pathOrder), std::begin(bestPathOrder));
        }
        currentGreedyPathLength -= fabs(distanceMatrix[pathOrder[depth - 1]][pathOrder[depth]].distance);
    }
    else {
        // Try all possible vertices at the current depth
        for (int j = 1; j < orderedVertices.size() - 1; j++) {
            bool isInPath = false;
            // Check if the vertex is already in the current path
            for (int k = 1; k < depth && !isInPath; k++) {
                isInPath = (j == pathOrder[k]);
            }
            // If the vertex is not in the path, explore it
            if (!isInPath) {
                pathOrder[depth] = j;
                currentGreedyPathLength += fabs(distanceMatrix[pathOrder[depth - 1]][pathOrder[depth]].distance);
                greedyBacktrackingRec(depth + 1);
                currentGreedyPathLength -= fabs(distanceMatrix[pathOrder[depth - 1]][pathOrder[depth]].distance);
            }
        }
    }
}

// Greedy approach to find an approximate TSP solution
CTrack SalesmanTrackBacktrackingGreedy(CGraph& graph, CVisits& visits) {
    orderedVertices = visits.m_Vertices;
    unorderedVertices = orderedVertices;

    // Initialize the distance matrix using Dijkstra's algorithm for all pairs of vertices
    int i = 0;
    for (CVertex* startVertex : orderedVertices) {
        Dijkstra(graph, startVertex);
        int j = 0;
        for (CVertex* endVertex : unorderedVertices) {
            distanceMatrix[i][j].distance = endVertex->m_DijkstraDistance;

            // Store the edges forming the shortest path using Dijkstra's previous pointers
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
    pathOrder[0] = 0;  // Start vertex is fixed
    pathOrder[orderedVertices.size() - 1] = orderedVertices.size() - 1;  // End vertex is fixed
    greedyBacktrackingRec(1);

    // Reconstruct the best path found by the greedy approach
    CTrack result(&graph);
    for (int m = 0; m < orderedVertices.size() - 1; m++) {
        while (!distanceMatrix[bestPathOrder[m]][bestPathOrder[m + 1]].pathEdges.empty()) {
            result.m_Edges.push_back(distanceMatrix[bestPathOrder[m]][bestPathOrder[m + 1]].pathEdges.front());
            distanceMatrix[bestPathOrder[m]][bestPathOrder[m + 1]].pathEdges.pop_front();
        }
    }

    return result;  // Return the best greedy TSP solution found
}
