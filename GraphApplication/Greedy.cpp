#include "pch.h"
#include "Graph.h"
// SalesmanTrackGreedy =========================================================

// Function to construct the greedy path for the traveling salesman problem
CTrack SalesmanTrackGreedy(CGraph& graph, CVisits& visits) {
    // Check if there are visits to process
    if (visits.m_Vertices.empty()) return CTrack(&graph);

    CTrack result(&graph);

    // Initialize starting vertex
    CVertex* pStartVertex = visits.m_Vertices.front();
    CVertex* pEndVertex = visits.m_Vertices.back();

    // Initialize candidates (all vertices except the first and last one)
    list<CVertex*> candidates(visits.m_Vertices.begin(), visits.m_Vertices.end());
    candidates.pop_front(); // Remove the first (start) vertex
    candidates.pop_back();  // Remove the last (end) vertex

    // Current vertex to track the path
    CVertex* currentVertex = pStartVertex;

    while (!candidates.empty()) {
        // Perform Dijkstra from the current vertex to find shortest paths
        DijkstraQueue(graph, currentVertex);

        // Find the nearest vertex from the current position
        CVertex* nearestVertex = nullptr;
        double minDistance = numeric_limits<double>::infinity();

        for (auto candidate : candidates) {
            if (candidate->m_DijkstraDistance < minDistance) {
                nearestVertex = candidate;
                minDistance = candidate->m_DijkstraDistance;
            }
        }

        // If a nearest vertex is found, build the path to it
        if (nearestVertex) {
            // Trace back using Dijkstra's path pointers to construct the path
            CVertex* v = nearestVertex;
            list<CEdge*> tempEdges;

            while (v != currentVertex && v->m_pDijkstraPrevious != nullptr) {
                tempEdges.push_front(v->m_pDijkstraPrevious);
                v = v->m_pDijkstraPrevious->m_pOrigin;
            }

            // Add all the edges to the track result
            for (auto pEdge : tempEdges) {
                result.AddLast(pEdge);
            }

            // Update current vertex to nearest vertex and remove it from candidates
            currentVertex = nearestVertex;
            candidates.remove(nearestVertex);
        }
    }

    // Finally, go from the last selected candidate to the final vertex
    if (currentVertex != pEndVertex) {
        DijkstraQueue(graph, currentVertex);

        // Trace back using Dijkstra's path pointers to construct the path
        CVertex* v = pEndVertex;
        list<CEdge*> tempEdges;

        while (v != currentVertex && v->m_pDijkstraPrevious != nullptr) {
            tempEdges.push_front(v->m_pDijkstraPrevious);
            v = v->m_pDijkstraPrevious->m_pOrigin;
        }

        // Add all the edges to the track result
        for (auto pEdge : tempEdges) {
            result.AddLast(pEdge);
        }
    }

    return result;
}
