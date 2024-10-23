#include "pch.h"
#include "Graph.h"
#include <queue> // Required for using priority_queue
#include <limits> // Required for using numeric_limits

// =============================================================================
// Dijkstra ====================================================================
// =============================================================================

void Dijkstra(CGraph& graph, CVertex* pStart) {
    // Inicialitzar les distàncies dels vèrtexs a infinit excepte la del vèrtex pStart que serà 0
    for (CVertex& vertex : graph.m_Vertices) {
        vertex.m_DijkstraDistance = (vertex.m_Name == pStart->m_Name) ? 0.0 : std::numeric_limits<double>::max();
        vertex.m_DijkstraVisit = false;
        vertex.m_pDijkstraPrevious = nullptr;
    }

    CVertex* pActual = pStart;

    while (pActual != nullptr) {
        // Recorre tots els veïns v de pActual
        for (CEdge* pEdge : pActual->m_Edges) {
            CVertex* v = pEdge->m_pDestination;

            // Si la distància de v és més grossa que la distància del vèrtex actual més la longitud de l’aresta que els uneix, actualitzar la distancia de v
            if (v->m_DijkstraDistance == -1.0 || v->m_DijkstraDistance > pActual->m_DijkstraDistance + pEdge->m_Length) {
                v->m_DijkstraDistance = pActual->m_DijkstraDistance + pEdge->m_Length;
                v->m_pDijkstraPrevious = pEdge;
            }
        }

        // Marcar pActual com visitat
        pActual->m_DijkstraVisit = true;

        // pActual = vèrtex no visitat amb distancia més petita o NULL si no hi ha vèrtexs no visitats
        pActual = nullptr;
        double minDistance = -1.0;
        for (CVertex& vertex : graph.m_Vertices) {
            if (!vertex.m_DijkstraVisit && vertex.m_DijkstraDistance != -1.0) {
                if (pActual == nullptr || vertex.m_DijkstraDistance < minDistance) {
                    pActual = &vertex;
                    minDistance = vertex.m_DijkstraDistance;
                }
            }
        }
    }
}

// =============================================================================
// DijkstraQueue ===============================================================
// =============================================================================

#include <algorithm> // For sorting edges

// Comparator for the priority queue (min-heap)
struct VertexComparator {
    bool operator()(const CVertex* a, const CVertex* b) {
        return a->m_DijkstraDistance > b->m_DijkstraDistance;
    }
};

void DijkstraQueue(CGraph& graph, CVertex* pStart) {
    // Inicialitzar les distàncies dels vèrtexs a infinit excepte la del vèrtex start que serà 0
    for (CVertex& vertex : graph.m_Vertices) {
        vertex.m_DijkstraDistance = (vertex.m_Name == pStart->m_Name) ? 0.0 : std::numeric_limits<double>::max();
        vertex.m_DijkstraVisit = false;
        vertex.m_pDijkstraPrevious = nullptr;
    }

    // Declare and initialize the priority queue
    std::priority_queue<CVertex*, std::vector<CVertex*>, VertexComparator> queue;
    queue.push(pStart);

    while (!queue.empty()) {
        // Get the vertex with the smallest distance and remove it from the queue
        CVertex* va = queue.top();
        queue.pop();

        // If the vertex is not visited
        if (!va->m_DijkstraVisit) {
            // Sort edges by destination vertex name for consistent edge processing
            va->m_Edges.sort([](const CEdge* a, const CEdge* b) {
                return a->m_pDestination->m_Name < b->m_pDestination->m_Name;
                });

            // Iterate through neighbors of va
            for (CEdge* pEdge : va->m_Edges) {
                CVertex* v = pEdge->m_pDestination;

                // If the distance to v is greater than the distance to va plus the edge length, update the distance and add v to the queue
                if (v->m_DijkstraDistance > va->m_DijkstraDistance + pEdge->m_Length) {
                    v->m_DijkstraDistance = va->m_DijkstraDistance + pEdge->m_Length;
                    v->m_pDijkstraPrevious = pEdge;
                    queue.push(v);
                }
            }

            // Mark va as visited
            va->m_DijkstraVisit = true;
        }
    }
}