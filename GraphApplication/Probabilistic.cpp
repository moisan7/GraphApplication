#include "pch.h"
#include "Graph.h"
#include <queue>
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm> // For std::swap
#include "Matrix.h"

// SalesmanTrackProbabilistic ==================================================
struct PathWithDistance {
    CTrack track;
    double distance;
    PathWithDistance() : distance(0.0) {}
};

// Calculates an approximate Traveling Salesman Path using a probabilistic approach.
CTrack SalesmanTrackProbabilistic(CGraph& graph, CVisits& visits)
{
    int numVisits = visits.GetNVertices();
    Matrix<PathWithDistance> pathsBetweenVisits(numVisits, numVisits);
    CTrack bestTrack(&graph); // Best track found so far.
    double bestDistance = numeric_limits<double>::max(); // Distance of the best track.

    // Precalculates the shortest paths and their distances between all pairs of visits.
    int startVisitIndex = 0;
    for (CVertex* startVertex : visits.m_Vertices) {
        Dijkstra(graph, startVertex);
        int endVisitIndex = 0;
        for (CVertex* endVertex : visits.m_Vertices) {
            pathsBetweenVisits(startVisitIndex, endVisitIndex).distance = endVertex->m_DijkstraDistance;

            // Reconstructs the path using Dijkstra's previous pointers.
            pathsBetweenVisits(startVisitIndex, endVisitIndex).track = CTrack(&graph);
            CVertex* currentVertex = endVertex;
            while (currentVertex->m_pDijkstraPrevious) {
                pathsBetweenVisits(startVisitIndex, endVisitIndex).track.m_Edges.push_front(currentVertex->m_pDijkstraPrevious);
                currentVertex = currentVertex->m_pDijkstraPrevious->m_pOrigin;
            }
            endVisitIndex++;
        }
        startVisitIndex++;
    }

    // Performs multiple iterations to find a good solution.
    std::random_device rd;
    std::mt19937 generator(rd());

    for (int iteration = 0; iteration < 1500 * numVisits; ++iteration) {
        std::vector<int> visitOrderIndices;
        visitOrderIndices.push_back(0); // Always starts with the first visit.

        // Creates a list of visit indices for random permutation.
        std::vector<int> remainingVisitsIndices;
        int numRemainingVisits = numVisits - 2; // Excluding the first and the last.
        for (int i = 0; i < numRemainingVisits; ++i) {
            remainingVisitsIndices.push_back(i + 1);
        }

        // Randomly permutes the order of intermediate visits.
        std::shuffle(remainingVisitsIndices.begin(), remainingVisitsIndices.end(), generator);

        // Adds the permuted visits to the route order.
        for (int visitIndex : remainingVisitsIndices) {
            visitOrderIndices.push_back(visitIndex);
        }
        visitOrderIndices.push_back(visits.m_Vertices.size() - 1); // Ends with the last visit.

        double currentDistance = 0;
        // Calculates the distance of the current route.
        for (size_t i = 0; i < visitOrderIndices.size() - 1; ++i) {
            currentDistance += pathsBetweenVisits(visitOrderIndices[i], visitOrderIndices[i + 1]).distance;
        }

        // Applies a local search to try and improve the current route.
        bool improved = true;
        while (improved) {
            improved = false;
            double bestLocalDistance = currentDistance;
            int bestSwapRow = -1, bestSwapCol = -1;

            // Evaluates all possible swaps to find the best improvement.
            for (size_t i = 1; i < visitOrderIndices.size() - 1; ++i) {
                for (size_t j = i + 1; j < visitOrderIndices.size() - 1; ++j) { // Avoid redundant and self-swaps
                    std::swap(visitOrderIndices[i], visitOrderIndices[j]); // Apply the swap

                    double swappedDistance = 0;
                    for (size_t k = 0; k < visitOrderIndices.size() - 1; ++k) {
                        swappedDistance += pathsBetweenVisits(visitOrderIndices[k], visitOrderIndices[k + 1]).distance;
                    }

                    if (swappedDistance < bestLocalDistance) {
                        bestLocalDistance = swappedDistance;
                        bestSwapRow = i;
                        bestSwapCol = j;
                    }

                    std::swap(visitOrderIndices[i], visitOrderIndices[j]); // Revert the swap
                }
            }

            // If a better swap is found, apply it.
            if (bestSwapRow != -1) {
                std::swap(visitOrderIndices[bestSwapRow], visitOrderIndices[bestSwapCol]);
                currentDistance = bestLocalDistance;
                improved = true;
            }
        }

        // If the current route is better than the best found so far, update the best route.
        if (currentDistance < bestDistance) {
            bestTrack.m_Edges.clear();
            bestDistance = currentDistance;
            for (size_t i = 0; i < visitOrderIndices.size() - 1; ++i) {
                bestTrack.Append(pathsBetweenVisits(visitOrderIndices[i], visitOrderIndices[i + 1]).track);
            }
        }
    }

    return bestTrack;
}