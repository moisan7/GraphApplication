#include "pch.h"
#include "Graph.h"
#include <queue>
#include <iostream>
#include <iomanip>
#include <limits> // Required for numeric_limits
#include <memory> // For smart pointers

//=================================================================================
// Branch and Bound Node Structure (Version 1 - Dijkstra-based)
//=================================================================================

/// @brief Represents a node in the Branch and Bound search tree for the first version,
///        using Dijkstra's distance for node ordering.
struct CBBNode1 {
    double lower_bound;       ///< The lower bound cost to reach this node.
    bool visited[20];        ///< Array indicating which vertices have been visited.
    unsigned vertex_indices[20]; ///< Array storing indices of visited vertices in order.
    unsigned index;            ///< Current index in the vertex_indices array.
    CVertex* current_vertex;    ///< Pointer to the current vertex represented by this node.
};

/// @brief Custom comparator for the priority queue used in the first Branch and Bound version.
///        Orders nodes based on the Dijkstra distance of their current vertex.
struct comparator1 {
    /// @brief Compares two CBBNode1 pointers based on Dijkstra distance.
    /// @param s1 Pointer to the first CBBNode1.
    /// @param s2 Pointer to the second CBBNode1.
    /// @return True if the Dijkstra distance of s1's vertex is less than s2's, false otherwise.
    bool operator()(const CBBNode1* s1, const CBBNode1* s2) {
        return s1->current_vertex->m_DijkstraDistance < s2->current_vertex->m_DijkstraDistance;
    }
};

//=================================================================================
// Salesman Track Branch and Bound (Version 1)
//=================================================================================

/// @brief Finds the shortest salesman track using Branch and Bound with Dijkstra-based node ordering.
///
/// @param graph The graph representing the network.
/// @param visits The set of vertices that must be visited.
/// @return The shortest track visiting all specified vertices.
CTrack SalesmanTrackBranchAndBound1(CGraph& graph, CVisits& visits) {
    // Initialize variables
    unsigned num_vertices = 0;
    double distance_matrix[20][20];
    CVertex* vertices_array[20];
    double new_lower_bound, lower_bound;
    unsigned column, current_vertex_index;
    double shortest_path_length = std::numeric_limits<double>::max();
    CBBNode1 initial_node, * current_node, * best_path_node = (CBBNode1*)malloc(sizeof(CBBNode1));
    CTrack shortest_track(&graph);
    CVertex* start_vertex, * end_vertex;

    // Initialize the initial node with the starting vertex
    initial_node.lower_bound = 0;
    initial_node.vertex_indices[0] = 0; // Index of the starting vertex
    initial_node.visited[0] = true;
    initial_node.index = 1;
    initial_node.current_vertex = visits.m_Vertices.front();

    // Create a priority queue for nodes, ordered by comparator1 (Dijkstra distance)
    std::priority_queue<CBBNode1*, std::vector<CBBNode1*>, comparator1> node_queue;
    node_queue.push(&initial_node);

    // Copy the vertices from visits.m_Vertices to vertices_array for indexed access
    for (CVertex* vertex : visits.m_Vertices) {
        vertices_array[num_vertices] = vertex;
        num_vertices++;
    }

    // Build the distance matrix using Dijkstra's algorithm
    for (unsigned row = 0; row < num_vertices - 1; row++) {
        Dijkstra(graph, vertices_array[row]);
        column = 0;
        for (CVertex* vertex : visits.m_Vertices) {
            distance_matrix[row][column] = vertex->m_DijkstraDistance;
            column++;
        }
    }

    // Handle the special case of only two vertices to visit
    if (visits.m_Vertices.size() == 2) {
        start_vertex = visits.m_Vertices.front();
        end_vertex = visits.m_Vertices.back();
        for (CEdge* edge : start_vertex->m_Edges) {
            if (edge->m_pDestination == end_vertex) {
                shortest_track.AddLast(edge);
                return shortest_track;
            }
        }
    }

    // Calculate the minimum and maximum distances for each column in the matrix
    column = 0;
    while (column < num_vertices) {
        distance_matrix[18][column] = shortest_path_length; // Initialize minimum distance to max
        distance_matrix[19][column] = 0;                  // Initialize maximum distance to 0
        for (unsigned row = 0; row < num_vertices - 1; row++) {
            double matrix_value = distance_matrix[row][column];
            if (matrix_value < distance_matrix[18][column] && row != column) {
                distance_matrix[18][column] = matrix_value; // Update minimum distance
            }
            if (matrix_value > distance_matrix[19][column]) {
                distance_matrix[19][column] = matrix_value; // Update maximum distance
            }
        }
        column++;
    }

    // Calculate the initial lower bound
    column = 0;
    while (column < num_vertices) {
        initial_node.lower_bound += distance_matrix[18][column];
        column++;
    }
    for (unsigned i = 1; i < num_vertices - 1; i++) initial_node.visited[i] = false;

    // Main Branch and Bound loop
    while (!node_queue.empty()) {
        current_node = node_queue.top();
        node_queue.pop();
        current_vertex_index = current_node->vertex_indices[current_node->index - 1];
        lower_bound = current_node->lower_bound;

        // Check if the current node represents a complete path
        if (current_node->index == num_vertices - 1) {
            if (lower_bound + distance_matrix[current_vertex_index][num_vertices - 1] - distance_matrix[18][num_vertices - 1] < shortest_path_length) {
                *best_path_node = *current_node; // Store the best path found so far
                shortest_path_length = lower_bound + distance_matrix[current_vertex_index][num_vertices - 1] - distance_matrix[18][num_vertices - 1];
            }
        }
        else {
            // Expand the current node by exploring unvisited vertices
            for (unsigned i = 1; i < num_vertices - 1; i++) {
                new_lower_bound = current_node->lower_bound + distance_matrix[current_vertex_index][i] - distance_matrix[18][i];
                if (!(current_node->visited[i]) && new_lower_bound < shortest_path_length) {
                    CBBNode1* new_node = (CBBNode1*)malloc(sizeof(CBBNode1));
                    *new_node = *current_node;
                    new_node->vertex_indices[new_node->index] = i;
                    new_node->visited[i] = true;
                    new_node->lower_bound = new_lower_bound;
                    new_node->index++;
                    new_node->current_vertex = vertices_array[i];
                    node_queue.push(new_node);
                }
            }
        }
    }

    best_path_node->vertex_indices[num_vertices - 1] = num_vertices - 1;

    // Construct the shortest track from the best path found
    for (unsigned j = num_vertices - 1; j > 0; j--) {
        start_vertex = vertices_array[best_path_node->vertex_indices[j - 1]];
        end_vertex = vertices_array[best_path_node->vertex_indices[j]];
        Dijkstra(graph, start_vertex);

        while (end_vertex != start_vertex) {
            shortest_track.AddFirst(end_vertex->m_pDijkstraPrevious);
            end_vertex = end_vertex->m_pDijkstraPrevious->m_pOrigin;
        }
    }
    return shortest_track;
}

//============================================================================
// Branch and Bound Node Structure (Version 2 - Lower and Upper Bound-based)
//============================================================================

/// @brief Represents a node in the Branch and Bound search tree for the second version,
///        using both lower and upper bounds for pruning.
struct CBBNode2 {
    double upper_bound;      ///< The upper bound cost to reach the destination from this node.
    double lower_bound;      ///< The lower bound cost to reach this node.
    bool visited[20];       ///< Array indicating which vertices have been visited.
    unsigned vertex_indices[20]; ///< Array storing indices of visited vertices in order.
    unsigned index;            ///< Current index in the vertex_indices array.
};

/// @brief Custom comparator for the priority queue used in the second Branch and Bound version.
///        Orders nodes based on their lower bound cost.
struct comparator2 {
    /// @brief Compares two CBBNode2 pointers based on their lower bounds.
    /// @param s1 Pointer to the first CBBNode2.
    /// @param s2 Pointer to the second CBBNode2.
    /// @return True if the lower bound of s1 is less than s2's, false otherwise.
    bool operator()(const CBBNode2* s1, const CBBNode2* s2) {
        return s1->lower_bound < s2->lower_bound;
    }
};

//============================================================================
// Salesman Track Branch and Bound (Version 2)
//============================================================================

/// @brief Finds the shortest salesman track using Branch and Bound with lower and upper bounds.
///
/// @param graph The graph representing the network.
/// @param visits The set of vertices that must be visited.
/// @return The shortest track visiting all specified vertices.
CTrack SalesmanTrackBranchAndBound2(CGraph& graph, CVisits& visits) {
    // Variables
    unsigned num_vertices = 0;
    double distance_matrix[20][20];
    CVertex* vertices_array[20];
    unsigned column;
    double shortest_path_length = std::numeric_limits<double>::max();
    CBBNode2 initial_node, * current_node, * best_path_node = (CBBNode2*)malloc(sizeof(CBBNode2));
    unsigned current_vertex_index;
    double new_lower_bound, lower_bound;
    CTrack shortest_track(&graph);
    CVertex* start_vertex, * end_vertex;
    constexpr double epsilon = 0.00001; // Small error margin for floating-point comparisons

    // Initialize values for the starting node
    initial_node.lower_bound = 0;
    initial_node.upper_bound = 0;
    initial_node.vertex_indices[0] = 0; // Index of the starting vertex
    initial_node.visited[0] = true;
    initial_node.index = 1;

    // Create a priority queue for nodes, ordered by comparator2 (lower bound)
    std::priority_queue<CBBNode2*, std::vector<CBBNode2*>, comparator2> node_queue;
    node_queue.push(&initial_node);

    // Copy the vertices from visits.m_Vertices to vertices_array for indexed access
    for (CVertex* vertex : visits.m_Vertices) {
        vertices_array[num_vertices] = vertex;
        num_vertices++;
    }

    // Build the distance matrix using Dijkstra's algorithm
    for (unsigned row = 0; row < num_vertices - 1; row++) {
        Dijkstra(graph, vertices_array[row]);
        column = 0;
        for (CVertex* vertex : visits.m_Vertices) {
            distance_matrix[row][column] = vertex->m_DijkstraDistance;
            column++;
        }
    }

    // Handle the special case of only two vertices to visit
    if (visits.m_Vertices.size() == 2) {
        start_vertex = visits.m_Vertices.front();
        end_vertex = visits.m_Vertices.back();
        for (CEdge* edge : start_vertex->m_Edges) {
            if (edge->m_pDestination == end_vertex) {
                shortest_track.AddLast(edge);
                return shortest_track;
            }
        }
    }

    // Calculate the minimum (lower bound) and maximum (upper bound) distances for each column in the matrix
    // Store the minimums in row 18 and maximums in row 19
    column = 0;
    while (column < num_vertices) {
        distance_matrix[18][column] = shortest_path_length; // Initialize minimum distance to max
        distance_matrix[19][column] = 0;                  // Initialize maximum distance to 0
        for (unsigned row = 0; row < num_vertices - 1; row++) {
            double matrix_value = distance_matrix[row][column];
            if (matrix_value < distance_matrix[18][column] && row != column) {
                distance_matrix[18][column] = matrix_value; // Update minimum distance (lower bound)
            }
            if (matrix_value > distance_matrix[19][column]) {
                distance_matrix[19][column] = matrix_value; // Update maximum distance (upper bound)
            }
        }
        column++;
    }

    // Calculate the initial lower and upper bounds
    column = 0;
    while (column < num_vertices) {
        initial_node.upper_bound += distance_matrix[19][column];  // Total upper bound
        initial_node.lower_bound += distance_matrix[18][column];  // Total lower bound
        column++;
    }

    for (unsigned i = 1; i < num_vertices - 1; i++) initial_node.visited[i] = false;

    // Main Branch and Bound loop
    while (!node_queue.empty()) {
        current_node = node_queue.top();
        node_queue.pop();
        current_vertex_index = current_node->vertex_indices[current_node->index - 1];
        lower_bound = current_node->lower_bound;

        // Check if the current node represents a complete path
        if (current_node->index == num_vertices - 1) {
            // If the path's cost is less than the current upper bound (plus a small error margin), update the best path
            if (lower_bound + distance_matrix[current_vertex_index][num_vertices - 1] - distance_matrix[18][num_vertices - 1] < shortest_path_length + epsilon) {
                *best_path_node = *current_node; // Store the best path found so far
                shortest_path_length = lower_bound + distance_matrix[current_vertex_index][num_vertices - 1] - distance_matrix[18][num_vertices - 1]; // Update the upper bound
            }
        }
        else {
            // Expand the current node by exploring unvisited vertices
            for (unsigned i = 1; i < num_vertices - 1; i++) {
                new_lower_bound = current_node->lower_bound + distance_matrix[current_vertex_index][i] - distance_matrix[18][i];
                if (!(current_node->visited[i]) && new_lower_bound < shortest_path_length + epsilon) {
                    CBBNode2* new_node = (CBBNode2*)malloc(sizeof(CBBNode2));
                    *new_node = *current_node;
                    new_node->vertex_indices[new_node->index] = i;
                    new_node->visited[i] = true;
                    new_node->lower_bound = new_lower_bound;
                    new_node->upper_bound += distance_matrix[current_vertex_index][i] - distance_matrix[19][i];
                    new_node->index++;

                    // Update the upper bound if the new node's upper bound is lower
                    if (new_node->upper_bound < shortest_path_length + epsilon) {
                        shortest_path_length = new_node->upper_bound;
                    }

                    node_queue.push(new_node);
                }
            }
        }
    }

    best_path_node->vertex_indices[num_vertices - 1] = num_vertices - 1;

    // Construct the shortest track from the best path found
    for (unsigned j = num_vertices - 1; j > 0; j--) {
        start_vertex = vertices_array[best_path_node->vertex_indices[j - 1]];
        end_vertex = vertices_array[best_path_node->vertex_indices[j]];
        Dijkstra(graph, start_vertex);

        while (end_vertex != start_vertex) {
            shortest_track.AddFirst(end_vertex->m_pDijkstraPrevious);
            end_vertex = end_vertex->m_pDijkstraPrevious->m_pOrigin;
        }
    }
    return shortest_track;
}

// SalesmanTrackBranchAndBound3 ===================================================


CTrack SalesmanTrackBranchAndBound3(CGraph& graph, CVisits& visits)
{
	return CTrack(&graph);
}

// SalesmanTrackBranchAndBound4 ===================================================


CTrack SalesmanTrackBranchAndBound4(CGraph& graph, CVisits& visits)
{
	return CTrack(&graph);
}
