#include "pch.h"
#include "Graph.h"
#include "GraphApplicationDlg.h"
#include <vector>
#include <limits>
#include <utility>
#include <cmath>
#include <algorithm>

// Function to find the centroid of the vertices (inlined for performance)
inline CVertex find_centroid(const std::vector<CVertex*>& vertices) {
    double x = 0.0;
    double y = 0.0;
    for (const auto& vertex : vertices) {
        x += vertex->m_Point.m_X;
        y += vertex->m_Point.m_Y;
    }
    x /= vertices.size();
    y /= vertices.size();
    return CVertex("centroid", x, y); // Centroid vertex
}

// Function to find the angle between two points (inlined for performance)
inline double find_angle(const CVertex& centroid, const CVertex* point) {
    return atan2(point->m_Point.m_Y - centroid.m_Point.m_Y, point->m_Point.m_X - centroid.m_Point.m_X);
}

// Function to sort vertices in counter-clockwise order around the centroid
std::vector<CVertex*> sort_vertices_counter_clockwise(std::vector<CVertex*>& vertices) {
    CVertex centroid = find_centroid(vertices);
    std::sort(vertices.begin(), vertices.end(), [&centroid](CVertex* a, CVertex* b) {
        return find_angle(centroid, a) > find_angle(centroid, b); // Notice the '>' for counter-clockwise
        });
    return vertices;
}

// Find the distance between vertices
inline double find_distance_between_vertices(CVertex* v1, CVertex* v2, CVertex* v3) {
    double a = v1->m_Point.m_Y - v2->m_Point.m_Y;
    double b = v2->m_Point.m_X - v1->m_Point.m_X;
    double c = v1->m_Point.m_X * v2->m_Point.m_Y - v2->m_Point.m_X * v1->m_Point.m_Y;
    return std::abs(a * v3->m_Point.m_X + b * v3->m_Point.m_Y + c) / std::sqrt(a * a + b * b);
}

// Divide the vertices of the graph into two segments based on a line from vertex_A to vertex_B
std::pair<std::vector<CVertex*>, std::vector<CVertex*>> create_segment(CVertex* vertex_A, CVertex* vertex_B,
    const std::vector<CVertex*>& vertices) {

    std::vector<CVertex*> above, below;
    if (vertex_B->m_Point.m_X - vertex_A->m_Point.m_X == 0) {
        return { above, below };  // Return empty segments if vertical
    }

    double m = (vertex_B->m_Point.m_Y - vertex_A->m_Point.m_Y) / (vertex_B->m_Point.m_X - vertex_A->m_Point.m_X);
    double c = vertex_A->m_Point.m_Y - m * vertex_A->m_Point.m_X;

    for (auto vertex : vertices) {
        if (vertex->m_Point.m_Y > m * vertex->m_Point.m_X + c) {
            above.push_back(vertex);
        }
        else if (vertex->m_Point.m_Y < m * vertex->m_Point.m_X + c) {
            below.push_back(vertex);
        }
    }

    return { above, below };
}

// Recursive function of the QuickHull Algorithm
std::vector<CVertex*> QuickHullRec(CVertex* vertex_A, CVertex* vertex_B, const std::vector<CVertex*>& segment, char flag) {
    std::vector<CVertex*> designated_half_Hull_list;

    if (segment.empty()) {
        return designated_half_Hull_list;
    }

    double farthest_distance = std::numeric_limits<double>::min();
    CVertex* farthest_point = nullptr;
    for (CVertex* vertex : segment) {
        double distance = find_distance_between_vertices(vertex_A, vertex_B, vertex);
        if (distance > farthest_distance) {
            farthest_distance = distance;
            farthest_point = vertex;
        }
    }
    designated_half_Hull_list.push_back(farthest_point);

    // Point is now in the convex hull, so we can create new segments excluding this point
    std::vector<CVertex*> remaining_vertices = segment;
    remaining_vertices.erase(std::remove(remaining_vertices.begin(), remaining_vertices.end(), farthest_point), remaining_vertices.end());

    auto segment_A = create_segment(vertex_A, farthest_point, remaining_vertices);
    auto segment_B = create_segment(vertex_B, farthest_point, remaining_vertices);

    if (flag == 'N') {
        auto hull_A = QuickHullRec(vertex_A, farthest_point, segment_A.first, flag);
        designated_half_Hull_list.insert(designated_half_Hull_list.end(), hull_A.begin(), hull_A.end());

        auto hull_B = QuickHullRec(farthest_point, vertex_B, segment_B.first, flag);
        designated_half_Hull_list.insert(designated_half_Hull_list.end(), hull_B.begin(), hull_B.end());
    }
    else if (flag == 'S') {
        auto hull_A = QuickHullRec(vertex_A, farthest_point, segment_A.second, flag);
        designated_half_Hull_list.insert(designated_half_Hull_list.end(), hull_A.begin(), hull_A.end());

        auto hull_B = QuickHullRec(farthest_point, vertex_B, segment_B.second, flag);
        designated_half_Hull_list.insert(designated_half_Hull_list.end(), hull_B.begin(), hull_B.end());
    }

    return designated_half_Hull_list;
}

// Main function of the QuickHull Algorithm
CConvexHull QuickHull(CGraph& g) {
    // Constants to use on the recursive split of the algorythm
    const char ABOVE = 'N';
    const char BELOW = 'S';

    // Create the template for our Hull
    CConvexHull Hull(&g);

    // In case that the graph has less than 3 vertices, there is no Hull
    if (g.m_Vertices.empty()) {
        return Hull;
    }
    if (g.GetNVertices() < 3) {
        Hull.m_Vertices.push_back(g.GetVertex(0));
        if (g.GetNVertices() == 2 &&
            (g.GetVertex(0)->m_Point.m_X != g.GetVertex(1)->m_Point.m_X || g.GetVertex(0)->m_Point.m_Y != g.GetVertex(1)->m_Point.m_Y)) {
            Hull.m_Vertices.push_back(g.GetVertex(1));
        }
        return Hull;
    }

    // Copy the vertix pointers of the graph into a custom list
    std::vector<CVertex*> horizontally_sorted_vertices;
    horizontally_sorted_vertices.reserve(g.GetNVertices());
    for (auto point : g.m_Vertices) {
        horizontally_sorted_vertices.push_back(new CVertex(point));
    }

    // Find the minimum and maximum points on the x-axis
    std::sort(horizontally_sorted_vertices.begin(), horizontally_sorted_vertices.end(), [](CVertex* a, CVertex* b) {
        return a->m_Point.m_X < b->m_Point.m_X;
        });

    CVertex* vertex_A = horizontally_sorted_vertices.front();
    CVertex* vertex_B = horizontally_sorted_vertices.back();

    // Remove from the list as they are now in the convex hull
    horizontally_sorted_vertices.erase(horizontally_sorted_vertices.begin());
    horizontally_sorted_vertices.pop_back();

    // Determine points above and below the line (first = above, second = below)
    auto segments = create_segment(vertex_A, vertex_B, horizontally_sorted_vertices);

    // Perform the recursive QuickHull on above and below segments
    auto upper_hull = QuickHullRec(vertex_A, vertex_B, segments.first, ABOVE);
    auto lower_hull = QuickHullRec(vertex_A, vertex_B, segments.second, BELOW);

    // Append the results to the Hull
    std::vector<CVertex*> full_hull;
    full_hull.push_back(vertex_A);
    full_hull.insert(full_hull.end(), upper_hull.begin(), upper_hull.end());
    full_hull.push_back(vertex_B);
    full_hull.insert(full_hull.end(), lower_hull.begin(), lower_hull.end());

    full_hull = sort_vertices_counter_clockwise(full_hull);
    Hull.m_Vertices.insert(Hull.m_Vertices.end(), full_hull.begin(), full_hull.end());

    return Hull;
}