#include "pch.h"
#include "Graph.h"
#include "GraphApplicationDlg.h"
#include <set>
#include <limits>
#include <utility>

// Find the distance between vertices
double find_distance_between_vertices(CVertex* v1, CVertex* v2, CVertex* v3) {
	// using the formula ax + bx + c = 0
	double a = v1->m_Point.m_Y - v2->m_Point.m_Y;
	double b = v2->m_Point.m_X - v1->m_Point.m_X;
	double c = v1->m_Point.m_X * v2->m_Point.m_Y - v2->m_Point.m_X * v1->m_Point.m_Y;

	// use dot product to find the distance between a line and a point
	return abs(a * v3->m_Point.m_X + b * v3->m_Point.m_Y + c) / sqrt(a * a + b * b);
}

// Divide the vertices of the graph in two segments based on 
// a line from vertex_A to vertex_B
pair<list<CVertex*>, list<CVertex*>> create_segment(CVertex* vertex_A, CVertex* vertex_B,
	list<CVertex*> vertices) {
	// We prepare the lists for both segments of vertices
	list<CVertex*> above;
	list<CVertex*> below;

	// If the line isn't vertical, the vertex group can be split
	if (vertex_B->m_Point.m_X - vertex_A->m_Point.m_X == 0) {
		// Return the pair of lists, above and below segments EMPTY
		return pair<list<CVertex*>, list<CVertex*>>(above, below);
	}

	// Calculate m and c from y = mx + c
	double m = (vertex_B->m_Point.m_Y - vertex_A->m_Point.m_Y) / (vertex_B->m_Point.m_X - vertex_A->m_Point.m_X);
	double c = -m * vertex_A->m_Point.m_X + vertex_A->m_Point.m_Y;

	// Loop through each coordinate and place it into the correct list
	for (auto vertex : vertices) {
		// y > mx + c means it is above the line
		if (vertex->m_Point.m_Y > m * vertex->m_Point.m_X + c) {
			above.push_back(vertex);
		}
		// y < mx + c means it is below the line
		else if (vertex->m_Point.m_Y < m * vertex->m_Point.m_X + c) {
			below.push_back(vertex);
		}
	}

	// Return the pair of lists, above and below segments
	return pair<list<CVertex*>, list<CVertex*>>(above, below); // first = above, second = below
}

// Recursive function of the QuickHull Agorythm
list<CVertex*> QuickHullRec(CVertex* vertex_A, CVertex* vertex_B, list<CVertex*> segment, char flag) {
	// Prepare the list
	list<CVertex*> designated_half_Hull_list;

	// Exit case for the recursion
	if (segment.empty() || !vertex_A || !vertex_B) {
		// Return half of the Hull
		return designated_half_Hull_list;
	}

	// Calculate the distance of every point from the line to find the farthest point
	double farthest_distance = numeric_limits<double>::min();
	CVertex* farthest_point = nullptr;
	for (CVertex* vertex : segment) {
		double distance = find_distance_between_vertices(vertex_A, vertex_B, vertex);
		if (distance > farthest_distance) {
			farthest_distance = distance;
			farthest_point = vertex;
		}
	}
	designated_half_Hull_list.push_back(farthest_point);

	// Point is now in the convex hull so remove it from the segment
	segment.remove(farthest_point);

	// Determine the segments formed two lines vA-farthest_point and vB-farthest_point
	// first = above, second = below
	auto segment_A = create_segment(vertex_A, farthest_point, segment);
	auto segment_B = create_segment(vertex_B, farthest_point, segment);

	// Only use the segments in the same direction, the opposite direction is contained in the convex hull
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

	// Return half of the Hull
	return designated_half_Hull_list;
}

// Main function of the QuickHull Algorythm
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
			(g.GetVertex(0)->m_Point.m_X != g.GetVertex(1)->m_Point.m_X ||
			g.GetVertex(0)->m_Point.m_Y != g.GetVertex(1)->m_Point.m_Y)) {
			Hull.m_Vertices.push_back(g.GetVertex(1));
		}
		return Hull;
	}

	// Copy the vertix pointers of the graph into a custom list
	std::list<CVertex*> horizontally_sorted_vertices;
	for (auto point : g.m_Vertices) {
		horizontally_sorted_vertices.push_back(new CVertex(point));
	}

	// Find the minimum and maximum points on the x-axis
	horizontally_sorted_vertices.sort([](CVertex* a, CVertex* b) {
		return a->m_Point.m_X < b->m_Point.m_X;
		});
	CVertex* vertex_A = horizontally_sorted_vertices.front();
	CVertex* vertex_B = horizontally_sorted_vertices.back();
	Hull.m_Vertices.push_back(vertex_A);
	Hull.m_Vertices.push_back(vertex_B);

	// Remove from the list as they are now in the convex hull
	horizontally_sorted_vertices.pop_front();
	horizontally_sorted_vertices.pop_back();

	// Determine points above and below the line (first = above, second = below)
	std::pair<std::list<CVertex*>, std::list<CVertex*>> segments = create_segment(vertex_A, vertex_B, horizontally_sorted_vertices);

	// Perform the recursive QuickHull on above and below segments
	auto upper_hull = QuickHullRec(vertex_A, vertex_B, segments.first, ABOVE);
	auto lower_hull = QuickHullRec(vertex_A, vertex_B, segments.second, BELOW);

	// Append the results to the Hull
	Hull.m_Vertices.insert(Hull.m_Vertices.end(), upper_hull.begin(), upper_hull.end());
	Hull.m_Vertices.insert(Hull.m_Vertices.end(), lower_hull.begin(), lower_hull.end());

	return Hull;
}