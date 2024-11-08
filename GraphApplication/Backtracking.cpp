#include "pch.h"
#include "Graph.h"
#include <set>
#include <queue>

CTrack CamiMesCurt(NULL);
double LongitudCamiMesCurt;

CVertex* pDesti;
CTrack CamiActual(NULL);
double LongitudCamiActual;

list<CVertex*> visites;

struct NodeCami {
	CEdge* m_pEdge;
	NodeCami* m_pAnterior;
};

struct comparator {
	bool operator()(CVertex* pEdge1, CVertex* pEdge2) {
		return pEdge1->m_DijkstraDistance > pEdge2->m_DijkstraDistance;
	}
};
	
struct newComparator {
	bool operator()(pair<int, int> pEdge1, pair<int, int> pEdge2) {
		return pEdge1.second > pEdge2.second;
	}
};

vector<priority_queue<CVertex*, std::vector<CVertex*>, comparator>> dijkstraDistancesMatrix;

struct index {
	index(CVertex* indexFirst){
		indexPrevi->first = indexFirst;
		m_pAnterior = NULL;
	}
	pair<CVertex*, CVertex*>* indexPrevi;
	index* m_pAnterior;
};

class Cell {
public:
	double dijkstraDistance;
	list<CEdge*> tram;
};

vector<vector<Cell>> matriu;
list<int> indexCamiMesCurt;

void SalesmanTrackBacktrackingRec(NodeCami* pAnterior, CVertex* pActual) 
{
	if (LongitudCamiActual < LongitudCamiMesCurt) {
		if (pActual == pDesti) {
			bool totsVisitats = true;
			
			for (CVertex* pVertex : visites) {
				if (pVertex->m_visitCount == 0) {
					totsVisitats = false;
					break;
				}
			}

			if (totsVisitats) {
				CamiMesCurt.Clear();
				while (pAnterior) {
					CamiMesCurt.m_Edges.push_front(pAnterior->m_pEdge);
					pAnterior = pAnterior->m_pAnterior;
				}
				LongitudCamiMesCurt = LongitudCamiActual;
			}
			else {
				NodeCami node;
				node.m_pAnterior = pAnterior;
				for (CEdge* pEdge : pActual->m_Edges) {
					if (!pEdge->m_Used) {
						pEdge->m_Used = true;
						node.m_pEdge = pEdge;
						LongitudCamiActual += pEdge->m_Length;
						SalesmanTrackBacktrackingRec(&node, pEdge->m_pDestination);
						LongitudCamiActual -= pEdge->m_Length;
						pEdge->m_Used = false;
					}
				}
			}
		}
		else {
			pActual->m_visitCount++;

			NodeCami node;
			node.m_pAnterior = pAnterior;
			for (CEdge* pEdge : pActual->m_Edges) {
				if (!pEdge->m_Used) {
					pEdge->m_Used = true;
					node.m_pEdge = pEdge;
					LongitudCamiActual += pEdge->m_Length;
					SalesmanTrackBacktrackingRec(&node, pEdge->m_pDestination);
					LongitudCamiActual -= pEdge->m_Length;
					pEdge->m_Used = false;
				}
			}

			pActual->m_visitCount--;
		}

	}
}

CTrack SalesmanTrackBacktracking(CGraph &graph, CVisits &visits){
	if (!graph.m_Edges.empty()) {
		pDesti = visits.m_Vertices.back();
		visites = visits.m_Vertices;
		visites.pop_back();
		visites.pop_front();

		LongitudCamiActual = numeric_limits<double>::max();

		SalesmanTrackBacktrackingRec(NULL, visits.m_Vertices.front());

		return CamiMesCurt;
	}

	return NULL;
}

void SalesmanTrackBacktrackingGreedyRec(list<int> possibleVisits, int currentIndexNode, list<int> cami) {
	if (LongitudCamiActual < LongitudCamiMesCurt) {
		if (currentIndexNode == pDesti->m_indexMatrix) {
			if (possibleVisits.empty()) {
				LongitudCamiMesCurt = LongitudCamiActual;
				indexCamiMesCurt = cami;
			}
		}
		else {
			priority_queue<pair<int, int>, std::vector<pair<int, int>>, newComparator> nearNodes;

			for (int indexNewNode : possibleVisits){
				nearNodes.push(std::make_pair(indexNewNode, matriu[currentIndexNode][indexNewNode].dijkstraDistance));
			}

			while (!nearNodes.empty()) {

			}
		}
	}
}

CTrack SalesmanTrackBacktrackingGreedy(CGraph& graph, CVisits& visits){
	return CTrack(&graph);
}
