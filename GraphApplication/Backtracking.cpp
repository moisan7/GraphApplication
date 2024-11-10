#include "pch.h"
#include "Graph.h"
#include <set>
#include <queue>


// =============================================================================
// SalesmanTrackBacktracking ===================================================
// =============================================================================


//we create these global variables that we don't modify in each call of the recursive function
CTrack CamiMesCurt2(NULL);
double LongitudCamiMesCurt2;

CVertex* pDesti2;
CTrack CamiActual2(NULL);
double LongitudCamiActual2;


list<CVertex*> VISITES2;

//we have created global variable to take advantage of the local memory in the functions, saving time accesing to the values
struct NodeCami2 {
	CEdge* m_pEdge;
	NodeCami2* m_pAnterior;
};


//versio de SALESMAN backtracking pur amb CTRACK per guarda cami actual


void SalesmanTrackBacktrackingRec(CVertex* pActual)
{
	//we are not going to explore a vertex if the current lengh is already longer that the shortest way
	if (LongitudCamiActual2 < LongitudCamiMesCurt2)
	{
		if (pActual == pDesti2) {


			bool correct(true);
			//int correct(0);

			//we check that we passed all the vertex that we have to visit
			for (CVertex* v1 : VISITES2)
			{
				//if we didn't pass through one vertex of visits, the way is no a possible answer

				if (v1->m_count == 0)
				{
					correct = false;
					break;
				}


				//solucio proposta per el profe, comprobar que surt una aresta dels vertex que hem de visitar
				//for (CEdge* edge:v1->m_Edges)
				//{
					//if (edge->m_Used) {

						//correct++;
						//break;
					//}
				//}


			}

			//if the way is a possible answer
			if (correct)
			{

				//CamiMesCurt2.Clear();
				//while (pAnterior) {
					//CamiMesCurt2.m_Edges.push_front(pAnterior->m_pEdge);
					//pAnterior = pAnterior->m_pAnterior;
				//}


				CamiMesCurt2 = CamiActual2;

				LongitudCamiMesCurt2 = LongitudCamiActual2;


			}
			else
			{
				//node to  save actual way
				//NodeCami2 node;
				//node.m_pAnterior = pAnterior;

				for (CEdge* pE : pActual->m_Edges) {
					//if (!pE->m_pDestination->m_JaHePassat) {

					if (!pE->m_Used)
					{
						pE->m_Used = true;
						//node.m_pEdge = pE;

						CamiActual2.m_Edges.push_back(pE);

						LongitudCamiActual2 += pE->m_Length;


						//cout << "WE ARE AT NODE: " << pActual->m_Name << endl;
						//cout << "Edge: " << pE->m_Name << endl;
						SalesmanTrackBacktrackingRec(pE->m_pDestination);
						//cout << "WE ARE BACK TO THE NODE: " << pActual->m_Name << endl;

						CamiActual2.m_Edges.pop_back();
						LongitudCamiActual2 -= pE->m_Length;
						pE->m_Used = false;
					}
				}

			}

		}
		else
		{

			//pActual->m_JaHePassat = true;

			//we pass +1 time trhough the vertex
			pActual->m_count++;

			//node to  save actual way
			//NodeCami2 node;
			//node.m_pAnterior = pAnterior;
			for (CEdge* pE : pActual->m_Edges) {

				if (!pE->m_Used)
				{
					pE->m_Used = true;
					//node.m_pEdge = pE;

					CamiActual2.m_Edges.push_back(pE);

					LongitudCamiActual2 += pE->m_Length;


					//cout << "WE ARE AT NODE: " << pActual->m_Name << endl;
					//cout << "Edge: " << pE->m_Name << endl;
					SalesmanTrackBacktrackingRec(pE->m_pDestination);
					//cout << "WE ARE BACK TO THE NODE: " << pActual->m_Name << endl;

					CamiActual2.m_Edges.pop_back();
					LongitudCamiActual2 -= pE->m_Length;
					pE->m_Used = false;
				}

			}

			//we get to the previous state
			pActual->m_count--;

		}
	}


	//return CamiMesCurt2;
}


//versio per fer crida de recursiu amb Ctrack per guardar cami actual

CTrack SalesmanTrackBacktracking(CGraph& graph, CVisits& visits)
{
	if (!graph.m_Edges.empty())
	{
		pDesti2 = visits.m_Vertices.back();
		//we take out the inicial node and the destination
		VISITES2 = visits.m_Vertices;
		VISITES2.pop_back();
		VISITES2.pop_front();

		//we initialize the shortest way
		LongitudCamiMesCurt2 = numeric_limits<double>::max();
		//NodeCami2 first;
		//we add the first
		//first.m_pEdge = visits.m_Vertices.front()->m_Edges.front();
		//we took the first neighbour that is the destination of the first edge of the beginner

		//actual = CTrack(&graph);

		SalesmanTrackBacktrackingRec(visits.m_Vertices.front());

		return CamiMesCurt2;
	}


	return NULL;
}

/*
//versio versio de SALESMAN backtracking puramb amb llista enlla�ada NodeCami2 per guarda cami actual

void SalesmanTrackBacktrackingRec(NodeCami2* pAnterior, CVertex* pActual)
{
	//we are not going to explore a vertex if the current lengh is already longer that the shortest way
	if (LongitudCamiActual2 < LongitudCamiMesCurt2)
	{

		if (pActual == pDesti2) {


			//this variable is to control that we have visited all the vertex that we must
			bool correct = true;

			//we check that we passed all the vertex that we have to visit
			for (CVertex* v1 : VISITES2)
			{
				//if we didn't pass through one vertex of visits, the way is no a possible answer
				if (v1->m_count == 0)
				{
					correct = false;
					break;
				}
			}

			//if the way is a possible answer
			if (correct)
			{
				//we update the optimal way with the way that we have found
				CamiMesCurt2.Clear();
				while (pAnterior) {
					CamiMesCurt2.m_Edges.push_front(pAnterior->m_pEdge);
					pAnterior = pAnterior->m_pAnterior;
				}
				LongitudCamiMesCurt2 = LongitudCamiActual2;



			}
			else
			{
				//node to  save actual way
				NodeCami2 node;
				node.m_pAnterior = pAnterior;
				for (CEdge* pE : pActual->m_Edges) {

					if (!pE->m_Used)
					{
						pE->m_Used = true;
						node.m_pEdge = pE;
						LongitudCamiActual2 += pE->m_Length;
						//cout << "Edge: " << pE->m_Name << endl;
						SalesmanTrackBacktrackingRec(&node, pE->m_pDestination);
						//cout << "WE ARE BACK TO THE NODE: " << pActual->m_Name << endl;
						LongitudCamiActual2 -= pE->m_Length;
						pE->m_Used = false;
					}
				}

			}

		}
		else
		{

			pActual->m_count++;

			//node to  save actual way
			NodeCami2 node;
			node.m_pAnterior = pAnterior;
			for (CEdge* pE : pActual->m_Edges) {


				if (!pE->m_Used)
				{
					pE->m_Used = true;
					node.m_pEdge = pE;
					LongitudCamiActual2 += pE->m_Length;
					//cout << "WE ARE AT NODE: " << pActual->m_Name << endl;
					//cout << "Edge: " << pE->m_Name << endl;
					SalesmanTrackBacktrackingRec(&node, pE->m_pDestination);
					//cout << "WE ARE BACK TO THE NODE: " << pActual->m_Name << endl;
					LongitudCamiActual2 -= pE->m_Length;
					pE->m_Used = false;
				}

			}

			//we get to the previous state
			pActual->m_count--;

		}
	}


}


CTrack SalesmanTrackBacktracking(CGraph& graph, CVisits& visits)
{
	if (!graph.m_Edges.empty())
	{
		pDesti2 = visits.m_Vertices.back();
		//we take out the inicial node and the destination
		VISITES2 = visits.m_Vertices;
		VISITES2.pop_back();
		VISITES2.pop_front();

		//we initialize the shortest way 
		LongitudCamiMesCurt2 = numeric_limits<double>::max();

		SalesmanTrackBacktrackingRec(NULL, visits.m_Vertices.front());

		return CamiMesCurt2;
	}


	return NULL;
}
*/


// =============================================================================
// SalesmanTrackBacktrackingGreedy =============================================
// =============================================================================


struct comparator {
	bool operator()(CVertex* pE1, CVertex* pE2) {
		return pE1->m_DijkstraDistance > pE2->m_DijkstraDistance;
	}
};


//we are going to use a priority queue because we are going to explore the vertex in order of the Dijkstra Distance. 
	//from the cheapest to the most expensive
vector<priority_queue<CVertex*, std::vector<CVertex*>, comparator>> dijkstraDistancesMatrix;



struct index {

	index(CVertex* indexFirst)
	{
		indexPrevi->first = indexFirst;
	}

	pair<CVertex*, CVertex*>* indexPrevi;

	index* m_pAnterior;
};


CVisits* VISITES;
struct newComparator {
	bool operator()(pair<int, int> pE1, pair<int, int> pE2) {
		return pE1.second > pE2.second;
	}
};


class Cell {
public:
	double dijkstraDistance;
	list<CEdge*> tram;
};


vector<vector<Cell>> matriu;
list<int> indexCamiMesCurt;

void SalesmanTrackBacktrackingGreedyRec(list<int> possibleVisits, int currentIndexNode, list<int> cami) {

	//if the current lenght is bigger we are not to continue
	if (LongitudCamiActual2 < LongitudCamiMesCurt2) {

		if (currentIndexNode == pDesti2->m_indexMatrix)
		{
			//if we don't have possible nodes to visit, it means that we have visit all of them
			if (possibleVisits.empty())
			{
				//we update the optimal way with the current
				LongitudCamiMesCurt2 = LongitudCamiActual2;

				//we update the optimal way with the way that we have found
				indexCamiMesCurt = cami;
			}


		}
		else {

			//we create priority queue to explore the near node from the current one
			priority_queue<pair<int, int>, std::vector<pair<int, int>>, newComparator> nearNodes;

			//we order the nodes in order of near to current node
			for (int indexNewNode : possibleVisits)
			{
				nearNodes.push(std::make_pair(indexNewNode, matriu[currentIndexNode][indexNewNode].dijkstraDistance));
			}

			while (!nearNodes.empty())
			{


				//we add the weight of the optimal way section to the total weight of the way
				LongitudCamiActual2 += nearNodes.top().second;


				//we update the current node with the new one
				cami.push_back(nearNodes.top().first);


				//we erase the new visited node from the possible visits and we added to the cami
				possibleVisits.erase(find(possibleVisits.begin(), possibleVisits.end(), nearNodes.top().first));


				//we sent the nearest vertex to the current vertex to the recursive function to explore it
				SalesmanTrackBacktrackingGreedyRec(possibleVisits, nearNodes.top().first, cami);//the nearest vertex is on the top of the priority queue

				//we take out the weight of the optimal way section to do a step back
				LongitudCamiActual2 -= nearNodes.top().second;
				cami.pop_back();
				possibleVisits.push_back(nearNodes.top().first);

				//we take it out and we try with the next near vertex
				nearNodes.pop();

			}

		}
	}


}




CTrack SalesmanTrackBacktrackingGreedy(CGraph& graph, CVisits& visits)
{
	//we control that the graf has edges
	if (!graph.m_Edges.empty())
	{
		//MATRIX CREATION

		//we create ways matrix with n vextors = n-1 vistis because we dont want the dijkstra distance from final vertex to the others
		matriu.resize(visits.m_Vertices.size() - 1);
		int index = 0;

		for (auto vOriginDijks = visits.m_Vertices.begin(); vOriginDijks != next(visits.m_Vertices.begin(), visits.m_Vertices.size() - 1); vOriginDijks++)
		{
			DijkstraQueue(graph, *vOriginDijks);

			//we save in each vertex its index from the matrix
			(*vOriginDijks)->m_indexMatrix = index;

			for (CVertex* vDestinationDijks : visits.m_Vertices) {

				Cell matrixCell;
				//we have to save in each vell of the matrix the DIJKSTRA way and de distance of this way
				matrixCell.dijkstraDistance = vDestinationDijks->m_DijkstraDistance;

				//we dont have dijkstraPrevios in the vertex from we explore because is th origin of Dijkstra
				if (*vOriginDijks != vDestinationDijks)
				{
					while (vDestinationDijks != *vOriginDijks)
					{
						matrixCell.tram.push_front(vDestinationDijks->m_pDijkstraPrevious);
						vDestinationDijks = vDestinationDijks->m_pDijkstraPrevious->m_pOrigin;
					}
				}

				matriu[index].push_back(matrixCell);

			}


			index++;
		}

		//we put the index of the column that represents the destination vertex in the matrix
		visits.m_Vertices.back()->m_indexMatrix = index;

		//------------------------------------------------------------------

		//we get the destination and we select all the visits that can be or not visitied with the current way
		pDesti2 = visits.m_Vertices.back();

		//we put the vertex that we need to visit as global because we need them in case we get to a possible solution
		//VISITES = &visits;

		//at the beggining we can visitid all the vertex
		list<int> possibleVisits;
		for (CVertex* v : visits.m_Vertices)
		{
			possibleVisits.push_back(v->m_indexMatrix);
		}

		//we take out the first one and we will start there
		int currentNode = possibleVisits.front();
		possibleVisits.pop_front();

		//we initialize the shortest way 
		LongitudCamiMesCurt2 = numeric_limits<double>::max();



		//we need to save the order in which we are going to visit the nodes
		list<int> cami;
		cami.push_back(currentNode);

		//we look for the best way to visit all the nodes
		SalesmanTrackBacktrackingGreedyRec(possibleVisits, currentNode, cami);

		//create the cTrack solucion with the index cami solution
		//we create the CTrack with th first section already addes
		//auto itAux = next(indexCamiMesCurt.begin(), 1);
		CamiMesCurt2 = CTrack(&graph, matriu[indexCamiMesCurt.front()][*next(indexCamiMesCurt.begin(), 1)].tram);

		//merge al de indexes trams that we created with Dijkstra and we had saved in the matrix
		for (auto it = next(indexCamiMesCurt.begin(), 1); it != --indexCamiMesCurt.end(); it++) {
			auto itNext = it;
			itNext = next(itNext, 1);
			for (CEdge* e : matriu[*it][*itNext].tram)
				CamiMesCurt2.m_Edges.push_back(e);
			//solution.m_Edges.merge(matriu[topNode->m_indexes[element]][topNode->m_indexes[element] + 1].tram);
		}


	}

	return CamiMesCurt2;
}