#include "pch.h"
#include "Graph.h"
#include <queue>
#include <iostream>
#include <iomanip> 

struct CBBNode1 {
	double CotaInferior; //cota_inf
	bool passat[20]; //llista_observados
	unsigned indexos[20]; //llista_indices
	unsigned index;
	CVertex* v;
};

struct comparator1 {
	bool operator()(const CBBNode1* s1, const CBBNode1* s2) {
		return s1->v->m_DijkstraDistance < s2->v->m_DijkstraDistance; //ordena la cua per dijkstra enlloc de cota 
	}
};

CTrack SalesmanTrackBranchAndBound1(CGraph& graph, CVisits& visits)
{
	//Inicialització de les variables
	unsigned num_vertex2 = 0;
	double matriz[20][20];
	CVertex* llistaAux[20];
	double nova_cota_inf, cota_inf;
	unsigned columna, actual;
	double Longitud_Max = numeric_limits<double>::max();
	CBBNode1 inicial, * auxiliar, * cami_final = (CBBNode1*)malloc(sizeof(CBBNode1));
	CTrack Final = CTrack(&graph);
	CVertex* v_Inici, * v_Final;

	//Inicialitzem la struct CBBNode1 amb el node inicial
	inicial.CotaInferior = 0;
	inicial.indexos[0] = 0; //lo del vertex inicial ho podem posar directament
	inicial.passat[0] = true;
	inicial.index = 1;
	inicial.v = visits.m_Vertices.front();


	//Creació de la cua amb prioritat, ja podem afegir la 1a visita --> pel while
	priority_queue<CBBNode1*, std::vector<CBBNode1*>, comparator1> cua;
	cua.push(&inicial);


	//Fem una còpia de la llista visits.m_Vertices a llistaAux --> per poder accedir als vertex per l'index
	for (CVertex* v1 : visits.m_Vertices) {
		llistaAux[num_vertex2] = v1;
		num_vertex2++;
	}

	//Creació de la matriu
	for (unsigned fila = 0; fila < num_vertex2 - 1; fila++) {
		Dijkstra(graph, llistaAux[fila]);
		columna = 0;
		for (CVertex* v : visits.m_Vertices) {
			matriz[fila][columna] = v->m_DijkstraDistance;
			columna++;
		}
	}

	//Cas particular del simple 
	if (visits.m_Vertices.size() == 2) {
		v_Inici = visits.m_Vertices.front();
		v_Final = visits.m_Vertices.back();
		for (CEdge* edgeAux : v_Inici->m_Edges) {
			if (edgeAux->m_pDestination == v_Final) {
				CEdge* EdgeAdded = edgeAux;
				Final.AddLast(EdgeAdded);
				return Final;
			}
		}
	}

	columna = 0;
	while (columna < num_vertex2) {
		matriz[18][columna] = Longitud_Max;
		matriz[19][columna] = 0;
		for (unsigned fila = 0; fila < num_vertex2 - 1; fila++) {
			double valorMatriu = matriz[fila][columna];
			if (valorMatriu < matriz[18][columna] && fila != columna) {
				matriz[18][columna] = valorMatriu;
			}
			if (valorMatriu > matriz[19][columna]) {
				matriz[19][columna] = valorMatriu;
			}
		}
		columna++;
	}

	columna = 0;
	while (columna < num_vertex2) {
		inicial.CotaInferior += matriz[18][columna];  //Total cota inferior --> Cota inferior global
		columna++;
	}
	for (unsigned i = 1; i < num_vertex2 - 1; i++) inicial.passat[i] = false;

	while (!cua.empty()) {
		auxiliar = cua.top(); //cua té l'inicial, que té tots els valors guardats a dalt
		cua.pop();
		actual = auxiliar->indexos[auxiliar->index - 1]; //Del vertex volem saber la posicio (fila)
		cota_inf = auxiliar->CotaInferior;
		if (auxiliar->index == num_vertex2 - 1) { //Si l'auxiliar és la última visita
			if (cota_inf + matriz[actual][num_vertex2 - 1] - matriz[18][num_vertex2 - 1] < Longitud_Max) {
				*cami_final = *auxiliar; //Per construir el track --> anirem tirant cap a les arestes anteriors
				Longitud_Max = cota_inf + matriz[actual][num_vertex2 - 1] - matriz[18][num_vertex2 - 1];
			}
		}
		else {
			for (unsigned i = 1; i < num_vertex2 - 1; i++) {
				nova_cota_inf = auxiliar->CotaInferior + matriz[actual][i] - matriz[18][i]; //Nova cota inferior
				if (!(auxiliar->passat[i]) && nova_cota_inf < Longitud_Max) {

					CBBNode1* nou_inicial = (CBBNode1*)malloc(sizeof(CBBNode1));
					*nou_inicial = *auxiliar; //serà el nou vèrtex inicial
					//Inicialitzem tot:
					nou_inicial->indexos[nou_inicial->index] = i;
					nou_inicial->passat[i] = true;
					nou_inicial->CotaInferior = nova_cota_inf;
					nou_inicial->index++;
					nou_inicial->v = llistaAux[i];
					cua.push(nou_inicial);
				}
			}
		}
	}

	cami_final->indexos[num_vertex2 - 1] = num_vertex2 - 1; //L'últim sabem que és això

	//Passem a Track --> igual que a Backtracking!!!

	for (unsigned j = num_vertex2 - 1; j > 0; j--) {
		v_Inici = llistaAux[cami_final->indexos[j - 1]];
		v_Final = llistaAux[cami_final->indexos[j]];
		Dijkstra(graph, v_Inici);

		while (v_Final != v_Inici) {
			Final.AddFirst(v_Final->m_pDijkstraPrevious);
			v_Final = v_Final->m_pDijkstraPrevious->m_pOrigin;
		}
	}
	return Final;
}

struct CBBNode2 {
	double CotaSuperior; //cota_sup
	double CotaInferior; //cota_inf
	bool passat[20]; //llista_observados
	unsigned indexos[20]; //llista_indices
	unsigned index;
};

struct comparator2 {
	bool operator()(const CBBNode2* s1, const CBBNode2* s2) {
		return s1->CotaInferior < s2->CotaInferior;
	}
};

CTrack SalesmanTrackBranchAndBound2(CGraph& graph, CVisits& visits)
{
	//variables
	unsigned num_vertex2 = 0;
	double matriz[20][20];
	CVertex* llistaAux[20];
	unsigned columna;
	double Longitud_Max = numeric_limits<double>::max();
	CBBNode2 inicial, * auxiliar, * cami_final = (CBBNode2*)malloc(sizeof(CBBNode2));
	unsigned actual;
	double nova_cota_inf, cota_inf;
	CTrack Final = CTrack(&graph);
	CVertex* v_Inici, * v_Final;

	//donem valors
	inicial.CotaInferior = 0;
	inicial.CotaSuperior = 0;
	inicial.indexos[0] = 0; //lo del vertex inicial ho podem posar directament
	inicial.passat[0] = true;
	inicial.index = 1;
	double error = 0.00001;

	//creació de la cua + ja podem afegir la 1a visita --> pel while
	priority_queue<CBBNode2*, std::vector<CBBNode2*>, comparator2> cua;
	cua.push(&inicial);

	//fem una còpia de la llista visits.m_Vertices a llistaAux --> per poder accedir als vertex per l'index
	for (CVertex* v1 : visits.m_Vertices) {
		llistaAux[num_vertex2] = v1;
		num_vertex2++;
	}

	//creació de la matriu
	for (unsigned fila = 0; fila < num_vertex2 - 1; fila++) {
		Dijkstra(graph, llistaAux[fila]);
		columna = 0;
		for (CVertex* v : visits.m_Vertices) {
			matriz[fila][columna] = v->m_DijkstraDistance;
			columna++;
		}
	}

	//Cas particular del simple 
	if (visits.m_Vertices.size() == 2) {
		v_Inici = visits.m_Vertices.front();
		v_Final = visits.m_Vertices.back();
		for (CEdge* edgeAux : v_Inici->m_Edges) {
			if (edgeAux->m_pDestination == v_Final) {
				CEdge* EdgeAdded = edgeAux;
				Final.AddLast(EdgeAdded);
				return Final;
			}
		}
	}

	//a la posicio 18 guarda la cota inferior i a la 19 la superior --> sabem que als tests com a màxim 17 visites
	columna = 0;
	while (columna < num_vertex2) {
		matriz[18][columna] = Longitud_Max;
		matriz[19][columna] = 0;
		for (unsigned fila = 0; fila < num_vertex2 - 1; fila++) {
			double valorMatriu = matriz[fila][columna];
			if (valorMatriu < matriz[18][columna] && fila != columna) { // si la dijkstraDistance < cota inferior
				matriz[18][columna] = valorMatriu; // passa a ser la cota inferior
			}
			if (valorMatriu > matriz[19][columna]) { //si la dijstraDistance > cota superior
				matriz[19][columna] = valorMatriu; //passa a ser la cota superior
			}
		}
		columna++;
	}

	columna = 0; //--> per estalviar-nos declarar variables per fors
	while (columna < num_vertex2) {
		inicial.CotaSuperior += matriz[19][columna];  //Total cota superior --> Cota superior global
		inicial.CotaInferior += matriz[18][columna];  //Total cota inferior --> Cota inferior global
		columna++;
	}

	for (unsigned i = 1; i < num_vertex2 - 1; i++) inicial.passat[i] = false;

	while (!cua.empty()) {
		auxiliar = cua.top(); //cua té l'inicial, que té tots els valors guardats a dalt
		cua.pop();
		actual = auxiliar->indexos[auxiliar->index - 1]; //del vertex volem saber la posicio (fila)
		cota_inf = auxiliar->CotaInferior;
		if (auxiliar->index == num_vertex2 - 1) { //si l'auxiliar és la última visita
			if (cota_inf + matriz[actual][num_vertex2 - 1] - matriz[18][num_vertex2 - 1] < Longitud_Max + error) { //si és més petit que la cota Superior
				*cami_final = *auxiliar; //per construir el track --> anirem tirant cap a les arestes anteriors
				Longitud_Max = cota_inf + matriz[actual][num_vertex2 - 1] - matriz[18][num_vertex2 - 1]; //passa a ser la cota Superior
			}
		}
		else {
			for (unsigned i = 1; i < num_vertex2 - 1; i++) {
				nova_cota_inf = auxiliar->CotaInferior + matriz[actual][i] - matriz[18][i]; //nova cota inferior
				if (!(auxiliar->passat[i]) && nova_cota_inf < Longitud_Max + error) { //si la cota inferior < cota superior

					CBBNode2* nou_inicial = (CBBNode2*)malloc(sizeof(CBBNode2));
					*nou_inicial = *auxiliar; //serà el nou vèrtex inicial
					//inicialitzem tot:
					nou_inicial->indexos[nou_inicial->index] = i;
					nou_inicial->passat[i] = true;
					nou_inicial->CotaInferior = nova_cota_inf;
					nou_inicial->CotaSuperior += matriz[actual][i] - matriz[19][i];
					nou_inicial->index++;

					if (nou_inicial->CotaSuperior < Longitud_Max + error) {
						Longitud_Max = nou_inicial->CotaSuperior;
					}

					cua.push(nou_inicial);
				}
			}
		}
	}

	cami_final->indexos[num_vertex2 - 1] = num_vertex2 - 1; //l'últim sabem que és això

	//Passem a Track --> igual que a Backtracking!!!

	for (unsigned j = num_vertex2 - 1; j > 0; j--) {
		v_Inici = llistaAux[cami_final->indexos[j - 1]];
		v_Final = llistaAux[cami_final->indexos[j]];
		Dijkstra(graph, v_Inici);

		while (v_Final != v_Inici) {
			Final.AddFirst(v_Final->m_pDijkstraPrevious);
			v_Final = v_Final->m_pDijkstraPrevious->m_pOrigin;
		}
	}
	return Final;
}

// SalesmanTrackBranchAndBound3 ===================================================

CTrack SalesmanTrackBranchAndBound3(CGraph& graph, CVisits& visits)
{
	//variables
	unsigned num_vertex2 = 0;
	double matriz[20][20];
	CVertex* llistaAux[20];
	unsigned columna;
	double Longitud_Max = numeric_limits<double>::max();
	CBBNode2 inicial, * auxiliar, * cami_final = (CBBNode2*)malloc(sizeof(CBBNode2));
	unsigned actual;
	double nova_cota_inf, cota_inf;
	CTrack Final = CTrack(&graph);
	CVertex* v_Inici, * v_Final;

	//donem valors
	inicial.CotaInferior = 0;
	inicial.CotaSuperior = 0;
	inicial.indexos[0] = 0; //lo del vertex inicial ho podem posar directament
	inicial.passat[0] = true;
	inicial.index = 1;
	double error = 0.00001;

	//creació de la cua + ja podem afegir la 1a visita --> pel while
	priority_queue<CBBNode2*, std::vector<CBBNode2*>, comparator2> cua;
	cua.push(&inicial);

	//Cas particular del simple 
	if (visits.m_Vertices.size() == 2) {
		v_Inici = visits.m_Vertices.front();
		v_Final = visits.m_Vertices.back();
		for (CEdge* edgeAux : v_Inici->m_Edges) {
			if (edgeAux->m_pDestination == v_Final) {
				CEdge* EdgeAdded = edgeAux;
				Final.AddLast(EdgeAdded);
				return Final;
			}
		}
	}


	//fem una còpia de la llista visits.m_Vertices a llistaAux --> per poder accedir als vertex per l'index
	for (CVertex* v1 : visits.m_Vertices) {
		llistaAux[num_vertex2] = v1;
		num_vertex2++;
	}

	//creació de la matriu
	for (unsigned fila = 0; fila < num_vertex2 - 1; fila++) {
		DijkstraQueue(graph, llistaAux[fila]);
		columna = 0;
		for (CVertex* v : visits.m_Vertices) {
			matriz[fila][columna] = v->m_DijkstraDistance;
			columna++;
		}
	}

	//a la posicio 18 guarda la cota inferior i a la 19 la superior --> sabem que als tests com a màxim 17 visites
	columna = 0;
	while (columna < num_vertex2) {
		matriz[18][columna] = Longitud_Max;
		matriz[19][columna] = 0;
		for (unsigned fila = 0; fila < num_vertex2 - 1; fila++) {
			double valorMatriu = matriz[fila][columna];
			if (valorMatriu < matriz[18][columna] && fila != columna) { // si la dijkstraDistance < cota inferior
				matriz[18][columna] = valorMatriu; // passa a ser la cota inferior
			}
			if (valorMatriu > matriz[19][columna]) { //si la dijstraDistance > cota superior
				matriz[19][columna] = valorMatriu; //passa a ser la cota superior
			}
		}
		columna++;
	}

	columna = 0; //--> per estalviar-nos declarar variables per fors
	while (columna < num_vertex2) {
		inicial.CotaSuperior += matriz[19][columna];  //Total cota superior --> Cota superior global
		inicial.CotaInferior += matriz[18][columna];  //Total cota inferior --> Cota inferior global
		columna++;
	}

	list<unsigned> indexJaPassat;
	indexJaPassat.push_back(inicial.index - 1);

	for (unsigned i = 1; i < num_vertex2 - 1; i++) inicial.passat[i] = false;

	while (!cua.empty()) {
		auxiliar = cua.top(); //cua té l'inicial, que té tots els valors guardats a dalt
		cua.pop();
		actual = auxiliar->indexos[auxiliar->index - 1]; //del vertex volem saber la posicio (fila)
		cota_inf = auxiliar->CotaInferior;
		if (auxiliar->index == num_vertex2 - 1) { //si l'auxiliar és la última visita
			if (cota_inf + matriz[actual][num_vertex2 - 1] - matriz[18][num_vertex2 - 1] < Longitud_Max + error) { //si és més petit que la cota Superior
				*cami_final = *auxiliar; //per construir el track --> anirem tirant cap a les arestes anteriors
				Longitud_Max = cota_inf + matriz[actual][num_vertex2 - 1] - matriz[18][num_vertex2 - 1]; //passa a ser la cota Superior
			}
		}
		else {
			bool trobat = false;
			for (unsigned indx : indexJaPassat) {
				if (auxiliar->index - 1 == indx - 1) trobat = true; break;
			}
			if (trobat != true) {
				for (unsigned i = 1; i < num_vertex2 - 1; i++) {
					nova_cota_inf = auxiliar->CotaInferior + matriz[actual][i] - matriz[18][i]; //nova cota inferior
					if (!(auxiliar->passat[i]) && nova_cota_inf < Longitud_Max + error) { //si la cota inferior < cota superior
						indexJaPassat.push_back(auxiliar->index - 1);
						CBBNode2* nou_inicial = (CBBNode2*)malloc(sizeof(CBBNode2));
						*nou_inicial = *auxiliar; //serà el nou vèrtex inicial
						//inicialitzem tot:
						nou_inicial->indexos[nou_inicial->index] = i;
						nou_inicial->passat[i] = true;
						nou_inicial->CotaInferior = nova_cota_inf;
						nou_inicial->CotaSuperior += matriz[actual][i] - matriz[19][i];
						nou_inicial->index++;

						if (nou_inicial->CotaSuperior < Longitud_Max + error) {
							Longitud_Max = nou_inicial->CotaSuperior;
						}

						cua.push(nou_inicial);
					}
				}
			}
		}
	}

	cami_final->indexos[num_vertex2 - 1] = num_vertex2 - 1; //l'últim sabem que és això

	//Passem a Track --> igual que a Backtracking!!!

	for (unsigned j = num_vertex2 - 1; j > 0; j--) {
		v_Inici = llistaAux[cami_final->indexos[j - 1]];
		v_Final = llistaAux[cami_final->indexos[j]];
		Dijkstra(graph, v_Inici);

		while (v_Final != v_Inici) {
			Final.AddFirst(v_Final->m_pDijkstraPrevious);
			v_Final = v_Final->m_pDijkstraPrevious->m_pOrigin;
		}
	}
	return Final;
}

// SalesmanTrackBranchAndBound4 ===================================================


CTrack SalesmanTrackBranchAndBound4(CGraph& graph, CVisits& visits)
{
	return CTrack(&graph);
}
