
// GraphApplication.cpp: define los comportamientos de las clases para la aplicación.
//

#include "pch.h"
#include "framework.h"
#include "GraphApplication.h"
#include "GraphApplicationDlg.h"
#include <io.h>
#include <fcntl.h>
#include <iomanip>
#include "Logs.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// =============================================================================
// IDENTIFICACION DE LOS ALUMNOS ===============================================
// =============================================================================

CString NombreAlumno1 = "Moisés";
CString ApellidosAlumno1 = "Sánchez Pin";
CString NIUAlumno1 = "1603611"; // NIU alumno1

bool CheckNIU(const char* niu)
{
	if (strlen(niu) != 7) return false;
	if (strcmp(niu, "0000000") == 0) return false;
	return true;
}

// ReadCommandArguments ========================================================

void ReadCommandArguments(int &nArgs, CString args[])
{
	nArgs = 0;
	char* pCommand = GetCommandLineA();
	//cout << pCommand << endl;
	char *pArg = pCommand;
	// Saltar el nombre del ejecutable
	while (*pArg <= ' ') if (*pArg) ++pArg; else return;
	if (*pArg == '\"') {
		// Leer hasta "
		++pArg;
		while (*pArg != '\"') {
			if (*pArg) ++pArg;
			else return;
		}
		++pArg;
	}
	else {
		while (*pArg > ' ') {
			if (*pArg) ++pArg;
			else return;
		}
	}
	// Saltar espacios
	while (*pArg == ' ') ++pArg;
	if (*pArg == '\0') return;
	//cout << "Linea de comandos " << (int) *pArg << endl;
	//AttachConsole(ATTACH_PARENT_PROCESS);
	//ofstream *pout = new ofstream("CONOUT$", ofstream::out | ofstream::app);
	//cout.rdbuf(pout->rdbuf());

	while (*pArg != '\0') {
		// Saltar espacios
		while (*pArg == ' ') ++pArg;
		if (*pArg == '\0') break;
		else if (*pArg == '\"') {
			++pArg;
			for (;;) {
				if (*pArg == '\0') break;
				else if (*pArg == '\"') break;
				args[nArgs].AppendChar(*pArg);
				++pArg;
			}
		}
		else {
			while (*pArg > ' ') {
				args[nArgs].AppendChar(*pArg);
				++pArg;
			}
		}
		++nArgs;
		if (nArgs >= 10) {
			return;
		}
	}
}

// EjecutarLineaDeComandos =====================================================

bool EjecutarLineaDeComandos(int nArgs, CString args[])
{
	cout << "NARGS " << nArgs << endl;
	if (nArgs > 1 && args[0].CompareNoCase("Display") == 0) {
		return false;
	}
	if (nArgs == 1 && args[0].CompareNoCase("Authors") == 0) {
		ofstream file("Authors.txt");
		file << NIUAlumno1 << endl;
		file.close();
		return true;
	}
	if (nArgs < 3) {
		cout << "Número de argumentos erroneo" << endl;
		return true;
	}
	cout << "PARAMETROS:" << endl;

	if (args[1].CompareNoCase("Prim") == 0) {
		if (nArgs != 4) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo....: " << args[1] << endl;
		cout << "Grafo........: " << args[2] << endl;
		cout << "SpanningTree.: " << args[3] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			double t0 = Clock();
			CSpanningTree tree=SpanningTreePrim(graph);
			double t1 = Clock();
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << t1 - t0 << endl;
			tree.Save(args[3]);
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("KrusKal") == 0) {
		if (nArgs != 4) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo....: " << args[1] << endl;
		cout << "Grafo........: " << args[2] << endl;
		cout << "SpanningTree.: " << args[3] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			double t0 = Clock();
			CSpanningTree tree = SpanningTreeKruskal(graph);
			double t1 = Clock();
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << t1 - t0 << endl;
			tree.Save(args[3]);
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("dijkstra") == 0) {
		// dijkstra grafo distancias
		if (nArgs != 4) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Distancias.: " << args[3] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVertex* pStart = graph.FindVertex("Start");
			if (!pStart) {
				if (graph.m_Vertices.empty()) throw MyException("No es poden calcular distancies en un graf sense vertex");
				pStart = &graph.m_Vertices.front();
			}
			double t0 = Clock();
			Dijkstra(graph, pStart);
			double t1 = Clock();
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << t1 - t0 << endl;
			graph.SaveDistances(args[3]);
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("dijkstraqueue") == 0) {
		// dijkstraqueue grafo distancias
		if (nArgs != 4) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Distancias.: " << args[3] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVertex* pStart = graph.FindVertex("Start");
			if (!pStart) {
				if (graph.m_Vertices.empty()) throw MyException("No es poden calcular distancies en un graf sense vertex");
				pStart = &graph.m_Vertices.front();
			}
			double t0 = Clock();
			DijkstraQueue(graph, pStart);
			double t1 = Clock();
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << t1 - t0 << endl;
			graph.SaveDistances(args[3]);
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("greedy") == 0) {
		if (nArgs != 5) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Visitas....: " << args[3] << endl;
		cout << "Track......: " << args[4] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVisits visits(&graph);
			visits.Load(args[3]);
			double t0 = Clock();
			CTrack track = SalesmanTrackGreedy(graph, visits);
			double t1 = Clock();
			track.Save(args[4]);
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << t1 - t0 << endl;
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("backtracking") == 0) {
		if (nArgs != 5) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Visitas....: " << args[3] << endl;
		cout << "Track......: " << args[4] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVisits visits(&graph);
			visits.Load(args[3]);
			double t0 = Clock();
			CTrack track = SalesmanTrackBacktracking(graph, visits);
			double t1 = Clock();
			track.Save(args[4]);
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << t1 - t0 << endl;
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("backtrackinggreedy") == 0) {
		if (nArgs != 5) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Visitas....: " << args[3] << endl;
		cout << "Track......: " << args[4] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVisits visits(&graph);
			visits.Load(args[3]);
			double t0 = Clock();
			CTrack track = SalesmanTrackBacktrackingGreedy(graph, visits);
			double t1 = Clock();
			track.Save(args[4]);
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << t1 - t0 << endl;
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("branchandbound1") == 0) {
		if (nArgs != 5) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Visitas....: " << args[3] << endl;
		cout << "Track......: " << args[4] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVisits visits(&graph);
			visits.Load(args[3]);
			double t0 = Clock();
			CTrack track = SalesmanTrackBranchAndBound1(graph, visits);
			double t1 = Clock();
			track.Save(args[4]);
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << t1 - t0 << endl;
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("branchandbound2") == 0) {
		if (nArgs != 5) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Visitas....: " << args[3] << endl;
		cout << "Track......: " << args[4] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVisits visits(&graph);
			visits.Load(args[3]);
			double t0 = Clock();
			CTrack track = SalesmanTrackBranchAndBound2(graph, visits);
			double t1 = Clock();
			track.Save(args[4]);
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << t1 - t0 << endl;
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("branchandbound3") == 0) {
		if (nArgs != 5) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Visitas....: " << args[3] << endl;
		cout << "Track......: " << args[4] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVisits visits(&graph);
			visits.Load(args[3]);
			double t0 = Clock();
			CTrack track = SalesmanTrackBranchAndBound3(graph, visits);
			double t1 = Clock();
			track.Save(args[4]);
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << setprecision(9) << t1 - t0 << endl;
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("branchandbound4") == 0) {
		if (nArgs != 5) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Visitas....: " << args[3] << endl;
		cout << "Track......: " << args[4] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVisits visits(&graph);
			visits.Load(args[3]);
			double t0 = Clock();
			CTrack track = SalesmanTrackBranchAndBound4(graph, visits);
			double t1 = Clock();
			track.Save(args[4]);
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << setprecision(9) << t1 - t0 << endl;
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
		}
	else if (args[1].CompareNoCase("DynamicProgramming1") == 0) {
		if (nArgs != 5) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Visitas....: " << args[3] << endl;
		cout << "Track......: " << args[4] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVisits visits(&graph);
			visits.Load(args[3]);
			double t0 = Clock();
			CTrack track = SalesmanTrackDynamicProgramming1(graph, visits);
			double t1 = Clock();
			track.Save(args[4]);
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << setprecision(9) << t1 - t0 << endl;
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("DynamicProgramming2") == 0) {
		if (nArgs != 5) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Visitas....: " << args[3] << endl;
		cout << "Track......: " << args[4] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVisits visits(&graph);
			visits.Load(args[3]);
			double t0 = Clock();
			CTrack track = SalesmanTrackDynamicProgramming2(graph, visits);
			double t1 = Clock();
			track.Save(args[4]);
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << setprecision(9) << t1 - t0 << endl;
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
		}
	else if (args[1].CompareNoCase("probabilistic") == 0) {
		if (nArgs != 5) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo..: " << args[1] << endl;
		cout << "Grafo......: " << args[2] << endl;
		cout << "Visitas....: " << args[3] << endl;
		cout << "Track......: " << args[4] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			graph.SetDistancesToEdgeLength();
			CVisits visits(&graph);
			visits.Load(args[3]);
			double t0 = Clock();
			CTrack track = SalesmanTrackProbabilistic(graph, visits);
			double t1 = Clock();
			track.Save(args[4]);
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << setprecision(9) << t1 - t0 << endl;
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else if (args[1].CompareNoCase("QuickHull") == 0) {
		if (nArgs != 4) {
			cout << "Número de argumentos erroneo" << endl;
			return true;
		}
		cout << "Algoritmo....: " << args[1] << endl;
		cout << "Grafo........: " << args[2] << endl;
		cout << "Convex Hull..: " << args[3] << endl;
		try {
			CGraph graph(false);
			graph.Load(args[2]);
			double t0 = Clock();
			CConvexHull convexHull = QuickHull(graph);
			double t1 = Clock();
			cout << endl << "TIEMPO DE EJECUCION: " << fixed << t1 - t0 << endl;
			convexHull.Save(args[3]);
		}
		catch (exception& ex) {
			cout << "EXCEPTION: " << ex.what() << endl;
		}
	}
	else {
		cout << "Error en el nombre del algoritmo" << endl;
	}
	cout << "FIN EJECUCION COMANDOS" << endl;
	return true;
}

// CGraphApplicationApp

BEGIN_MESSAGE_MAP(CGraphApplicationApp, CWinApp)
	ON_COMMAND(ID_HELP, &CWinApp::OnHelp)
END_MESSAGE_MAP()


// Construcción de CGraphApplicationApp

CGraphApplicationApp::CGraphApplicationApp()
{
	// TODO: agregar aquí el código de construcción,
	// Colocar toda la inicialización importante en InitInstance
}


// Único objeto CGraphApplicationApp

CGraphApplicationApp theApp;


// Inicialización de CGraphApplicationApp

BOOL CGraphApplicationApp::InitInstance()
{
	ofstream sortida;
	int nArgs;
	CString args[10];
	ReadCommandArguments(nArgs, args);
	CString outputName = "Output.txt";
	if (nArgs > 1) outputName = args[0];
	bool lineaDeComandos = false;
	if (AttachConsole(ATTACH_PARENT_PROCESS)) {
		ofstream *pout = new ofstream(outputName);
		cout.rdbuf(pout->rdbuf());
		lineaDeComandos = true;
	}
	else {
		AllocConsole();
		sortida.open("CONOUT$");
		cout.rdbuf(sortida.rdbuf());
		//ofstream *pout = new ofstream("CONOUT$");
		//cout.rdbuf(pout->rdbuf());
	}
	bool error = false;
	cout << "Nom Alumne 1.....: " << NombreAlumno1 << endl;
	cout << "Cognoms Alumne 1.: " << ApellidosAlumno1 << endl;
	cout << "NIU Alumne 1.....: " << NIUAlumno1 << endl;
	if (!CheckNIU(NIUAlumno1)) {
		cout << " NIU ERRONEO!!!";
		error = true;
	}
	cout << endl;
	if (error) {
		if (!lineaDeComandos) AfxMessageBox("ERROR EN LA IDENTIFICACION DE LOS ALUMNOS", MB_OK | MB_ICONERROR);
		return FALSE;
	}
	/*
	LARGE_INTEGER cnt;
	LARGE_INTEGER fre;

	if (QueryPerformanceFrequency(&fre) && QueryPerformanceCounter(&cnt)) {
		//cout << 1.0 / (double)fre.QuadPart << endl;
		cout << "PRECISION CLOCK " << 1.0 / (double)fre.QuadPart << endl;
	}
	else {
		//return 0.0;
		cout << "Error clock" << endl;
	}
	*/

	// Funcionamiento de la aplicación sólo en modo consola si tiene parámetros
	if (nArgs > 0 && EjecutarLineaDeComandos(nArgs, args)) return FALSE;
	OpenLog();
	// Windows XP requiere InitCommonControlsEx() si un manifiesto de
	// aplicación especifica el uso de ComCtl32.dll versión 6 o posterior para habilitar
	// estilos visuales.  De lo contrario, se generará un error al crear ventanas.
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// Establecer para incluir todas las clases de control comunes que desee utilizar
	// en la aplicación.
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();


	AfxEnableControlContainer();

	// Crear el administrador de shell, en caso de que el cuadro de diálogo contenga
	// alguna vista de árbol de shell o controles de vista de lista de shell.
	CShellManager *pShellManager = new CShellManager;

	// Activar el administrador visual "nativo de Windows" para habilitar temas en controles MFC
	CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerWindows));

	// Inicialización estándar
	// Si no utiliza estas funcionalidades y desea reducir el tamaño
	// del archivo ejecutable final, debe quitar
	// las rutinas de inicialización específicas que no necesite
	// Cambie la clave del Registro en la que se almacena la configuración
	// TODO: debe modificar esta cadena para que contenga información correcta
	// como el nombre de su compañía u organización
	SetRegistryKey(_T("Aplicaciones generadas con el Asistente para aplicaciones local"));

	CGraphApplicationDlg dlg;
	if (nArgs > 1 && args[0].CompareNoCase("Display") == 0) {
		CString tmp("Display");
		for (int i = 1; i < nArgs; ++i) {
			tmp += " ";
			tmp += args[1];
		}
		PrintLog("%s", tmp);
		for (int i = 1; i < nArgs; ++i) {
			//AfxMessageBox(args[i]);
			char ext[10];
			_splitpath_s(args[i], NULL, 0, NULL, 0, NULL, 0, ext, 10);
			//AfxMessageBox(ext);
			if (_stricmp(ext, ".GR") == 0) dlg.m_Graph.Load(args[i]);
			else if (_stricmp(ext, ".DIS") == 0) {
				dlg.m_Graph.LoadDistances(args[i]);
				dlg.SetViewDistances(true);
			}
			else if (_stricmp(ext, ".VIS") == 0) dlg.m_Visits.Load(args[i]);
			else if (_stricmp(ext, ".TRK") == 0) dlg.m_Track.Load(args[i]);
			else if (_stricmp(ext, ".CH") == 0) dlg.m_ConvexHull.Load(args[i]);
			else cout << "No se puede leer " << args[i] << " por extension desconocida" << endl;
		}
	}
	m_pMainWnd = &dlg;
	INT_PTR nResponse = dlg.DoModal();
	if (nResponse == IDOK)
	{
		// TODO: insertar aquí el código para controlar
		//  cuándo se descarta el cuadro de diálogo con Aceptar
	}
	else if (nResponse == IDCANCEL)
	{
		// TODO: insertar aquí el código para controlar
		//  cuándo se descarta el cuadro de diálogo con Cancelar
	}
	else if (nResponse == -1)
	{
		TRACE(traceAppMsg, 0, "Advertencia: la aplicación está finalizando porque hubo un error al crear el cuadro de diálogo.\n");
		TRACE(traceAppMsg, 0, "Advertencia: si usa controles MFC en el cuadro de diálogo, no puede usar #define _AFX_NO_MFC_CONTROLS_IN_DIALOGS.\n");
	}

	// Eliminar el administrador de shell creado anteriormente.
	if (pShellManager != nullptr)
	{
		delete pShellManager;
	}

#if !defined(_AFXDLL) && !defined(_AFX_NO_MFC_CONTROLS_IN_DIALOGS)
	ControlBarCleanUp();
#endif

	// Dado que el cuadro de diálogo se ha cerrado, devolver FALSE para salir
	//  de la aplicación en vez de iniciar el suministro de mensajes de dicha aplicación.
	CloseLog();
	return FALSE;
}

