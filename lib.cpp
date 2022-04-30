#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <fstream>
#include <time.h>
#include <cstring>
#include <chrono>
#include <queue>
#include <cmath>
#include <vector>
#include <algorithm>
#include "libgraph.h"
#include "lib.h"
#include "math.h"
#include "sort_ways.hpp"
using namespace std;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
#pragma warning(disable : 4996)
Graph *graph = new Graph();
vector<Pthread_args> *bspvec = NULL;
Algorithm *algorithm = new Algorithm();
Network *network = new Network();
double budget = 8.70E+03;
double overhead = 0.0;
int Seed = 0;
double alphaT = 0;
double L = 0;
double betaC = 0;
double ginger_cost = 0;
int bsplock = 0;
long mvpercents = 0;
double t1 = 0;
double t3 = 0;
double t4 = 0;
int bsp = 0;
int F1 = 0;
int F2 = 0;
double global_mvcost = 0;
double current_mvcost = 0;
double moverand = 0;
int batchsize = 0;
float raw_data_size = 10; // MB;
double *ug = NULL;
double *dg = NULL;
double *ua = NULL;
double *da = NULL;
int optitime = 1;
char *FN = new char[100];
int MAX_ITER = 10;
int MAX_THREADS_NUM = 50;
/* Newly add start */
int test_character_num = 0;		  // 1 for in-going degree, 2 for mirror, 3 for DC_speed, 4 for DC_price, 0 for default, means no action
int pthread_counter = 0;		  // a counter for multi-thread
double ratio_of_nodes = -1;		  // By default, ratio is meaningless
bool deleted_zero_or_not = false; // By default should not delete zero
/* Newly add end */
double timeAll = 0;
double priceAll = 0;
double timeOld = 0;
double priceOld = 0;
double dataunit = 0.000008;
int theta = 0;
double t5 = 0;
double t2 = 0;
double t6 = 0;
ofstream outfile;
map<int, int> trained = map<int, int>();
map<int, vector<int>> PSS = map<int, vector<int>>();
double stepbudget = 0;
double addbudget = 0;
double randpro = 1;
int *optimize = NULL;
char *outputfile = new char[200];
pthread_barrier_t barrier;
pthread_mutex_t mutex;
pthread_mutex_t mutex1;
struct A
{
	int key;
	int value;
};
struct cmp1
{
	bool operator()(A &a, A &b)
	{
		return a.value > b.value; //最小值优先
	}
};
// NOTE: eachThreadNum means the node can be modified in each thread
void multiprocessing_pool(int *eachThreadsNum)
{
	Vertex *v = graph->getVertexs();
	int k = graph->getNum_Nodes() % MAX_THREADS_NUM;
	int num = floor((double)(graph->getNum_Nodes() / MAX_THREADS_NUM));
	for (int i = 0; i < MAX_THREADS_NUM; i++)
		if (i < k)
			eachThreadsNum[i] = num + 1;
		else
			eachThreadsNum[i] = num;

	int sum = 0;
	for (int i = 0; i < MAX_THREADS_NUM; i++)
	{
		int low = 0, high = 0;
		for (int j = sum; j < sum + eachThreadsNum[i]; j++)
		{
			if (v[j].getType() == 0)
				high++;
			else
				low++;
		}
		// cout<<"Thread "<<i<<" : low, "<<low<<" ,high "<<high<<endl;
		sum += eachThreadsNum[i];
	}
	/*for (int i = 0; i < MAX_THREADS_NUM; i++)
		cout << "Chunks " << i << " size : " << eachThreadsNum[i] << endl;*/
	vector<int> *vec = new vector<int>[MAX_THREADS_NUM];
	int *veclow = new int[MAX_THREADS_NUM];
	int *vechigh = new int[MAX_THREADS_NUM];
	for (int i = 0; i < MAX_THREADS_NUM; i++)
	{
		veclow[i] = 0;
		vechigh[i] = 0;
	}
	time_t start = time(NULL);
	for (int i = 0; i < graph->getNum_Nodes(); i++)
	{
		if (v[i].getType() == 0)
		{
			int min_index = min_value_index(vechigh, MAX_THREADS_NUM);
			vechigh[min_index]++;
			vec[min_index].push_back(i);
		}
		else
		{
			int min_index = min_value_index(veclow, MAX_THREADS_NUM);
			veclow[min_index]++;
			vec[min_index].push_back(i);
		}
	}
	double end = time(NULL) - start;
	cout << "overhead on addressing straggle is: " << end * MAX_ITER << endl;
	for (int i = 0; i < MAX_THREADS_NUM; i++)
	{
		// MARK: PSS is an important variable, but for now have no idea about it
		PSS.insert(pair<int, vector<int>>(i, vec[i]));
	}

	for (int i = 0; i < MAX_THREADS_NUM; i++)
	{
		// cout<<"Thread "<<i<<" : low, "<<veclow[i]<<" ,high "<<vechigh[i]<<endl;
	}

	delete[] vec;
	delete[] veclow;
	delete[] vechigh;
}
char *input_format_specifier(char *filename, char *delimiter)
{
	string line;
	char *data;
	int count = 0;
	ifstream in(filename);
	getline(in, line);
	data = (char *)line.data();
	// cout << line << endl;
	char *token = strtok(data, delimiter);
	while (token)
	{
		token = strtok(NULL, delimiter);
		count++;
	}
	if (count == 1)
	{
		// cout << "The file format is incorrect." << endl;
		return NULL;
	}
	char *format = new char[15];
	switch (count)
	{
	case 2:
		strcpy(format, "%u %u");
		format[5] = '\0';
		break;
	case 3:
		strcpy(format, "%u %u %f");
		format[8] = '\0';
		break;
	}
	// cout << count << endl;
	// cout << format << endl;
	return format;
}
void initialNode(FILE *file_descriptor, char *format) //��������ļ�ֻ������
{
	int offset = ftell(file_descriptor);
	int maxVertexID = 0;
	unsigned int source = 0;
	unsigned int max_source = 0;
	unsigned int destination = 0;
	unsigned int max_destination = 0;
	while (fscanf(file_descriptor, format, &source, &destination) == 2)
	{
		if (source == destination)
		{ // printf("source: %d,destination: %d\n",source,destination);
		  // getchar();
		}
		graph->addEdge();
		if (source > max_source)
			max_source = source;
		if (destination > max_destination)
			max_destination = destination;
	}
	if (max_source > max_destination)
		maxVertexID = max_source;
	else
		maxVertexID = max_destination;
	cout << "max: " << maxVertexID << endl;
	// getchar();

	fseek(file_descriptor, offset, SEEK_SET);

	int *active = new int[maxVertexID + 1];
	for (int i = 0; i <= maxVertexID; i++)
		active[i] = 0;

	while (fscanf(file_descriptor, format, &source, &destination) == 2)
	{
		if (active[source] == 0)
		{
			graph->addNode();
			active[source] = 1;
		}
		if (active[destination] == 0)
		{
			graph->addNode();
			active[destination] = 1;
		}
	}

	fseek(file_descriptor, offset, SEEK_SET);
	int num_nodes = graph->getNum_Nodes();
	Vertex *v = new Vertex[num_nodes];
	int count = 0;
	long long *mapped = new long long[maxVertexID + 1];
	while (fscanf(file_descriptor, format, &source, &destination) == 2)
	{
		if (active[source] == 1)
		{
			v[count].setSignal(0);
			v[count++].setVertexID(source);
			mapped[source] = count - 1;
			active[source] = 0;
		}
		if (active[destination] == 1)
		{
			v[count].setSignal(0);
			v[count++].setVertexID(destination);
			mapped[destination] = count - 1;
			active[destination] = 0;
		}
	}
	graph->setMapped(mapped);
	fseek(file_descriptor, offset, SEEK_SET);
	graph->setNodes(v);

	cout << "#num_nodes: " << graph->getNum_Nodes() << endl;
	cout << "#num_edges: " << graph->getNum_Edges() << endl;
	delete[] active;
}
void initialNode(FILE *file_descriptor, char *format, double edgeValue)
{
	int offset = ftell(file_descriptor);
	int maxVertexID = 0;
	unsigned int source = 0;
	unsigned int max_source = 0;
	unsigned int destination = 0;
	unsigned int max_destination = 0;
	double edge_value = 0;
	while (fscanf(file_descriptor, format, &source, &destination, &edgeValue) == 3)
	{
		graph->addEdge();
		if (source > max_source)
			max_source = source;
		if (destination > max_destination)
			max_destination = destination;
	}
	if (max_source > max_destination)
		maxVertexID = max_source;
	else
		maxVertexID = max_destination;

	fseek(file_descriptor, offset, SEEK_SET);

	int *active = new int[maxVertexID + 1];
	for (int i = 0; i <= maxVertexID; i++)
		active[i] = 0;

	while (fscanf(file_descriptor, format, &source, &destination, &edgeValue) == 3)
	{
		if (active[source] == 0)
		{
			graph->addNode();
			active[source] = 1;
		}
		if (active[destination] == 0)
		{
			graph->addNode();
			active[destination] = 1;
		}
	}

	fseek(file_descriptor, offset, SEEK_SET);
	int num_nodes = graph->getNum_Nodes();
	Vertex *v = new Vertex[num_nodes];
	int count = 0;
	long long *mapped = new long long[maxVertexID + 1];
	while (fscanf(file_descriptor, format, &source, &destination, &edgeValue) == 3)
	{
		if (active[source] == 1)
		{
			// v[count].initSignal_();
			v[count].setSignal(0);
			v[count++].setVertexID(source);
			mapped[source] = count - 1;
			active[source] = 0;
		}
		if (active[destination] == 1)
		{
			// v[count].initSignal_();
			v[count].setSignal(0);
			v[count++].setVertexID(destination);
			mapped[destination] = count - 1;
			active[destination] = 0;
		}
	}
	graph->setMapped(mapped);
	fseek(file_descriptor, offset, SEEK_SET);
	graph->setNodes(v);

	cout << "#num_nodes: " << graph->getNum_Nodes() << endl;
	cout << "#num_edges: " << graph->getNum_Edges() << endl;
	delete[] active;
}
// 0 means not directed
// 1 means directed
void initialEdge(FILE *file_descriptor, char *format, bool directed)
{
	int offset = ftell(file_descriptor);
	unsigned int source = 0;
	unsigned int destination = 0;
	long long *mapped = graph->getMapped();
	Vertex *v = graph->getVertexs();
	double edgeValue = 0;
	srand((unsigned)time(NULL));
	int a = 1;
	int b = 100;
	double random = 0;
	while (fscanf(file_descriptor, format, &source, &destination) == 2)
	{
		if (source == destination)
		{
			// cout<<"两个相同的点构成一条边"<<endl;
			//  getchar();
		}
		random = (rand() % (b - a + 1)) + a;
		Edge *e = new Edge(source, destination, random);
		v[mapped[source]].addOutgoingEdge(e);
		v[mapped[destination]].addIngoingEdge(e);

		if (!directed)
		{
			e = new Edge(destination, source, random);
			v[mapped[source]].addIngoingEdge(e);
			v[mapped[destination]].addOutgoingEdge(e);
		}
	}
	fseek(file_descriptor, offset, SEEK_SET);
}
void initialEdge(FILE *file_descriptor, char *format, double edgeWeight, bool directed)
{
	int offset = ftell(file_descriptor);
	unsigned int source = 0;
	unsigned int destination = 0;
	long long *mapped = graph->getMapped();
	Vertex *v = graph->getVertexs();
	double edgeValue = 0;
	int a = 1;
	int b = 100;
	while (fscanf(file_descriptor, format, &source, &destination, &edgeValue) == 3)
	{
		Edge *e = new Edge(source, destination, edgeValue);
		v[mapped[source]].addOutgoingEdge(e);
		v[mapped[destination]].addIngoingEdge(e);

		if (!directed)
		{
			e = new Edge(destination, source, edgeValue);
			v[mapped[source]].addIngoingEdge(e);
			v[mapped[destination]].addOutgoingEdge(e);
		}
	}
	fseek(file_descriptor, offset, SEEK_SET);
}
int read_input_network(char *filename, double *upload, double *download, double *upprice)
{
	FILE *file_descriptor = fopen(filename, "r");
	if (!file_descriptor)
	{
		cout << "Error : can't open the input file!" << endl;
		return -1;
	}
	float uload = 0;
	float dload = 0;
	float uprice = 0;
	char format[] = "%f %f %f";
	int count = 0;
	while (fscanf(file_descriptor, format, &uload, &dload, &uprice) == 3)
	{
		upload[count] = uload;
		download[count] = dload;
		upprice[count++] = uprice;
		// cout << upload[count - 1] << " " << download[count - 1] << " " << upprice[count - 1] << endl;
	}
	return 0;
}

// NOTE: Initialize the nodes and edges
int read_input_file(char *filename, bool directed)
{
	FILE *file_descriptor = fopen(filename, "r");
	double edgeValue = 0;
	if (!file_descriptor)
	{
		cout << "Error : can't open the input file!" << endl;
		return -1;
	}

	char delimiter[] = " ,'\t'";
	char *format = input_format_specifier(filename, delimiter);
	// cout << fscanf(file_descriptor, "%u %u", &source, &destination) << endl;
	// cout << format << endl;
	if (!format)
	{
		cout << "input file format is error��" << endl;
		return -1;
	}
	if (strcmp(format, "%u %u") == 0)
	{
		initialNode(file_descriptor, format);
	}
	else if (strcmp(format, "%u %u %f") == 0)
	{
		initialNode(file_descriptor, format, edgeValue);
	}

	// NOTE:
	int offset = ftell(file_descriptor);
	if (strcmp(format, "%u %u") == 0)
	{
		initialEdge(file_descriptor, format, directed);
	}
	else if (strcmp(format, "%u %u %f") == 0)
	{
		initialEdge(file_descriptor, format, edgeValue, directed);
	}
	return 0;
}

// NOTE: Initialize the algorithm and the DC which every nodes belong to
void initAlgorithm(char *algorithmName, int numDCs)
{
	algorithm->InitAlgorithm(algorithmName, numDCs); //����ָ��ռ�
	algorithm->InitLabel();							 //�������
	algorithm->InitProbability();					 //��ʼ������
}

void initPartitioning()
{
	list<Edge *> oedges;
	list<Edge *> iedges;
	list<Edge *>::iterator it;
	long long *mapped = graph->getMapped();
	Vertex *v = graph->getVertexs();
	double *Upprice = network->getUpprice();

	// NOTE: No practical idea about whether "mirrors" is changed or not
	// NOTE: Suppose "mirrors" is not changed
	Mirror **mirrors = graph->getMirrors();
	double max_price_dc = max_value_index(Upprice, algorithm->DCnums); // NOTE: Get the dc which upload price is highest
	for (int i = 0; i < algorithm->getDCNums(); i++)
		for (int j = 0; j < graph->getNum_Nodes(); j++)
		{
			mirrors[i][j].reset();
		}
	// global_mvcost+=raw_data_size*Upprice[v[j].iniLabel];

	// NOTE: Identifying high degree node and low degree node
	// NOTE: 0 for high degree in-ging edge node, 1 for low degree in-going edge node
	for (int i = 0; i < graph->getNum_Nodes(); i++)
	{
		if (v[i].getIngoingEdges().size() >= theta)
			v[i].setType(0);
		else
			v[i].setType(1);

		if (v[i].iniLabel != max_price_dc)
			// MARK: Move what to the max price dc?, or why should the data move to max price dc
			// MARK: What is the raw data size?
			// NOTE: setting the wrost condition
			global_mvcost += raw_data_size * Upprice[v[i].iniLabel];
	}

	// NOTE: Using hybrid cut to part the nodes and edges
	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		for (int j = 0; j < graph->getNum_Nodes(); j++)
		{
			// cout<<j<<endl;
			// NOTE: if the node is in the DC and low degree, using edge-cut, which means marking
			if (v[j].getLabel() == i && v[j].getType() == 1)
			{
				// NOTE: for in-going edges
				iedges = v[j].getIngoingEdges();
				for (it = iedges.begin(); it != iedges.end(); it++)
				{
					// NOTE: if the source node is not in the DC
					if (v[mapped[(*it)->getsourceID()]].getLabel() != i)
					{
						// NOTE: then the source node is a mirror, and add it to the DC's mirror array
						mirrors[i][mapped[(*it)->getsourceID()]].setID(mapped[(*it)->getsourceID()]);
						mirrors[i][mapped[(*it)->getsourceID()]].addOutD();
					}
				}
				// NOTE: for out-going edges
				oedges = v[j].getOutgoingEdges();
				for (it = oedges.begin(); it != oedges.end(); it++)
				{
					if (v[mapped[(*it)->getdestID()]].getLabel() != i && v[mapped[(*it)->getdestID()]].getType() == 0)
					{
						mirrors[i][mapped[(*it)->getdestID()]].setID(mapped[(*it)->getdestID()]);
						mirrors[i][mapped[(*it)->getdestID()]].addInD();
					}
				}
			}
			// NOTE: if the node is high degree, using vertex cut
			// MARK: why only modifying out-going edge in this part?

			else if (v[j].getLabel() == i && v[j].getType() == 0)
			{
				oedges = v[j].getOutgoingEdges();
				for (it = oedges.begin(); it != oedges.end(); it++)
				{
					if (v[mapped[(*it)->getdestID()]].getLabel() != i && v[mapped[(*it)->getdestID()]].getType() == 0)
					{
						mirrors[i][mapped[(*it)->getdestID()]].setID(mapped[(*it)->getdestID()]);
						mirrors[i][mapped[(*it)->getdestID()]].addInD();
					}
				}
			}
		}
	}
	for (int j = 0; j < graph->getNum_Nodes(); j++)
	{
		if (v[j].type == 1) // low-degree
		{
			v[j].idegreeinlabel = v[j].ingoingEdges.size();
			v[j].odegreeinlabel = 0;

			oedges = v[j].outgoingEdges;
			for (it = oedges.begin(); it != oedges.end(); it++)
			{
				// NOTE: if the low-degree note connected to a low-degree node, and the source node is in the same DC with destination node
				if (v[mapped[(*it)->destID]].type == 1 && v[mapped[(*it)->destID]].label == v[j].label)
					v[j].odegreeinlabel++;
				// mirrors[originL][ver].addOutD();
				// NOTE: if the low-degree node connected to a high-degree node, then store the edge in the DC to release the stress of high-degree node
				else if (v[mapped[(*it)->destID]].type == 0)
					v[j].odegreeinlabel++;
				// mirrors[originL][ver].addOutD();
			}
		}
		else
		{
			iedges = v[j].ingoingEdges;
			v[j].idegreeinlabel = 0;
			v[j].odegreeinlabel = 0;
			for (it = iedges.begin(); it != iedges.end(); it++)
			{
				long long sourceId = mapped[(*it)->sourceID];
				// NOTE: if the source node is in the same DC, then store the edge
				if (v[sourceId].label == v[j].label)
					v[j].idegreeinlabel++;
				// mirrors[originL][ver].addInD();
			}
			oedges = v[j].outgoingEdges;
			for (it = oedges.begin(); it != oedges.end(); it++)
			{
				long long destId = mapped[(*it)->destID];
				// NOTE: if the dest node is high-degree node, or the dest node is low-degree node but in the same DC with it
				// NOTE: store the outgoing edge in the DC
				if (v[destId].type == 0 || (v[destId].type == 1 && v[destId].label == v[j].label))
					v[j].odegreeinlabel++;
				// mirrors[originL][ver].addOutD();
			}
		}
	}
}
void *test(void *arguments)
{
	Pthread_args *args = (Pthread_args *)arguments;
	cout << "id: " << args->id << endl;
	delete args;
}
void Mirrortest()
{
	Mirror **m = graph->getMirrors();
	list<Edge *> oedges;
	list<Edge *> iedges;
	list<Edge *>::iterator it;
	long long *mapped = graph->getMapped();
	Vertex *v = graph->getVertexs();
	Mirror **mirrors = new Mirror *[algorithm->getDCNums()];

	for (int i = 0; i < algorithm->getDCNums(); i++)
		mirrors[i] = new Mirror[graph->getNum_Nodes()];

	for (int i = 0; i < algorithm->getDCNums(); i++)
		for (int j = 0; j < graph->getNum_Nodes(); j++)
			mirrors[i][j].reset();

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		for (int j = 0; j < graph->getNum_Nodes(); j++)
		{
			// cout<<j<<endl;
			if (v[j].getLabel() == i && v[j].getType() == 1)
			{
				iedges = v[j].getIngoingEdges();
				for (it = iedges.begin(); it != iedges.end(); it++)
				{
					if (v[mapped[(*it)->getsourceID()]].getLabel() != i)
					{
						mirrors[i][mapped[(*it)->getsourceID()]].setID(mapped[(*it)->getsourceID()]);
						mirrors[i][mapped[(*it)->getsourceID()]].addOutD();
					}
				}
				oedges = v[j].getOutgoingEdges();
				for (it = oedges.begin(); it != oedges.end(); it++)
				{
					if (v[mapped[(*it)->getdestID()]].getLabel() != i && v[mapped[(*it)->getdestID()]].getType() == 0)
					{
						mirrors[i][mapped[(*it)->getdestID()]].setID(mapped[(*it)->getdestID()]);
						mirrors[i][mapped[(*it)->getdestID()]].addInD();
					}
				}
			}
			else if (v[j].getLabel() == i && v[j].getType() == 0)
			{
				oedges = v[j].getOutgoingEdges();
				for (it = oedges.begin(); it != oedges.end(); it++)
				{
					if (v[mapped[(*it)->getdestID()]].getLabel() != i && v[mapped[(*it)->getdestID()]].getType() == 0)
					{
						mirrors[i][mapped[(*it)->getdestID()]].setID(mapped[(*it)->getdestID()]);
						mirrors[i][mapped[(*it)->getdestID()]].addInD();
					}
				}
			}
		}
	}
	int k = 0;
	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		for (int j = 0; j < graph->getNum_Nodes(); j++)
		{
			if (mirrors[i][j].getID() != m[i][j].getID() || mirrors[i][j].getInD() != m[i][j].getInD() || mirrors[i][j].getOutD() != m[i][j].getOutD())
			{
				cout << "Mirrors前后不一致！" << endl;
				if (mirrors[i][j].getID() != m[i][j].getID())
				{
					cout << "顶点 " << j << " ID不一致！" << endl;
					cout << "ID: " << mirrors[i][j].getID() << "  " << m[i][j].getID() << endl;
				}
				if (mirrors[i][j].getInD() != m[i][j].getInD())
				{
					cout << "顶点 " << j << " 入度不一致！" << endl;
					cout << "InD: " << mirrors[i][j].getInD() << "  " << m[i][j].getInD() << endl;
				}
				if (mirrors[i][j].getOutD() != m[i][j].getOutD())
				{
					cout << "顶点 " << j << " 出度不一致！" << endl;
					cout << "OutD: " << mirrors[i][j].getOutD() << "  " << m[i][j].getOutD() << endl;
				}
				k = 1;
				break;
			}
		}
		if (k == 1)
			break;
	}
	if (k == 0)
		cout << "Mirrors前后一致！" << endl;
	for (int i = 0; i < algorithm->getDCNums(); i++)
		delete[] mirrors[i];
	delete[] mirrors;
}
void Simulate()
{
	list<Edge *> oedges;
	list<Edge *> iedges;
	list<Edge *>::iterator it;
	long long *mapped = graph->getMapped();
	double *Upprice = network->getUpprice();
	Vertex *v = graph->getVertexs();
	int *numOfVIM = new int[algorithm->getDCNums()]; //每个DC的顶点数，包括mirror
	int *numOfV = new int[algorithm->getDCNums()];	 //每个DC的顶点数，不包括mirror
	int *numOfE = new int[algorithm->getDCNums()];	 //每个DC的边的数目
	double *up_gather = new double[algorithm->getDCNums()];
	double *down_gather = new double[algorithm->getDCNums()];
	double *up_apply = new double[algorithm->getDCNums()];
	double *down_apply = new double[algorithm->getDCNums()];
	double movecost = 0;

	Mirror **mirrors = graph->getMirrors();
	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		numOfVIM[i] = 0;
		numOfV[i] = 0;
		numOfE[i] = 0;
		up_gather[i] = 0;
		up_apply[i] = 0;
		down_gather[i] = 0;
		down_apply[i] = 0;
	}
	double total_move_cost = 0;
	double max_price_dc = max_value_index(Upprice, algorithm->DCnums);

	for (int i = 0; i < graph->getNum_Nodes(); i++) //统计每个DC的master以及边的数目
	{
		numOfVIM[v[i].getLabel()]++;
		numOfV[v[i].getLabel()]++;
		if (v[i].getType() == 1) //处理low-degree的边
		{
			numOfE[v[i].getLabel()] += v[i].getIngoingEdges().size();
			oedges = v[i].getOutgoingEdges();
			for (it = oedges.begin(); it != oedges.end(); it++)
			{
				long long destId = mapped[(*it)->getdestID()];
				if (v[destId].getType() == 0)
					numOfE[v[i].getLabel()]++;
			}
		}
		else //处理high-degree的点
		{
			oedges = v[i].getOutgoingEdges();
			for (it = oedges.begin(); it != oedges.end(); it++)
			{
				long long destId = mapped[(*it)->getdestID()];
				if (v[destId].getType() == 0)
					numOfE[v[i].getLabel()]++;
			}
		}
		if (v[i].iniLabel != v[i].label)
			movecost += raw_data_size * Upprice[v[i].iniLabel];

		if (v[i].iniLabel != max_price_dc)
			total_move_cost += raw_data_size * Upprice[v[i].iniLabel];
	}

	for (int i = 0; i < algorithm->getDCNums(); i++)
		for (int j = 0; j < graph->getNum_Nodes(); j++)
		{
			if (mirrors[i][j].getID() != -1)
				numOfVIM[i]++;
		}
	long numofver = 0;
	for (int i = 0; i < algorithm->getDCNums(); i++)
		numofver += numOfVIM[i];
	// cout<<"-----------------------------------------------"<<endl;
	outfile << "-----------------------------------------------" << endl;
	cout << "#vertex replication: " << (double)numofver / graph->getNum_Nodes() << endl;
	cout << "#edge replication: " << 1 << endl;

	outfile << "#vertex replication: " << (double)numofver / graph->getNum_Nodes() << endl;
	outfile << "#edge replication: " << 1 << endl;

	long maxV = numOfV[0], minV = numOfV[0];
	long maxE = numOfE[0], minE = numOfE[0];

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		if (maxV < numOfV[i])
			maxV = numOfV[i];
		if (minV > numOfV[i])
			minV = numOfV[i];
		if (maxE < numOfE[i])
			maxE = numOfE[i];
		if (minE > numOfE[i])
			minE = numOfE[i];
	}

	cout << "#vertex diff: " << maxV - minV << endl;
	cout << "#edge diff: " << maxE - minE << endl;

	outfile << "#vertex diff: " << maxV - minV << endl;
	outfile << "#edge diff: " << maxE - minE << endl;

	double wan_usage = 0;
	if (budget < 0) // for powerlyra
	{
		/**************Gather 阶段**************************/
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			for (int j = 0; j < graph->getNum_Nodes(); j++)
				if (mirrors[i][j].getID() != -1 && v[j].getType() == 0)
				{
					wan_usage += dataunit;
					wan_usage += dataunit;
					up_gather[i] += dataunit;
					down_gather[v[j].getLabel()] += dataunit;
				}
		}
		/**************Apply 阶段**************************/
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			for (int j = 0; j < graph->getNum_Nodes(); j++)
			{
				if (mirrors[i][j].getID() != -1)
				{
					wan_usage += dataunit;
					wan_usage += dataunit;
					down_apply[i] += dataunit;
					up_apply[v[j].getLabel()] += dataunit;
				}
			}
		}
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			// cout<<"DC  "<<i<<" : up_gather: "<< up_gather[i]<<"; down_gather: "<<down_gather[i]<<"; up_apply: "<<up_apply[i]<<"; down_apply: "<<down_apply[i]<<endl;
		}
	}
	else // for LA-cut
	{
		/**************Gather 阶段**************************/
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			for (int j = 0; j < graph->getNum_Nodes(); j++)
				if (mirrors[i][j].getID() != -1 && v[j].getType() == 0 && mirrors[i][j].getInD() != 0)
				{
					wan_usage += dataunit;
					wan_usage += dataunit;
					up_gather[i] += dataunit;
					down_gather[v[j].getLabel()] += dataunit;
				}
		}
		/**************Apply 阶段**************************/
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			for (int j = 0; j < graph->getNum_Nodes(); j++)
			{
				if (mirrors[i][j].getID() != -1 && mirrors[i][j].getOutD() != 0)
				{
					wan_usage += dataunit;
					wan_usage += dataunit;
					down_apply[i] += dataunit;
					up_apply[v[j].getLabel()] += dataunit;
				}
			}
		}
	}
	cout << "-----------------------------------------------" << endl;
	outfile << "-----------------------------------------------" << endl;
	if (budget >= 0) // for LA-cut
	{
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			cout << "DC  " << i << " : uploadG: " << up_gather[i] << endl;
			outfile << "DC  " << i << " : uploadG: " << up_gather[i] << endl;
		}
		cout << "-----------------------------------------------" << endl;
		outfile << "-----------------------------------------------" << endl;
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			cout << "DC  " << i << " : downloadG: " << down_gather[i] << endl;
			outfile << "DC  " << i << " : downloadG: " << down_gather[i] << endl;
		}
		cout << "-----------------------------------------------" << endl;
		outfile << "-----------------------------------------------" << endl;
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			cout << "DC  " << i << " : uploadA: " << up_apply[i] << endl;
			outfile << "DC  " << i << " : uploadA: " << up_apply[i] << endl;
		}
		cout << "-----------------------------------------------" << endl;
		outfile << "-----------------------------------------------" << endl;
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			cout << "DC  " << i << " : downloadA: " << down_apply[i] << endl;
			outfile << "DC  " << i << " : downloadA: " << down_apply[i] << endl;
		}
	}
	else // for powerlyra
	{
		int *mapped = new int[algorithm->DCnums];
		int *map = new int[algorithm->DCnums];
		for (int i = 0; i < algorithm->DCnums; i++)
		{
			mapped[i] = -1;
			map[i] = -1;
		}
		int max_index_upG = max_value_index(up_gather, algorithm->DCnums);
		int min_index_upload = min_value_index_double(network->upload, algorithm->DCnums);
		int max_index_downG = max_value_index(down_gather, algorithm->DCnums);
		int min_index_download = min_value_index_double(network->download, algorithm->DCnums);
		int loadnum, dcnum;
		int loadnumA, dcnumA;
		if (up_gather[max_index_upG] / network->upload[min_index_upload] > down_gather[max_index_downG] / network->download[min_index_download])
		{
			mapped[max_index_upG] = min_index_upload;
			map[min_index_upload] = max_index_upG;
			loadnum = max_index_upG;
			dcnum = min_index_upload;
		}
		else
		{
			mapped[max_index_downG] = min_index_download;
			map[min_index_download] = max_index_downG;
			loadnum = max_index_downG;
			dcnum = min_index_download;
		}

		double max = up_apply[loadnum] / network->upload[dcnum] > down_apply[loadnum] / network->download[dcnum] ? up_apply[loadnum] / network->upload[dcnum] : down_apply[loadnum] / network->download[dcnum];
		int max_index_upA = max_value_index_besides(up_apply, algorithm->DCnums, loadnum);
		min_index_upload = min_value_index_double_besides(network->upload, algorithm->DCnums, dcnum);
		if (max > up_apply[max_index_upA] / network->upload[min_index_upload])
		{
			loadnumA = loadnum;
			dcnumA = dcnum;
		}
		else
		{
			loadnumA = max_index_upA;
			dcnumA = min_index_upload;
			max = up_apply[max_index_upA] / network->upload[min_index_upload];
		}
		int max_index_downA = max_value_index_besides(down_apply, algorithm->DCnums, loadnum);
		min_index_download = min_value_index_double_besides(network->download, algorithm->DCnums, dcnum);
		if (max < down_apply[max_index_downA] / network->download[min_index_download])
		{
			loadnumA = max_index_downA;
			dcnumA = min_index_download;
		}
		mapped[loadnumA] = dcnumA;
		map[dcnumA] = loadnumA;
		int start = 0;
		for (int m = 0; m < algorithm->DCnums; m++)
		{
			if (map[m] == -1)
			{
				for (int n = start; n < algorithm->DCnums; n++)
				{
					if (mapped[n] == -1)
					{
						map[m] = n;
						mapped[n] = m;
						start = n + 1;
						break;
					}
				}
			}
		}

		double maxArray[algorithm->DCnums];
		double maxArray2[algorithm->DCnums];
		double uploadnumsum[algorithm->DCnums];

		for (int j = 0; j < algorithm->DCnums; j++)
		{
			double divide = up_gather[map[j]] / network->upload[j];
			double divided = down_gather[map[j]] / network->download[j];
			maxArray[j] = divide > divided ? divide : divided;

			double divide2 = up_apply[map[j]] / network->upload[j];
			double divided2 = down_apply[map[j]] / network->download[j];
			maxArray2[j] = divide2 > divided2 ? divide2 : divided2;

			uploadnumsum[j] = up_gather[map[j]] + up_apply[map[j]];
		}

		timeAll = max_value(maxArray) + max_value(maxArray2);
		priceAll = sumWeights(uploadnumsum, network->upprice);

		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			cout << "DC  " << i << " : uploadG: " << up_gather[map[i]] << endl;
			outfile << "DC  " << i << " : uploadG: " << up_gather[map[i]] << endl;
		}
		cout << "-----------------------------------------------" << endl;
		outfile << "-----------------------------------------------" << endl;
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			cout << "DC  " << i << " : downloadG: " << down_gather[map[i]] << endl;
			outfile << "DC  " << i << " : downloadG: " << down_gather[map[i]] << endl;
		}
		cout << "-----------------------------------------------" << endl;
		outfile << "-----------------------------------------------" << endl;
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			cout << "DC  " << i << " : uploadA: " << up_apply[map[i]] << endl;
			outfile << "DC  " << i << " : uploadA: " << up_apply[map[i]] << endl;
		}
		cout << "-----------------------------------------------" << endl;
		outfile << "-----------------------------------------------" << endl;
		for (int i = 0; i < algorithm->getDCNums(); i++)
		{
			cout << "DC  " << i << " : downloadA: " << down_apply[map[i]] << endl;
			outfile << "DC  " << i << " : downloadA: " << down_apply[map[i]] << endl;
		}
		delete[] mapped;
		delete[] map;
	}

	cout << "-----------------------------------------------" << endl;
	outfile << "-----------------------------------------------" << endl;
	cout << "#wan for communication: " << wan_usage << endl;
	cout << "#time: " << timeAll << endl;
	cout << "#cost: " << priceAll << endl;
	cout << "#budget: " << budget << endl;
	cout << "#movecost: " << movecost << endl;
	cout << "#movecost is " << movecost / total_move_cost << " of the total_move_cost" << endl;
	cout << "#totalcost: " << movecost + priceAll * mvpercents << endl;
	cout << "#overhead: " << overhead << endl;
	// cout<<"#current mvcost: "<<current_mvcost<<endl;

	outfile << "#wan for communication: " << wan_usage << endl;
	outfile << "#time: " << timeAll << endl;
	outfile << "#cost: " << priceAll << endl;
	outfile << "#budget: " << budget << endl;
	outfile << "#movecost: " << movecost << endl;
	outfile << "# movecost is " << movecost / total_move_cost << " of the total_move_cost" << endl;
	outfile << "#totalcost: " << movecost + priceAll * mvpercents << endl;
	outfile << "#overhead: " << overhead << endl;
	outfile.close();
	delete[] numOfE;
	delete[] numOfV;
	delete[] numOfVIM;
	delete[] up_gather;
	delete[] up_apply;
	delete[] down_gather;
	delete[] down_apply;
}
void *keepMirrorsBackUpStable(void *arguments) //	mirrors_backup[thread_id][i]=mirrors[thread_id][i];
{
	/*Pthread_args * args=(Pthread_args* )arguments;
	unsigned int low=args->getLow();
	unsigned int high=args->getHigh();
	int dc=args->getDc();
	int thread_id=args->getId();
	Mirror** mirrors=graph->getMirrors();
	Mirror** mirrors_backup=graph->getBackUpMirrors();
	for(int i=low;i<=high;i++)
	{
		mirrors_backup[dc][i]=mirrors[dc][i];
	}


	delete args;*/
	// cout<<"Thread： "<<thread_id<<" finished!"<<endl;
}
void *keepMirrorsStable(void *arguments) // mirrors[thread_id][i]=mirrors_backup[thread_id][i];
{
	/*Pthread_args * args=(Pthread_args* )arguments;
	unsigned int low=args->getLow();
	unsigned int high=args->getHigh();
	int dc=args->getDc();
	int thread_id=args->getId();
	Mirror** mirrors=graph->getMirrors();
	Mirror** mirrors_backup=graph->getBackUpMirrors();
	for(int i=low;i<=high;i++)
	{
		mirrors[dc][i]=mirrors_backup[dc][i];
	}


	delete args;*/
	// cout<<"Thread： "<<thread_id<<" finished!"<<endl;
}
/*
1.����һ�������ƶ���cost�����������еĶ��㶼�Ƶ�ͬһ��DC

2.�޸ľ��ߣ�cost<=budget

3.����ʱ��̫�������㷨�ĳɲ���

*/
void topK()
{
	Mirror **mirrors = graph->mirrors;
	Vertex *v = graph->nodes;
	int K = randpro * graph->num_nodes;
	int *mirrornum = new int[graph->num_nodes];
	priority_queue<int, vector<A>, cmp1> minHeap;
	int count = 0;
	for (int i = 0; i < graph->num_nodes; i++)
	{
		mirrornum[i] = 0;
		for (int j = 0; j < algorithm->DCnums; j++)
		{
			if (mirrors[j][i].id != -1)
				mirrornum[i]++;
		}
		if (count < K)
		{
			A t;
			t.key = i;
			t.value = mirrornum[i];
			minHeap.push(t);
			trained.insert(pair<int, int>(i, 1));
		}
		else
		{
			if (mirrornum[i] <= minHeap.top().value)
				continue;
			else
			{
				A tmp = minHeap.top();
				A t;
				t.key = i;
				t.value = mirrornum[i];
				minHeap.pop();
				minHeap.push(t);
				trained.erase(tmp.key);
				trained.insert(pair<int, int>(i, 1));
			}
		}
		count++;
	}
	cout << "K :" << K << endl;
	cout << "Trained.size= " << trained.size() << endl;
	delete[] mirrornum;
}

int compare(const void *a, const void *b)
{
	return (*(int *)a - *(int *)b);
}

void revolver(int TUI_OR_SCRIPT) //�������û��ƽ��budget
{
	int iter = 0;
	int flag = 0;
	double *signal_ = NULL;
	double timeOld = 0;
	// NOTE: vertex is not sorted by number
	// NOTE: need to use with mapped array
	Vertex *v = graph->getVertexs();
	double *uploadnumG = new double[algorithm->getDCNums()];
	double *downloadnumG = new double[algorithm->getDCNums()];
	double *uploadnumA = new double[algorithm->getDCNums()];
	double *downloadnumA = new double[algorithm->getDCNums()];
	double *maxArray = new double[algorithm->getDCNums()];
	double *maxArray2 = new double[algorithm->getDCNums()];
	double *uploadnumsum = new double[algorithm->getDCNums()];
	double *Upload = network->getUpload();
	double *Download = network->getDownload();
	double *Upprice = network->getUpprice();

	// double ratio = -1;
	// ratio = -1;
	long percentage_degree = -1;
	// delete_zero_or_not = false; // false means not delete, true means delete
	double *SR = new double[MAX_ITER + 1];	 // SR for what
	double *TOEI = new double[MAX_ITER + 1]; // TOEI for what

	/* Newly add code start */
	// 获得不同百分比对应的度数
	// int *in_going_degree_percentage;
	if (test_character_num == 1)
	{
		// 以10%作为分隔，所以数组长度为10
		// 默认用10%，可以自定义，但是还没有实现
		// double ratio;
		long max_in_going_degree = -1;

		if (TUI_OR_SCRIPT == 1)
		{
			char yes_or_no;

			cout << "Please input the ratio need to test: ";
			cin >> ratio_of_nodes;
			cout << "Delete 0 or not(Y/n): ";
			cin >> yes_or_no;
			if (yes_or_no == 'y' || yes_or_no == 'Y')
			{
				deleted_zero_or_not = true;
			}
		}
		// Should have a else here, ignoring it now

		// in_going_degree_percentage = new int[11]; //release: abandon
		long long nodes_num = graph->getNum_Nodes();
		cout << "The number of nodes: " << nodes_num << endl;
		long *every_vertex_in_going_degree = new long[nodes_num];
		long long ratio_node_num = ratio_of_nodes * nodes_num;
		cout << "Ratio node num: " << ratio_node_num << endl;
		long *ratio_vertex_in_going_degree = new long[ratio_node_num];
		long long deleted_zero_size = 0;
		if (!deleted_zero_or_not)
		{
			for (long long i = 0; i < nodes_num; i++)
			{
				every_vertex_in_going_degree[i] = v[i].ingoingEdges.size();
				if (every_vertex_in_going_degree[i] > max_in_going_degree)
				{
					max_in_going_degree = every_vertex_in_going_degree[i];
				}
				if (i < ratio_node_num)
				{
					ratio_vertex_in_going_degree[i] = v[i].ingoingEdges.size();
				}
			}
		}
		else // using deleted_zero_size to control the array
		{
			for (long long i = 0; i < nodes_num; i++)
			{
				if (v[i].ingoingEdges.size() > 0)
				{
					every_vertex_in_going_degree[deleted_zero_size] = v[i].ingoingEdges.size();
					if (every_vertex_in_going_degree[deleted_zero_size] > max_in_going_degree)
					{
						max_in_going_degree = every_vertex_in_going_degree[deleted_zero_size];
					}
					if (deleted_zero_size < ratio_node_num)
					{
						ratio_vertex_in_going_degree[deleted_zero_size] = v[i].ingoingEdges.size();
					}
					deleted_zero_size++;
				}
			}
		}
		// dev: check the percentage of every degree number // release: maybe abandon
		long min_non_tree_index = ratio_node_num / 2 - 1;
		// int *the_number_of_degree = new int[nodes_num];
		// // long *degree_number = new long[nodes_num];
		// int *degree_map = new int[max_in_going_degree + 1];
		// for (int i = 0; i <= max_in_going_degree; i++)
		// {
		// 	degree_map[i] = -1;
		// }
		// for (int i = 0; i < nodes_num; i++)
		// {
		// 	the_number_of_degree[i] = 0;
		// 	// degree_number[i] = -1;
		// }
		// int null_position = 0;
		// for (int i = 0; i < nodes_num; i++)
		// {
		// 	if (degree_map[every_vertex_in_going_degree[i]] == -1)
		// 	{
		// 		the_number_of_degree[null_position]++;
		// 		degree_map[every_vertex_in_going_degree[i]] = null_position++;
		// 	}
		// 	else
		// 	{
		// 		the_number_of_degree[degree_map[every_vertex_in_going_degree[i]]]++;
		// 	}
		// }

		// dev: looking up the character of graph // release: abandon
		// cout << "Max in-going degree: " << max_in_going_degree << endl;
		// int total_num = 0;
		// for (int i = 0; i <= max_in_going_degree; i++)
		// {
		// 	if (degree_map[i] != -1)
		// 	{
		// 		double percentage = (double)the_number_of_degree[degree_map[i]] / nodes_num;
		// 		total_num += the_number_of_degree[degree_map[i]];
		// 		if (percentage > 0.1)
		// 		{
		// 			cout << "Degree " << i << ": " << percentage * 100 << "%" << endl;
		// 		}
		// 	}
		// }
		// dev: check the statistic method // release: abandon
		// if (total_num == nodes_num)
		// {
		// 	cout << "The result is correct" << endl;
		// }
		// else
		// {
		// 	cout <<	"Total num: " << total_num << endl;
		// 	cout << "Nodes num: " << nodes_num << endl;
		// 	cout << "Statistic method has problem" << endl;
		// }

		for (long i = min_non_tree_index; i >= 0; i--)
		{
			max_heaptify_long(ratio_vertex_in_going_degree, i, ratio_node_num);
		}
		if (deleted_zero_or_not == false)
		{
			for (long i = ratio_node_num; i < nodes_num; i++)
			{
				if (ratio_vertex_in_going_degree[0] > every_vertex_in_going_degree[i])
				{
					ratio_vertex_in_going_degree[0] = every_vertex_in_going_degree[i];
					max_heaptify_long(ratio_vertex_in_going_degree, 0, ratio_node_num);
				}
			}
			percentage_degree = ratio_vertex_in_going_degree[0];
		}
		else
		{
			for (long i = ratio_node_num; i < deleted_zero_size; i++)
			{
				if (ratio_vertex_in_going_degree[0] > every_vertex_in_going_degree[i])
				{
					ratio_vertex_in_going_degree[0] = every_vertex_in_going_degree[i];
					max_heaptify_long(ratio_vertex_in_going_degree, 0, ratio_node_num);
				}
			}
			percentage_degree = ratio_vertex_in_going_degree[0];
		}

		// dev: Find out the percentage of 0 degree nodes // release: abandon
		// long zero_degree_node_num = 0;
		// for (long long i = 0; i < nodes_num; i++)
		// {
		// 	if (every_vertex_in_going_degree[i] == 0)
		// 	{
		// 		zero_degree_node_num++;
		// 	}
		// }
		// cout << "zero degree nodes percentage: " << (double)zero_degree_node_num / nodes_num * 100 << "%" << endl;
		//  dev: Check out the heapify method // release: too slow to not use
		// if (!delete_zero_or_not)
		// {
		// 	// qsort(every_vertex_in_going_degree, nodes_num, sizeof(long), compare);
		// }
		// else
		// {
		// 	// qsort(every_vertex_in_going_degree, deleted_zero_size, sizeof(long), compare);
		// }

		cout << "max degree of in going degree: " << every_vertex_in_going_degree[nodes_num - 1] << endl;
		cout << "Heaptify " << ratio_of_nodes * 100 << "% in_going_degree: " << percentage_degree << endl;
		// cout << "Quick sort " << ratio * 100 << "% in going degree: " << every_vertex_in_going_degree[ratio_node_num] << endl;
		delete[] every_vertex_in_going_degree;
		delete[] ratio_vertex_in_going_degree;
		// delete[] the_number_of_degree;
		// delete[] degree_map;
	}
	else if (test_character_num == 5)
	{
		if (TUI_OR_SCRIPT == 1)
		{
			cout << "Please input the ratio need to test: ";
			cin >> ratio_of_nodes;
		}
		// in_going_degree_percentage = new int[11]; //release: abandon

		long long nodes_num = graph->getNum_Nodes();
		cout << "The number of nodes: " << nodes_num << endl;
		long *every_vertex_degree = new long[nodes_num];
		long long ratio_node_num = ratio_of_nodes * nodes_num;
		cout << "Ratio node num: " << ratio_node_num << endl;
		long *ratio_vertex_degree = new long[ratio_node_num];
		long max_degree = -1;
		for (long long i = 0; i < nodes_num; i++)
		{
			every_vertex_degree[i] = v[i].ingoingEdges.size() + v[i].outgoingEdges.size();
			if (every_vertex_degree[i] > max_degree)
			{
				max_degree = every_vertex_degree[i];
			}

			if (i < ratio_node_num)
			{
				ratio_vertex_degree[i] = v[i].ingoingEdges.size() + v[i].outgoingEdges.size();
			}
		}
		long min_non_tree_index = ratio_node_num / 2 - 1;
		for (long i = min_non_tree_index; i >= 0; i--)
		{
			max_heaptify_long(ratio_vertex_degree, i, ratio_node_num);
		}
		for (long i = ratio_node_num; i < nodes_num; i++)
		{
			if (ratio_vertex_degree[0] > every_vertex_degree[i])
			{
				ratio_vertex_degree[0] = every_vertex_degree[i];
				max_heaptify_long(ratio_vertex_degree, 0, ratio_node_num);
			}
		}

		percentage_degree = ratio_vertex_degree[0];

		// dev: Find out the percentage of 0 degree nodes release: abandon
		// long one_degree_node_num = 0;
		// for (long long i = 0; i < nodes_num; i++)
		// {
		// 	if (every_vertex_degree[i] == 0)
		// 	{
		// 		one_degree_node_num++;
		// 	}
		// }

		// dev: check the percentage of every degree number // release: maybe abandon
		/* Initialize start */
		// int *the_number_of_degree = new int[nodes_num];
		// int *degree_map = new int[max_degree + 1];
		// for (int i = 0; i <= max_degree; i++)
		// {
		// 	degree_map[i] = -1;
		// }
		// for (int i = 0; i < nodes_num; i++)
		// {
		// 	the_number_of_degree[i] = 0;
		// 	// degree_number[i] = -1;
		// }
		// /* Initialize end */
		// int null_position = 0;
		// for (int i = 0; i < nodes_num; i++)
		// {
		// 	if (degree_map[every_vertex_degree[i]] == -1)
		// 	{
		// 		the_number_of_degree[null_position]++;
		// 		degree_map[every_vertex_degree[i]] = null_position++;
		// 	}
		// 	else
		// 	{
		// 		the_number_of_degree[degree_map[every_vertex_degree[i]]]++;
		// 	}
		// }

		// dev: Seeing the character of the graph // release: abandon
		// cout << "Max in-going degree: " << max_degree << endl;
		// int total_num = 0;
		// for (int i = 0; i <= max_degree; i++)
		// {
		// 	if (degree_map[i] != -1)
		// 	{
		// 		double percentage = (double)the_number_of_degree[degree_map[i]] / nodes_num;
		// 		total_num += the_number_of_degree[degree_map[i]];
		// 		if (percentage > 0.1)
		// 		{
		// 			cout << "Degree " << i << ": " << percentage * 100 << "%" << endl;
		// 		}
		// 	}
		// }

		// dev: Check out the correction of heapify method // release: too slow to not use
		// qsort(every_vertex_degree, nodes_num, sizeof(long), compare);
		// cout << "zero degree nodes percentage: " << (double)one_degree_node_num / nodes_num * 100 << "%" << endl;

		cout << "max degree of degree: " << every_vertex_degree[nodes_num - 1] << endl;
		cout << "Heaptify " << ratio_of_nodes * 100 << "% degree: " << percentage_degree << endl;

		// cout << "Quick sort " << ratio * 100 << "% in going degree: " << every_vertex_degree[ratio_node_num] << endl;
		delete[] every_vertex_degree;
		delete[] ratio_vertex_degree;
		// delete[] the_number_of_degree;
		// delete[] degree_map;
	}

	else if (test_character_num > 0)
	{
		cout << "Working..." << endl;
	}

	/* Newly add code end */

	initPartitioning();
	calculateTimeAndPrice(timeAll, priceAll);
	Mirror **mirrors = graph->getMirrors();
	Partition *partition = network->getPartition();
	pthread_t tid[MAX_THREADS_NUM];
	pthread_t tidd[MAX_THREADS_NUM];
	pthread_t tid_bsp[batchsize];
	int *eachThreadsNum = new int[MAX_THREADS_NUM];
	cout << "batchsize: " << batchsize << endl;
	cout << "before executing the algorithm:" << endl;
	cout << "#timeAll: " << timeAll << endl;
	cout << "#priceAll: " << priceAll << endl;
	cout << "The seed is: " << Seed << endl;
	// ofstream outfile(outputfile);
	// NOTE: outfile open the file named outputfile
	// NOTE: outputfile is declared in the argument 8th
	outfile = ofstream(outputfile);
	// outfile<<"asdf"<<endl;
	ofstream outtrain(FN);
	if (priceAll <= budget)
	{
		F1 = 0;
		F2 = 1;
		flag = 0;
	}
	else
	{
		F1 = 1;
		F2 = 0;
		flag = 1;
	}
	// thb
	// outfile<<"The seed is: "<<Seed<<endl;
	for (int i = 0; i < graph->getNum_Nodes(); i++)
		v[i].initSignal_();
	srand(time(0));
	// NOTE: initialize the multiprocessing
	multiprocessing_pool(eachThreadsNum); //并行处理顶点
										  //  trained.clear();
	// topK();

	for (int i = 0; i < graph->getNum_Nodes(); i++)
		v[i].action = v[i].label;

	budget = global_mvcost * budget + mvpercents * ginger_cost;
	cout << "The budget is " << budget << endl;
	while (iter <= MAX_ITER)
	{
		/********************************计算当前Gather、Apply阶段各个DC的数据传传输量*************************************/
		for (int i = 0; i < algorithm->DCnums; i++)
		{
			uploadnumG[i] = ug[i];
			downloadnumG[i] = dg[i];
			uploadnumA[i] = ua[i];
			downloadnumA[i] = da[i];
			maxArray[i] = 0;
			// cout<<ug[i]<<" "<<dg[i]<<" "<<ua[i]<<" "<<da[i]<<endl;
		}

		time_t start = time(NULL);
		// double t1= 0;  //action select时间
		// double t2= 0;  //vertex move时间
		// double t3= 0;  //sampling时间
		// double t4= 0;  //probability update时间

		if (iter > 0)
		{
			// ActionSelect();  //��ÿ��vertex����action����
			/************************************并行执行MyActionSelect***************************************/
			// cout<<"--------------------------"<<"Action Select Starts!"<<"-----------------------------"<<endl;
			// NOTE: high1 - low1 = eachThreadNum
			// NOTE: status1 get the multithread status
			int low1 = 0;
			int high1 = eachThreadsNum[0] - 1;
			int status1 = 0;
			start = time(NULL);
			high_resolution_clock::time_point beginTime = high_resolution_clock::now();
			for (int j = 0; j < MAX_THREADS_NUM; j++)
			{
				Pthread_args *args = new Pthread_args();
				args->setLow(low1);
				args->setHigh(high1);
				args->setIter(iter);
				status1 = pthread_create(&tid[j], NULL, MyActionSelect, (void *)args);
				low1 = high1 + 1;
				if (j < MAX_THREADS_NUM - 1)
					high1 = low1 + eachThreadsNum[j + 1] - 1;
			}

			for (int j = 0; j < MAX_THREADS_NUM; j++)
				pthread_join(tid[j], NULL);

			// t1=time(NULL)-start;
			high_resolution_clock::time_point endTime = high_resolution_clock::now();
			milliseconds timeInterval = std::chrono::duration_cast<milliseconds>(endTime - beginTime);
			t1 = timeInterval.count();
			// cout<<"--------------------------"<<"Action Select Finished!"<<"-----------------------------"<<endl;
			// MyActionSelect();
			/********************************************************************************************************/
			// moverand=rand()%100/(double)99;
			/**********************执行顶点迁移***************************/
			timeOld = timeAll;
			priceOld = priceAll;
			if (bsp == 0)
			{
				// cout<<"--------------------------"<<"Vertex Migrate Starts--no bsp!"<<"-----------------------------"<<endl;
				int percent = 0;
				int decidedtomove = 0;
				start = time(NULL);
				map<int, int>::iterator itora = trained.begin();
				while (itora != trained.end())
				{
					int i = itora->first;
					int ori = v[i].label;
					// start=time(NULL);
					if (v[i].label != v[i].action)
					{
						decidedtomove++;
						vertexMigrate(i, downloadnumG, uploadnumG, downloadnumA, uploadnumA, v[i].action);

						/************************Gather阶段Time***********************************/

						for (int j = 0; j < algorithm->DCnums; j++)
						{
							double divide = uploadnumG[j] / Upload[j];
							double divided = downloadnumG[j] / Download[j];
							maxArray[j] = divide > divided ? divide : divided;
						}

						double time = max_value(maxArray);

						/************************Apply阶段Time***********************************/
						for (int j = 0; j < algorithm->DCnums; j++)
						{
							double divide = uploadnumA[j] / Upload[j];
							double divided = downloadnumA[j] / Download[j];
							maxArray[j] = divide > divided ? divide : divided;
						}
						time += max_value(maxArray);

						for (int j = 0; j < algorithm->DCnums; j++)
						{
							uploadnumsum[j] = uploadnumG[j] + uploadnumA[j];
						}
						double cost = sumWeights(uploadnumsum, Upprice);

						if (time <= timeAll)
						{
							timeAll = time;
							priceAll = cost;
						}
						else //并行地维护两个DC的顶点i的邻居的mirrors
						{
							percent++;
							vertexMigrate(i, downloadnumG, uploadnumG, downloadnumA, uploadnumA, ori);
						}
					}
					itora++;
				}

				for (int i = 0; i < algorithm->DCnums; i++)
				{
					ug[i] = uploadnumG[i];
					dg[i] = downloadnumG[i];
					ua[i] = uploadnumA[i];
					da[i] = downloadnumA[i];
				}

				t2 = time(NULL) - start;
				// cout<<"有 "<<(double)percent/decidedtomove<<" 的顶点是没有真正移动的"<<endl;

				// cout<<"--------------------------"<<"Vertex Migrate Finished!--no bsp"<<"-----------------------------"<<endl;
			}
			else if (bsp == 1)
			{
				// cout<<"--------------------------"<<"Vertex Migrate Starts-- bsp!"<<"-----------------------------"<<endl;
				for (int j = 0; j < batchsize; j++)
				{
					Pthread_args *args = new Pthread_args();
					args->id = j;
					status1 = pthread_create(&tid_bsp[j], NULL, moveVertexBsp, (void *)args);
				}
				map<int, int>::iterator itora = trained.begin();
				int *verID = new int[batchsize];
				int *sig = new int[batchsize]; // 0不移，1移
				double **ulg = new double *[batchsize];
				double **dlg = new double *[batchsize];
				double **ula = new double *[batchsize];
				double **dla = new double *[batchsize];
				for (int i = 0; i < batchsize; i++)
				{
					ulg[i] = new double[algorithm->DCnums];
					dlg[i] = new double[algorithm->DCnums];
					ula[i] = new double[algorithm->DCnums];
					dla[i] = new double[algorithm->DCnums];
				}
				start = time(NULL);
				high_resolution_clock::time_point beginTime = high_resolution_clock::now();
				double o = 0;
				time_t s = time(NULL);
				int vertype = 0;
				while (itora != trained.end())
				{
					int numver = 0;
					while (itora != trained.end() && numver < batchsize)
					{
						if (v[itora->first].label != v[itora->first].action && v[itora->first].type == vertype)
						{
							verID[numver] = itora->first;
							sig[numver] = 0;
							numver++;
						}
						itora++;
					}
					if (itora == trained.end() && vertype == 0)
					{
						itora = trained.begin();
						vertype = 1;
					}
					for (int i = 0; i < numver; i++)
					{
						for (int j = 0; j < algorithm->DCnums; j++)
						{
							ulg[i][j] = ug[j];
							dlg[i][j] = dg[j];
							ula[i][j] = ua[j];
							dla[i][j] = da[j];
						}
					}

					/*-------------------并行部分***************************/
					s = time(NULL);
					for (int j = 0; j < numver; j++)
					{
						// cout<<"j: "<<j<<endl;
						Pthread_args args = Pthread_args();
						args.id = j;
						args.iteration = iter;
						args.verindex = j;
						args.ver = verID[j];
						args.downloadnumG = dlg[j];
						args.uploadnumG = ulg[j];
						args.downloadnumA = dla[j];
						args.uploadnumA = ula[j];
						args.destdc = v[verID[j]].action;
						args.s = sig;
						bspvec[j].push_back(args);
						// moveVertex(verID[j],dlg[j],ulg[j],dla[j],dla[j],v[verID[j]].action);
						// status1 = pthread_create(&tid_bsp[j], NULL, moveVertexBsp , (void*)args);
						// status1 = pthread_create(&tid_bsp[j], NULL, test , (void*)args);
					}
					while (bsplock != numver)
					{
						// cout<<"bsplock: "<<bsplock<<" ; numver: "<<numver<<endl;
						//  getchar();
					}
					// cout<<"hello"<<endl;
					bsplock = 0;
					/*  for(int j=0;j<numver;j++)
			pthread_join(tid_bsp[j],NULL);*/

					o += time(NULL) - s;
					/*-------------------并行部分结束***************************/

					for (int i = 0; i < numver; i++)
					{

						if (sig[i] == 1)
						{
							if (iter == 1 && v[verID[i]].iniLabel != v[verID[i]].action)
							{
								current_mvcost += raw_data_size * Upprice[v[verID[i]].iniLabel];
							}
							else if (iter > 1)
							{
								if (v[verID[i]].iniLabel == v[verID[i]].label && v[verID[i]].label != v[verID[i]].action)
									current_mvcost += raw_data_size * Upprice[v[verID[i]].iniLabel];
								else if (v[verID[i]].iniLabel != v[verID[i]].label && v[verID[i]].label != v[verID[i]].action)
								{
									if (v[verID[i]].iniLabel == v[verID[i]].action)
										current_mvcost -= raw_data_size * Upprice[v[verID[i]].iniLabel];
								}
							}
							vertexMigrate(verID[i], dg, ug, da, ua, v[verID[i]].action);
							/************************Gather阶段Time***********************************/

							for (int j = 0; j < algorithm->DCnums; j++)
							{
								double divide = ug[j] / Upload[j];
								double divided = dg[j] / Download[j];
								maxArray[j] = divide > divided ? divide : divided;

								double divide2 = ua[j] / Upload[j];
								double divided2 = da[j] / Download[j];
								maxArray2[j] = divide2 > divided2 ? divide2 : divided2;

								uploadnumsum[j] = ug[j] + ua[j];
							}

							double time = max_value(maxArray) + max_value(maxArray2);

							/************************Apply阶段Time***********************************/
							/*   for (int j = 0; j < algorithm->DCnums; j++)
			{
				double divide=ua[j] / Upload[j];
				double divided=da[j] / Download[j];
				maxArray[j] =divide > divided ? divide : divided;
			}
			time+=max_value(maxArray);*/

							/* for (int j = 0; j < algorithm->DCnums; j++)
			{
				uploadnumsum[j] =ug[j]+ua[j];
			}*/
							double cost = sumWeights(uploadnumsum, Upprice);
							timeAll = time;
							priceAll = cost;
						}
					}
				}
				for (int j = 0; j < batchsize; j++)
					pthread_cancel(tid_bsp[j]);

				for (int j = 0; j < batchsize; j++)
					pthread_join(tid_bsp[j], NULL);

				// cout<<"judge: "<<o<<endl;
				//  t2=time(NULL)-start;
				high_resolution_clock::time_point endTime = high_resolution_clock::now();
				milliseconds timeInterval = std::chrono::duration_cast<milliseconds>(endTime - beginTime);
				t2 = timeInterval.count();
				delete[] verID;
				delete[] sig;
				for (int i = 0; i < batchsize; i++)
				{
					delete[] ulg[i];
					delete[] dlg[i];
					delete[] ula[i];
					delete[] dla[i];
				}
				delete[] ulg;
				delete[] dlg;
				delete[] ula;
				delete[] dla;
				//	cout<<"--------------------------"<<"Vertex Migrate Finished!-- bsp"<<"-----------------------------"<<endl;
			}
			/*if(flag==1)    //大于budget
		{
		   F1=1;
		   F2=0;
		}
		else
		{
			double impro=(timeOld-timeAll)/timeOld;
			if(F1==1)
			{
				F1=1;
				F2=0;
			}
			else if(priceAll>budget||impro<=0.05)
			{
				F1=1;
				F2=0;
			}
		}*/
			TOEI[iter - 1] = (t1 + t2 + t3 + t4 + t5) / 1000;
			SR[iter - 1] = randpro;
			double usedtime = 0;
			for (int p = 0; p < iter; p++)
				usedtime += TOEI[p];
			// MARK: leftteach for what, L for what
			double lefteach = (L - usedtime) / (MAX_ITER - iter);
			cout << "--------lefteach:" << lefteach << endl;
			int recent = 3;
			if (lefteach <= 0)
				randpro = 0;
			else
			{
				if (iter >= recent)
				{
					double total_sr = 0;
					for (int p = iter - 1; p >= iter - recent; p--)
						total_sr += lefteach / TOEI[p] * SR[p];
					randpro = total_sr / recent;
					cout << "--------total_sr:" << total_sr << endl;
				}
				else
				{
					double total_sr = 0;
					for (int p = 0; p <= iter - 1; p++)
						total_sr += lefteach / TOEI[p] * SR[p];
					randpro = total_sr / recent;
				}
			}
			// randpro=0.01;
			cout << "--------randpro:" << randpro << endl;
			if (randpro > 1)
				randpro = 1;
		} //(iter>0)
		/********************************并行重新选择被训练的点***************************************/
		// cout<<"--------------------------"<<"Sampling Starts!"<<"-----------------------------"<<endl;

		if (true)
		{
			trained.clear();
			// srand(time(NULL));
			/*	for(int i=0;i<graph->getNum_Nodes();i++)
	   {
		double random=rand()%100/(double)99;
		if(random<=randpro)
		{
			trained.insert(pair<int,int>(i,1));
		}
	  }*/
			int low1 = 0;
			int high1 = eachThreadsNum[0] - 1;
			int status1 = 0;
			start = time(NULL);
			high_resolution_clock::time_point beginTime = high_resolution_clock::now();
			for (int j = 0; j < MAX_THREADS_NUM; j++)
			{
				Pthread_args *args = new Pthread_args();
				args->setLow(low1);
				args->setHigh(high1);
				args->setIter(iter);
				args->setPercentage(ratio_of_nodes);
				args->setTestCharacter(test_character_num);
				args->setThreshold(percentage_degree);
				// MARK: PLACE NEED TO CHANGE
				if (test_character_num < 1 || test_character_num > 5)
				{
					status1 = pthread_create(&tid[j], NULL, Sampling, (void *)args);
				}
				else
				{
					if (test_character_num == 1 || test_character_num == 5)
					{
						if (deleted_zero_or_not)
						{
							args->setDeleteZero(1); // means the zero is deleted
						}
						else
						{
							args->setDeleteZero(0); // means the zero is reserve
						}
					}
					// cout << "Thread " << j << ":" << endl;
					// cout << "Threshold: " << percentage_degree << endl;
					// cout << "Test_character_num: " << test_character_num << endl;
					// cout << "Ratio: " << ratio << endl;
					status1 = pthread_create(&tid[j], NULL, Sampling_adv, (void *)args);
				}
				low1 = high1 + 1;
				if (j < MAX_THREADS_NUM - 1)
					high1 = low1 + eachThreadsNum[j + 1] - 1;
			}

			pthread_counter = 0;
			for (int j = 0; j < MAX_THREADS_NUM; j++)
				pthread_join(tid[j], NULL);

			high_resolution_clock::time_point endTime = high_resolution_clock::now();
			milliseconds timeInterval = std::chrono::duration_cast<milliseconds>(endTime - beginTime);
			t3 = timeInterval.count();
			// t3=time(NULL)-start;
			t3 = 0;
		}
		else
		{
			trained.clear();
			srand(time(NULL));
			topK();
		}
		// cout<<"--------------------------"<<"Sampling Finished!"<<"-----------------------------"<<endl;

		//---之前的if(iter>0)
		cout << "trained size: " << trained.size() << endl;
		/**************************************Debug**********************************************/
		/*for (int i = 0; i < algorithm->getDCNums(); i++)
			{
			   cout<<"DC "<<i<<": "<<"uploadG: "<<uploadnumG[i]<<" ，downloadG: "<<downloadnumG[i]<<" , uploadA: "<<uploadnumA[i]<<" , downloadA: "<<downloadnumA[i]<<endl;
			}
			for (int j = 0; j < algorithm->getDCNums(); j++)
			{
				maxArray[j] =uploadnumG[j] / Upload[j] > downloadnumG[j] / Download[j] ? uploadnumG[j] / Upload[j] : downloadnumG[j] / Download[j];
			}

			double time = max_value(maxArray);
			int la=max_value_index(maxArray,algorithm->getDCNums());
			cout<<"Gather阶段最大传输时间的DC是： "<<la<<" : "<<time<<endl;
			  for (int j = 0; j < algorithm->getDCNums(); j++)
			{
				maxArray[j] =uploadnumA[j] / Upload[j] > downloadnumA[j] / Download[j] ? uploadnumA[j] / Upload[j] : downloadnumA[j] / Download[j];
			}
			time=max_value(maxArray);
			la=max_value_index(maxArray,algorithm->getDCNums());
			cout<<"Apply阶段最大传输时间的DC是： "<<la<<" : "<<time<<endl;*/
		/*****************************************************************************************/
		// initPartitioning();
		// reCalTimeAndCost(timeAll,priceAll);

		cout << "the " << iter << " iteration :" << endl;
		cout << " #the sampling rate: " << randpro << endl;
		// cout << "stepbudget: " << stepbudget << endl;
		cout << "#timeAll: " << timeAll << ";#priceAll: " << priceAll << endl;
		cout << "#mvcost: " << current_mvcost << " ;#totalcost: " << current_mvcost + mvpercents * priceAll << endl;
		// thb
		/*
			outfile << "the " << iter  << " iteration :" << endl;
			outfile << "#timeAll: " << timeAll << ";#priceAll: " << priceAll << endl;
			*/
		// NOTE: outfile for arg[8]
		// MARK: Reset the format
		if (iter == 0)
		{
			outfile << "Iteration,Time_cost,Price_cost,Overhead" << endl;
		}

		outfile << iter << ',' << timeAll << ',' << priceAll << ',';
		if (iter == 0)
			outfile << endl;

		// pthread_barrier_init(&barrier, NULL, MAX_THREADS_NUM);
		if (iter % 4 == 0)
		{
			betaC = (iter / 4) * 0.25;
			alphaT = 1 - betaC;
		}
		// cout<<"alpha: "<<alphaT<<" ,beta: "<<betaC<<endl;
		moverand = rand() % 100 / (double)99;
		cout << "moverand: " << moverand << endl;
		if (moverand <= 0.5)
			cout << "------------------------obey-----------------------" << endl;
		/**********************************并行执行score计算以及objective function**********************************/
		//	cout<<"--------------------------"<<"Score And Objective Function Starts!"<<"-----------------------------"<<endl;
		// cout<<"F1: "<<F1<<" , F2: "<<F2<<endl;
		int low = 0;
		int high = eachThreadsNum[0] - 1;
		int status = 0;
		high_resolution_clock::time_point beginTime = high_resolution_clock::now();
		start = time(NULL);
		for (int j = 0; j < MAX_THREADS_NUM; j++)
		{
			// cout<<"low: "<<low<<", high: "<<high<<endl;
			Pthread_args *args = new Pthread_args();
			args->setId(j);
			args->setLow(low);
			args->setHigh(high);
			args->setIter(iter);
			// STOPLINE
			status = pthread_create(&tidd[j], NULL, Parallel_Score_Signal_function, (void *)args);
			low = high + 1;
			if (j < MAX_THREADS_NUM - 1)
				high = low + eachThreadsNum[j + 1] - 1;
		}
		// cout<<"Main Thread finished!"<<endl;
		for (int j = 0; j < MAX_THREADS_NUM; j++)
			pthread_join(tidd[j], NULL);
		high_resolution_clock::time_point endTime = high_resolution_clock::now();
		milliseconds timeInterval = std::chrono::duration_cast<milliseconds>(endTime - beginTime);
		t5 = timeInterval.count();
		// t5=time(NULL)-start;
		// cout<<"All threads finished!"<<endl;
		// getchar();
		//	  cout<<"--------------------------"<<"Score And Objective Function Finished!"<<"-----------------------------"<<endl;
		/************************************************************************************************************/

		/*******************************************并行执行概率更新************************************************/
		//	cout<<"--------------------------"<<"Probability Update Starts!"<<"-----------------------------"<<endl;
		low = 0;
		high = eachThreadsNum[0] - 1;
		status = 0;
		start = time(NULL);
		beginTime = high_resolution_clock::now();
		for (int j = 0; j < MAX_THREADS_NUM; j++)
		{
			// cout<<"low: "<<low<<", high: "<<high<<endl;
			Pthread_args *args = new Pthread_args();
			args->setLow(low);
			args->setHigh(high);
			args->setIter(iter);
			status = pthread_create(&tid[j], NULL, probability_Update, (void *)args);
			low = high + 1;
			if (j < MAX_THREADS_NUM - 1)
				high = low + eachThreadsNum[j + 1] - 1;
		}
		for (int j = 0; j < MAX_THREADS_NUM; j++)
			pthread_join(tid[j], NULL);
		endTime = high_resolution_clock::now();
		timeInterval = std::chrono::duration_cast<milliseconds>(endTime - beginTime);
		t4 = timeInterval.count();
		// t4=time(NULL)-start;
		//     cout<<"--------------------------"<<"Probability Update Finished!"<<"-----------------------------"<<endl;
		// probability_Update();
		/******************************************************************************************************/
		/*-----------------------------------------输出概率测试--------------------------------------------*/
		/*  double** probability = algorithm->getProbability();
		for(int w=0;w<algorithm->getDCNums();w++)
		cout<<probability[0][w]<<"     ";
		cout<<endl;
		getchar();*/

		// probability_Update();
		if (iter >= 1)
		{
			cout << "Action Select: " << t1 << endl;
			cout << "Vertex Move: " << t2 << endl;
			cout << "Sampling: " << t3 << endl;
			cout << "Probability Update: " << t4 << endl;
			cout << "Score Function: " << t5 << endl;
			double T = t1 + t2 + t3 + t4 + t5;
			overhead += T;
			// thb
			// outfile << "The overhead is: "<< T  << endl;
			// outfile << T << endl;
			outfile << TOEI[iter - 1] << endl;
			cout << "******************************************************************" << endl;
			cout << "The overhead is: " << TOEI[iter - 1] << endl;
			cout << "The T is: " << T << endl;
			cout << "******************************************************************" << endl;
		}
		outtrain << endl;
		outtrain << "---------------------------"
				 << "第 " << iter << " 次迭代"
				 << "----------------------------------" << endl;
		outtrain << endl;
		map<int, int>::iterator itertemp;
		itertemp = trained.begin();
		int debug = 1;
		while (itertemp != trained.end())
		{
			if (debug % 25 == 0)
				outtrain << itertemp->first << "   " << endl;
			else
				outtrain << itertemp->first << "   ";
			debug++;
			itertemp++;
		}
		outtrain << endl;
		iter++;

		// delete args;
	}
	overhead = 0;
	for (int p = 0; p < MAX_ITER; p++)
		overhead += TOEI[p];
	// outfile.close();
	outtrain.close();
	delete[] eachThreadsNum;
	delete[] uploadnumA;
	delete[] uploadnumG;
	delete[] downloadnumA;
	delete[] downloadnumG;
	delete[] maxArray;
	delete[] maxArray2;
	delete[] uploadnumsum;
}
void *Parallel_Score_Signal_function(void *arguments) //可并行
{
	Pthread_args *args = (Pthread_args *)arguments;
	unsigned int low = args->low;
	unsigned int high = args->high;
	// srand(time(NULL) - low);
	int thread_id = args->id;
	int iter = args->iteration;
	double *uploadnumG = new double[algorithm->DCnums];
	double *downloadnumG = new double[algorithm->DCnums];
	double *uploadnumA = new double[algorithm->DCnums];
	double *downloadnumA = new double[algorithm->DCnums];
	double *uploadnumsum = new double[algorithm->DCnums];
	long long *mapped = graph->mapped;
	double *maxArray = new double[algorithm->DCnums];
	double *maxArray2 = new double[algorithm->DCnums];
	double *Upload = network->upload;
	double *Download = network->download;
	double *Upprice = network->upprice;
	Vertex *v = graph->nodes;
	list<Edge *>::iterator it;
	// int * labels=new int[graph->getNum_Nodes()];
	Mirror **mirrors = graph->mirrors;
	list<Edge *> oedges;
	list<Edge *> iedges;
	double *temp_score_array = new double[algorithm->DCnums];

	for (int i = 0; i < algorithm->DCnums; i++)
	{
		uploadnumG[i] = ug[i];
		downloadnumG[i] = dg[i];
		uploadnumA[i] = ua[i];
		downloadnumA[i] = da[i];
	}

	double totalwan = 0;
	for (int i = 0; i < algorithm->DCnums; i++)
	{
		totalwan = totalwan + ug[i] + dg[i] + ua[i] + da[i];
	}

	vector<int> vec = PSS[thread_id];
	int vecsize = vec.size();

	for (int p = 0; p < vecsize; p++)
	{
		int ver = vec[p];
		// for(int ver=low;ver<=high;ver++)
		//{
		map<int, int>::iterator itora = trained.find(ver);
		if (itora != trained.end())
		{
			// int originL=v[ver].getLabel();
			// vector<int> vec;
			for (int i = 0; i < algorithm->DCnums; i++)
			{
				double rr = 0;
				if (v[ver].iniLabel == v[ver].label && v[ver].label != i)
					rr = raw_data_size * Upprice[v[ver].iniLabel];
				else if (v[ver].iniLabel != v[ver].label && v[ver].label != i && v[ver].iniLabel == i)
					rr = -raw_data_size * Upprice[v[ver].iniLabel];

				if (i == v[ver].label)
				{
					// double bw=(2.0 * rand() / RAND_MAX - 1.0)*0.000000001/10;
					//  double aw=bw<0?1+bw:1-bw;
					// if(priceAll<=budget&&!(iter>=MAX_ITER-10&&priceAll>budget))
					// temp_score_array[i] = bw*(rand()%100/(double)99);
					// else
					temp_score_array[i] = 0;
					// if(betaC<=0.75)
					// temp_score_array[i] = betaC*(budget-priceAll)/budget;
					/*else
				{
					temp_score_array[i]=priceAll>budget?-1000:0;
				}*/
					// temp_score_array[i] = priceAll;
					//  migration[i]=priceAll;
				}
				else
				{
					//---------------------------------将ver移至DC i---------------------------------
					// moveVertex(ver,downloadnumG,uploadnumG,downloadnumA,uploadnumA,i);
					moveVertex(ver, downloadnumG, uploadnumG, downloadnumA, uploadnumA, i);
					//计算score
					/************************Gather阶段Time***********************************/
					for (int j = 0; j < algorithm->DCnums; j++)
					{
						double divide = uploadnumG[j] / Upload[j];
						double divided = downloadnumG[j] / Download[j];
						maxArray[j] = divide > divided ? divide : divided;

						double divide2 = uploadnumA[j] / Upload[j];
						double divided2 = downloadnumA[j] / Download[j];
						maxArray2[j] = divide2 > divided2 ? divide2 : divided2;

						uploadnumsum[j] = uploadnumA[j] + uploadnumG[j];
					}

					double time = max_value(maxArray) + max_value(maxArray2);
					// double minT=min_value(maxArray);

					/************************Apply阶段Time***********************************/
					/*   for (int j = 0; j < algorithm->DCnums; j++)
			{
				double divide=uploadnumA[j] / Upload[j];
				double divided=downloadnumA[j] / Download[j] ;
				maxArray[j] =divide > divided ? divide : divided;
			}
			time+=max_value(maxArray);*/
					// minT+=min_value(maxArray);
					/*  for (int j = 0; j < algorithm->DCnums; j++)
			{
				uploadnumsum[j] =uploadnumA[j] +uploadnumG[j];
			}*/
					double cost = sumWeights(uploadnumsum, Upprice);
					double cost_A = priceAll * mvpercents + current_mvcost;
					double cost_B = cost * mvpercents + current_mvcost + rr;
					double wan = 0;
					for (int j = 0; j < algorithm->DCnums; j++)
					{
						wan = wan + uploadnumG[j] + downloadnumG[j] + uploadnumA[j] + downloadnumA[j];
					}
					// temp_score_array[i] = timeAll - time;
					int overbudget = 0;
					if (priceAll > budget)
						overbudget = 1;
					double downbudget = overbudget == 1 ? 0 : 1;
					double costdown = cost > budget ? cost : budget;
					double over = cost > budget ? 1 : 0;
					// temp_score_array[i]=(timeAll-time)/timeAll+2*overbudget*(priceAll-costdown)/priceAll+downbudget*over*-1000;
					// if(betaC<=0.75)
					// temp_score_array[i]=alphaT*((timeAll-time)/timeAll)+betaC*((budget-cost)/budget);
					// cout<<"budegt: "<<budget<<endl;
					// double belta=(priceAll-budget)/priceAll;
					double belta = (cost_A - budget) / cost_A;
					double alpha = 1 - belta;
					if (belta <= 0)
					{
						belta = 0;
						alpha = 1;
					}
					// double firstw,secondw;
					if (alpha == 1)
					{
						/*if(moverand<=0.5)
			   {
				   alpha=0;
				   belta=1;
			   }*/
						// secondw=(2.0 * rand() / RAND_MAX - 1.0)*0.000000001/10;
						// firstw=secondw<0?1+secondw:1-secondw;
						/*  belta=iter*(1.0/(MAX_ITER-1));
			alpha=1-belta;
			if(belta>0.5)
			{
				belta=belta-0.5;
				alpha=1-belta;
			}*/
						// alpha=alpha>=belta?alpha:belta;
						// belta=1-alpha;
						// belta=(budget-priceAll)/budget;
						// alpha=1-belta;
						// belta=belta;
						// alpha=1-belta;
						// belta=-0.01;
						// alpha=0.99;
					}
					if (iter >= MAX_ITER - MAX_ITER / 2 && cost_A > budget)
					{
						alpha = 0;
						belta = 1;
					}
					// double newscore=alpha*(timeAll-time)/timeAll+belta*(priceAll-cost)/priceAll;
					if (alpha != 1)
						temp_score_array[i] = alpha * (timeAll - time) / timeAll + belta * (cost_A - cost_B) / cost_A;
					else
					{

						/*if((iter/4)%2==0)
				{
					a=1;
					b=0;
				}
				else
				{
					a=0;
					b=1;
				}*/
						double a, b;
						if (iter == 0)
						{
							a = 1;
							b = 0;
							optitime = 1;
						}
						else if (optitime == 0)
						{
							a = 1;
							b = 0;
							optitime = 1;
						}
						else if (optitime == 1 && (timeOld - timeAll) / timeOld >= 0.3)
						{
							a = 1;
							b = 0;
							optitime = 1;
						}
						else
						{
							a = 0;
							b = 1;
							optitime = 0;
						}
						temp_score_array[i] = a * (timeAll - time) / timeAll + b * (totalwan - wan) / totalwan;
					}
					//--------------------------测试------------------------------
					/*else
			{
				if(cost>budget)
				temp_score_array[i]=-10000;
				else
				temp_score_array[i]=(timeAll-time)/timeAll;
			}*/

					// temp_score_array[i]=(timeAll-time)/timeAll+F1*((priceAll-cost)/priceAll)-F2*((priceAll-cost)/priceAll);
					/* if(iter>=MAX_ITER/2)
		   temp_score_array[i]=(timeAll-time)/timeAll-((priceAll-cost)/priceAll);
		   else
		   temp_score_array[i]=(timeAll-time)/timeAll+((priceAll-cost)/priceAll);*/

					//----------------------------------将ver移回原来的dc
					// moveVertex(ver,downloadnumG,uploadnumG,downloadnumA,uploadnumA,originL);
				}

				for (int k = 0; k < algorithm->DCnums; k++)
				{
					uploadnumG[k] = ug[k];
					downloadnumG[k] = dg[k];
					uploadnumA[k] = ua[k];
					downloadnumA[k] = da[k];
				}
			}

			/*--------------------------objective function-----------------------------*/
			// int max_index=vec[random_at_most(vec.size()-1)];
			int max_index = max_value_index(temp_score_array, algorithm->DCnums);
			double *signal_ = v[ver].signal_;
			// pthread_mutex_lock(&mutex);
			signal_[max_index] += 1;
			// pthread_mutex_unlock(&mutex);
		}
	}
	// cout<<"Thread "<<thread_id<<" finished!"<<endl;
	delete args;
	delete[] uploadnumsum;
	delete[] uploadnumG;
	delete[] downloadnumG;
	delete[] uploadnumA;
	delete[] downloadnumA;
	delete[] maxArray;
	delete[] maxArray2;
	delete[] temp_score_array;
}
void vertexMigrate(int ver, double *downloadnumG, double *uploadnumG, double *downloadnumA, double *uploadnumA, int destdc) //改变loadnum以及mirrors
{
	long long *mapped = graph->mapped;
	Vertex *v = graph->nodes;
	list<Edge *>::iterator it;
	// int * labels=new int[graph->getNum_Nodes()];
	Mirror **mirrors = graph->mirrors;
	list<Edge *> oedges;
	list<Edge *> iedges;
	int originL = v[ver].label;

	// v[ver].setLabel(destdc);//第一件事
	v[ver].label = destdc;

	if (v[ver].type == 1) //对于low-degree
	{
		//第一步，移点
		for (int i = 0; i < algorithm->DCnums; i++)
			if (mirrors[i][ver].outD > 0)
				uploadnumA[originL] -= dataunit;

		mirrors[originL][ver].inD = v[ver].idegreeinlabel;
		mirrors[originL][ver].outD = v[ver].odegreeinlabel;

		if (mirrors[originL][ver].outD > 0 || mirrors[originL][ver].inD > 0)
			mirrors[originL][ver].id = ver;
		// mirrors[originL][ver].setID(ver);

		if (mirrors[originL][ver].inD > 0)
			uploadnumG[originL] += dataunit;
		if (mirrors[originL][ver].outD > 0)
			downloadnumA[originL] += dataunit;

		if (mirrors[destdc][ver].outD > 0)
			downloadnumA[destdc] -= dataunit;

		if (mirrors[originL][ver].inD > 0)
			downloadnumG[destdc] += dataunit;

		for (int i = 0; i < algorithm->DCnums; i++)
		{
			if (i != destdc && mirrors[i][ver].outD > 0)
				uploadnumA[destdc] += dataunit;
		}

		if (mirrors[destdc][ver].id != -1)
		{
			v[ver].idegreeinlabel = mirrors[destdc][ver].inD;
			v[ver].odegreeinlabel = mirrors[destdc][ver].outD;
		}
		else
		{
			v[ver].idegreeinlabel = 0;
			v[ver].odegreeinlabel = 0;
		}

		//第二步，移动入边
		iedges = v[ver].ingoingEdges;
		for (it = iedges.begin(); it != iedges.end(); it++)
		{
			v[ver].idegreeinlabel++;
			long long sourceId = mapped[(*it)->sourceID];
			if (v[sourceId].type == 1 && v[sourceId].label == originL) // low-degree master
			{
				// mirrors[originL][ver].subInD();
				mirrors[originL][ver].inD--;
				v[sourceId].odegreeinlabel--;
				if (mirrors[originL][ver].inD == 0)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[destdc] -= dataunit;
					if (mirrors[originL][ver].outD == 0)
					{
						mirrors[originL][ver].inD = 0;
						mirrors[originL][ver].outD = 0;
						mirrors[originL][ver].id = -1;
					}
					// mirrors[originL][ver].reset();
				}
				if (mirrors[destdc][sourceId].id == -1)
				{
					mirrors[destdc][sourceId].id = sourceId;
					mirrors[destdc][sourceId].outD++;
					uploadnumA[originL] += dataunit;
					downloadnumA[destdc] += dataunit;
				}
				else
				{
					mirrors[destdc][sourceId].outD++;
				}
			}
			else if (v[sourceId].type == 1 && v[sourceId].label != originL) // low-degree Mirror
			{
				mirrors[originL][ver].inD--;
				mirrors[originL][sourceId].outD--;
				if (mirrors[originL][ver].inD == 0)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[destdc] -= dataunit;
					if (mirrors[originL][ver].outD == 0)
					{
						mirrors[originL][ver].inD = 0;
						mirrors[originL][ver].outD = 0;
						mirrors[originL][ver].id = -1;
					}
				}
				if (mirrors[originL][sourceId].outD == 0)
				{
					downloadnumA[originL] -= dataunit;
					uploadnumA[v[sourceId].label] -= dataunit;
					// mirrors[originL][sourceId].reset();
					mirrors[originL][sourceId].inD = 0;
					mirrors[originL][sourceId].outD = 0;
					mirrors[originL][sourceId].id = -1;
				}
				if (v[sourceId].label != destdc)
				{
					if (mirrors[destdc][sourceId].id == -1)
					{
						mirrors[destdc][sourceId].id = sourceId;
						mirrors[destdc][sourceId].outD++;
						downloadnumA[destdc] += dataunit;
						uploadnumA[v[sourceId].label] += dataunit;
					}
					else
					{
						mirrors[destdc][sourceId].outD++; //当初找了好久的bug
					}
				}
				else
				{
					v[sourceId].odegreeinlabel++;
				}
			}
			else if (v[sourceId].type == 0 && v[sourceId].label == originL) // high-degree master
			{
				mirrors[originL][ver].inD--;
				v[sourceId].odegreeinlabel--;
				if (mirrors[originL][ver].inD == 0)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[destdc] -= dataunit;
					if (mirrors[originL][ver].outD == 0)
					{
						mirrors[originL][ver].inD = 0;
						mirrors[originL][ver].outD = 0;
						mirrors[originL][ver].id = -1;
					}
				}
				if (mirrors[destdc][sourceId].id != -1)
				{
					mirrors[destdc][sourceId].outD++;
					if (mirrors[destdc][sourceId].outD == 1)
					{
						downloadnumA[destdc] += dataunit;
						uploadnumA[originL] += dataunit;
					}
				}
				else
				{
					mirrors[destdc][sourceId].id = sourceId;
					mirrors[destdc][sourceId].outD++;
					downloadnumA[destdc] += dataunit;
					uploadnumA[originL] += dataunit;
				}
			}
			else if (v[sourceId].type == 0 && v[sourceId].label != originL) // high-degree mirror
			{
				mirrors[originL][ver].inD--;
				mirrors[originL][sourceId].outD--;
				if (mirrors[originL][ver].inD == 0)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[destdc] -= dataunit;
					if (mirrors[originL][ver].outD == 0)
					{
						mirrors[originL][ver].inD = 0;
						mirrors[originL][ver].outD = 0;
						mirrors[originL][ver].id = -1;
					}
				}
				if (mirrors[originL][sourceId].outD == 0)
				{
					downloadnumA[originL] -= dataunit;
					uploadnumA[v[sourceId].label] -= dataunit;
					if (mirrors[originL][sourceId].inD == 0)
					{
						mirrors[originL][sourceId].inD = 0;
						mirrors[originL][sourceId].outD = 0;
						mirrors[originL][sourceId].id = -1;
					}
				}
				if (v[sourceId].label != destdc)
				{
					if (mirrors[destdc][sourceId].id == -1)
					{
						mirrors[destdc][sourceId].id = sourceId;
						mirrors[destdc][sourceId].outD++;
						uploadnumA[v[sourceId].label] += dataunit;
						downloadnumA[destdc] += dataunit;
					}
					else
					{
						mirrors[destdc][sourceId].outD++;
						if (mirrors[destdc][sourceId].outD == 1)
						{
							uploadnumA[v[sourceId].label] += dataunit;
							downloadnumA[destdc] += dataunit;
						}
					}
				}
				else
				{
					v[sourceId].odegreeinlabel++;
				}
			}
		}
		//第三步，移动出边------------
		oedges = v[ver].outgoingEdges;
		for (it = oedges.begin(); it != oedges.end(); it++)
		{
			long long destId = mapped[(*it)->destID];
			if (v[destId].type == 0 && v[destId].label == originL) // high-degree master
			{
				v[destId].idegreeinlabel--;
				v[ver].odegreeinlabel++;
				mirrors[originL][ver].outD--;
				if (mirrors[originL][ver].outD == 0)
				{
					mirrors[originL][ver].inD = 0;
					mirrors[originL][ver].outD = 0;
					mirrors[originL][ver].id = -1;
					downloadnumA[originL] -= dataunit;
					uploadnumA[destdc] -= dataunit;
				}
				if (mirrors[destdc][destId].id == -1)
				{
					mirrors[destdc][destId].id = destId;
					mirrors[destdc][destId].inD++;
					uploadnumG[destdc] += dataunit;
					downloadnumG[originL] += dataunit;
				}
				else
				{
					mirrors[destdc][destId].inD++;
					if (mirrors[destdc][destId].inD == 1)
					{
						uploadnumG[destdc] += dataunit;
						downloadnumG[originL] += dataunit;
					}
				}
			}
			else if (v[destId].type == 0 && v[destId].label != originL) // high-degree mirror
			{
				v[ver].odegreeinlabel++;
				mirrors[originL][ver].outD--;
				mirrors[originL][destId].inD--;
				if (mirrors[originL][ver].outD == 0)
				{
					mirrors[originL][ver].inD = 0;
					mirrors[originL][ver].outD = 0;
					mirrors[originL][ver].id = -1;
					downloadnumA[originL] -= dataunit;
					uploadnumA[destdc] -= dataunit;
				}
				if (mirrors[originL][destId].inD == 0)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[v[destId].label] -= dataunit;
					if (mirrors[originL][destId].outD == 0)
					{
						mirrors[originL][destId].inD = 0;
						mirrors[originL][destId].outD = 0;
						mirrors[originL][destId].id = -1;
					}
				}
				if (v[destId].label != destdc && mirrors[destdc][destId].id == -1)
				{
					mirrors[destdc][destId].id = destId;
					mirrors[destdc][destId].inD++;
					uploadnumG[destdc] += dataunit;
					downloadnumG[v[destId].label] += dataunit;
				}
				else if (v[destId].label != destdc && mirrors[destdc][destId].id != -1)
				{
					mirrors[destdc][destId].inD++;
					if (mirrors[destdc][destId].inD == 1)
					{
						uploadnumG[destdc] += dataunit;
						downloadnumG[v[destId].label] += dataunit;
					}
				}
				if (v[destId].label == destdc)
					v[destId].idegreeinlabel++;
			}
		}
		// mirrors[destdc][ver].reset();
		mirrors[destdc][ver].inD = 0;
		mirrors[destdc][ver].outD = 0;
		mirrors[destdc][ver].id = -1;
	}
	else if (v[ver].type == 0) //对于high-degree顶点
	{
		//第一步，移点
		for (int i = 0; i < algorithm->DCnums; i++)
		{
			if (mirrors[i][ver].id != -1 && mirrors[i][ver].inD > 0)
				downloadnumG[originL] -= dataunit;
			if (mirrors[i][ver].id != -1 && mirrors[i][ver].outD > 0)
				uploadnumA[originL] -= dataunit;
		}

		mirrors[originL][ver].inD = v[ver].idegreeinlabel;
		mirrors[originL][ver].outD = v[ver].odegreeinlabel;

		if (mirrors[originL][ver].inD != 0 || mirrors[originL][ver].outD != 0)
			mirrors[originL][ver].id = ver;

		if (mirrors[originL][ver].inD > 0)
			uploadnumG[originL] += dataunit;
		if (mirrors[originL][ver].outD > 0)
			downloadnumA[originL] += dataunit;

		if (mirrors[destdc][ver].inD > 0)
			uploadnumG[destdc] -= dataunit;
		if (mirrors[destdc][ver].outD > 0)
			downloadnumA[destdc] -= dataunit;
		for (int i = 0; i < algorithm->DCnums; i++)
		{
			if (i != destdc && mirrors[i][ver].inD > 0)
				downloadnumG[destdc] += dataunit;
			if (i != destdc && mirrors[i][ver].outD > 0)
				uploadnumA[destdc] += dataunit;
		}

		if (mirrors[destdc][ver].id != -1)
		{
			v[ver].idegreeinlabel = mirrors[destdc][ver].inD;
			v[ver].odegreeinlabel = mirrors[destdc][ver].outD;
		}
		else
		{
			v[ver].idegreeinlabel = 0;
			v[ver].odegreeinlabel = 0;
		}

		//第二步，移动入边
		/*iedges=v[ver].getIngoingEdges();
		for(it=iedges.begin();it!=iedges.end();it++)
		{
			long long sourceId=mapped[(*it)->getsourceID()];
			if(v[sourceId].getType()==1&&v[sourceId].getLabel()==originL)    //low-degree master
			{

			}
			if(v[sourceId].getType()==0&&v[sourceId].getLabel()==originL)    //high-degree master
			{

			}
		}*/

		//第三步，移动出边
		oedges = v[ver].outgoingEdges;
		for (it = oedges.begin(); it != oedges.end(); it++)
		{
			long long destId = mapped[(*it)->destID];
			if (v[destId].type == 0 && v[destId].label == originL) // high-degree master
			{
				v[ver].odegreeinlabel++;
				v[destId].idegreeinlabel--;
				mirrors[originL][ver].outD--;
				if (mirrors[originL][ver].outD == 0)
				{
					downloadnumA[originL] -= dataunit;
					uploadnumA[destdc] -= dataunit;
					if (mirrors[originL][ver].inD == 0)
					{
						mirrors[originL][ver].inD = 0;
						mirrors[originL][ver].outD = 0;
						mirrors[originL][ver].id = -1;
					}
				}
				if (mirrors[destdc][destId].id == -1) //---------
				{
					mirrors[destdc][destId].id = destId;
					mirrors[destdc][destId].inD++;
					uploadnumG[destdc] += dataunit;
					downloadnumG[originL] += dataunit;
				}
				else
				{
					mirrors[destdc][destId].inD++;
					if (mirrors[destdc][destId].inD == 1)
					{
						uploadnumG[destdc] += dataunit;
						downloadnumG[originL] += dataunit;
					}
				}
			}
			else if (v[destId].type == 0 && v[destId].label != originL) // high-degree mirror
			{
				v[ver].odegreeinlabel++;
				mirrors[originL][ver].outD--;
				mirrors[originL][destId].inD--;
				if (mirrors[originL][ver].outD == 0)
				{
					downloadnumA[originL] -= dataunit;
					uploadnumA[destdc] -= dataunit;
					if (mirrors[originL][ver].inD == 0)
					{
						mirrors[originL][ver].inD = 0;
						mirrors[originL][ver].outD = 0;
						mirrors[originL][ver].id = -1;
					}
				}
				if (mirrors[originL][destId].inD == 0)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[v[destId].label] -= dataunit;
					if (mirrors[originL][destId].outD == 0)
					{
						mirrors[originL][destId].inD = 0;
						mirrors[originL][destId].outD = 0;
						mirrors[originL][destId].id = -1;
					}
				}
				if (v[destId].label != destdc && mirrors[destdc][destId].id == -1)
				{
					mirrors[destdc][destId].id = destId;
					mirrors[destdc][destId].inD++;
					uploadnumG[destdc] += dataunit;
					downloadnumG[v[destId].label] += dataunit;
				}
				else if (v[destId].label != destdc && mirrors[destdc][destId].id != -1)
				{
					mirrors[destdc][destId].inD++;
					if (mirrors[destdc][destId].inD == 1)
					{
						uploadnumG[destdc] += dataunit;
						downloadnumG[v[destId].label] += dataunit;
					}
				}
				if (v[destId].label == destdc)
					v[destId].idegreeinlabel++;
			}
		}
		// mirrors[destdc][ver].reset();
		mirrors[destdc][ver].inD = 0;
		mirrors[destdc][ver].outD = 0;
		mirrors[destdc][ver].id = -1;
	}
}
void moveVertex(int ver, double *downloadnumG, double *uploadnumG, double *downloadnumA, double *uploadnumA, int destdc) //只改变loadnum不改变mirrors
{
	long long *mapped = graph->mapped;
	Vertex *v = graph->nodes;
	list<Edge *>::iterator it;
	Mirror **mirrors = graph->mirrors;
	list<Edge *> oedges;
	list<Edge *> iedges;
	int originL = v[ver].label;
	int masterin = 0;
	int mirrorin = 0;
	int masterout = 0;
	int mirrorout = 0;

	// v[ver].setLabel(destdc);//第一件事
	if (v[ver].type == 1) // low-degree
	{
		//第一步，移点
		for (int i = 0; i < algorithm->DCnums; i++)
			if (mirrors[i][ver].outD > 0)
				uploadnumA[originL] -= dataunit;

		masterin = v[ver].idegreeinlabel;
		masterout = v[ver].odegreeinlabel;

		// if(mirrors[originL][ver].getInD()>0)
		if (masterin > 0)
			uploadnumG[originL] += dataunit;
		// if(mirrors[originL][ver].getOutD()>0)
		if (masterout > 0)
			downloadnumA[originL] += dataunit;

		mirrorin = mirrors[destdc][ver].inD;
		mirrorout = mirrors[destdc][ver].outD;

		// if(mirrors[destdc][ver].getOutD()>0)
		if (mirrorout > 0)
			downloadnumA[destdc] -= dataunit;
		/*if(mirrors[originL][ver].getInD()>0)
		  downloadnum[destdc]+=dataunit;*/
		if (masterin > 0)
			downloadnumG[destdc] += dataunit;

		for (int i = 0; i < algorithm->DCnums; i++)
		{
			if (i != destdc && i != originL && mirrors[i][ver].outD > 0)
				uploadnumA[destdc] += dataunit;
		}

		if (masterout > 0)
			uploadnumA[destdc] += dataunit;

		//第二步，移动入边
		iedges = v[ver].ingoingEdges;
		for (it = iedges.begin(); it != iedges.end(); it++)
		{
			long long sourceId = mapped[(*it)->sourceID];
			if (v[sourceId].type == 1 && v[sourceId].label == originL) // low-degree master
			{
				// mirrors[originL][ver].subInD();
				masterin--;
				if (masterin == 0)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[destdc] -= dataunit;
					// if(mirrors[originL][ver].getOutD()==0)
					//  mirrors[originL][ver].reset();
				}
				if (mirrors[destdc][sourceId].id == -1)
				{
					// mirrors[destdc][sourceId].setID(sourceId);
					// mirrors[destdc][sourceId].addOutD();
					uploadnumA[originL] += dataunit;
					downloadnumA[destdc] += dataunit;
				}
				else
				{
					// mirrors[destdc][sourceId].addOutD();
				}
			}
			else if (v[sourceId].type == 1 && v[sourceId].label != originL) // low-degree Mirror
			{
				// mirrors[originL][ver].subInD();
				masterin--;
				// mirrors[originL][sourceId].subOutD();
				if (masterin == 0)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[destdc] -= dataunit;
					// if(mirrors[originL][ver].getOutD()==0)
					// mirrors[originL][ver].reset();
				}
				if (mirrors[originL][sourceId].outD == 1)
				{
					downloadnumA[originL] -= dataunit;
					uploadnumA[v[sourceId].label] -= dataunit;
					// mirrors[originL][sourceId].reset();
				}
				if (v[sourceId].label != destdc)
				{
					if (mirrors[destdc][sourceId].id == -1)
					{
						// irrors[destdc][sourceId].setID(sourceId);
						// mirrors[destdc][sourceId].addOutD();
						downloadnumA[destdc] += dataunit;
						uploadnumA[v[sourceId].label] += dataunit;
					}
				}
			}
			else if (v[sourceId].type == 0 && v[sourceId].label == originL) // high-degree master
			{
				// mirrors[originL][ver].subInD();
				masterin--;
				if (masterin == 0)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[destdc] -= dataunit;
					// if(mirrors[originL][ver].getOutD()==0)
					// mirrors[originL][ver].reset();
				}
				if (mirrors[destdc][sourceId].id != -1)
				{
					// mirrors[destdc][sourceId].addOutD();
					if (mirrors[destdc][sourceId].outD == 0)
					{
						downloadnumA[destdc] += dataunit;
						uploadnumA[originL] += dataunit;
					}
				}
				else
				{
					// mirrors[destdc][sourceId].setID(sourceId);
					// mirrors[destdc][sourceId].addOutD();
					downloadnumA[destdc] += dataunit;
					uploadnumA[originL] += dataunit;
				}
			}
			else if (v[sourceId].type == 0 && v[sourceId].label != originL) // high-degree mirror
			{
				// mirrors[originL][ver].subInD();
				masterin--;
				// mirrors[originL][sourceId].subOutD();
				if (masterin == 0)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[destdc] -= dataunit;
					// if(mirrors[originL][ver].getOutD()==0)
					//  mirrors[originL][ver].reset();
				}
				if (mirrors[originL][sourceId].outD == 1)
				{
					downloadnumA[originL] -= dataunit;
					uploadnumA[v[sourceId].label] -= dataunit;
					// if(mirrors[originL][sourceId].getInD()==0)
					//  mirrors[originL][sourceId].reset();
				}
				if (v[sourceId].label != destdc)
				{
					if (mirrors[destdc][sourceId].id == -1)
					{
						// mirrors[destdc][sourceId].setID(sourceId);
						// mirrors[destdc][sourceId].addOutD();
						uploadnumA[v[sourceId].label] += dataunit;
						downloadnumA[destdc] += dataunit;
					}
					else
					{
						// mirrors[destdc][sourceId].addOutD();
						if (mirrors[destdc][sourceId].outD == 0)
						{
							uploadnumA[v[sourceId].label] += dataunit;
							downloadnumA[destdc] += dataunit;
						}
					}
				}
			}
		}
		//第三步，移动出边
		oedges = v[ver].outgoingEdges;
		for (it = oedges.begin(); it != oedges.end(); it++)
		{
			long long destId = mapped[(*it)->destID];
			if (v[destId].type == 0 && v[destId].label == originL) // high-degree master
			{
				// mirrors[originL][ver].subOutD();
				masterout--;
				if (masterout == 0)
				{
					// mirrors[originL][ver].reset();
					downloadnumA[originL] -= dataunit;
					uploadnumA[destdc] -= dataunit;
				}
				if (mirrors[destdc][destId].id == -1)
				{
					// mirrors[destdc][destId].setID(destId);
					// mirrors[destdc][destId].addInD();
					uploadnumG[destdc] += dataunit;
					downloadnumG[originL] += dataunit;
				}
				else
				{
					// mirrors[destdc][destId].addInD();
					if (mirrors[destdc][destId].inD == 0)
					{
						uploadnumG[destdc] += dataunit;
						downloadnumG[originL] += dataunit;
					}
				}
			}
			else if (v[destId].type == 0 && v[destId].label != originL) // high-degree mirror
			{
				// mirrors[originL][ver].subOutD();
				masterout--;
				// mirrors[originL][destId].subInD();
				if (masterout == 0)
				{
					// mirrors[originL][ver].reset();
					downloadnumA[originL] -= dataunit;
					uploadnumA[destdc] -= dataunit;
				}
				if (mirrors[originL][destId].inD == 1)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[v[destId].label] -= dataunit;
					// if(mirrors[originL][destId].getOutD()==0)
					// mirrors[originL][destId].reset();
				}
				if (v[destId].label != destdc && mirrors[destdc][destId].id == -1)
				{
					// mirrors[destdc][destId].setID(destId);
					// mirrors[destdc][destId].addInD();
					uploadnumG[destdc] += dataunit;
					downloadnumG[v[destId].label] += dataunit;
				}
				else if (v[destId].label != destdc && mirrors[destdc][destId].id != -1)
				{
					// mirrors[destdc][destId].addInD();
					if (mirrors[destdc][destId].inD == 0)
					{
						uploadnumG[destdc] += dataunit;
						downloadnumG[v[destId].label] += dataunit;
					}
				}
			}
		}
		// mirrors[destdc][ver].reset();
	}
	else if (v[ver].type == 0) //对于high-degree顶点
	{
		//第一步，移点
		for (int i = 0; i < algorithm->DCnums; i++)
		{
			if (mirrors[i][ver].id != -1 && mirrors[i][ver].inD > 0)
				downloadnumG[originL] -= dataunit;
			if (mirrors[i][ver].id != -1 && mirrors[i][ver].outD > 0)
				uploadnumA[originL] -= dataunit;
		}
		masterin = v[ver].idegreeinlabel;
		masterout = v[ver].odegreeinlabel;

		if (masterin > 0)
			uploadnumG[originL] += dataunit;
		if (masterout > 0)
			downloadnumA[originL] += dataunit;

		if (mirrors[destdc][ver].inD > 0)
			uploadnumG[destdc] -= dataunit;
		if (mirrors[destdc][ver].outD > 0)
			downloadnumA[destdc] -= dataunit;
		for (int i = 0; i < algorithm->DCnums; i++)
		{
			if (i != destdc && i != originL && mirrors[i][ver].inD > 0)
				downloadnumG[destdc] += dataunit;
			if (i != destdc && i != originL && mirrors[i][ver].outD > 0)
				uploadnumA[destdc] += dataunit;
		}

		if (masterin > 0)
			downloadnumG[destdc] += dataunit;
		if (masterout > 0)
			uploadnumA[destdc] += dataunit;

		//第三步，移动出边
		oedges = v[ver].outgoingEdges;
		for (it = oedges.begin(); it != oedges.end(); it++)
		{
			long long destId = mapped[(*it)->destID];
			if (v[destId].type == 0 && v[destId].label == originL) // high-degree master
			{
				// mirrors[originL][ver].subOutD();
				masterout--;
				if (masterout == 0)
				{
					downloadnumA[originL] -= dataunit;
					uploadnumA[destdc] -= dataunit;
					/* if(mirrors[originL][ver].getInD()==0)
					 mirrors[originL][ver].reset();*/
				}
				if (mirrors[destdc][destId].id == -1)
				{
					// mirrors[destdc][destId].setID(destId);
					// mirrors[destdc][destId].addInD();
					uploadnumG[destdc] += dataunit;
					downloadnumG[originL] += dataunit;
				}
				else
				{
					// mirrors[destdc][destId].addInD();
					if (mirrors[destdc][destId].inD == 0)
					{
						uploadnumG[destdc] += dataunit;
						downloadnumG[originL] += dataunit;
					}
				}
			}
			else if (v[destId].type == 0 && v[destId].label != originL) // high-degree mirror
			{
				// mirrors[originL][ver].subOutD();
				masterout--;
				// mirrors[originL][destId].subInD();
				if (masterout == 0)
				{
					downloadnumA[originL] -= dataunit;
					uploadnumA[destdc] -= dataunit;
					/*if(mirrors[originL][ver].getInD()==0)
					 mirrors[originL][ver].reset();*/
				}
				if (mirrors[originL][destId].inD == 1)
				{
					uploadnumG[originL] -= dataunit;
					downloadnumG[v[destId].label] -= dataunit;
					/*if(mirrors[originL][destId].getOutD()==0)
					 mirrors[originL][destId].reset();*/
				}
				if (v[destId].label != destdc && mirrors[destdc][destId].id == -1)
				{
					// mirrors[destdc][destId].setID(destId);
					// mirrors[destdc][destId].addInD();
					uploadnumG[destdc] += dataunit;
					downloadnumG[v[destId].label] += dataunit;
				}
				else if (v[destId].label != destdc && mirrors[destdc][destId].id != -1)
				{
					// mirrors[destdc][destId].addInD();
					if (mirrors[destdc][destId].inD == 0)
					{
						uploadnumG[destdc] += dataunit;
						downloadnumG[v[destId].label] += dataunit;
					}
				}
			}
		}

		// mirrors[destdc][ver].reset();
	}
}
void *moveVertexBsp(void *arguments)
{
	// pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
	Pthread_args *argss = (Pthread_args *)arguments;
	int thread_id = argss->id;
	delete argss;
	// srand(time(NULL) + thread_id);

	/*	int verindex=args->verindex;
	int ver=args->ver;
	double* downloadnumG=args->downloadnumG;
	double* uploadnumG=args->uploadnumG;
	double* downloadnumA=args->downloadnumA;
	double* uploadnumA=args->uploadnumA;
	int destdc=args->destdc;
	int* s=args->s;*/
	double *Upload = network->upload;
	double *Download = network->download;
	double *Upprice = network->upprice;
	long long *mapped = graph->mapped;
	double uploadnumsum[algorithm->DCnums];
	Vertex *v = graph->nodes;
	list<Edge *>::iterator it;
	Mirror **mirrors = graph->mirrors;
	list<Edge *> oedges;
	list<Edge *> iedges;
	// int originL=v[ver].label;
	int masterin = 0;
	int mirrorin = 0;
	int masterout = 0;
	int mirrorout = 0;
	// double* maxArray=new double[algorithm->DCnums];
	// double* maxArray2=new double[algorithm->DCnums];
	// cout<<"thread: "<<thread_id<<endl;
	while (true)
	{
		// printf("123\n");
		// cout<<"thread: "<<thread_id<<endl;
		pthread_testcancel();
		if (!bspvec[thread_id].empty())
		{
			// cout<<"thread: "<<thread_id<<endl;
			Pthread_args args = bspvec[thread_id][0];
			bspvec[thread_id].pop_back();
			int verindex = args.verindex;
			int ver = args.ver;
			double *downloadnumG = args.downloadnumG;
			double *uploadnumG = args.uploadnumG;
			double *downloadnumA = args.downloadnumA;
			double *uploadnumA = args.uploadnumA;
			int destdc = args.destdc;
			int iter = args.iteration;
			int *s = args.s;
			int originL = v[ver].label;
			double *maxArray = new double[algorithm->DCnums];
			double *maxArray2 = new double[algorithm->DCnums];
			// v[ver].setLabel(destdc);//第一件事
			if (v[ver].type == 1) // low-degree
			{
				//第一步，移点
				for (int i = 0; i < algorithm->DCnums; i++)
					if (mirrors[i][ver].outD > 0)
						uploadnumA[originL] -= dataunit;

				masterin = v[ver].idegreeinlabel;
				masterout = v[ver].odegreeinlabel;

				// if(mirrors[originL][ver].getInD()>0)
				if (masterin > 0)
					uploadnumG[originL] += dataunit;
				// if(mirrors[originL][ver].getOutD()>0)
				if (masterout > 0)
					downloadnumA[originL] += dataunit;

				mirrorin = mirrors[destdc][ver].inD;
				mirrorout = mirrors[destdc][ver].outD;

				// if(mirrors[destdc][ver].getOutD()>0)
				if (mirrorout > 0)
					downloadnumA[destdc] -= dataunit;
				/*if(mirrors[originL][ver].getInD()>0)
		  downloadnum[destdc]+=dataunit;*/
				if (masterin > 0)
					downloadnumG[destdc] += dataunit;

				for (int i = 0; i < algorithm->DCnums; i++)
				{
					if (i != destdc && i != originL && mirrors[i][ver].outD > 0)
						uploadnumA[destdc] += dataunit;
				}

				if (masterout > 0)
					uploadnumA[destdc] += dataunit;

				//第二步，移动入边
				iedges = v[ver].ingoingEdges;
				for (it = iedges.begin(); it != iedges.end(); it++)
				{
					long long sourceId = mapped[(*it)->sourceID];
					if (v[sourceId].type == 1 && v[sourceId].label == originL) // low-degree master
					{
						// mirrors[originL][ver].subInD();
						masterin--;
						if (masterin == 0)
						{
							uploadnumG[originL] -= dataunit;
							downloadnumG[destdc] -= dataunit;
							// if(mirrors[originL][ver].getOutD()==0)
							//  mirrors[originL][ver].reset();
						}
						if (mirrors[destdc][sourceId].id == -1)
						{
							// mirrors[destdc][sourceId].setID(sourceId);
							// mirrors[destdc][sourceId].addOutD();
							uploadnumA[originL] += dataunit;
							downloadnumA[destdc] += dataunit;
						}
						else
						{
							// mirrors[destdc][sourceId].addOutD();
						}
					}
					else if (v[sourceId].type == 1 && v[sourceId].label != originL) // low-degree Mirror
					{
						// mirrors[originL][ver].subInD();
						masterin--;
						// mirrors[originL][sourceId].subOutD();
						if (masterin == 0)
						{
							uploadnumG[originL] -= dataunit;
							downloadnumG[destdc] -= dataunit;
							// if(mirrors[originL][ver].getOutD()==0)
							// mirrors[originL][ver].reset();
						}
						if (mirrors[originL][sourceId].outD == 1)
						{
							downloadnumA[originL] -= dataunit;
							uploadnumA[v[sourceId].label] -= dataunit;
							// mirrors[originL][sourceId].reset();
						}
						if (v[sourceId].label != destdc)
						{
							if (mirrors[destdc][sourceId].id == -1)
							{
								// irrors[destdc][sourceId].setID(sourceId);
								// mirrors[destdc][sourceId].addOutD();
								downloadnumA[destdc] += dataunit;
								uploadnumA[v[sourceId].label] += dataunit;
							}
						}
					}
					else if (v[sourceId].type == 0 && v[sourceId].label == originL) // high-degree master
					{
						// mirrors[originL][ver].subInD();
						masterin--;
						if (masterin == 0)
						{
							uploadnumG[originL] -= dataunit;
							downloadnumG[destdc] -= dataunit;
							// if(mirrors[originL][ver].getOutD()==0)
							// mirrors[originL][ver].reset();
						}
						if (mirrors[destdc][sourceId].id != -1)
						{
							// mirrors[destdc][sourceId].addOutD();
							if (mirrors[destdc][sourceId].outD == 0)
							{
								downloadnumA[destdc] += dataunit;
								uploadnumA[originL] += dataunit;
							}
						}
						else
						{
							// mirrors[destdc][sourceId].setID(sourceId);
							// mirrors[destdc][sourceId].addOutD();
							downloadnumA[destdc] += dataunit;
							uploadnumA[originL] += dataunit;
						}
					}
					else if (v[sourceId].type == 0 && v[sourceId].label != originL) // high-degree mirror
					{
						// mirrors[originL][ver].subInD();
						masterin--;
						// mirrors[originL][sourceId].subOutD();
						if (masterin == 0)
						{
							uploadnumG[originL] -= dataunit;
							downloadnumG[destdc] -= dataunit;
							// if(mirrors[originL][ver].getOutD()==0)
							//  mirrors[originL][ver].reset();
						}
						if (mirrors[originL][sourceId].outD == 1)
						{
							downloadnumA[originL] -= dataunit;
							uploadnumA[v[sourceId].label] -= dataunit;
							// if(mirrors[originL][sourceId].getInD()==0)
							//  mirrors[originL][sourceId].reset();
						}
						if (v[sourceId].label != destdc)
						{
							if (mirrors[destdc][sourceId].id == -1)
							{
								// mirrors[destdc][sourceId].setID(sourceId);
								// mirrors[destdc][sourceId].addOutD();
								uploadnumA[v[sourceId].label] += dataunit;
								downloadnumA[destdc] += dataunit;
							}
							else
							{
								// mirrors[destdc][sourceId].addOutD();
								if (mirrors[destdc][sourceId].outD == 0)
								{
									uploadnumA[v[sourceId].label] += dataunit;
									downloadnumA[destdc] += dataunit;
								}
							}
						}
					}
				}
				//第三步，移动出边
				oedges = v[ver].outgoingEdges;
				for (it = oedges.begin(); it != oedges.end(); it++)
				{
					long long destId = mapped[(*it)->destID];
					if (v[destId].type == 0 && v[destId].label == originL) // high-degree master
					{
						// mirrors[originL][ver].subOutD();
						masterout--;
						if (masterout == 0)
						{
							// mirrors[originL][ver].reset();
							downloadnumA[originL] -= dataunit;
							uploadnumA[destdc] -= dataunit;
						}
						if (mirrors[destdc][destId].id == -1)
						{
							// mirrors[destdc][destId].setID(destId);
							// mirrors[destdc][destId].addInD();
							uploadnumG[destdc] += dataunit;
							downloadnumG[originL] += dataunit;
						}
						else
						{
							// mirrors[destdc][destId].addInD();
							if (mirrors[destdc][destId].inD == 0)
							{
								uploadnumG[destdc] += dataunit;
								downloadnumG[originL] += dataunit;
							}
						}
					}
					else if (v[destId].type == 0 && v[destId].label != originL) // high-degree mirror
					{
						// mirrors[originL][ver].subOutD();
						masterout--;
						// mirrors[originL][destId].subInD();
						if (masterout == 0)
						{
							// mirrors[originL][ver].reset();
							downloadnumA[originL] -= dataunit;
							uploadnumA[destdc] -= dataunit;
						}
						if (mirrors[originL][destId].inD == 1)
						{
							uploadnumG[originL] -= dataunit;
							downloadnumG[v[destId].label] -= dataunit;
							// if(mirrors[originL][destId].getOutD()==0)
							// mirrors[originL][destId].reset();
						}
						if (v[destId].label != destdc && mirrors[destdc][destId].id == -1)
						{
							// mirrors[destdc][destId].setID(destId);
							// mirrors[destdc][destId].addInD();
							uploadnumG[destdc] += dataunit;
							downloadnumG[v[destId].label] += dataunit;
						}
						else if (v[destId].label != destdc && mirrors[destdc][destId].id != -1)
						{
							// mirrors[destdc][destId].addInD();
							if (mirrors[destdc][destId].inD == 0)
							{
								uploadnumG[destdc] += dataunit;
								downloadnumG[v[destId].label] += dataunit;
							}
						}
					}
				}
				// mirrors[destdc][ver].reset();
			}
			else if (v[ver].type == 0) //对于high-degree顶点
			{
				//第一步，移点
				for (int i = 0; i < algorithm->DCnums; i++)
				{
					if (mirrors[i][ver].id != -1 && mirrors[i][ver].inD > 0)
						downloadnumG[originL] -= dataunit;
					if (mirrors[i][ver].id != -1 && mirrors[i][ver].outD > 0)
						uploadnumA[originL] -= dataunit;
				}
				masterin = v[ver].idegreeinlabel;
				masterout = v[ver].odegreeinlabel;

				if (masterin > 0)
					uploadnumG[originL] += dataunit;
				if (masterout > 0)
					downloadnumA[originL] += dataunit;

				if (mirrors[destdc][ver].inD > 0)
					uploadnumG[destdc] -= dataunit;
				if (mirrors[destdc][ver].outD > 0)
					downloadnumA[destdc] -= dataunit;
				for (int i = 0; i < algorithm->DCnums; i++)
				{
					if (i != destdc && i != originL && mirrors[i][ver].inD > 0)
						downloadnumG[destdc] += dataunit;
					if (i != destdc && i != originL && mirrors[i][ver].outD > 0)
						uploadnumA[destdc] += dataunit;
				}

				if (masterin > 0)
					downloadnumG[destdc] += dataunit;
				if (masterout > 0)
					uploadnumA[destdc] += dataunit;

				//第三步，移动出边
				oedges = v[ver].outgoingEdges;
				for (it = oedges.begin(); it != oedges.end(); it++)
				{
					long long destId = mapped[(*it)->destID];
					if (v[destId].type == 0 && v[destId].label == originL) // high-degree master
					{
						// mirrors[originL][ver].subOutD();
						masterout--;
						if (masterout == 0)
						{
							downloadnumA[originL] -= dataunit;
							uploadnumA[destdc] -= dataunit;
							/* if(mirrors[originL][ver].getInD()==0)
					 mirrors[originL][ver].reset();*/
						}
						if (mirrors[destdc][destId].id == -1)
						{
							// mirrors[destdc][destId].setID(destId);
							// mirrors[destdc][destId].addInD();
							uploadnumG[destdc] += dataunit;
							downloadnumG[originL] += dataunit;
						}
						else
						{
							// mirrors[destdc][destId].addInD();
							if (mirrors[destdc][destId].inD == 0)
							{
								uploadnumG[destdc] += dataunit;
								downloadnumG[originL] += dataunit;
							}
						}
					}
					else if (v[destId].type == 0 && v[destId].label != originL) // high-degree mirror
					{
						// mirrors[originL][ver].subOutD();
						masterout--;
						// mirrors[originL][destId].subInD();
						if (masterout == 0)
						{
							downloadnumA[originL] -= dataunit;
							uploadnumA[destdc] -= dataunit;
							/*if(mirrors[originL][ver].getInD()==0)
					 mirrors[originL][ver].reset();*/
						}
						if (mirrors[originL][destId].inD == 1)
						{
							uploadnumG[originL] -= dataunit;
							downloadnumG[v[destId].label] -= dataunit;
							/*if(mirrors[originL][destId].getOutD()==0)
					 mirrors[originL][destId].reset();*/
						}
						if (v[destId].label != destdc && mirrors[destdc][destId].id == -1)
						{
							// mirrors[destdc][destId].setID(destId);
							// mirrors[destdc][destId].addInD();
							uploadnumG[destdc] += dataunit;
							downloadnumG[v[destId].label] += dataunit;
						}
						else if (v[destId].label != destdc && mirrors[destdc][destId].id != -1)
						{
							// mirrors[destdc][destId].addInD();
							if (mirrors[destdc][destId].inD == 0)
							{
								uploadnumG[destdc] += dataunit;
								downloadnumG[v[destId].label] += dataunit;
							}
						}
					}
				}

				// mirrors[destdc][ver].reset();
			}
			/************************Gather阶段Time***********************************/

			for (int j = 0; j < algorithm->DCnums; j++)
			{
				double divide = uploadnumG[j] / Upload[j];
				double divided = downloadnumG[j] / Download[j];
				maxArray[j] = divide > divided ? divide : divided;

				double divide2 = uploadnumA[j] / Upload[j];
				double divided2 = downloadnumA[j] / Download[j];
				maxArray2[j] = divide2 > divided2 ? divide2 : divided2;

				uploadnumsum[j] = uploadnumG[j] + uploadnumA[j];
			}

			double time = max_value(maxArray) + max_value(maxArray2);
			double cost = sumWeights(uploadnumsum, Upprice);

			double wan = 0;
			double totalwan = 0;
			for (int j = 0; j < algorithm->DCnums; j++)
			{
				wan = wan + uploadnumG[j] + downloadnumG[j] + uploadnumA[j] + downloadnumA[j];
				totalwan = totalwan + ug[j] + dg[j] + ua[j] + da[j];
			}
			/************************Apply阶段Time***********************************/
			/*   for (int j = 0; j < algorithm->DCnums; j++)
			{
				double divide=uploadnumA[j] / Upload[j];
				double divided=downloadnumA[j] / Download[j];
				maxArray[j] =divide > divided ? divide : divided;
			}
			time+=max_value(maxArray);*/
			/*int overbudget=0;
			if(priceAll>budget)
			  overbudget=1;
			double downbudget=overbudget==1?0:1;
			double costdown=cost>budget?cost:budget;
			//double over=cost>budget?1:0;
			double score=(timeAll-time)/timeAll+4*overbudget*(priceAll-costdown)/priceAll;
			if(score>=0)
			s[verindex]=1;*/
			double costdown = cost > budget ? cost : budget;
			double score = (timeAll - time) / timeAll + (priceAll - costdown) / priceAll;
			double score2 = (timeAll - time) / timeAll - (priceAll - cost) / priceAll;

			/*	 if(betaC<=0.5&&time<=timeAll)
			s[verindex]=1;
		else if(betaC>0.5)
		{
			 if(priceAll<=budget&&cost<=budget&&time<=timeAll)
			 {
				 s[verindex]=1;
			 }
			 else if(priceAll>budget&&cost<priceAll)
				 s[verindex]=1;
		}*/
			// double belta=-(budget-cost)/budget;
			double rr = 0;
			if (v[ver].iniLabel == v[ver].label && v[ver].label != v[ver].action)
				rr = raw_data_size * Upprice[v[ver].iniLabel];
			else if (v[ver].iniLabel != v[ver].label && v[ver].label != v[ver].action && v[ver].iniLabel == v[ver].action)
				rr = -raw_data_size * Upprice[v[ver].iniLabel];

			double cost_A = priceAll * mvpercents + current_mvcost;
			double cost_B = cost * mvpercents + current_mvcost + rr;

			double belta = (cost_A - budget) / cost_A;
			double alpha = 1 - belta;
			if (belta < 0)
			{
				belta = 0;
				alpha = 1;
			}
			// double firstw,secondw;

			if (alpha == 1)
			{
			}
			if ((iter >= MAX_ITER - (MAX_ITER / 2 - 1) && cost_A > budget) || (iter == MAX_ITER && ((budget - cost_A) / budget) < 1e-5))
			{
				alpha = 0;
				belta = 1;
			}
			double newscore = 0;
			if (alpha != 1)
				newscore = alpha * (timeAll - time) / timeAll + belta * (cost_A - cost_B) / cost_A;
			else
			{
				double b = ((iter - 1) / 4) * 0.2;
				double a = 1 - b;
				if (optitime == 1)
				{
					a = 1;
					b = 0;
				}
				else
				{
					a = 0;
					b = 1;
				}
				newscore = a * (timeAll - time) / timeAll + b * (totalwan - wan) / totalwan;
			}

			if (newscore >= 0)
				s[verindex] = 1;

			// else if(betaC>0.75&&time<=timeAll&&((cost<=budget)||(priceAll>budget&&cost<priceAll)))
			//  s[verindex]=1;

			delete[] maxArray;
			delete[] maxArray2;
			pthread_mutex_lock(&mutex);
			bsplock++;
			pthread_mutex_unlock(&mutex);
		}
	}
	// delete args;

	// cout<<"Thread: "<<thread_id<<" finished!"<<endl;
}
void *Parallel_Score_Migration_function(void *arguments)
{
	/*Pthread_args * args=(Pthread_args* )arguments;
	unsigned int low=args->getLow();
	unsigned int high=args->getHigh();
	double* uploadnum = new double[algorithm->getDCNums()];
	double* downloadnum = new double[algorithm->getDCNums()];
	long long *mapped = graph->getMapped();
	double *maxArray = new double[algorithm->getDCNums()];
	double* Upload = network->getUpload();
	double* Download = network->getDownload();
	double* Upprice = network->getUpprice();
	Vertex* v = graph->getVertexs();
	list<Edge*>::iterator it;
	list<int>::iterator iter;
	int * labels=new int[graph->getNum_Nodes()];
	list<int>* mirrorsbk=new list<int>[algorithm->getDCNums()];
	list<int>* mirrors=graph->getMirrors();
	list<Edge*>oedges;
	list<Edge*>iedges;
	double *temp_score_array=new double[algorithm->getDCNums()];
	double *migration=new double[algorithm->getDCNums()];
	//cout<<high<<endl;
	for(int ver=low;ver<=high;ver++)
	{
		cout<<ver<<endl;
	 for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		if (i == v[ver].getLabel())
		{
			if(priceAll<=budget)
			 temp_score_array[i] = 0;
			 else
			 temp_score_array[i]=-100;
		}
		else
		{
			for(int j=0;j<algorithm->getDCNums();j++)
			   {
				   uploadnum[j]=0;
				   downloadnum[j]=0;
			   }
		  /*  for(int j=0;j<graph->getNum_Nodes();j++)
			   labels[j]=v[j].getLabel();
			labels[ver]=i;*/
	/*	for(int j=0;j<algorithm->getDCNums();j++)
			  mirrorsbk[j]=mirrors[j];
		 /*  mirrorsbk[v[ver].getLabel()].clear();
			mirrorsbk[i].clear();

		for(int j=0;j < graph->getNum_Nodes(); j++)
		{
			if((labels[j]==i||labels[j]==v[ver].getLabel())&&v[j].getType()==1)
			{
			   iedges=v[j].getIngoingEdges();
			for (it = iedges.begin(); it != iedges.end(); it++)
			{
				if(labels[mapped[(*it)->getsourceID()]]!=labels[j])
				  mirrorsbk[labels[j]].push_back(mapped[(*it)->getsourceID()]);
			}
			 oedges=v[j].getOutgoingEdges();
			for (it = oedges.begin(); it != oedges.end(); it++)
			{
				if(labels[mapped[(*it)->getdestID()]]!=labels[j]&&v[mapped[(*it)->getdestID()]].getType()==0)
				  mirrorsbk[labels[j]].push_back(mapped[(*it)->getdestID()]);
			}

			}
			else if((labels[j]==i||labels[j]==v[ver].getLabel())&&v[j].getType()==0)
			{
				oedges=v[j].getOutgoingEdges();
				for (it = oedges.begin(); it != oedges.end(); it++)
			{
				if(labels[mapped[(*it)->getdestID()]]!=i&&v[mapped[(*it)->getdestID()]].getType()==0)
				  mirrorsbk[i].push_back(mapped[(*it)->getdestID()]);
			}

			}
		}
		mirrorsbk[i].unique();
		mirrorsbk[v[ver].getLabel()].unique();*/
	//-------------------------------Gather Stage----------------------------------
	/*  for(int j=0;j<algorithm->getDCNums();j++)
	{
	for (iter = mirrorsbk[j].begin(); iter != mirrorsbk[j].end(); iter++)
	   if(v[(*iter)].getType()==0)
	   {
		   uploadnum[j]+=dataunit;
		   downloadnum[v[(*iter)].getLabel()]+=dataunit;
	   }
	}
	 //--------------------------------Scatter Stage----------------------------------
	for(int j=0;j<algorithm->getDCNums();j++)
	{
	 for (iter = mirrorsbk[j].begin(); iter != mirrorsbk[j].end(); iter++)
	{
		   downloadnum[j]+=dataunit;
		   uploadnum[v[(*iter)].getLabel()]+=dataunit;
	 }
	}


			for (int j = 0; j < algorithm->getDCNums(); j++)
			{
				maxArray[j] =uploadnum[j] / Upload[j] > downloadnum[j] / Download[j] ? uploadnum[j] / Upload[j] : downloadnum[j] / Download[j];
			}
			double time = max_value(maxArray);
			migration[i] = sumWeights(uploadnum, Upprice);
			temp_score_array[i] = timeAll - time;
			if(migration[i]<=budget||migration[i]<priceAll)
			  temp_score_array[i]=2*temp_score_array[i];
			else
			  temp_score_array[i]=-100;
		}
	}
	   double *signal_ = NULL;
	   int max_index = max_value_index(temp_score_array,algorithm->getDCNums());
	   list<Edge*>e = v[ver].getIngoingEdges();
		for (it = e.begin(); it != e.end(); it++)
		{

					signal_ = v[mapped[(*it)->getsourceID()]].getSignal_();
					if(v[mapped[(*it)->getsourceID()]].getAction()==max_index)
					   pthread_mutex_lock(&mutex);
						 signal_[max_index]+=1;
					   pthread_mutex_unlock(&mutex);

		}
		  signal_ = v[ver].getSignal_();
		  pthread_mutex_lock(&mutex);
		  signal_[max_index]+=1;
		  pthread_mutex_unlock(&mutex);

 }
	delete[] labels;
	delete[] mirrorsbk;
	delete[] uploadnum;
	delete[] downloadnum;
	delete[] maxArray;
	delete[] migration;
	delete[] temp_score_array;
	//cout<<"222222222222222222222222222222"<<endl;*/
}
void Score_Migration_function(double *temp_score_array, double *migration, int index, double *temp_uploadnum, double *temp_downloadnum) //Ϊÿһ��vertex��ÿһ��partition��֣�����ڸö���������ڵķ�������ֻ������time
{
	/*	double* uploadnum = new double[algorithm->getDCNums()];
	double* downloadnum = new double[algorithm->getDCNums()];
	Partition* partition = network->getPartition();
	long long *mapped = graph->getMapped();
	double *maxArray = new double[algorithm->getDCNums()];
	double* Upload = network->getUpload();
	double* Download = network->getDownload();
	double* Upprice = network->getUpprice();
	Vertex* v = graph->getVertexs();
	list<Edge*>::iterator it;
	list<int>::iterator iter;
	int * labels=new int[graph->getNum_Nodes()];
	list<int>* mirrorsbk=new list<int>[algorithm->getDCNums()];
	list<int>* mirrors=graph->getMirrors();
	list<Edge*>oedges;
	list<Edge*>iedges;

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{

		if (i == v[index].getLabel())
		{
			if(priceAll<=budget)
			 temp_score_array[i] = 0;
			 else
			 temp_score_array[i]=-100;
			migration[i] = 0;
			for (int j = 0; j < algorithm->getDCNums(); j++)
			   {
								uploadnum[j] = partition[j].getUploadNum();
								downloadnum[j] = partition[j].getDownloadNum();
							  //  cout<<"calcu: "<<"partition "<<j<<":"<<uploadnum[j]<<"   "<<downloadnum[j]<<endl;
				}
		}
		else
		{
			for(int j=0;j<algorithm->getDCNums();j++)
			   {
				   uploadnum[j]=0;
				   downloadnum[j]=0;
			   }
			for(int j=0;j<graph->getNum_Nodes();j++)
			   labels[j]=v[j].getLabel();
			labels[index]=i;
			for(int j=0;j<algorithm->getDCNums();j++)
			  mirrorsbk[j]=mirrors[j];
			mirrorsbk[v[index].getLabel()].clear();
			mirrorsbk[i].clear();

		for(int j=0;j < graph->getNum_Nodes(); j++)
		{
			if((labels[j]==i||labels[j]==v[index].getLabel())&&v[j].getType()==1)
			{
			   iedges=v[j].getIngoingEdges();
			for (it = iedges.begin(); it != iedges.end(); it++)
			{
				if(labels[mapped[(*it)->getsourceID()]]!=labels[j])
				  mirrorsbk[labels[j]].push_back(mapped[(*it)->getsourceID()]);
			}
			 oedges=v[j].getOutgoingEdges();
			for (it = oedges.begin(); it != oedges.end(); it++)
			{
				if(labels[mapped[(*it)->getdestID()]]!=labels[j]&&v[mapped[(*it)->getdestID()]].getType()==0)
				  mirrorsbk[labels[j]].push_back(mapped[(*it)->getdestID()]);
			}

			}
		}
		mirrorsbk[i].unique();
		mirrorsbk[v[index].getLabel()].unique();
		 /*--------------------------------Gather Stage----------------------------------*/
	/* for(int j=0;j<algorithm->getDCNums();j++)
	{
	for (iter = mirrorsbk[j].begin(); iter != mirrorsbk[j].end(); iter++)
	   if(v[(*iter)].getType()==0)
	   {
		   uploadnum[j]+=dataunit;
		   downloadnum[v[(*iter)].getLabel()]+=dataunit;
	   }
	}
	 /*--------------------------------Scatter Stage----------------------------------*/
	/*	for(int j=0;j<algorithm->getDCNums();j++)
	{
	 for (iter = mirrorsbk[j].begin(); iter != mirrorsbk[j].end(); iter++)
	{
		   downloadnum[j]+=dataunit;
		   uploadnum[v[(*iter)].getLabel()]+=dataunit;
	 }
	}

	/*--------------------------------------------------------------------------------*/

	/*	for (int j = 0; j < algorithm->getDCNums(); j++)
			{
				maxArray[j] =uploadnum[j] / Upload[j] > downloadnum[j] / Download[j] ? uploadnum[j] / Upload[j] : downloadnum[j] / Download[j];
			}
			double time = max_value(maxArray);
			migration[i] = sumWeights(uploadnum, Upprice);
			temp_score_array[i] = timeAll - time;
			if(migration[i]<=budget||migration[i]<priceAll)
			  temp_score_array[i]=2*temp_score_array[i];
			else
			  temp_score_array[i]=-100;
		}

  }
	delete[] labels;
	delete[] mirrorsbk;
	delete[] uploadnum;
	delete[] downloadnum;
	delete[] maxArray;*/
}
void calculateTimeAndPrice(double &timeAll, double &priceAll) //������ʱ����ܻ���
{
	// cout<<"helllo"<<endl;
	Vertex *v = graph->getVertexs();
	list<int>::iterator it;
	Mirror **mirrors = graph->getMirrors();
	long long *mapped = graph->getMapped();
	double *Upload = network->getUpload();
	double *Download = network->getDownload();
	double *Upprice = network->getUpprice();
	double *Uptime = new double[algorithm->getDCNums()];   //ÿ������upload����
	double *Downtime = new double[algorithm->getDCNums()]; //ÿ������download����

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		Uptime[i] = 0;
		Downtime[i] = 0;
	}
	// cout<<"hhhhhhhhhhhhhhhhhhhhhhhhhh"<<endl;
	/*--------------------------------Gather Stage----------------------------------*/
	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		for (int j = 0; j < graph->getNum_Nodes(); j++)
			if (mirrors[i][j].getID() != -1 && v[j].getType() == 0 && mirrors[i][j].getInD() != 0)
			{
				Uptime[i] += dataunit;
				Downtime[v[j].getLabel()] += dataunit;
			}
	}

	for (int w = 0; w < algorithm->DCnums; w++)
	{
		ug[w] = Uptime[w];
		dg[w] = Downtime[w];
		// cout<<Uptime[w]<<" "<<Downtime[w]<<" " <<network->uploadnumG[w]<<" "<<network->downloadnumG[w]<<endl;
	}

	double *UptimeS = new double[algorithm->getDCNums()];
	double *DowntimeS = new double[algorithm->getDCNums()];
	double *maxArray = new double[algorithm->getDCNums()];

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{

		maxArray[i] = Uptime[i] / Upload[i] > Downtime[i] / Download[i] ? Uptime[i] / Upload[i] : Downtime[i] / Download[i];
		UptimeS[i] = 0;
		DowntimeS[i] = 0;
	}

	timeAll = max_value(maxArray);

	// cout<<"ggggggggggggggggggggggggggggggggg"<<endl;
	/*--------------------------------Apply Stage----------------------------------*/
	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		for (int j = 0; j < graph->getNum_Nodes(); j++)
		{
			if (mirrors[i][j].getID() != -1 && mirrors[i][j].getOutD() != 0)
			{
				Downtime[i] += dataunit;
				Uptime[v[j].getLabel()] += dataunit;
				DowntimeS[i] += dataunit;
				UptimeS[v[j].getLabel()] += dataunit;
			}
		}
	}
	// cout<<"2222222222222222222222"<<endl;
	for (int w = 0; w < algorithm->DCnums; w++)
	{
		ua[w] = UptimeS[w];
		da[w] = DowntimeS[w];
	}
	/*--------------------------------------------------------------------------------*/
	priceAll = sumWeights(Uptime, Upprice);

	/*for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		network->uploadnumG[i]=Uptime[i]-UptimeS[i];
		network->downloadnumG[i]=Downtime[i]-DowntimeS[i];
		network->uploadnumA[i]=UptimeS[i];
		network->downloadnumA[i]=DowntimeS[i];
		//cout<<ug[i]<<" "<<dg[i]<<" "<<ua[i]<<" "<<da[i]<<endl;
		//cout<<Uptime[i]-UptimeS[i]<<" "<<Downtime[i]-DowntimeS[i]<<" "<<UptimeS[i]<<" "<<DowntimeS[i]<<endl;
	}*/
	Partition *partitions = new Partition[algorithm->getDCNums()];

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		partitions[i].setUploadNum(Uptime[i]);
		partitions[i].setDownloadNum(Downtime[i]);
		maxArray[i] = UptimeS[i] / Upload[i] > DowntimeS[i] / Download[i] ? UptimeS[i] / Upload[i] : DowntimeS[i] / Download[i];
	}

	network->setPartition(partitions);
	timeAll += max_value(maxArray);

	delete[] Uptime;
	delete[] Downtime;
	delete[] maxArray;
	delete[] UptimeS;
	delete[] DowntimeS;
	// cout<<"111111111111111111111111"<<endl;
}
void powerlyrareCalTimeAndCost(double &timeAll, double &priceAll)
{
	Vertex *v = graph->getVertexs();
	list<int>::iterator it;
	Mirror **mirrors = graph->getMirrors();
	long long *mapped = graph->getMapped();
	double *Upload = network->getUpload();
	double *Download = network->getDownload();
	double *Upprice = network->getUpprice();
	double *maxArray = new double[algorithm->getDCNums()];
	double *Uptime = new double[algorithm->getDCNums()];	//ÿ������upload����
	double *Downtime = new double[algorithm->getDCNums()];	//ÿ������download����
	double *UptimeS = new double[algorithm->getDCNums()];	//ÿ������upload����
	double *DowntimeS = new double[algorithm->getDCNums()]; //ÿ������download����

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		Uptime[i] = 0;
		Downtime[i] = 0;
		UptimeS[i] = 0;
		DowntimeS[i] = 0;
	}
	/*--------------------------------Gather Stage----------------------------------*/
	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		for (int j = 0; j < graph->getNum_Nodes(); j++)
			if (mirrors[i][j].getID() != -1 && v[j].getType() == 0)
			{
				Uptime[i] += dataunit;
				Downtime[v[j].getLabel()] += dataunit;
			}
	}

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		maxArray[i] = Uptime[i] / Upload[i] > Downtime[i] / Download[i] ? Uptime[i] / Upload[i] : Downtime[i] / Download[i];
	}

	timeAll = max_value(maxArray);

	/*--------------------------------Apply Stage----------------------------------*/
	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		for (int j = 0; j < graph->getNum_Nodes(); j++)
		{
			if (mirrors[i][j].getID() != -1)
			{
				Downtime[i] += dataunit;
				Uptime[v[j].getLabel()] += dataunit;
				DowntimeS[i] += dataunit;
				UptimeS[v[j].getLabel()] += dataunit;
			}
		}
	}

	/*--------------------------------------------------------------------------------*/
	priceAll = sumWeights(Uptime, Upprice);

	/*for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		//cout << Uptime[i] << "and"<<Downtime[i] << endl;
	}*/

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		maxArray[i] = UptimeS[i] / Upload[i] > DowntimeS[i] / Download[i] ? UptimeS[i] / Upload[i] : DowntimeS[i] / Download[i];
	}

	timeAll += max_value(maxArray);

	delete[] Uptime;
	delete[] Downtime;
	delete[] UptimeS;
	delete[] DowntimeS;
	delete[] maxArray;
}
void reCalTimeAndCost(double &timeAll, double &priceAll)
{
	Vertex *v = graph->getVertexs();
	list<int>::iterator it;
	Mirror **mirrors = graph->getMirrors();
	long long *mapped = graph->getMapped();
	double *Upload = network->getUpload();
	double *Download = network->getDownload();
	double *Upprice = network->getUpprice();
	double *maxArray = new double[algorithm->getDCNums()];
	double *Uptime = new double[algorithm->getDCNums()];	//ÿ������upload����
	double *Downtime = new double[algorithm->getDCNums()];	//ÿ������download����
	double *UptimeS = new double[algorithm->getDCNums()];	//ÿ������upload����
	double *DowntimeS = new double[algorithm->getDCNums()]; //ÿ������download����

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		Uptime[i] = 0;
		Downtime[i] = 0;
		UptimeS[i] = 0;
		DowntimeS[i] = 0;
	}
	/*--------------------------------Gather Stage----------------------------------*/
	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		for (int j = 0; j < graph->getNum_Nodes(); j++)
			if (mirrors[i][j].getID() != -1 && v[j].getType() == 0 && mirrors[i][j].getInD() != 0)
			{
				Uptime[i] += dataunit;
				Downtime[v[j].getLabel()] += dataunit;
			}
	}

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		maxArray[i] = Uptime[i] / Upload[i] > Downtime[i] / Download[i] ? Uptime[i] / Upload[i] : Downtime[i] / Download[i];
	}

	timeAll = max_value(maxArray);

	/*--------------------------------Scatter Stage----------------------------------*/
	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		for (int j = 0; j < graph->getNum_Nodes(); j++)
		{
			if (mirrors[i][j].getID() != -1 && mirrors[i][j].getOutD() != 0)
			{
				Downtime[i] += dataunit;
				Uptime[v[j].getLabel()] += dataunit;
				DowntimeS[i] += dataunit;
				UptimeS[v[j].getLabel()] += dataunit;
			}
		}
	}

	/*--------------------------------------------------------------------------------*/
	priceAll = sumWeights(Uptime, Upprice);

	/*for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		//cout << Uptime[i] << "and"<<Downtime[i] << endl;
	}*/

	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		maxArray[i] = UptimeS[i] / Upload[i] > DowntimeS[i] / Download[i] ? UptimeS[i] / Upload[i] : DowntimeS[i] / Download[i];
	}

	timeAll += max_value(maxArray);

	delete[] Uptime;
	delete[] Downtime;
	delete[] UptimeS;
	delete[] DowntimeS;
	delete[] maxArray;
}
void objective_function(double *score, int index)
{
	list<Edge *>::iterator it;
	long long *mapped = graph->getMapped();
	double *signal_ = NULL;
	Vertex *v = graph->getVertexs();
	int max_index = max_value_index(score, algorithm->getDCNums());
	list<Edge *> ei = v[index].getIngoingEdges();
	list<Edge *> eo = v[index].getOutgoingEdges();
	if (v[index].getType() == 1) // low-degree
	{
		for (it = ei.begin(); it != ei.end(); it++) //处理入边
		{
			if (v[mapped[(*it)->getsourceID()]].getType() == 1) //处理low-degree
			{
			}
		}
	}
	/*
	   list<Edge*>e = v[index].getIngoingEdges();
		for (it = e.begin(); it != e.end(); it++)
		{
					signal_ = v[mapped[(*it)->getsourceID()]].getSignal_();

					if(v[mapped[(*it)->getsourceID()]].getAction()==max_index)
						 signal_[max_index]+=1;

		}

		  signal_ = v[index].getSignal_();
		  signal_[max_index]+=1;*/
}
double Benifit(int index, int label)
{
	Partition *partition = network->getPartition();
	Vertex *v = graph->getVertexs();
	long long *mapped = graph->getMapped();
	double *Upprice = network->getUpprice();
	double *uploadnum = new double[algorithm->getDCNums()];
	double *downloadnum = new double[algorithm->getDCNums()];
	list<Edge *>::iterator it;
	if (label == v[index].getLabel())
	{
		delete[] uploadnum;
		delete[] downloadnum;
		return priceAll;
	}
	else
	{
		for (int j = 0; j < algorithm->getDCNums(); j++)
		{
			uploadnum[j] = partition[j].getUploadNum();
			downloadnum[j] = partition[j].getDownloadNum();
		}

		list<Edge *> outEdges = v[index].getOutgoingEdges();
		for (it = outEdges.begin(); it != outEdges.end(); it++)
		{
			if (v[mapped[(*it)->getdestID()]].getLabel() == v[index].getLabel())
			{
				downloadnum[v[index].getLabel()] += dataunit;
				uploadnum[label] += dataunit;
			}
			else if (v[mapped[(*it)->getdestID()]].getLabel() == label)
			{
				uploadnum[v[index].getLabel()] -= dataunit;
				downloadnum[label] -= dataunit;
			}
			else
			{
				uploadnum[v[index].getLabel()] -= dataunit;
				uploadnum[label] += dataunit;
			}
		}
		list<Edge *> inEdges = v[index].getIngoingEdges();
		for (it = inEdges.begin(); it != inEdges.end(); it++)
		{
			if (v[mapped[(*it)->getsourceID()]].getLabel() == v[index].getLabel())
			{
				uploadnum[v[index].getLabel()] += dataunit;
				downloadnum[label] += dataunit;
			}
			else if (v[mapped[(*it)->getsourceID()]].getLabel() == label)
			{
				downloadnum[v[index].getLabel()] -= dataunit;
				uploadnum[label] -= dataunit;
			}
			else
			{
				downloadnum[v[index].getLabel()] -= dataunit;
				downloadnum[label] += dataunit;
			}
		}
	}
	double newprice = sumWeights(uploadnum, Upprice);
	delete[] uploadnum;
	delete[] downloadnum;
	return newprice;
}
void weighted_probability_update(int iter)
{
	double alpha = 1; // LA reward signal
	double beta = .1;
	int num = algorithm->getDCNums();
	// cout << num << endl;
	unsigned int i = 0;
	unsigned int j = 0;
	int k = 0;
	int k_ = 0;
	int kk = 0;
	int kkk = 0;
	double sum_signal = 0;
	int positive_num = 0;
	int negative_num = 0;
	double seprator = 0;
	double error_sum = 0;
	double threshold = (double)1 / (num * sqrt(num));
	double w0 = 0.9;
	double w1 = 0.4;
	double weighgt = ((w0 - w1) * iter * sqrt(num)) / MAX_ITER;
	int weight_index = 0;
	int idx = -1;
	double *alpha_ = new double[num];
	double *beta_ = new double[num];
	double *sig_ = new double[num];
	int *indices = new int[num];
	char *sign = new char[num + 1];
	double *temp_prob = new double[num];
	Vertex *v = graph->getVertexs();
	double **probability = algorithm->getProbability();
	double *signal_ = NULL;

	for (i = 0; i < num; i++)
	{
		alpha_[i] = 0;
		beta_[i] = 0;
		sig_[i] = 0;
		indices[i] = 0;
		sign[i] = 0;
		temp_prob[i] = 0;
	}
	sign[num] = '\0';
	for (i = 0; i < graph->getNum_Nodes(); i++)
	{
		// cout << i << endl;
		signal_ = v[i].getSignal_();
		positive_num = 0;
		negative_num = 0;
		seprator = sum_value(v[i].getSignal_(), num) / num;
		weight_index = max_value_index(v[i].getSignal_(), num);
		signal_[weight_index] = signal_[weight_index] + signal_[weight_index] * weighgt;
		for (k = 0; k < num; k++)
		{
			if (signal_[k] >= seprator)
			{
				sign[k] = 0;
				positive_num++;
			}
			else
			{
				sign[k] = 1;
				negative_num++;
			}
		}
		if (i == idx)
		{
			for (int w = 0; w < num; w++)
			{
				temp_prob[w] = probability[i][w];
			}
		}
		double *positive_part = new double[positive_num];
		int *positive_ind = new int[positive_num];

		for (int w = 0; w < positive_num; w++)
		{
			positive_part[w] = 0;
			positive_ind[w] = 0;
		}
		double *negative_part = new double[negative_num];
		int *negative_ind = new int[negative_num];
		for (int w = 0; w < negative_num; w++)
		{
			negative_part[w] = 0;
			negative_ind[w] = 0;
		}
		positive_num = 0;
		negative_num = 0;
		for (k = 0; k < num; k++)
		{
			if (signal_[k] >= seprator)
			{
				positive_part[positive_num] = signal_[k];
				positive_ind[positive_num] = k;
				positive_num++;
			}
			else
			{
				negative_part[negative_num] = signal_[k];
				negative_ind[negative_num] = k;
				negative_num++;
			}
		}

		sum_signal = sum_value(positive_part, positive_num);
		if (sum_signal > 0)
		{
			for (k = 0; k < positive_num; k++)
			{
				if (positive_part[k])
					positive_part[k] = positive_part[k] / sum_signal;
			}
		}
		for (k = 0; k < negative_num; k++)
		{
			if (negative_part[k] < 0)
				negative_part[k] = -negative_part[k];
		}

		sum_signal = sum_value(negative_part, negative_num);
		if (sum_signal > 0)
		{
			for (k = 0; k < negative_num; k++)
				negative_part[k] = negative_part[k] / sum_signal;
		}
		else
		{
			for (k = 0; k < negative_num; k++)
			{
				negative_part[k] = (double)1 / negative_num;
			}
		}
		positive_num = 0;
		negative_num = 0;
		for (k = 0; k < num; k++)
		{
			if (!sign[k])
			{
				sig_[k] = positive_part[positive_num];
				positive_num++;
			}
			else
			{
				sig_[k] = negative_part[negative_num];
				negative_num++;
			}
		}

		bubble_sort(positive_part, positive_ind, positive_num);
		bubble_sort(negative_part, negative_ind, negative_num);

		for (k = 0; k < negative_num; k++)
		{
			sig_[negative_ind[k]] = 0;
		}
		memcpy(indices, negative_ind, sizeof(negative_ind));
		memcpy(indices + negative_num, positive_ind, sizeof(positive_ind));

		for (k = 0; k < num; k++)
		{
			kk = indices[k];
			if (!sign[kk] && sig_[kk])
			{
				alpha_[kk] = sig_[kk] * alpha;
				probability[i][kk] = probability[i][kk] + alpha_[kk] * (1 - probability[i][kk]);
				for (k_ = 0; k_ < num; k_++)
				{
					if (k_ != kk)
						probability[i][k_] = (1 - alpha_[kk]) * probability[i][k_];
				}
			}
			else if (sign[kk] && sig_[kk])
			{
				beta_[kk] = sig_[kk] * beta;
				probability[i][kk] = (1 - beta_[kk]) * probability[i][kk];
				for (k_ = 0; k_ < num; k_++)
				{
					if (k_ != kk)
						probability[i][k_] = (beta_[kk] / (num - 1)) + (1 - beta_[kk]) * probability[i][k_];
				}
			}
		}

		kk = num - 1;
		for (k = num - 1; k >= 0; k--)
		{
			kkk = indices[k];
			if (probability[i][kkk] < threshold)
			{
				if (probability[i][kkk] > 0)
				{
					probability[i][kk] += probability[i][kkk];
					probability[i][kkk] = 0;
					error_sum = sum_value(probability[i], num);
					if ((1 - error_sum) > 0.0001)
					{
						probability[i][kk] += 1 - error_sum;
					}
					kk--;
				}
				else if (probability[i][kkk] == 0)
					break;
			}
		}

		memset(signal_, 0, num * sizeof(double));
		error_sum = sum_value(probability[i], num);
		if ((1 - error_sum) > 0.0001)
		{
			exit(0);
		}
		delete[] positive_part;
		delete[] positive_ind;
		delete[] negative_part;
		delete[] negative_ind;
	}
	delete[] alpha_;
	delete[] beta_;
	delete[] sig_;
	delete[] indices;
	delete[] sign;
	delete[] temp_prob;
}
void allocBudget()
{
	Partition *partitions = network->getPartition();
	double *upload = network->getUpload();
	double *price = network->getUpprice();
	double *multiple = new double[algorithm->getDCNums()];
	double sum = 0;
	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		multiple[i] = upload[i] * price[i];
		sum += multiple[i];
	}
	for (int i = 0; i < algorithm->getDCNums(); i++)
	{
		partitions[i].setBudget(budget * (multiple[i] / sum));
	}
}
void ActionSelect()
{
	double **probability = algorithm->getProbability();
	Vertex *v = graph->getVertexs();
	for (int i = 0; i < graph->getNum_Nodes(); i++)
	{
		// int index = max_value_index(probability[i],algorithm->getDCNums());
		v[i].setAction(trimmer_tool(probability[i], algorithm->getDCNums()));
	}
}
// NOTE: Select the action node DC
// MARK: need to change
void *MyActionSelect(void *arguments)
{
	Pthread_args *args = (Pthread_args *)arguments;
	unsigned int low = args->getLow();
	unsigned int high = args->getHigh();
	int iter = args->getIter();
	double **probability = algorithm->getProbability();
	Vertex *v = graph->getVertexs();
	// MARK: what is epison for, what's the use of it
	double epison = 0.8;
	// srand(time(NULL) + low);
	/***************************选择概率最大的DC**************************************************/
	for (int i = low; i <= high; i++)
	{
		map<int, int>::iterator itora = trained.find(i);
		if (itora != trained.end())
		{
			double random = rand() % 100 / (double)99;
			int index = max_value_index(probability[i], algorithm->getDCNums());
			// cout<<index<<endl;
			// cout<<random_at_most(algorithm->getDCNums()-1)<<endl;

			// NOTE: if random less than epison, then choose the max possibility DC, else random choose another one
			if (random <= epison)
				v[i].setAction(index);
			else
				v[i].setAction(random_at_most(algorithm->getDCNums() - 1));
		}
	}
	delete args;
	/***************************轮盘赌算法**************************************************/
	/*double *array=new double[algorithm->getDCNums()];
	double randpara;
	for(int i=0;i < graph->getNum_Nodes(); i++)
	{
		array[0]=probability[i][0];
		for(int j=1;j<algorithm->getDCNums();j++)
	  {
		array[j]=array[j-1]+probability[i][j];
	   }
	   randpara=rand()%1001/(double)1000;
	   if(randpara<=array[0])
		v[i].setAction(0);
		else
		{
		  for(int k=1;k<algorithm->getDCNums();k++)
		  {
			  if(randpara>array[k-1]&&randpara<=array[k])
				{
					v[i].setAction(k);
					break;
				}
		  }
		}
	}
	delete[] array;*/
}
// TODO: making the sampling not random
void *Sampling(void *arguments)
{
	Pthread_args *args = (Pthread_args *)arguments;
	unsigned int low = args->low;
	unsigned int high = args->high;
	// srand(time(NULL) + high);
	for (int i = low; i <= high; i++)
	{
		double random = rand() % 100 / (double)99;
		if (random <= randpro)
		{
			pthread_mutex_lock(&mutex);
			trained.insert(pair<int, int>(i, 1));
			pthread_mutex_unlock(&mutex);
		}
	}
	delete args;
}

void *Sampling_adv(void *arguments)
{
	Pthread_args *args = (Pthread_args *)arguments;
	unsigned int low = args->low;
	unsigned int high = args->high;

	int threshold = args->getThreshold();
	double percent = args->getPercentage();
	long percent_nodes_num = percent * graph->getNum_Nodes();
	Vertex *temp_v = graph->getVertexs();
	int test_character = args->getTestCharacter();
	int delete_zero_or_not = args->getDeleteZero();
	// srand(time(NULL) + high);
	for (int i = low; i <= high; i++)
	{
		// double random = rand() % 100 / (double)99;
		if (test_character == 1)
		{
			if (delete_zero_or_not == 0)
			{
				if (temp_v[i].getIngoingEdges().size() <= threshold)
				{
					pthread_mutex_lock(&mutex);
					if (pthread_counter <= percent_nodes_num)
					{
						trained.insert(pair<int, int>(i, 1));
						pthread_counter++;
					}
					pthread_mutex_unlock(&mutex);
				}
			}
			// MARK: need to check
			else if (delete_zero_or_not == 1)
			{
				if (temp_v[i].getIngoingEdges().size() <= threshold && temp_v[i].getIngoingEdges().size() != 0)
				{
					pthread_mutex_lock(&mutex);
					if (pthread_counter <= percent_nodes_num)
					{
						trained.insert(pair<int, int>(i, 1));
						pthread_counter++;
					}
					pthread_mutex_unlock(&mutex);
				}
			}
		}
		else if (test_character == 5)
		{
			// cout << "Working..." << endl;
			if ((temp_v[i].getIngoingEdges().size() + temp_v[i].getOutgoingEdges().size()) <= threshold)
			{
				pthread_mutex_lock(&mutex);
				if (pthread_counter <= percent_nodes_num)
				{
					trained.insert(pair<int, int>(i, 1));
					pthread_counter++;
				}
				pthread_mutex_unlock(&mutex);
			}
		}
	}
	delete args;
}

int trimmer_tool(double *probability, int num_event)
{
	int returned_index = -1;
	double max_probality = max_value(probability);
	double random_number = 0;
	if ((1 - max_probality) < 0.001)
	{
		returned_index = max_value_index(probability, algorithm->getDCNums());
		return returned_index;
	}
	int factor = 2;
	int num_events = num_event;
	double *probabilities = new double[num_events];
	int *indices = new int[num_events];
	for (int i = 0; i < num_events; i++)
	{
		probabilities[i] = probability[i];
		indices[i] = i;
	}
	while (num_events > factor)
	{
		random_number = (double)rand() / RAND_MAX;
		trim_probability(&probabilities, &indices, &num_events, random_number);
	}
	if (num_events == 1)
		returned_index = indices[0];
	else if (num_events == factor)
	{
		random_number = (double)rand() / RAND_MAX;
		returned_index = (random_number < probabilities[0]) ? indices[0] : indices[1];
	}
	else
	{
		cout << "Error: trimmer_tool()" << endl;
		exit(-1);
	}
	delete[] probabilities;
	delete[] indices;
	return returned_index;
}
void trim_probability(double **probabilities, int **indices, int *num_events, double random_number)
{
	double probability = 0.0;
	double sum_probability = 1.0;
	int factor = 2;
	double seprator = (double)sum_probability / factor;
	double epsilon = 1e-6;
	int left_max = 0;
	int left_num = 0;
	int right_min = 0;
	int right_num = 0;
	int min_idx = (*indices)[0];
	int k = 0;
	int j = 0;
	do
	{
		probability += (*probabilities)[k];
		k++;
	} while (abs_double(probability - seprator) > epsilon && (probability < seprator));

	right_min = k;
	left_max = right_min - 1;
	left_num = left_max + 1;
	right_num = *num_events - left_num;
	if (abs_double(probability - seprator) > epsilon)
	{
		right_min--;
		right_num = right_num + 1;
	}
	if (random_number < seprator)
	{
		int *left_indices = new int[left_num];
		double *left_probabilities = new double[left_num];
		for (k = left_max; k >= left_num - left_max - 1; k--)
		{
			j = left_num - k - 1;
			left_indices[j] = j + min_idx;
			left_probabilities[j] = (*probabilities)[j];
		}
		if (abs_double(probability - seprator) > epsilon)
			left_probabilities[left_num - 1] = left_probabilities[left_num - 1] - (sum_value(left_probabilities, left_num) - seprator);
		for (k = 0; k < left_num; k++)
		{
			left_probabilities[k] *= factor;
		}
		*num_events = left_num;
		delete[] * probabilities;
		delete[] * indices;
		*probabilities = new double[*num_events];
		*indices = new int[*num_events];
		for (int i = 0; i < *num_events; i++)
		{
			(*probabilities)[i] = left_probabilities[i];
			(*indices)[i] = left_indices[i];
		}
		delete[] left_probabilities;
		delete[] left_indices;
	}
	else
	{
		int *right_indices = new int[right_num];
		double *right_probabilities = new double[right_num];
		for (k = right_min; k < right_num + right_min; k++)
		{
			j = k - right_min;
			right_indices[j] = k + min_idx;
			right_probabilities[j] = (*probabilities)[k];
		}
		if (abs_double(probability - seprator) > epsilon)
			right_probabilities[0] = right_probabilities[0] - (sum_value(right_probabilities, right_num) - seprator);
		for (k = 0; k < right_num; k++)
		{
			right_probabilities[k] *= factor;
		}
		*num_events = right_num;
		delete[] * probabilities;
		delete[] * indices;
		*probabilities = new double[*num_events];
		*indices = new int[*num_events];
		for (int i = 0; i < *num_events; i++)
		{
			(*probabilities)[i] = right_probabilities[i];
			(*indices)[i] = right_indices[i];
		}
		delete[] right_probabilities;
		delete[] right_indices;
	}
}

void update_Download_Upload(int index)
{
	double *uploadnum = new double[algorithm->getDCNums()];
	double *downloadnum = new double[algorithm->getDCNums()];
	Vertex *v = graph->getVertexs();
	list<Edge *>::iterator it;
	long long *mapped = graph->getMapped();
	Partition *partition = network->getPartition();

	for (int j = 0; j < algorithm->getDCNums(); j++)
	{
		uploadnum[j] = partition[j].getUploadNum();
		downloadnum[j] = partition[j].getDownloadNum();
	}
	list<Edge *> outEdges = v[index].getOutgoingEdges();
	for (it = outEdges.begin(); it != outEdges.end(); it++)
	{
		if (v[mapped[(*it)->getdestID()]].getLabel() == v[index].getLabel())
		{
			downloadnum[v[index].getLabel()] += sizeof(v[index].getValue());
			uploadnum[v[index].getAction()] += sizeof(v[index].getValue());
		}
		else if (v[mapped[(*it)->getdestID()]].getLabel() == v[index].getAction())
		{
			uploadnum[v[index].getLabel()] -= sizeof(v[index].getValue());
			downloadnum[v[index].getAction()] -= sizeof(v[index].getValue());
		}
		else
		{
			uploadnum[v[index].getLabel()] -= sizeof(v[index].getValue());
			uploadnum[v[index].getAction()] += sizeof(v[index].getValue());
		}
	}
	list<Edge *> inEdges = v[index].getIngoingEdges();
	for (it = inEdges.begin(); it != inEdges.end(); it++)
	{
		if (v[mapped[(*it)->getsourceID()]].getLabel() == v[index].getLabel())
		{
			uploadnum[v[index].getLabel()] += sizeof(v[mapped[(*it)->getsourceID()]].getValue());
			downloadnum[v[index].getAction()] += sizeof(v[mapped[(*it)->getsourceID()]].getValue());
		}
		else if (v[mapped[(*it)->getsourceID()]].getLabel() == v[index].getAction())
		{
			downloadnum[v[index].getLabel()] -= sizeof(v[mapped[(*it)->getsourceID()]].getValue());
			uploadnum[v[index].getAction()] -= sizeof(v[mapped[(*it)->getsourceID()]].getValue());
		}
		else
		{
			downloadnum[v[index].getLabel()] -= sizeof(v[mapped[(*it)->getsourceID()]].getValue());
			downloadnum[v[index].getAction()] += sizeof(v[mapped[(*it)->getsourceID()]].getValue());
		}
	}
	partition[v[index].getLabel()].setDownloadNum(downloadnum[v[index].getLabel()]);
	partition[v[index].getLabel()].setUploadNum(uploadnum[v[index].getLabel()]);
	partition[v[index].getAction()].setDownloadNum(downloadnum[v[index].getAction()]);
	partition[v[index].getAction()].setUploadNum(uploadnum[v[index].getAction()]);
	delete[] uploadnum;
	delete[] downloadnum;
}
void *probability_Update(void *arguments)
{
	Pthread_args *args = (Pthread_args *)arguments;
	unsigned int low = args->getLow();
	unsigned int high = args->getHigh();
	int iter = args->getIter();
	double **probability = algorithm->getProbability();
	double alpha = 0.6;
	double belta = 0.2;
	Vertex *v = graph->getVertexs();
	double *signal = NULL;
	for (int i = low; i <= high; i++)
	{
		map<int, int>::iterator itora = trained.find(i);
		signal = v[i].getSignal_();
		if (itora != trained.end())
		{
			int index = max_value_index(signal, algorithm->getDCNums());
			for (int j = 0; j < algorithm->getDCNums(); j++)
			{
				/*if (v[i].getSignal() == 0)
			{
				if (v[i].getAction() == j)
					probability[i][j] += alpha * (1 - probability[i][j]);
				else
					probability[i][j] = probability[i][j] * (1 - alpha);
			}
			else
			{
				if (v[i].getAction() == j)
					probability[i][j] = probability[i][j]*(1-beta);
				else
					probability[i][j] = probability[i][j] * (1 - beta)+beta/(algorithm->getDCNums()-1);
			}
					   */
				if (index == v[i].action || true)
				{
					if (index == j)
						probability[i][j] += alpha * (1 - probability[i][j]);
					else
						probability[i][j] = probability[i][j] * (1 - alpha);
				}
				else
				{
					if (index == j)
						probability[i][j] = probability[i][j] * (1 - belta);
					else
						probability[i][j] = probability[i][j] * (1 - belta) + belta / (algorithm->DCnums - 1);
				}
			}
		}
		for (int k = 0; k < algorithm->getDCNums(); k++)
		{
			signal[k] = 0;
		}
	}
	delete args;
}
