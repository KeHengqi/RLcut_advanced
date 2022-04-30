//#define CRTDBG_MAP_ALLOC
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
//#include<crtdbg.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include "libgraph.h"
#include "lib.h"
#include "math.h"
#include "sort_ways.hpp"
using namespace std;
#define N 999
#pragma warning(disable : 4996)
char algorithmName[] = "revolver-hybrid";

int Initialize_the_parameter(int argc, char *argv[])
{
	if (!graph)
		return (-1);
	// MARK: CHANGE HERE ./main TUI_or_script app input network dcnum budget theta threads outputfilename randpro trainedfile bsp batchsize L n(app_iter*called) ginger_cost /*X trained_new_file X*/ test_character
	/*//NOTE:
		1. ./main : means the programme
		2. TUI_or_script : using tui or running by script
	// All the parameters behind are provided to run the script
		3. app: application
		4. input: the graph need to input
		5. network: the network need to input
		6. dcnum: distributed data center number, depends on network
		7. budget: budget can be used
		8. theta: threshold for identifying high and low degree node
		9. threads: threads num want to use, depends on machine
		10. outputfilename: Result of the seperation
		11. randpro: a standard number in a random process
		12. trainedfile: results of every iteration
		13. bsp: 0 for not using bsp mode, and 1 for using
		14. batchsize: variable should be defined if using bsp mode
		15. L: left time
		16. n: no idea, not important
		17. ginger cost: no idea, not important
		18. Test character: in-going, mirror, DC_speed, DC_price, degree

		// Following option is seperated by the testing character
		I. In-going
		19. ratio_of_nodes: the percentage of nodes need to train
		20. deleted zero or not: not counting zero in-going degree nodes
		II. Degree
		19. ratio_of_nodes: the percentage of nodes need to train
		// STOPLINE April 27th: Haven't done the code
		III. mirror
		IV. DC_speed
		V. degree

		19. trained new file: my newly add output file	// Abandon
		khq note
	*/
	// NOTE: trained_new file is the new argument
	char application[100];
	char filename_graph[100];
	char filename_network[100];
	char test_character[100];
	bool directed = true;
	int numDCs;
	// int TUI_OR_SCRIPT = 2;	// 1 for tui, and 2 for script
	if (!strcmp(argv[1], "TUI"))
	{
		// TUI_OR_SCRIPT = 1;
		cout << "Welcome to graph seperation simulation and evaluation platform!" << endl;
		cout << "Now input the application of graph: ";
		cin >> application;
		cout << "Now input the graph file you want to test: ";
		cin >> filename_graph;
		cout << "Now input the network you want to simulate: ";
		cin >> filename_network;
		cout << "Now input the Data Center's number: ";
		cin >> numDCs;
		cout << "Now input the budget: ";
		cin >> budget;
		cout << "Now input the threshold indentifying high and low degree: ";
		cin >> theta;
		cout << "Now set the thread num: ";
		cin >> MAX_THREADS_NUM;
		cout << "Now set the path you want to put the result: ";
		cin >> outputfile;
		cout << "Now set the stability of the algorithm: ";
		cin >> randpro;
		cout << "Now input the file you want to store the training data: ";
		cin >> FN;

		cout << "Running in batch(Y or n): ";
		char yes_or_no;
		cin >> yes_or_no;
		if (yes_or_no == 'y' || yes_or_no == 'Y')
		{
			bsp = 1;
			cout << "Now set the batch size: ";
			cin >> batchsize;
		}
		else
		{
			bsp = 0;
			batchsize = -1;
		}
		cout << "Now input the time limit: ";
		cin >> L;
		cout << "Now set the move percent: ";
		cin >> mvpercents;
		getchar();
		mvpercents = 1;
		ginger_cost = 0;
		cout << "Here are the strategy can choose: " << endl;
		cout << "1. in-going" << endl;
		cout << "2. mirror" << endl;
		cout << "3. DC_speed" << endl;
		cout << "4. DC_price" << endl;
		cout << "5. degree" << endl;
		cout << "Choose the strategy you want to test: ";
		cin >> test_character;
		cout << "Now loading..." << endl;
	}
	else
	{
		// Warning: should be updated as the condition is not static
		if (argc < 18)
		{
			cout << "The number of parameter has problem, please input the parameters later." << endl;
			return (-1);
		}

		strcpy(application, argv[2]);
		cout << application << endl;
		strcpy(filename_graph, argv[3]);
		cout << filename_graph << endl;
		strcpy(filename_network, argv[4]);
		cout << filename_network << endl;

		numDCs = atoi(argv[5]);
		cout << numDCs << endl;
		budget = atof(argv[6]);
		cout << budget << endl;
		theta = atoi(argv[7]);
		cout << theta << endl;
		MAX_THREADS_NUM = atoi(argv[8]);
		cout << MAX_THREADS_NUM << endl;
		// Set up the output file name
		strcpy(outputfile, "");
		strcat(outputfile, argv[9]);
		// Fix the randpro, maybe should change in the future
		randpro = atof(argv[10]);
		cout << randpro << endl;
		strcpy(FN, argv[11]);
		// randpro = 0.01; //从1%开始
		bsp = atoi(argv[12]);
		batchsize = atoi(argv[13]);
		L = atof(argv[14]);
		mvpercents = atoi(argv[15]);
		ginger_cost = atof(argv[16]);
		strcpy(test_character, argv[17]);

		/* Initialize test character start */
		if (!strcmp(test_character, "in-going"))
		{
			cout << "Now test the in-going effect" << endl;
			test_character_num = 1;
			ratio_of_nodes = atof(argv[18]);
			cout << ratio_of_nodes << endl;
			cout << argv[19] << endl;
			if (!strcmp(argv[19], "y") || !strcmp(argv[19], "Y"))
			{
				deleted_zero_or_not = true;
			}
			else
			{
				deleted_zero_or_not = false;
			}
		}
		else if (!strcmp(test_character, "mirror"))
		{
			cout << "Now test the mirror effect" << endl;
			test_character_num = 2;
		}
		else if (!strcmp(test_character, "DC_speed"))
		{
			cout << "Now test the DC_speed effect" << endl;
			test_character_num = 3;
		}
		else if (!strcmp(test_character, "DC_price"))
		{
			cout << "Now test the DC_price effect" << endl;
			test_character_num = 4;
		}
		else if (!strcmp(test_character, "degree"))
		{
			cout << "Now test the degree effect" << endl;
			test_character_num = 5;
			ratio_of_nodes = atof(argv[18]);
		}
		else if (!strcmp(test_character, "N") || !strcmp(test_character, "n"))
		{
			cout << "Now running original algorithm" << endl;
			test_character_num = -1;
		}
		else
		{
			cout << "Error input" << endl;
			exit(EXIT_FAILURE);
		}
		/* Initialize test character end */
	}

	/* Initialize application start */
	if (strcmp(application, "pagerank") == 0)
		dataunit = 0.000008;
	else if (strcmp(application, "sssp") == 0)
		dataunit = 0.000004;
	else if (strcmp(application, "subgraph") == 0)
	{
		dataunit = 0.001;
		directed = false;
	}
	else
	{
		cout << "application error!" << endl;
	}
	/* Initialize application end */

	/* Reading graph start */
	int status = read_input_file(filename_graph, directed);
	if (status == -1)
	{
		cout << "Graph initial is error!" << endl;
		exit(-1);
	}
	/* Reading graph end */

	initAlgorithm(algorithmName, numDCs);
	
	/* Reading network start */
	double *upload = new double[algorithm->getDCNums()];
	double *download = new double[algorithm->getDCNums()];
	double *upprice = new double[algorithm->getDCNums()];
	status = read_input_network(filename_network, upload, download, upprice);
	if (status)
	{
		cout << "network definition is error!" << endl;
		exit(-1);
	}
	network->initNetwork(upload, download, upprice);
	/* Reading network end */

	/* Initialize budget start */
	if (budget == 0)
	{
		cout << "Sorry, the budget is too low!" << endl;
		return -1;
	}
	stepbudget = 0;
	addbudget = budget / MAX_ITER;
	/* Initialize budget end */

	/* Initialize batch start */
	bspvec = new vector<Pthread_args>[batchsize];
	/* Initialize batch end */

	/* Initialize DC setting start */
	ug = new double[algorithm->DCnums];
	dg = new double[algorithm->DCnums];
	ua = new double[algorithm->DCnums];
	da = new double[algorithm->DCnums];
	/* Initialize DC setting end */

	/* Initialize graph seperation algorithm start */
	// MARK: the algorithmName variable is fixed, and this should be flexible with the user input
	graph->initMirrors(numDCs, graph->getNum_Nodes());
	/* Initialize graph seperation algorithm end */
	return 0;
	/*Algorithm and Network*/
	// theta=stoi(argv[6]);
	// char algorithmName[] = "revolver-hybrid";
}

int main(int argc, char *argv[])
{
	int TUI_OR_SCRIPT = 2;
	if (!strcmp(argv[1], "TUI"))
	{
		TUI_OR_SCRIPT = 1;
	}

	int initialize_result = Initialize_the_parameter(argc, argv);

	if (initialize_result == -1)
	{
		exit(EXIT_FAILURE);
	}
	optimize = new int[MAX_ITER];
	for (int i = 0; i < MAX_ITER; i++)
		optimize[i] = 1;
	/****************************************测试moveVertex函数***********************************************/
	/*test();
	return 0;*/
	Vertex *v = graph->getVertexs();
	long long *mapped = graph->getMapped();
	/*--------------------------------------Powerlyra算法--------------------------------------------------*/
	if (budget < 0)
	{
		// clock_t start=clock();
		time_t start = time(NULL);
		/*	cout<<"--------------------Hash powerlyra算法---------------------"<<endl;
			for(int i=0;i<graph->getNum_Nodes();i++)
			{
				int dcnum=algorithm->getDCNums();
				int label=i%dcnum;
				v[i].setLabel(label);
			}
			overhead=time(NULL)-start;
			initPartitioning();
			powerlyrareCalTimeAndCost(timeAll,priceAll);
			//cout<< "#timeAll: " << timeAll << ";  #priceAll: " << priceAll << endl;

			//Simulate();
			//overhead=(clock() - start) * 1.0 / CLOCKS_PER_SEC;
			//cout<<"Overhead: "<<overhead<<" seconds!"<<endl;*/

		cout << "--------------------Heuristic powerlyra算法---------------------" << endl;
		algorithm->InitLabel();
		double *score = new double[algorithm->getDCNums()];
		double *load = new double[algorithm->getDCNums()];
		double *eload = new double[algorithm->getDCNums()];
		double *locality = new double[algorithm->getDCNums()];
		double *proc_balance = new double[algorithm->getDCNums()];
		list<Edge *> oedges;
		list<Edge *> iedges;
		list<Edge *>::iterator it;
		double vtoe = (double)graph->getNum_Nodes() / (double)graph->getNum_Edges();
		int low = 0, high = 0;
		double gamma = 1.5;
		double alpha = sqrt(algorithm->getDCNums()) * double(graph->getNum_Edges()) / pow(graph->getNum_Nodes(), gamma);
		for (int i = 0; i < graph->getNum_Nodes(); i++)
		{
			if (v[i].getIngoingEdges().size() >= theta)
			{
				v[i].setType(0);
				high++;
			}
			else
			{
				v[i].setType(1);
				low++;
			}
		}
		// cout<<"low: "<<low<<endl;
		// cout<<"high: "<<high<<endl;
		start = time(NULL);
		for (int j = 0; j < algorithm->getDCNums(); j++)
		{
			load[j] = 0;
			eload[j] = 0;
			proc_balance[j] = 0;
		}
		for (int i = 0; i < (graph->getNum_Nodes()) / 1 && (rand() % (N + 1) / (float)(N + 1) <= 1); i++) // lalalal
		{
			int followedges = 0;
			for (int j = 0; j < algorithm->getDCNums(); j++)
			{
				score[j] = 0;
				locality[j] = 0;
			}
			if (v[i].getType() == 1) //先处理low-degree
			{
				iedges = v[i].getIngoingEdges();
				for (it = iedges.begin(); it != iedges.end(); it++)
				{
					followedges++;
					long long sourceId = mapped[(*it)->getsourceID()];
					if (sourceId < i)
						locality[v[sourceId].getLabel()]++;
				}
				oedges = v[i].getOutgoingEdges();
				for (it = oedges.begin(); it != oedges.end(); it++)
				{
					long long destId = mapped[(*it)->getdestID()];
					if (v[destId].getType() == 0)
						followedges++;
					if (v[destId].getType() == 0 && destId < i)
						locality[v[destId].getLabel()]++;
				}
			}
			else //处理high-degree
			{
				oedges = v[i].getOutgoingEdges();
				for (it = oedges.begin(); it != oedges.end(); it++)
				{
					long long destId = mapped[(*it)->getdestID()];
					if (v[destId].getType() == 0)
						followedges++;
					if (v[destId].getType() == 0 && destId < i)
						locality[v[destId].getLabel()]++;
				}
			}
			double localitysum = 0, balancesum = 0;
			for (int j = 0; j < algorithm->getDCNums(); j++)
			{
				localitysum += locality[j];
				balancesum += (load[j] + vtoe * eload[j]) / 2;
			}
			for (int j = 0; j < algorithm->getDCNums(); j++)
			{
				double p1 = localitysum == 0 ? 0 : locality[j] / localitysum;
				double p2 = balancesum == 0 ? 0 : ((load[j] + vtoe * eload[j]) / 2) / balancesum;
				score[j] = locality[j] - alpha * gamma * pow(proc_balance[j], (gamma - 1));
			}
			int index = max_value_index(score, algorithm->getDCNums());
			v[i].setLabel(index);
			load[index]++;
			eload[index] += followedges;
			proc_balance[index]++;
			proc_balance[index] += followedges * vtoe;
		}
		overhead = time(NULL) - start;
		for (int j = 0; j < algorithm->getDCNums(); j++)
		{
			// cout<<"DC "<<j<<" vload: "<<load[j]<<", eload: "<<eload[j]<<endl;
		}
		initPartitioning();
		powerlyrareCalTimeAndCost(timeAll, priceAll);
		Simulate();

		// cout<< "#timeAll: " << timeAll << ";  #priceAll: " << priceAll << endl;

		delete[] score;
		delete[] load;
		delete[] eload;
		delete[] locality;
		return 0;
	}

	/**********************选择训练的点*********************************************/
	// NOTE: Using rand() to select the nodes need to train
	//  char* FN=new char[100];

	// TODO: getting the conclusion, then can change the code here
	// strcpy(FN, argv[10]);
	//  ofstream out(FN);
	trained.clear();
	srand(time(NULL));
	Seed = rand();
	srand(Seed);
	// out<<"The seed is: "<<Seed<<endl;
	for (int i = 0; i < graph->getNum_Nodes(); i++)
	{
		double random = rand() % 100 / (double)99;
		if (random <= randpro)
		{
			trained.insert(pair<int, int>(i, 1));
			// out<<i<<endl;
		}
	}
	/*----------------------------------------------LA-cut算法----------------------------------------------------------------------*/

	if (strcmp(algorithmName, "revolver-hybrid") == 0)
		revolver(TUI_OR_SCRIPT);

	// Vertex* v = graph->getVertexs();
	/*  ActionSelect();
	  for (int i = 0; i < graph->getNum_Nodes(); i++)
	  {
		  if (v[i].getLabel()!=v[i].getAction())
		  {
			  v[i].setLabel(v[i].getAction());
		  }
	  }
		  initPartitioning();
		  reCalTimeAndCost(timeAll,priceAll);
  cout << "after executing the algorithm:" << endl;
  //Vertex* v = graph->getVertexs();
  /*for (int i = 0; i < graph->getNum_Nodes(); i++)
	  cout << v[i].getVertexID() << " : " << v[i].getLabel() << endl;*/
	cout << "---------------------------------------------------" << endl;
	/*cout << "after executing the algorithm:" << endl;
	cout << "#timeAll: " << timeAll << endl;
	cout << "#priceAll: " << priceAll << endl;
	cout<<"Overhead: "<<overhead<<" seconds!"<<endl;*/
	//	Simulate();
	// delete[] eachThreadsNum;

	delete[] ug;
	delete[] dg;
	delete[] ua;
	delete[] da;
	return 0;
}
