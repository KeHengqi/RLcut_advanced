#include"libgraph.h"
#include"math.h"
#include"lib.h"
#include<iostream>
#include<cstring>
#include<stdlib.h>    
#include<time.h> 
using namespace std;
#pragma warning(disable:4996)
/*libgraph.h�ļ�����Ķ���*/
Mirror::Mirror()
{
	this->id=-1;
	this->outD=0;
	this->inD=0;
}
void Mirror::reset()
{
	this->id=-1;
	this->outD=0;
	this->inD=0;
}
void Mirror::subOutD()
{
(this->outD)--;
}
void Mirror::subInD()
{
	(this->inD)--;
}
void Mirror::setID(int id)
{
	this->id=id;
}
void Mirror::setInD(int ind)
{
    this->inD=ind;
}
void Mirror::setOutD(int outd)
{
    this->outD=outd;
}
int Mirror::getID()
{
	return this->id;
}
void Mirror::addOutD()
{
	(this->outD)++;
}
void Mirror::addInD()
{
	(this->inD)++;
}
void Mirror::addInDMulti(int num)
{
	this->inD=(this->inD)+num;
}
int Mirror::getOutD()
{
	return this->outD;
}
int Mirror::getInD()
{
	return this->inD;
}
Message::Message(int source, int dest, double msgValue)
{
	this->source = source;
	this->dest = dest;
	this->value = msgValue;
}
void Partition::recMsg(Message msg)
{
	megs.push_front(msg);
}
void Partition::setUploadNum(double uploadnum)
{
	this->UploadNum = uploadnum;
}
double Partition::getUploadNum()
{
	return UploadNum;
}
void Partition::setDownloadNum(double downloadnum)
{
	this->DownloadNum = downloadnum;
}
double Partition::getDownloadNum()
{
	return DownloadNum;
}
double Partition::getBudget()
{
	return budget;
}
void Partition::setBudget(double budget)
{
	this->budget = budget;
}
Edge::Edge(int source, int dest, double value)
{
	this->sourceID = source;
	this->destID = dest;
	this->value = value;
}
int Edge::getsourceID()
{
	return sourceID;
}
void Edge::setsourceID(int source)
{
	this->sourceID = source;
}
int Edge::getdestID()
{
	return destID;
}
void Edge::setdestID(int dest)
{
	this->destID = dest;
}
double Edge::getValue()
{
	return this->value;
}
void Edge::setValue(double value)
{
	this->value = value;
}
Vertex::Vertex()
{
	type=-1;
}
Vertex::~Vertex()
{
	list<Edge*>::iterator it;
	for (it = outgoingEdges.begin(); it != outgoingEdges.end(); it++) 
		delete *it;
	for (it = ingoingEdges.begin(); it != ingoingEdges.end(); it++)
		delete *it;
	delete[] signal_;
}
void Vertex::setType(int type)
{
	this->type=type;
}
int Vertex::getType()
{
    return type;
}
int Vertex::getVertexID()
{
	return vertexID;
}
void Vertex::setVertexID(int vertexID)
{
	this->vertexID = vertexID;
}
double Vertex::getValue()
{
	return value;
}
void Vertex::setValue(double value)
{
	this->value = value;
}
void Vertex::setAction(int action)
{
	this->action = action;
}
int Vertex::getAction()
{
	return action;
}
void Vertex::resetSignal_()
{
        //signal_ = new double[algorithm->getDCNums()];
        //cout << algorithm->getDCNums() << endl;
        for (int i = 0; i < algorithm->getDCNums(); i++)
                signal_[i] = 0;
}
void Vertex::initSignal_()
{
	signal_ = new double[algorithm->getDCNums()];
	//cout << algorithm->getDCNums() << endl;
	for (int i = 0; i < algorithm->getDCNums(); i++)
		signal_[i] = 0;
}
double* Vertex::getSignal_()
{
	return this->signal_;
}
int Vertex::getSignal()
{
	return signal;
}
void Vertex::setSignal(int signal)
{
	this->signal = signal;
}
int Vertex::getLabel()
{
	return label;
}
void Vertex::setLabel(int label)
{
	this->label = label;
}
list<Edge*> Vertex::getOutgoingEdges()
{
	return this->outgoingEdges;
}
list<Edge*> Vertex::getIngoingEdges()
{
	return this->ingoingEdges;
}
void Vertex::addOutgoingEdge(Edge *e)
{
	this->outgoingEdges.push_back(e);
}
void Vertex::addIngoingEdge(Edge* e)
{
	this->ingoingEdges.push_back(e);
}
void Vertex::sendMessage(int destination, double msgValue)
{
	Message msg(vertexID, destination, msgValue);
	partition.recMsg(msg);

}
void Vertex::recMessage(Message mes)
{
	this->msgs.push_front(mes);
}
list<Message> Vertex::getMessage()
{
	return msgs;
}
void Vertex::clearMessage()
{
	msgs.clear();
}
void Vertex::showVertex()
{
	cout << "Vertex " << vertexID << " outgoing edges :" << endl;
	list<Edge*>::iterator it;
	for (it = outgoingEdges.begin(); it != outgoingEdges.end(); it++) {
		cout << (*it)->getsourceID() << "---->" << (*it)->getdestID() << " weights: " << (*it)->getValue() << endl;
	}
	cout << "Vertex " << vertexID << " ingoing edges :" << endl;
	for (it = ingoingEdges.begin(); it != ingoingEdges.end(); it++) {
		cout << (*it)->getsourceID() << "---->" << (*it)->getdestID() << " weights: " << (*it)->getValue() << endl;
	}
}
Graph::~Graph()
{
	delete[] nodes;
	delete[] mapped;
	for(int i=0;i<algorithm->getDCNums();i++)
	   delete[] this->mirrors[i];
	delete[] this->mirrors;
}
Graph::Graph()
{
	num_edges = 0;
	num_nodes = 0;
}

// mirror : A two-dimention array, the first dimention means dc, and the second means which node
// the second dimention is the same with "v"
void Graph::initMirrors(int num,int nodesnum)
{
	this->mirrors=new Mirror*[num];
	for(int i=0;i<num;i++)
	    this->mirrors[i]=new Mirror[nodesnum];

}

Mirror** Graph::getMirrors()
{
	return this->mirrors;
}
void Graph::setMapped(long long* mapped)
{
	this->mapped = mapped;
}
long long* Graph::getMapped()
{
	return mapped;
}
void Graph::addNode()
{
	num_nodes++;
}
void Graph::addEdge()
{
	num_edges++;
}
long long Graph::getNum_Nodes()
{
	return num_nodes;
}
long long Graph::getNum_Edges()
{
	return num_edges;
}
void Graph::setNodes(Vertex* v)
{
	nodes = v;
}
Vertex* Graph::getVertexs()
{
	return nodes;
}
void Graph::showGraph()
{
	for (int i = 0; i < num_nodes; i++)
	{
		nodes[i].showVertex();
	}
}

// Initialize the algorithm
void Algorithm::InitAlgorithm(char* name,int numDCs)
{
	this->name = new char[50];   //
	strcpy(this->name, name);
	this->name[strlen(name)] = '\0';
	// Probability that the node will go to which DC
	// In every choice, the node will choose to go to the highest probability dc
	probability = new double*[graph->getNum_Nodes()];
	for (int i = 0; i < graph->getNum_Nodes(); i++)
		probability[i] = new double[numDCs];

	this->DCnums = numDCs;
}
Algorithm::~Algorithm()
{
	delete[] name;
	for (int i = 0; i < graph->getNum_Nodes(); i++)
		delete[] probability[i];
	delete[] probability;

}

// Initialize the node in which DC
// Important!
// Maybe need to change
void Algorithm::InitLabel()
{
	Vertex* v = graph->getVertexs();
	int num = graph->getNum_Nodes() / DCnums;
	int label = 0;
	// srand(time(NULL));
	for (int i = 0; i < graph->getNum_Nodes(); i++)
	{
		if (i%num == 0&&label<DCnums-1&&i!=0)
			label++;
		v[i].setLabel(label);	// Set the node to the DC
		v[i].iniLabel=label;	// Annouce the initialize label
    
	   /*  int rl=random_at_most(DCnums-1);
         v[i].setLabel(rl);
		 v[i].iniLabel=rl;*/
		 
		/*  int l=rand()%DCnums;
		  v[i].setLabel(l);
		  v[i].iniLabel=l;*/
	}


}
int Algorithm::getDCNums()
{
	return DCnums;
}

// Set the possibility of every edge using average model
void Algorithm::InitProbability()
{
	for (int i = 0; i < graph->getNum_Nodes(); i++)
		for (int j = 0; j < algorithm->getDCNums(); j++)
			probability[i][j] = (double)1 / algorithm->getDCNums();	// Using average model
     /* srand(time(0));
      double* random=new double[DCnums-1];
      for (int i = 0; i < graph->getNum_Nodes(); i++)
        {
          int count=0;
          while(count<DCnums-1)
         { 
          random[count++]=rand() / double(RAND_MAX); 
         }
         bubble_sort(random,DCnums-1);
         
         for (int j = 0; j < algorithm->getDCNums(); j++)
         {
            if(j==0)
              probability[i][j]=random[j]-0;
            else if(j==DCnums-1)
              probability[i][j]=1-random[j-1];
            else
              probability[i][j]=random[j]-random[j-1];
         }
       /*  for(int w=0;w<algorithm->getDCNums();w++)
            cout<<probability[i][w]<<"  ";
            cout<<endl;
        }      */
}
double** Algorithm::getProbability()
{
	return probability;
}
Network::Network()
{
	upload = new double[algorithm->getDCNums()];
	download = new double[algorithm->getDCNums()];
	upprice = new double[algorithm->getDCNums()];
	uploadnumG=new double[algorithm->getDCNums()];
	uploadnumA=new double[algorithm->getDCNums()];
	downloadnumG=new double[algorithm->getDCNums()];
	downloadnumA=new double[algorithm->getDCNums()];
}
Network::~Network()
{
	delete[] uploadnumA;
	delete[] uploadnumG;
	delete[] downloadnumA;
	delete[] downloadnumG;
	delete[] upload;
	delete[] download;
	delete[] upprice;
	delete[] partitions;
}
void Network::initNetwork(double* upload, double* download, double* upprice)
{
	this->upload = upload;
	this->download = download;
	this->upprice = upprice;

}
double* Network::getUpload()
{
	return upload;
}
double* Network::getDownload()
{
	return download;
}
double* Network::getUpprice()
{
	return upprice;
}
void Network::setPartition(Partition* partitions)
{
	this->partitions = partitions;
}
Partition* Network::getPartition()
{
	return partitions ;
}
void Pthread_args::setIter(int iter)
{
	iteration = iter;
}
int  Pthread_args::getIter()
{
	return iteration;
}
void Pthread_args::setLow(int low)
{
	this->low = low;
}
int  Pthread_args::getId()
{
	return id;
}
void Pthread_args::setId(int id)
{
	this->id = id;
}
int  Pthread_args::getLow()
{
	return low;
}
void Pthread_args::setHigh(int high)
{
	this->high = high;
}
int Pthread_args::getHigh()
{
	return high;
}
void Pthread_args::setDc(int dc)
{
	this->dc=dc;
}
int Pthread_args::getDc()
{
	return dc;
}

void Pthread_args::setPercentage(double percent)
{
	this->percentage = percent;
}

double Pthread_args::getPercentage()
{
	return percentage;
}
void Pthread_args::setTestCharacter(int test_chara)
{
	this->test_character = test_chara;
}
int Pthread_args::getTestCharacter()
{
	return test_character;
}
void Pthread_args::setThreshold(int threshold)
{
	this->threshold_number = threshold;
}
int Pthread_args::getThreshold()
{
	return threshold_number;
}
void Pthread_args::setDeleteZero(int delete_zero_or_not)
{
	this->delete_zero = delete_zero_or_not;
}

int Pthread_args::getDeleteZero()
{
	return delete_zero;
}