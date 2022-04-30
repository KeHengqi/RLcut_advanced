#ifndef LIB_H
#define LIB_H
#include<iostream>
#include<map>
#include<vector>
#include"libgraph.h"
extern double alphaT;
extern double betaC;
extern Graph* graph;
extern double* ug;
extern double* dg;
extern double* ua;
extern double* da;
extern long mvpercents;
extern char* FN;
extern double ginger_cost;
extern double global_mvcost;
extern double current_mvcost;
/* Newly add part start */
extern int test_character_num; // 1 for in-going degree, 2 for mirror, 3 for DC_speed, 4 for DC_price, 0 for default, means no action
extern int pthread_counter; // a counter for the multi-thread
extern double ratio_of_nodes; // If test_character_num has meaning, then ratio should be set
extern bool deleted_zero_or_not;
/* Newly add part end */
extern double L;
extern int optitime;
extern int F1;
extern int F2;
extern double t1;
extern double t3;
extern double t4;
extern double t5;  //score时间
extern double t2;   //vertex move时间
extern double t6;
extern double moverand;
extern ofstream outfile;
extern float raw_data_size;
extern vector<Pthread_args> *bspvec;
extern map<int,vector<int>>PSS;
extern double dataunit;
extern int bsp;
extern int bsplock;
extern int batchsize;
extern Algorithm* algorithm;
extern Network* network;
extern int MAX_ITER;
extern int MAX_THREADS_NUM;
extern double timeAll;
extern double priceAll;
extern double timeOld ;
extern double priceOld ;
extern int theta;
extern char* outputfile;
extern double budget;
extern double stepbudget;
extern double addbudget;
extern double randpro;
extern double overhead;
extern int Seed;
extern map<int,int>trained;
extern int* optimize;
using namespace std;
void Simulate();
void* test(void* arguments);
void Mirrortest();
void reCalTimeAndCost(double&,double&);
void powerlyrareCalTimeAndCost(double&,double&);
void* Parallel_Score_Migration_function(void* arguments);
void* keepMirrorsStable(void* arguments);
void* keepMirrorsBackUpStable(void* arguments);
void* Parallel_Score_Signal_function(void* arguments);
char* input_format_specifier(char *,char*);
int read_input_file(char*,bool);
int read_input_network(char*,double*upload,double*download,double* upprice);
void initialNode(FILE* , char*);
void initialNode(FILE* , char* , double);
void initialEdge(FILE*, char*,bool);
void initialEdge(FILE*, char*,double,bool);
void initAlgorithm(char*,int);
void moveVertex(int ver,double* downloadnumG,double* uploadnumG,double* downloadnumA,double* uploadnumA,int destdc);
void* moveVertexBsp(void* arguments);
void vertexMigrate(int ver,double* downloadnumG,double* uploadnumG,double* downloadnumA,double* uploadnumA,int destdc);
void revolver(int TUI_OR_SCRIPT);
void ActionSelect();
void* MyActionSelect(void*);
void* Sampling(void*);
void* Sampling_adv(void*);
void Score_Migration_function(double*,double*,int,double*,double*);
void calculateTimeAndPrice(double&,double&);
void initPartitioning();
void* probability_Update(void*);
void update_Download_Upload(int);
void allocBudget();
int trimmer_tool(double*,int);
void trim_probability(double**,int**,int*,double);
void weighted_probability_update(int);
void multiprocessing_pool(int*);
double Benifit(int,int);
void objective_function(double*,int);
#endif
