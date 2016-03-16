#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "MT_rand.h"
#include "Dimacs.h"
#include "settings.h"
#include "maxclique.h"
#include "AntsOps.h" //header for Ants operations s.a moving etc
#include "localOptimize.h"
#include "localOptimizePAR.h"




void runStageLoop(const int nVertices,const int nEdges); 

int main(int argc, char *argv[]){
  
  if (argc < 2){
    cerr << "Usage: " << argv[0]  <<  " infile";
    exit(1);}

  if (argc==3) seed_t = atoi(argv[2]);
  else seed_t = time(0);


  init_genrand(seed_t);//randomize seed according to time

  //BUG: make sure that nothing is randomize in the slaves (except ants movement),
  //all randomize , ie.e local optimization is to be done on the MASTER
  
#ifdef D1
#warning "D1: Seeding info on"
  int myRandN = genrand_int32()%1000;
  cout << "Seed: " << seed_t << ", rand num " << myRandN << endl;
#endif


  //read in graph
  char *sInFile=argv[1];   int nEdges, nVertices, nAnts;
  readDIMACSBinaryFormat(sInFile,nVertices,nEdges);

  nAnts=antMultipler*nVertices;
  MAX_ANTS=static_cast<int>(nAnts*antFactor);
  ANTS_PER_EDGE = (double)nAnts/nEdges;
  ANTS_PER_VERTEX = (double)nAnts/nVertices;

#ifdef _DEBUG
  printf("Graph %s, #vertices: %d, #nedges: %d, #ants: %d\n",
	 sInFile,nVertices,nEdges,nAnts);
#endif



#ifdef TESTEDGES//testing edges
#warning "TESTEDGES: turn on testing edges"
  cout<<"Seed ";
  printf("Who ");
  printf(" EdgeNum ");
  tEdges=genrand_int32()%nEdges;

  for (int i = 0 ; i < MAX_STAGE ;++i){
	printf(" Stage%d",i);
  }

  printf("\n");
  cout << seed_t << " ";
#ifdef RIZZO
  printf("Rizzo ");
#else
  printf("Tvn ");
#endif

  printf("%d ", tEdges);
#endif




#ifdef _OPENMP
  nthreads=omp_get_max_threads();
#else
  nthreads=1;
#endif 


  //initilizing 
  initVerticesAndEdges(nVertices,nEdges);  
  sv = new vertex* [nVertices];  sv2 = new vertex* [nVertices];
  vector<vertex *>theGreedyClique;

  setGreedyClique(theGreedyClique,nVertices);//find a greedy vertice

#ifdef _DEBUG
  printf("The Greedy Clique size is %d\n",theGreedyClique.size());
  assert(verifyClique(theGreedyClique));
#endif

  sol.push_back(theGreedyClique);//save the sol to vector, its index is 0


  //init ants to the greedy clique
  setAntsOnGreedyClique(theGreedyClique,nVertices,nAnts); 

  //printAntsOnVertices(antsOnVertices);





  //Begin the evolution
  
  runStageLoop(nVertices,nEdges);





  printBestSol(sol,nVertices,nEdges);  //print the best result


  cleanUp(nVertices,nEdges);//delete mem allocated



#ifdef D1
  cout << "Seed: " << seed_t << ", rand num " << myRandN << endl;
#endif


#ifdef USE_MPI
  MPI_Finalize();
#endif

  return 0;
}


void runStageLoop(const int nVertices,
				  const int nEdges){


  bool edgeToUpdate[nEdges]; //create a shared memory variable
  for (int iE=0; iE < nEdges; ++iE){edgeToUpdate[iE]=false;}

  VertexInfo mVertices[nVertices*nthreads]; //creating a (quite large) matrix info for data privating

  for (int iStage=0;iStage<MAX_STAGE;++iStage){

#ifdef _DEBUG
	  printf( "\n\nstage: %d\n\n",iStage);
#endif
	
	int iCycle = 0 ; //this is some weird thing, just to be consistent w/ his
	
	iCycle = runCycleLoop(iStage,nVertices,nEdges,edgeToUpdate);
	
	
	
#ifdef _DEBUG	
	printDHasP("after runCycle" , nEdges,false);
#endif
	


	//local optimization
	/*
	printf("stage %d\n",iStage);
	printDHasP("before Optimization" , nEdges,false);
	*/

	


	parallelLocalOpt(sol,mVertices,iStage,iCycle,nVertices,nEdges);
	
	

	//}


  }//end iStage
}
