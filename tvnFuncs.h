//COMMON UTILS SECTION
void printClique(const vector<vertex* >&clique){
  printf("Clique size: %zu => ",clique.size());
  for (int i = 0 ; i < clique.size() ; ++ i){
    printf("%d ",clique.at(i)->v_id);
  }
  cout << "\n";
}

bool verifyClique(const vector<vertex *> &clique){
  for (int i = 0 ;i <clique.size(); ++i){
    for (int j = 0 ; j < clique.size();++j){
      if (i!=j){
	if(!getDIMACSBinaryEdgeSwap(clique.at(i)->v_id,clique.at(j)->v_id)){
	  return false;
	}
      }
    }
  } 
  return true;
}

template<class T>
void printVector(const vector<T> &vec){
  int count = 0;
  for (int i = 0 ; i < vec.size() ; ++ i){
    count++;
    cout << vec[i];
    if (i<vec.size()-1){ cout << ",";}
    if ((count % 10)==0){ cout << "\n";	}
  }
  cout << "\n";
}


template<class T>
void printVectorPtr(const vector<T> &vec){
  for (int i = 0 ; i < vec.size() ; ++ i){
    cout << vec[i]->v_id ; 	
    if (i<vec.size()-1)  cout << ",";
  }
}



//INITALIZING SECTION

void initVerticesAndEdges(const int nVertices,const int nEdges){

  pVertices = new vertex* [nVertices];
  pEdges = new edge* [nEdges];

  int edgeCount=0;
  int vertexCount = 0; 
  for (int i = 0 ;  i < nVertices ; ++i){
    //vertex v ;	v.v_id = i;	v.score=0.0; v.status=-1; v.intDeg = 0;
    vertex *vPtr = new vertex;
    vPtr->v_id = i ; 
    vPtr->status = -100000 ; 
    vPtr->intDeg = -10000000;
    vPtr->score = -1000000.0 ; 
    vPtr->numAnts = 0;

    pVertices[i]=vPtr;
    vertexCount++;

    //todo: when done, change EdgeSwap to just Edge 
    //for performance
    for (int j = 0 ; j < i ;++j){
	  
#ifdef _DEBUG
      assert (!(i==j && getDIMACSBinaryEdgeSwap(i,j)));
#endif

      if (getDIMACSBinaryEdgeSwap(i,j)){
	edge *edgePtr = new edge; 
	edgePtr->e_id =edgeCount;
	edgePtr->pOnEdge=0;//init value
	edgePtr->lastUpdate=0;//never updated before
	edgePtr->pointA=i;edgePtr->pointB=j;
	pEdges[edgeCount]=edgePtr; //todo: might remove later, not if I decide to 
	//tie everything into a struct


	pVertices[i]->adj.push_back(j);
	pVertices[i]->edgeList.push_back(edgeCount);		

	pVertices[j]->adj.push_back(i);
	pVertices[j]->edgeList.push_back(edgeCount);				
				
	edgeCount++;
		
      }//end if (getDIMACSBinaryEdgeSwap(i,j))
    }//end for j 
  }//end for i

#ifdef _DEBUG  
  assert(edgeCount==nEdges && vertexCount==nVertices);
#endif
}


void setGreedyClique(vector<vertex *>&theGreedyClique,const int nVertices){
  vertex *pVertex ;

  //openmp , probably can parallize this loop
  for (int i = 0;  i < nVertices; i++) {
    pVertex = pVertices[i];
    pVertex->status = POSSIBLY_IN_OUTPUT;
    sv[i] = pVertex;
  }

  qsort((void*)sv, (size_t)nVertices, sizeof(vertex*), compare_degree);

  for (int i = 0; i < nVertices; i++){
    pVertex = sv[i];

    if (pVertex->status == POSSIBLY_IN_OUTPUT){
      pVertex->status = IN_OUTPUT;
      theGreedyClique.push_back(pVertex);

      //probably can parallel this loop
      for (int j = i+1 ; j < nVertices;++j){
	if (!getDIMACSBinaryEdgeSwap(pVertex->v_id,sv[j]->v_id)){
	  sv[j]->status=NOT_IN_OUTPUT;
	}
      }//end for j
    }//end     if (pVertex->status == POSSIBLY_IN_OUTPUT){
  }

  

#ifdef _DEBUG
  printClique(theGreedyClique);
#endif
}

void setAntsOnGreedyClique(const vector<vertex *>&theGreedyClique,
			   const int nVertices,
			   const int nAnts){

  for (int i = 0 ; i < nAnts; ++ i){
    ant *a = new ant();	//init ant
    a->a_id = i;  	a->age=0; 

    if (genrand_int32()%100 < probAntsOnGreedyClique){
      a->current = theGreedyClique.at(genrand_int32()%theGreedyClique.size());
    }
    else{
      //this also includes the vertices in the greedyClique
      int randIndex =genrand_int32()%nVertices;
#ifdef _DEBUG
      assert(pVertices[randIndex]->v_id==randIndex);
#endif
      a->current=pVertices[randIndex];
    }

    //undecided yet
    a->last = NULL;	a->next = NULL; a->edgeToNext = NULL; 
	
    pVertices[a->current->v_id]->numAnts++;
    vAnts.push_back(a);
  }
}







//Local Optimize Utilities
void updateVertex(vertex* pVertex){
  if(vAnts.size()<MAX_ANTS- nExtraAnts){
    while(pVertex->numAnts<nExtraAnts){
      ant *a = new ant(); 
      a->a_id = vAnts.size(); //not size()-1 since already has ant w/ that id	
      a->age = 0;  a->current=pVertex; //add this ant to this vertex
      a->next = NULL; a->last = NULL; a->edgeToNext = NULL;


      pVertex->numAnts++;
      vAnts.push_back(a);
    }//end while
    updatePheromoneAdjEdges(pVertex->edgeList,
			    REINFORCEMENT_PHEROMONE);
  }
}



void setScores(const int dTime,
	       const int nVertices,const int nEdges,
	       const double ALPHA,
	       const double BETA){

#ifdef PSCORE
#warning "PSCORE: OMP for setScore"
#pragma omp parallel for
#endif
  for (int i=0;i<nVertices;++i){
    // vector<int>adjE=pVertices[i]->edgeList; //edges to vertices adjacent to old(curr)Pos
    int pOnAdjEdges=0;
    for (int j = 0 ; j < pVertices[i]->edgeList.size() ; ++j){
      pOnAdjEdges+=pEdges[pVertices[i]->edgeList.at(j)]->pOnEdge;
    }


    pVertices[i]->score=(ALPHA*pVertices[i]->numAnts + BETA*pOnAdjEdges);

  }//end for i
}


// END SECTION
void printBestSol(const vector<vector<vertex *> >&sol,const int nVertices, const int nEdges){
  int bestIndex=0;  int bestSize = sol.at(bestIndex).size();
  for (int i = 1 ; i < sol.size();++i){
    //assert(verifyClique(sol.at(i)));
    if (bestSize<sol.at(i).size()){
      bestSize=sol.at(i).size();
      bestIndex=i;
    }
  }
  
  
  printf("nthreads: %d cliq_siz: %zu vertices: %d edges: %d ",
	 nthreads,sol.at(bestIndex).size(),nVertices,nEdges);

  printf("clique: (");
  printVectorPtr(sol.at(bestIndex));
  printf(") ");
  cout << "seed: " << seed_t << endl;
}


//DEBUG SECTION

void printDHasP(const char* c, const int nEdges,int printAlot){
  //dhasP code
  int testPheromoneAmount=0;	int dHasP=0;
  for (int i = 0 ; i < nEdges; ++ i){
    if (pEdges[i]->pOnEdge>0){
      testPheromoneAmount+=pEdges[i]->pOnEdge; dHasP++;
      if (printAlot==true) printf("dHasP edge %d has %d\n",i,pEdges[i]->pOnEdge);
      assert(i==pEdges[i]->e_id);
    }
    assert(pEdges[i]->pOnEdge>=0);
  }
  printf("DHasP: ");
  printf(c);
  printf(", p amount %d, dHasP %d, ants %zu \n",
	 testPheromoneAmount,dHasP,vAnts.size());
  //dHasPCode
}


void printAntsOnVertices(const vector<vector<int> > &antsOnVertices){
  for (int i = 0 ; i < antsOnVertices.size() ; ++i){
    printf("\nvertex %d has %zu ants:\n",i,antsOnVertices[i].size());
    printVector(antsOnVertices[i]);
  }

}
bool stillAClique(const vector<vertex *> &clique,
		  const vertex *v){

  for (int i = 0 ; i < clique.size() ;++i){
    if (!getDIMACSBinaryEdgeSwap(clique.at(i)->v_id,v->v_id)){
      return false;
    }
  }
  return true;
}


//testing vars
#ifdef _DEBUG
int antActivated;
int antUnactivated;
int antMove;
int antMoveR;
int antMoveH;
int antStay;
#endif



/// Utilities /////////
////////////////// Utilities /////////////////////////
/*
  void testNumAnts(const char* c, 
  const int nVertices){
  printf(c);
  printf("\n");
  for (int i = 0 ; i < nVertices; ++ i){
  vertex *pVertex = pVertices[i];
  if (antsOnVertices.at(pVertex->v_id).size()!=pVertex->numAnts){
  printf("v%d,  aov %d , numAnts %d\n",
  pVertex->v_id,
  antsOnVertices.at(pVertex->v_id).size(),
  pVertex->numAnts);
  assert(1==0);
  }
	
  }

  }
*/
