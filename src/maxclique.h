//COMMON UTILS SECTION
void printClique(const vector<vertex* >&clique){
  printf("Clique size: %zu => ",clique.size());
  for (auto &x : clique){
    printf("%d ", x->v_id);
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


// WAS IN SETTINGS
void evaporatePheromone(edge *theEdge,const int &theAmount){
  theEdge->pOnEdge -= theAmount;
  if (theEdge->pOnEdge < MIN_EDGE_PHEROMONE){
    theEdge->pOnEdge = MIN_EDGE_PHEROMONE;
  }
}

void addPheromone(edge *theEdge,const int &theAmount){
  theEdge->pOnEdge += theAmount;
  if (theEdge->pOnEdge > MAX_EDGE_PHEROMONE){
    theEdge->pOnEdge = MAX_EDGE_PHEROMONE;
  }
}

void updatePheromoneAdjEdges(const vector<int> &theAdjEdges,
			     const int &amount){
    for (auto &x : theAdjEdges){
    assert(x == pEdges[x]->e_id);
    addPheromone(pEdges[x], amount);
  }
}


void cleanUp(const int &nVertices, const int &nEdges){
  for (int i = 0 ; i < nVertices ; ++i){
    delete pVertices[i];
  }
  delete pVertices;
  delete sv;//dangling pointer but will be fine when go out of scope
  delete sv2;

  for (int i = 0 ; i < nEdges ; ++i){
    delete pEdges[i];
  }
  delete pEdges;

  for (int i = 0 ; i <vAnts.size();++i){
    delete vAnts[i];
  }
  vAnts.clear();
}


int compare_int_deg(const void *v1, const void *v2 ){
  vertex* a1 = *(vertex**)v1;  vertex* a2 = *(vertex**)v2;  
  return a2->intDeg - a1->intDeg;
}

int compare_score (const void *v1, const void *v2 ){
  vertex* a1 = *(vertex**)v1; vertex* a2 = *(vertex**)v2;
  if (a2->score>a1->score){	return 1; }
  if (a2->score==a1->score){return 0; }
  if (a2->score<a1->score){return -1; }
  return 0; //default
}

int compare_degree (const void *v1, const void *v2 ){
  vertex* a1 = *(vertex**)v1;  vertex* a2 = *(vertex**)v2;
  return a2->adj.size() - a1->adj.size();
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
  printf("%s",c);
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



// ANTS OPERATIONS
vertex *chooseNextMoveHeuristiclyTVN(const int iTime,const ant *pAnt,
				     int &adjIndex,//stores edge id to new pos
				     bool edgeToUpdate[],
				     const int &nVertices,
				     const int &nEdges){

  vertex* oldPos = pAnt->current;

  vertex* newPos = NULL;//dummy, will be changed later
  adjIndex=-1;

  vector<int>adjV=oldPos->adj;//vertices adjacent to old(curr)Pos
  vector<int>adjE=oldPos->edgeList; //edges to vertices adjacent to old(curr)Pos


  //num of adj vertices to consider, used to be = getNSampleCount()
  //(should be mostly same as adjE.size())
  int iSampleCount=adjE.size();

  double vertexScoreList[iSampleCount];

  double dPheromoneScore = 0.0, dPopulationScore = 0.0, iPheromoneScore = 0.0,iPopulationScore = 0.0;
  long ScoreSum = 0;

  

  for (int i = 0; i < iSampleCount; i++){

    vertex *pVertex = pVertices[adjV.at(i)];


    // Check tabu list to disallow backtracking, except when on a peninsula.
    //Rizzo does not have (pAnt->last==NULL) , probably equivalent to mine , since 
    //if last==NULL then randVertex address is not same 
    if(pAnt->last==NULL||pVertex->v_id!=pAnt->last->v_id||adjV.size()==1){

      edge *pEdge = pEdges[adjE.at(i)];

      //I don't think the iTeim > will ever be true, since this is move later.
      if (edgeToUpdate[pEdge->e_id]==false&&iTime > pEdge->lastUpdate){  
	edgeToUpdate[pEdge->e_id]=true;
      }
		


      dPheromoneScore = pEdge->pOnEdge/(ANTS_PER_EDGE*TOTAL_CLOCK_TICK*20.0);
      dPopulationScore = pVertex->numAnts/(ANTS_PER_VERTEX);



      iPheromoneScore = dPheromoneScore * MAX_SCORE;
      iPopulationScore = dPopulationScore * MAX_SCORE;

	  
      // Bound scores.
      iPheromoneScore = (iPheromoneScore < 0.0) ? 0.0 : iPheromoneScore;
      iPheromoneScore = (iPheromoneScore > MAX_SCORE) ? MAX_SCORE : iPheromoneScore;
      iPopulationScore = (iPopulationScore < 0.0) ? 0.0 : iPopulationScore;
      iPopulationScore = (iPopulationScore > MAX_SCORE) ? MAX_SCORE : iPopulationScore;

    }//(adjV.at(iVertexIndex)!=pAnt->lastVertex||adjV.size()==1)
	
	
	
    vertexScoreList[i] =  (iPheromoneScore * fPheromoneWeight)+
      (iPopulationScore * fPopulationWeight) ;
	
    ScoreSum += static_cast<int>(vertexScoreList[i]);
	
  }//isample count loop
  
  
  
  //Roulette Wheel
  long needle = genrand_int32()%ScoreSum;

  //printf("at iTime %d, ScoreSum=%d, need= %d\n",iTime,ScoreSum,needle);

  
  for (int i = 0 ; i < iSampleCount ;++i){
    if (needle-vertexScoreList[i] > 0.0){
      needle-=static_cast<int>(vertexScoreList[i]);
    }
    else{
      newPos = pVertices[adjV.at(i)];	  adjIndex=i;
      break;
    }
  }


  return newPos;

}


edge* findEdge(int &adjIndex,
	       const vertex *v1,
	       const vertex *v2){

  
  edge* theEdge=NULL;  adjIndex=-1;
  
  for (int i = 0 ; i <v1->adj.size() ;++i){
    if (v1->adj.at(i)==v2->v_id){
      adjIndex=i;
      theEdge=pEdges[v1->edgeList.at(i)];
      //printf("[%d %d],size %d\n",adjIndex,theEdge->e_id,v1->adj.size());
      return theEdge;
    }
  }
  return theEdge;
}


vertex *chooseNextMoveRandomly(const vertex *oldPos,
			       const int nVertices){


  vertex *newPos=NULL;//dummy, will be changed later
  do {
    newPos=pVertices[genrand_int32()%nVertices];
  }while(newPos==oldPos);//move somewhere new
  return newPos;
}

void chooseNextMove(const int iTime,ant* pAnt,
		    bool edgeToUpdate[],
		    const int nVertices,const int nEdges){


  vertex* oldPos = pAnt->current;
  vertex* newPos = NULL;//undecided yet



  edge* theEdgeToNewPos=NULL;

  if (genrand_int32()%100<probAntMoveRandomly){

    newPos = chooseNextMoveRandomly(oldPos,nVertices);




    pAnt->next=newPos;//update the next move

    int adjIndex;	
    //if there's an edge to it 
    theEdgeToNewPos =findEdge(adjIndex,oldPos,newPos);
    if(theEdgeToNewPos!=NULL){	  
      pAnt->edgeToNext=theEdgeToNewPos;
    }//end findEdge()
	
  }
  else{


    //adj vertices to the old (current) position
    int adjIndex;
    newPos = chooseNextMoveHeuristiclyTVN(iTime,pAnt,
					  adjIndex,//index of the chosen stored here
					  edgeToUpdate,
					  nVertices,nEdges);

	

    theEdgeToNewPos=pEdges[oldPos->edgeList.at(adjIndex)];


    pAnt->next=newPos;//update the next move
    pAnt->edgeToNext=theEdgeToNewPos;

  }
}


void activateAnt(const int iTime,ant *pAnt,bool edgeToUpdate[],
		 const int nVertices,const int nEdges){

  if(genrand_int32()%100<probAntActivate){
    int probAntMove=100-pAnt->age;//the younger, the more likely will move
    if (probAntMove > AS_MAX) probAntMove = AS_MAX;
    if (probAntMove < AS_MIN) probAntMove = AS_MIN;

    if(genrand_int32()%100<probAntMove){
      chooseNextMove(iTime,pAnt,edgeToUpdate,nVertices,nEdges);
      pAnt->age++;
    }
  }
  
}


void moveAnts(const int iTime, ant *pAnt){

  if (pAnt->next!=NULL){

    //there's an edge between current and next (i.e, not random move)
    if (pAnt->edgeToNext!=NULL){
      addPheromone(pAnt->edgeToNext,EXPLORATORY_PHEROMONE);
    }


    //keep a history of current 
    pAnt->last=pAnt->current;	


    pAnt->current=pAnt->next;	//move the Ant to new vertex

    //reset the next Vertex
    pAnt->next=NULL;
    pAnt->edgeToNext=NULL;

    //update ant numbers
    pAnt->current->numAnts++;	
    pAnt->last->numAnts--; 
  }

#ifdef D2
#warning D2: output ant movements
  if (pAnt->last!=NULL)
    printf("ant %d [%d -> %d]\n",pAnt->a_id,pAnt->last->v_id,pAnt->current->v_id);
#endif


}

//Local optimizations
void growClique(vector<vertex *>&clique,
		const int nVertices,
		const int nEdges){
  
  //IT SEEMS that RIZZO's FIRST TIME EVEN GET A BIG CLIQUE
  //THUS NOT TOO MANY MORE VERTEX CAN GET ADDED 
  //HENCE NOT SO MANY ADDITIONAL EDGES,
  //MINE FIRST CLIQUE IS SMALL
  //HIS SEEMS TO BE STEADY, E.G later stage not 
  //finding any thing new so not many new added
  //mine keeps finding new 


  //CHECK WHY THE FIRST STAGE WHERE PHEROMONE is 0
  //in all edges.  In that case how does his score
  //sorted different than mine ? 

  //printf("Before Grow, clique size %d\n",clique.size());	
  //printVector(clique);


  vertex *pVertex;
  int startIndex=genrand_int32()%nVertices;
  int vertexIndex;
  for (int i = 0 ; i < nVertices ; ++i){
    vertexIndex=(startIndex+i)%nVertices ; //wrapper around  
    pVertex=pVertices[vertexIndex];

    if (pVertex->status==NOT_IN_OUTPUT){//if not in the clique
      //check to see adding it still make a clique
      if (!stillAClique(clique,pVertex)){ 
	//do nothing
      }
      else{
		
	//		//tvnAdded++;
	//		//printf("added %d\n",i);

	pVertex->status=IN_OUTPUT;
	clique.push_back(pVertex);
		

		
		
	//now adds some ants and pheromone to this vertex 
		
	//Jul14 , try apply openmp here, i.e take this function out 
	//create a vector that holds the newly added vertices
	updateVertex(pVertex);

	//printf("PREGET HERE\n");		
		
	//printf("GET HERE\n");		
	//add pheromone to all edges adj to this vertex
	//updatePheromoneAdjEdges(vertices.at(vertexIndex).edgeList,//edges adj to the vertex
	//								pOnEdges,
	//								REINFORCEMENT_PHEROMONE);
	//printf("SHOULDNOT GET HERE\n");		
	/*
	  assert(vertices.at(vertexIndex).edgeList.size()==vertices.at(vertexIndex).adj.size());
	*/

      }//end if can add vertex to clique

    }//end if != IN_OUTPUT
  }


  /*
    int aVAddedP=0;
    int edgesAddedP = 0 ;

    for (int i = 0 ; i < pOnEdges.size(); ++i){
    if(TEMP_EDGE[i]==1){edgesAddedP++;}
    }
    for (int i = 0 ; i < nVertices; ++i){
    if(TEMP_V[i]==1){aVAddedP++;}
    }
  */

  //this seems to be too big , edgesAddedP on the first run (stage 0)
  //FIND OUT WHY

  /*
    printf("After Grow, clique size %d, tvnAdded %d, pHasDTemp %d, %d\n",
    clique.size(),tvnAdded,edgesAddedP,aVAddedP);
  */


  assert(verifyClique(clique));//assert that it's a clique

  //testNumAnts("After Grow Clique",nVertices);

}

void setMimicGreedyClique(vector<vertex *>&theGreedyClique,
			  const int sv2Count){



#ifdef _DEBUG
  printf("Beg mimicgreedy, set candidate is %d, passingreedyClique is %d\n",
	 sv2Count,theGreedyClique.size());
#endif

  int iChecked =  0 ; 

  //BUG : for some reason my int DEG is much less than his on 1st iteration and so on
  //0th iteration is fine.  Uncomment block code below to see
  //so it seems that the setCandidate I picked is not in findClique is not as good
  //BIGCHANGE: I'll just get the biggest one , much different from his 
  //has to check later if something's wrong

  //  for (int i = 0 ; i < sv2Count; ++i){
  //	printf("v %d in sc has %d int deg\n",setCandidate.at(i).v_id,setCandidate.at(i).intDeg);
  //  }
  //end BUG

  vertex *pVertex, *pVertex2, *pVertex3;

  while (iChecked < sv2Count){
    //printf("iChecked %d\n",iChecked);
    int iThisDeg = -1;
    int iMaxIndex = -1;
    int iMax = -1;

    //printf("test 1\n");

    //6/24 12:15
    //OMP: can parallel this part
    for (int i = 0; i < sv2Count; i++) {

      pVertex = sv2[i];
      iThisDeg = pVertex->intDeg;
      if (iThisDeg > iMax)
	{
	  iMax = iThisDeg;
	  iMaxIndex = i;
	}

      //	printf("test 1 , v %d in sc has %d int deg\n",
      //		   setCandidate.at(i).v_id,setCandidate.at(i).intDeg);

    }//end for

    pVertex = sv2[iMaxIndex];


    if (pVertex->status==POSSIBLY_IN_OUTPUT){

      pVertex->status=IN_OUTPUT;
      pVertex->intDeg= -1; //prevent it from being selected again
      theGreedyClique.push_back(pVertex);
      iChecked++;


      //12:17 6/24 
      //OMP: probably can OMP this part 

      for (int i = 0; i < sv2Count; i++){
	//vertex anotherV(i) = setCandidate.at(i);
	pVertex2 = sv2[i];

	if (pVertex2->status == POSSIBLY_IN_OUTPUT){
	  if (!getDIMACSBinaryEdgeSwap(pVertex->v_id,pVertex2->v_id)){
			
            pVertex2->status = NOT_IN_OUTPUT;
            pVertex2->intDeg = -1;
            iChecked++;
			
            // Decrement internal degrees within candidateCount.
            for (int j = 0; j < sv2Count; j++){
	      pVertex3 = sv2[j];
              if (getDIMACSBinaryEdgeSwap(pVertex2->v_id,pVertex3->v_id)) {
                pVertex3->intDeg--;
              }
            }
	  }
	  else{          
            // In candidate clique, so decrement internal degree.
            pVertex2->intDeg--;
	  }
	}//	 		if (pVertex2->status == POSSIBLY_IN_OUTPUT){

      }	//for (int i = 0; i < sv2Count; i++){

    } // if (pVertex->status==POSSIBLY_IN_OUTPUT){

	
  }//end while()


#ifdef _DEBUG
  assert(verifyClique(theGreedyClique)==true);
  printf("mimicGreedy, the size %d\n", theGreedyClique.size());  
  printf("end mimic greedy\n");
#endif

}



/* //this seems to the the most time consuming part */
/* take too long */
/* //find out WHY */
void expandCandidate(const int dTime,
		     int &sv2Count,
		     const int nVertices,
		     const int nEdges,
		     const double DELTA){


  vertex *pVertex;
  
  int howMany = static_cast<int>(sv2Count*DELTA); 
  
  
  if (howMany+sv2Count> nVertices){ 
#ifdef _DEBUG
    printf("howMany+setCandidate.size()> nVertices");
#endif
	  
    howMany = nVertices-sv2Count;//max it can do 
  } 
  
#ifdef _DEBUG
  printf("Before expand, candSize %d, DELTA (DELTA) %f, Trying to add %d\n", 
	 sv2Count,DELTA,howMany); 
#endif
  
  
  int addedSoFar=0;

  //OMP: can I parallelizing this ?  using a for loop
  //for (int addedSoFar=0; addedSoFar < howMany ; )
  //probably not ,  since of the status 
  //also needs to take care of sv2Count (reduction)

  while(addedSoFar<howMany){ 
    int iThisDeg = - 1;
    int iMaxIndex = -1;//these 2 should be 
    int iMax = -1;//reset everytime  

    //OMP: probably can parallel this part
    //use MAX reduction
    for (int i = 0 ; i < nVertices ; ++i){ 
      pVertex = pVertices[i];
      if (pVertex->status==NOT_IN_OUTPUT){//consider this vertex i */
	iThisDeg = pVertex->intDeg;
	if (iThisDeg>iMax){ 
	  iMaxIndex = i;  iMax=iThisDeg;
	} 
      }
    }//end linear scan 
	
    /* 	//could be -1 if 1)THERE's BUG in my code or  */
    /* 	//2)the setCandidate is of size S where S + howMany > nVertices */
	
	
    if (iMaxIndex!=-1) {
	  
      pVertex = pVertices[iMaxIndex];
      pVertex->status=POSSIBLY_IN_OUTPUT;

      //6/24 12:13
      //OMP: cannot parallel this part, 
      //since pVertices[adj[i]] might be 
      //accessed concurrently 
      //have to made it in crit section

      //tvn DEBUG 6/24 12:13
      //why do you increase intDeg , the pVertex[iMaxIndex]
      //does not neccessarily IN the clique
      //I think this is OK, since in MimicGreedy,
      //there are some intDeg--
      for(int j=0;j<pVertex->adj.size();++j){ 
	pVertices[pVertex->adj.at(j)]->intDeg++; 
      }
      sv2[sv2Count]=pVertex;
      addedSoFar++;
      sv2Count++;
    } 
    else{
      printf("debug: setCandidate size %d, addesSoFar %d, iMaxIndex = %d\n" ,
	     sv2Count,addedSoFar,iMaxIndex); 
      printf ("could be -1 if 1)THERE's BUG in my code\n"); 
      //in this case, probably need to break 
      printf("2)the setCandidate is of size S where S + howMany > nVertices\n");
      //tvn change 12:06 6/24 
      //this should not happens since iMaxDeg is -1 , even if no connections
      //then it still 0
      printf("3)No outside vertices connects to the ones in the candidate set"); 
      assert(iMaxIndex!=-1); 
      break; 
    }
  }//end while(addedSoFar<howMany){ 
}




void findClique(vector<vertex *> &clique,
		const int dTime,
		const int nVertices,
		const int nEdges){
  

  //TVN: Todo , made this 0 to be consistent with Rizzo's since
  //his results makes it always 0
  //  double ALPHA=0.0;
  
  double ALPHA=double(alphaSlope * dTime) + alphaInit;;
  double BETA=double(betaSlope * dTime) + betaInit;
  double GAMMA =double(gammaGrowBySlope * dTime) + gammaGrowByInit;
  double DELTA =double(deltaGrowBySlope * dTime) + deltaGrowByInit;

#ifdef _DEBUG
  printf("dTime %d, ALPHA %f, BETA %f, GAMMA %f, DELTA %f\n",
	 dTime,ALPHA,BETA,GAMMA,DELTA);
#endif




  setScores(dTime,nVertices,nEdges,ALPHA,BETA);


#ifdef _DEBUG
  printDHasP("after set score",nEdges,false);
#endif

    
  qsort((void*)sv, (size_t)nVertices, sizeof(vertex*), compare_score);

#ifdef _DEBUG
  printf("first score %f, last score %f\n",
	 sv[0]->score,sv[nVertices-1]->score);
#endif
  


  int howMany = (int)floor(nVertices*GAMMA);



  vertex *pVertex,*pVertex2;
  int sv2Count = 0;  // # of vertices in candidate set.
  for (int i = 0; i < nVertices; i++)  {

    pVertex = sv[i];//sorted vertices (sv) according to score, 0 is the largest
    pVertex->intDeg = 0;

    //	printf("i: %d, testing v %d\n",i,pVertex->m_iID);
    if (i < howMany) {//top something PERCENT

      // In candidate solution.
      pVertex->status = POSSIBLY_IN_OUTPUT;
      // Set this vertex's internal degree.
      for (int j = 0; j < sv2Count; j++) {
        pVertex2 = sv2[j];
        if (getDIMACSBinaryEdgeSwap(pVertex->v_id, pVertex2->v_id)){
	  //assert (pVertex->v_id != pVertex2->v_id);
          pVertex->intDeg++;
          pVertex2->intDeg++;
	  //printf("pair %d %d\n",pVertex->m_iID,pVertex2->m_iID);
        }
      }
      //	  printf("added %d\n",pVertex->m_iID);
      sv2[sv2Count] = pVertex;
      sv2Count++;
    }
    else
      {
	pVertex->status = NOT_IN_OUTPUT;
	// Get degree of connection to candidate set.
	int iEdgeCount = pVertex->edgeList.size();
	assert(iEdgeCount==pVertex->adj.size());
	//OMP : cannot parallelizing since it asks for 
	//vertices that adj of this one.
	//i.e proc 1 doing vertex 1 , its adj is vertex 10 handling by proc 2
	//vertex 10 is read POSSIBLY_IN_OUTPUT by proc 1 but 
	//proc2 might changed it to something else
	//YEP, that's the whole reason why the 'status' is here
	for (int j = 0; j < iEdgeCount; j++) {
	  pVertex2 = pVertices[pVertex->adj.at(j)];
	  if (pVertex2->status == POSSIBLY_IN_OUTPUT) {
	    pVertex->intDeg++;
	  }
	}
      }
  }


#ifdef _DEBUG
  printf("GAMMA stuffs: GAMMA %f, candSize %d\n",GAMMA,sv2Count);
  printDHasP("after adding candidates",nEdges,false);
#endif


  expandCandidate(dTime,sv2Count,nVertices,nEdges,DELTA);

#ifdef _DEBUG
  printf("After expand, candSize %d\n", sv2Count);
#endif

  setMimicGreedyClique(clique,sv2Count);  //mimic greedyclique
}




#ifdef TESTEDGES//testing edges
int tEdges;
#endif




void localOptimization(vector<vector<vertex *> >&sol,
		       const int iStage,
		       const int iCycle,
		       const int nVertices, 
		       const int nEdges){
  
  vector<vertex *>localBestClique;

  findClique(localBestClique,
	     iStage*MAX_CYCLE+iCycle,//Rizzo's dtime 
	     nVertices,
	     nEdges);
	
	
#ifdef _DEBUG
  int iCliqueSizeNotUsed=localBestClique.size();
  int iNAntsNotUsed = vAnts.size();
#endif
	

	
  //update pheromone etc 
  growClique(localBestClique,nVertices,nEdges);
	
#ifdef _DEBUG
  printf("\nStage: %d, before grow, %d, ants %d, after : %d, ants %d, init %d\n\n",
	 iStage,iCliqueSizeNotUsed,iNAntsNotUsed,
	 localBestClique.size(),vAnts.size(),
	 sol.at(0).size());
	
  printDHasP("after grow Clique",nEdges,false);
#endif
	
#ifdef _DEBUG	
  assert(verifyClique(localBestClique)==true);
#endif
	
	
	
  //shuffleAnts(vertices,ants,antsOnVertices,pOnEdges,nVertices);
	
	
	
  sol.push_back(localBestClique);//remove when done, this is seq part

#ifdef TESTEDGES
  //  printf("Rizzo -> Stage %d, Edge %d = %d , last update %d\n",
  //		 iStage,tEdges,pEdges[tEdges]->pOnEdge,pEdges[tEdges]->lastUpdate);
  printf("%d ",pEdges[tEdges]->pOnEdge);
#endif
}




//localOptimizePAR
void syncCP(VertexInfo mVertices[],const int sIndex,const int nVertices,const int SYNC_OPT){
  if (SYNC_OPT==SYNC_TO_mVERTICES){
    for (int i = 0 ; i < nVertices ; ++i){
      mVertices[sIndex+i].viStatus=pVertices[i]->status;
      mVertices[sIndex+i].viIntDeg=pVertices[i]->intDeg;
      mVertices[sIndex+i].viScore=pVertices[i]->score;
      mVertices[sIndex+i].viNumAnts=pVertices[i]->numAnts;
    }
  }
  else{
#ifdef PSYNC
#warning "_PSYNC: OMP on syncCP()"
#pragma omp parallel for
#endif
    for (int i = 0 ; i < nVertices ; ++i){
      pVertices[i]->status=mVertices[sIndex+i].viStatus;
      pVertices[i]->intDeg=mVertices[sIndex+i].viIntDeg;
      pVertices[i]->score=mVertices[sIndex+i].viScore;
      pVertices[i]->numAnts=mVertices[sIndex+i].viNumAnts;
    }
	
  }
}

void expandCandidatePAR(VertexInfo mVertices[],
			vertex** priv_sv,  //BUG: pass by pointer, make sure values changed
			vertex** priv_sv2,
			const int sIndex,
			const int dTime,
			int &sv2Count,
			const int nVertices,
			const int nEdges,
			const double DELTA){
  
  
  vertex *pVertex;
  
  int howMany = static_cast<int>(sv2Count*DELTA); 
  
  
  if (howMany+sv2Count> nVertices){  howMany = nVertices-sv2Count;} 
  
  int addedSoFar=0;

  while(addedSoFar<howMany){ 
    int iThisDeg = -1;
    int iMaxIndex = -1;//these 2 should be 
    int iMax = -1;//reset everytime  

    //OMP: probably can parallel this part
    //use MAX reduction
    for (int i = 0 ; i < nVertices ; ++i){ 
      pVertex = pVertices[i];
      // 	  if (pVertex->status==NOT_IN_OUTPUT){//consider this vertex i */
      if (mVertices[sIndex+pVertex->v_id].viStatus==NOT_IN_OUTPUT){//consider this vertex i */
	iThisDeg = mVertices[sIndex+pVertex->v_id].viIntDeg;	//iThisDeg = pVertex->intDeg;
	if (iThisDeg>iMax){ 
	  iMaxIndex = i;  iMax=iThisDeg;
	} 
      }
    }//end linear scan 
	
	
    if (iMaxIndex!=-1) {
	  
      pVertex = pVertices[iMaxIndex]; 
      mVertices[sIndex+pVertex->v_id].viStatus=POSSIBLY_IN_OUTPUT;//pVertex->status=POSSIBLY_IN_OUTPUT;

      //tvn DEBUG 6/24 12:13
      //why do you increase intDeg , the pVertex[iMaxIndex]
      //does not neccessarily IN the clique
      //I think this is OK, since in MimicGreedy,
      //there are some intDeg--
      for(int j=0;j<pVertex->adj.size();++j){ 
	//pVertices[pVertex->adj.at(j)]->intDeg++; 
	mVertices[sIndex+pVertex->adj.at(j)].viIntDeg++; //recheck this,  should be equiv with the line above
      }
      priv_sv2[sv2Count]=pVertex;
      addedSoFar++;	  sv2Count++;
    } 
    else{
      printf("debug: setCandidate size %d, addesSoFar %d, iMaxIndex = %d\n" ,
	     sv2Count,addedSoFar,iMaxIndex); 
      printf ("could be -1 if 1)THERE's BUG in my code\n"); 
      //in this case, probably need to break 
      printf("2)the setCandidate is of size S where S + howMany > nVertices\n");
      //tvn change 12:06 6/24 
      //this should not happens since iMaxDeg is -1 , even if no connections
      //then it still 0
      printf("3)No outside vertices connects to the ones in the candidate set"); 
      assert(iMaxIndex!=-1); 
      break; 
    }
  }//end while(addedSoFar<howMany){ 
}



void setMimicGreedyCliquePAR(vector<vertex *>&theGreedyClique,
			     VertexInfo mVertices[],							 
			     vertex** priv_sv2,
			     const int sIndex,
			     const int sv2Count){



  int iChecked =  0 ; 


  vertex *pVertex, *pVertex2, *pVertex3;


  while (iChecked < sv2Count){
    //printf("iChecked %d\n",iChecked);
    int iThisDeg = -1;	int iMaxIndex = -1;	int iMax = -1;


    for (int i = 0; i < sv2Count; i++) {

      pVertex = priv_sv2[i];
      iThisDeg = mVertices[sIndex+pVertex->v_id].viIntDeg; //iThisDeg = pVertex->intDeg;
      if (iThisDeg > iMax)  {
        iMax = iThisDeg; iMaxIndex = i;
      }

    }//end for


    pVertex = priv_sv2[iMaxIndex];

    //if (pVertex->status==POSSIBLY_IN_OUTPUT){
    if (mVertices[sIndex+pVertex->v_id].viStatus==POSSIBLY_IN_OUTPUT){

      mVertices[sIndex+pVertex->v_id].viStatus=IN_OUTPUT;//pVertex->status=IN_OUTPUT;
      mVertices[sIndex+pVertex->v_id].viIntDeg= -1; //pVertex->intDeg= -1; //prevent it from being selected again
      theGreedyClique.push_back(pVertex);
      iChecked++;


      for (int i = 0; i < sv2Count; i++){
	//vertex anotherV(i) = setCandidate.at(i);
	pVertex2 = priv_sv2[i];

	//		if (pVertex2->status == POSSIBLY_IN_OUTPUT){
	if (mVertices[sIndex+pVertex2->v_id].viStatus == POSSIBLY_IN_OUTPUT){
	  if (!getDIMACSBinaryEdgeSwap(pVertex->v_id,pVertex2->v_id)){
			
	    // pVertex2->status = NOT_IN_OUTPUT;
	    mVertices[sIndex+pVertex2->v_id].viStatus=NOT_IN_OUTPUT;
	    //            pVertex2->intDeg = -1;
	    mVertices[sIndex+pVertex2->v_id].viIntDeg=-1;

            iChecked++;
			
            // Decrement internal degrees within candidateCount.
            for (int j = 0; j < sv2Count; j++){
	      pVertex3 = priv_sv2[j];
              if (getDIMACSBinaryEdgeSwap(pVertex2->v_id,pVertex3->v_id)) {
                //pVertex3->intDeg--;
		mVertices[sIndex+pVertex3->v_id].viIntDeg--;
              }
            }
	  }
	  else{          
            // In candidate clique, so decrement internal degree.
            //pVertex2->intDeg--;
	    mVertices[sIndex+pVertex2->v_id].viIntDeg--;

	  }
	}//	 		if (pVertex2->status == POSSIBLY_IN_OUTPUT){

      }	//for (int i = 0; i < sv2Count; i++){
    } // if (pVertex->status==POSSIBLY_IN_OUTPUT){
  }//end while()

}


void findCliquePAR(vector<vertex *> &clique,
		   VertexInfo mVertices[],
		   vertex **priv_sv,
		   vertex **priv_sv2,
		   const int sIndex,
		   const int dTime,
		   const int nVertices,
		   const int nEdges){
  


  double GAMMA =double(gammaGrowBySlope * dTime) + gammaGrowByInit;
  double DELTA =double(deltaGrowBySlope * dTime) + deltaGrowByInit;

  int howMany = (int)floor(nVertices*GAMMA);



  vertex *pVertex,*pVertex2;
  int sv2Count = 0;  // # of vertices in candidate set.

  for (int i = 0; i < nVertices; i++)  {

    pVertex = priv_sv[i];//sorted vertices (sv) according to score, 0 is the largest

    mVertices[sIndex+pVertex->v_id].viIntDeg=0;  //pVertex->intDeg = 0;


    if (i < howMany) {//top something PERCENT
      // In candidate solution.
      mVertices[sIndex+pVertex->v_id].viStatus=POSSIBLY_IN_OUTPUT;	  //pVertex->status = POSSIBLY_IN_OUTPUT;
      // Set this vertex's internal degree.
      for (int j = 0; j < sv2Count; j++) {
        pVertex2 = priv_sv2[j];
        if (getDIMACSBinaryEdgeSwap(pVertex->v_id, pVertex2->v_id)){
	  mVertices[sIndex+pVertex->v_id].viIntDeg++;   //pVertex->intDeg++;    
	  mVertices[sIndex+pVertex2->v_id].viIntDeg++; //pVertex2->intDeg++;
        }
      }
      priv_sv2[sv2Count] = pVertex;
      sv2Count++;
    }
    else{
      mVertices[sIndex+pVertex->v_id].viStatus=NOT_IN_OUTPUT; //pVertex->status = NOT_IN_OUTPUT;

      // Get degree of connection to candidate set.
      for (int j = 0; j <  pVertex->edgeList.size(); j++) {
        pVertex2 = pVertices[pVertex->adj.at(j)];
	//        if (pVertex2->status == POSSIBLY_IN_OUTPUT) {
        if (  mVertices[sIndex+pVertex2->v_id].viStatus ==POSSIBLY_IN_OUTPUT) {

	  mVertices[sIndex+pVertex->v_id].viIntDeg++; //pVertex->intDeg++;
        }
      }
    }
  }




  expandCandidatePAR(mVertices,priv_sv,priv_sv2,sIndex,dTime,sv2Count,nVertices,nEdges,DELTA);
  setMimicGreedyCliquePAR(clique,mVertices,priv_sv2,sIndex,sv2Count);  //mimic greedyclique

}


void singleLocalOptimization(vector<vertex *> &localBestClique,
			     VertexInfo mVertices[],
			     const int threadID,
			     const int dTime,
			     const int nVertices, 
			     const int nEdges){
  
  


  vertex **priv_sv,   **priv_sv2 ;
  priv_sv = new vertex*[nVertices];
  priv_sv2 = new vertex*[nVertices];

  int sIndex = threadID*nVertices ;

  //TVN: Todo , made this 0 to be consistent with Rizzo's since
  //his results makes it always 0
  double ALPHA=0.0;  //double ALPHA=double(alphaSlope * dTime) + alphaInit;;

  double BETA=double(betaSlope * dTime) + betaInit;


  //todo :  have setscore Done outside
  //then sort to sv
  //then copy priv_sv[i]=sv[i];

  setScores(dTime,nVertices,nEdges,ALPHA,BETA);

  

  //to avoid this , make it a global var e.g,  vertex ** priv_sv[nthreads]
  for (int i = 0 ; i < nVertices ;++i){	priv_sv[i]=pVertices[i]; }

  qsort((void*)priv_sv, (size_t)nVertices, sizeof(vertex*), compare_score);

  vertex *tempVertex;

  //mangling the score position
  if (threadID>0){
    for (int i = 0 ; i < nVertices ;++i){
      if (i<nVertices-1){
	if(priv_sv[i]->score == priv_sv[i+1]->score&&genrand_int32()%100<=50){
	  //printf("thread id %d, switch priv_sv\n",threadID);
	  //the tvn BUBBLE switch style ...
	  tempVertex=priv_sv[i];
	  priv_sv[i]=priv_sv[i+1];
	  priv_sv[i+1]=tempVertex;
	  i++;
	}
      }
    }
  }
  


  /*  printf("thread %d ,vid %d , first score %f, vid %d ,last score %f\n",
      threadID,
      priv_sv[0]->v_id,
      priv_sv[0]->score,
      priv_sv[nVertices-1]->v_id,
      priv_sv[nVertices-1]->score);
  */



  //imitating the operation of privating VertexInfo
  syncCP(mVertices,sIndex,nVertices,SYNC_TO_mVERTICES);


  findCliquePAR(localBestClique,mVertices,priv_sv,priv_sv2,
		sIndex,dTime,nVertices,nEdges);

  


	

  //growCliquePAR(localBestClique,mVertices,sIndex,nVertices,nEdges);	//update pheromone etc 

  //shuffleAnts(vertices,ants,antsOnVertices,pOnEdges,nVertices);

}


void parallelLocalOpt(vector<vector<vertex *> >&sol,
		      VertexInfo mVertices[],
		      const int iStage,
		      const int iCycle,
		      const int nVertices, 
		      const int nEdges){


  //init values
  vector<vector<vertex *> > localSol;
  for (int threadID = 0 ; threadID < nthreads ; ++threadID){
    vector<vertex *>localBestClique;
    localSol.push_back(localBestClique);
  }


  int dTime = iStage*MAX_CYCLE+iCycle;

#ifdef _OPENMP
#warning "PLOCAL: OMP localoptimizing() on"
#pragma omp parallel for
#endif
  for (int threadID = 0 ; threadID < nthreads ; ++threadID){
    singleLocalOptimization(localSol.at(threadID),mVertices,
			    threadID,dTime,nVertices,nEdges);
  }


  //done with extract a (small) clique
  int bestLocalIndex;
  for (int threadID =0 ; threadID < nthreads; ++threadID){  //scan for best one
    //printf("thread %d clique size %d\n",threadID,localSol.at(threadID).size());
    if (threadID==0){
      bestLocalIndex=threadID;
    }
    else{
      if (localSol.at(threadID).size()>localSol.at(bestLocalIndex).size()){
	//printf("Local Opt: threadID %d found a clique with size %d, better than threadID%d with  size %d\n",
	//threadID,localSol.at(threadID).size(),bestLocalIndex,bestLocalSize);
	bestLocalIndex=threadID;
      }
    }
  }
  

  //sync the vertexInfo
  syncCP(mVertices,bestLocalIndex*nVertices,nVertices,SYNC_TO_pVERTICES);



  //update pheromone etc 
  //  printf("before grow %d",localSol.at(bestLocalIndex).size());
  growClique(localSol.at(bestLocalIndex),nVertices,nEdges);
  //  printf(", after grow %d\n",localSol.at(bestLocalIndex).size());

  assert(verifyClique(localSol.at(bestLocalIndex))==true);



  //push back best sol
  sol.push_back(localSol.at(bestLocalIndex));

}

