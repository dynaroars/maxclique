

vertex *chooseNextMoveHeuristiclyTVN(const int iTime,const ant *pAnt,
									 int &adjIndex,//stores edge id to new pos
									 bool edgeToUpdate[],
									 const int nVertices,const int nEdges){

  vertex* oldPos = pAnt->current;

  vertex* newPos=NULL;//dummy, will be changed later
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




int runCycleLoop(const int iStage,
				 const int nVertices,
				 const int nEdges,
				 bool edgeToUpdate[]){

  int iCycle = 0 ; 
  for (;iCycle<MAX_CYCLE;++iCycle){



	int iTime = iStage * MAX_CYCLE + iCycle ;

#ifdef _DEBUG
	printf("Cycle %d, iTime: %d\n",iCycle,iTime);
	antMove = 0 ;antActivated=0; antUnactivated=0; antMoveH=0;antMoveR=0;antStay=0;
#endif


	int antSize = vAnts.size();


	//activate and determine move for ants


#ifdef _OPENMP
#warning "_opemp: OMP activateAnt() on"
#pragma omp parallel for
#endif
	for (int iAnt=0; iAnt<antSize ; ++iAnt){
	  activateAnt(iTime,vAnts.at(iAnt),edgeToUpdate,nVertices,nEdges);
	}




#ifdef _OPENMP
#warning "_openmp: OMP for evaporatePheromone Edges on"
#pragma omp parallel for
#endif
	for (int iE=0; iE < nEdges; ++iE){
	  if (edgeToUpdate[iE]==true){
		evaporatePheromone(pEdges[iE],
						   ((iTime - pEdges[iE]->lastUpdate)*EVAPORATION_RATE));
		pEdges[iE]->lastUpdate = iTime;
		edgeToUpdate[iE]=false;//reset
	  }
	}//	for (int iE=0; iE < nEdges; ++iE){




		
#ifdef DHASP
	printf("stage %d, after ant decision, before running moveAnts , localOpt\n",iStage);
	printDHasP("stage loop" , nEdges,false);
#endif


	//move the ants, sequentially
	for (int iAnt=0;iAnt<vAnts.size();++iAnt){
	  moveAnts(iTime,vAnts.at(iAnt));
	}

	
  }//end of iCycle loop

  return iCycle;
}




