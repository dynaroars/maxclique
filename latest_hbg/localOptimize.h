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



