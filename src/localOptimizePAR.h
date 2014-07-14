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

