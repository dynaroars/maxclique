

// Param settings

const int MAX_STAGE = 25;
const int MAX_CYCLE = 10;
const int TOTAL_CLOCK_TICK = MAX_STAGE*MAX_CYCLE;

const int antMultipler = 10;  //nAnts = nvertices * antMultipler
const int nExtraAnts = 12 ; //# of ants to add to newly found vertex
const double antFactor = 1.25 ;


//this is a bit too much,  90% of the ants population concentrated here
const int probAntsOnGreedyClique = 90;  
const int probAntActivate = 75;
const int AS_MAX = 99;
const int AS_MIN = 40;
const int probAntMoveRandomly = 10;
const double fPheromoneWeight = 1.0;
const double fPopulationWeight = 5.0;
double factorC = 10; //will change (from something_larger_thanAB to AB)
double alphaInit=10.0;
double alphaFinal=5.0;
double alphaSlope = (alphaFinal-alphaInit)/TOTAL_CLOCK_TICK; //will change (from 10 to 5)

double betaInit=.1;
double betaFinal=1.0;
double betaSlope=(betaFinal-betaInit)/TOTAL_CLOCK_TICK;

double gammaGrowByInit = .10 ; //will change (from .10 to .25)
double gammaGrowByFinal = .23 ; //will change (from .10 to .23)
double gammaGrowBySlope =  (gammaGrowByFinal - gammaGrowByInit)/TOTAL_CLOCK_TICK;
double deltaGrowByInit = 0.0 ; 
double deltaGrowByFinal = 1.3 ; //DELTA
double deltaGrowBySlope = (deltaGrowByFinal - deltaGrowByInit)/TOTAL_CLOCK_TICK;

const int EXPLORATORY_PHEROMONE  = 25;       // Amount dropped in a usual ant move.
const int REINFORCEMENT_PHEROMONE  = EXPLORATORY_PHEROMONE * 10;  // Amount dropped when adding ant(s) to a vertex.
const int MIN_EDGE_PHEROMONE  = 0;
const int MAX_EDGE_PHEROMONE   = 30000;
const int EVAPORATION_RATE     = 5;       // Subtracted from all edges at end of each minor loop.


double  ANTS_PER_EDGE = 0.0;
double  ANTS_PER_VERTEX = 0.0;
const double MAX_SCORE = 1000.0;

/* sample paramters */
const double   VA_SAMPLE_PERCENTAGE          = 1.0;    // % of neighbors to score.
const int     VA_MIN_SAMPLE_SIZE            = 5;       // Min # of neighbors to score - overrules the %.
const int     VA_MAX_SAMPLE_SIZE            = 1000;     // Max # of neighbors to score - overrules the %.


/* Shuffling parameters */
const int probVerticesShuffled = 40 ; 
const int probPheromoneReductionOnEdges = 10; 
const int amountPheromoneReduction = 3 ; //percent

/* end Shuffling parameters */



enum V_STATUS{ NOT_IN_OUTPUT,  POSSIBLY_IN_OUTPUT,  IN_OUTPUT};



struct vertex{
  int v_id; 
  vector<int>adj;
  vector<int>edgeList;
  int status; 
  int intDeg; 
  double score;//should remove this ?  only use in findClique
  int numAnts; //seems like Rizzo's doesn't maintain a list of ants on vertice[i];
  //seems like that's okay, he doesn't need to 
};

struct edge{
  int e_id;
  int pOnEdge;
  int lastUpdate;
  int pointA;
  int pointB;
};

struct ant{
  int a_id;
  int age;

  vertex *current;
  vertex *next;
  vertex *last;
  edge *edgeToNext;
};


struct VertexInfo{
  int viStatus;  int viIntDeg;
  double viScore;  int viNumAnts;
};


const int SYNC_TO_mVERTICES = 0;
const int SYNC_TO_pVERTICES = 1;



//////////// GLOBAL PARAMS ///////////////
time_t seed_t;
vertex **pVertices , **sv, **sv2;
edge **pEdges;
vector <ant *> vAnts;  


vector<vector<vertex *> >sol;

int MAX_ANTS ; //of Ants * 1.25


int nthreads;

///////////////////////////////////////////

void evaporatePheromone(edge *theEdge,const int theAmount){
  //printf("original %d, deduct %d\n",pAmount,theAmount);
  theEdge->pOnEdge -= theAmount;
  if (theEdge->pOnEdge < MIN_EDGE_PHEROMONE){
	theEdge->pOnEdge = MIN_EDGE_PHEROMONE;//Todo: use the a= a>b?i:j  style
  }
}

void addPheromone(edge *theEdge,const int theAmount){
  
  theEdge->pOnEdge += theAmount;
  if (theEdge->pOnEdge > MAX_EDGE_PHEROMONE){
	theEdge->pOnEdge = MAX_EDGE_PHEROMONE;//Todo: use the a= a>b?i:j  style
  }
}

void updatePheromoneAdjEdges(const vector<int> &theAdjEdges,
							 const int amount){
  //  printf("adjsize %d\n", theAdjEdges.size());

  for (int i = 0 ; i < theAdjEdges.size() ;++i){
	//printf("amount %d\n", amount);
	assert(theAdjEdges.at(i)==pEdges[theAdjEdges.at(i)]->e_id);
	addPheromone(pEdges[theAdjEdges.at(i)], amount);
	//printf("Edge ID %d now contains %d\n",theAdjEdges.at(i),
	//pOnEdges.at(theAdjEdges.at(i)));
  }
}



void cleanUp(const int nVertices, const int nEdges){
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



