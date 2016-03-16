// Param settings
const auto MAX_STAGE = 25;
const auto MAX_CYCLE = 10;
const auto TOTAL_CLOCK_TICK = MAX_STAGE*MAX_CYCLE;

const auto antMultipler = 10;  //nAnts = nvertices * antMultipler
const auto nExtraAnts = 12 ; //# of ants to add to newly found vertex
const auto antFactor = 1.25 ;


//this is a bit too much,  90% of the ants population concentrated here
const auto probAntsOnGreedyClique = 90;  
const auto probAntActivate = 75;
const auto AS_MAX = 99;
const auto AS_MIN = 40;
const auto probAntMoveRandomly = 10;
const auto fPheromoneWeight = 1.0;
const auto fPopulationWeight = 5.0;
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

const auto EXPLORATORY_PHEROMONE  = 25;       // Amount dropped in a usual ant move.
const auto REINFORCEMENT_PHEROMONE  = EXPLORATORY_PHEROMONE * 10;  // Amount dropped when adding ant(s) to a vertex.
const auto MIN_EDGE_PHEROMONE  = 0;
const auto MAX_EDGE_PHEROMONE   = 30000;
const auto EVAPORATION_RATE     = 5;       // Subtracted from all edges at end of each minor loop.


double  ANTS_PER_EDGE = 0.0;
double  ANTS_PER_VERTEX = 0.0;
const auto MAX_SCORE = 1000.0;

/* sample paramters */
const auto VA_SAMPLE_PERCENTAGE          = 1.0;    // % of neighbors to score.
const auto VA_MIN_SAMPLE_SIZE            = 5;       // Min # of neighbors to score - overrules the %.
const auto VA_MAX_SAMPLE_SIZE            = 1000;     // Max # of neighbors to score - overrules the %.


/* Shuffling parameters */
const auto probVerticesShuffled = 40 ; 
const auto probPheromoneReductionOnEdges = 10; 
const auto amountPheromoneReduction = 3 ; //percent

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


const auto SYNC_TO_mVERTICES = 0;
const auto SYNC_TO_pVERTICES = 1;



//////////// GLOBAL PARAMS ///////////////
time_t seed_t;
vertex **pVertices , **sv, **sv2;
edge **pEdges;
vector <ant *> vAnts;  


vector<vector<vertex *> >sol;

int MAX_ANTS ; //of Ants * 1.25
int nthreads;
