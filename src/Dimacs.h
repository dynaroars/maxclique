// genbin.h code block
/* If you change MAX_NR_VERTICES, change MAX_NR_VERTICESdiv8 to be 
the 1/8th of it */
#define MAX_NR_VERTICES      5000  //tvn: probably need to be BIGGER (todo)
#define MAX_NR_VERTICESdiv8  625 
#define BOOL                 char
#define MAX_PREAMBLE         10000
static char Preamble[MAX_PREAMBLE];
// end genbin.h code block



BOOL Bitmap[MAX_NR_VERTICES][MAX_NR_VERTICESdiv8];

// Dimacs graph format translator to and from a binary, more efficient
// format. Written by Tamas Badics (badics@rutcor.rutgers.edu), 
// using the technique of Marcus Peinado, Boston University. 
// 
// Corrected to number nodes from 1 to n (not 0 to n-1)
// MT

// masks[i] = pow(2, i);
char masks[8];


void initDimacs(){
  masks[0] = 0x01;
  masks[1] = 0x02;
  masks[2] = 0x04;
  masks[3] = 0x08;
  masks[4] = 0x10;
  masks[5] = 0x20;
  masks[6] = 0x40;
  masks[7] = 0x80;
}



int getDIMACSBinaryParams(int &nVertices , int &nEdges)
    /* getting nVertices and Nr_edge from the preamble string, 
       containing Dimacs format "p ??? num num" */
{
  char  c, *tmp;                           // Holds problem description
  char* pp = Preamble;                     // The graph preamble
  int   stop = 0;                          // Flag used to stop processing
  tmp = (char *)calloc(100, sizeof(char));
  nVertices = nEdges = 0;
  
  // Loop through the first characters of each line in the preamble.
  while (!stop && (c = *pp++) != '\0'){
    switch (c)
      {
      case 'c':
        // Comment. Skip this line by finding the end.
        while ((c = *pp++) != '\n' && c != '\0');
        break;

      case 'p':
        // Found the problem description line!  Read the values and stop.
        sscanf(pp, "%s %d %d\n", tmp, &nVertices, &nEdges);
        stop = 1;
        break;
        
      default:
        break;
      }
  }
  free(tmp);
  if (nVertices == 0 || nEdges == 0)    return 0;  /* error */
  else return 1;
}


BOOL getDIMACSBinaryEdge(int i, int j){ //i>=j

  char masks[8];
  int byte, bit;
  char mask;
  
  bit  = 7-(j & 0x00000007);  // bit  = 7 - (j % 8);
  byte = j >> 3;              // byte = j / 8;
  
  mask = masks[bit];
  return( (Bitmap[i][byte] & mask)==mask );
}



/* Same as above except it accepts parameters in any order.  */
bool getDIMACSBinaryEdgeSwap(int i,  int j){
  int byte, bit;  char mask;

  if (i < j)  {
    int t = i;
    i = j;
    j = t;
  }

  bit  = 7-(j & 0x00000007);  // bit  = 7 - (j % 8);
  byte = j >> 3;              // byte = j / 8;
  
  mask = masks[bit];
  return( (Bitmap[i][byte] & mask)==mask );
}



/* The idea here is to see if the bitmap changes over time. */
void showBitmap(){
   for (int i = 0; i < 50; i++)   {
     for (int j = 0; j <= (i / 8); j++)     {
        printf("%c", Bitmap[i][j]);
     }
   }
   printf("\n\n");
}



void readDIMACSBinaryFormat(char *file, int &nVertices, int &nEdges){
  initDimacs();

  int i, length = 0;
  FILE *fp;
  
  // SUNS: fp=fopen(file, "r"))   The extra 'b' opens it in binary (untranslated)
  //                              mode as opposed to text (translated) mode.  See
  //                              MSDN help for fopen for details.

  if ( (fp=fopen(file,"rb"))==NULL )
    { printf("ERROR: Cannot open infile\n"); scanf("%d", &i); }

  if (!fscanf(fp, "%d\n", &length))
    { printf("ERROR: Corrupted preamble.\n"); scanf("%d", &i); }

  if(length >= MAX_PREAMBLE)
    { printf("ERROR: Preamble is too long.\n"); scanf("%d", &i); }
       
  fread(Preamble, 1, length, fp);
  Preamble[length] = '\0';
  

  if (!getDIMACSBinaryParams(nVertices,nEdges))
      { printf("ERROR: Corrupted preamble.\n"); scanf("%d", &i); }


  // Load the edge info into Bitmap.
  for ( i = 0 ; i < nVertices && fread(Bitmap[i], 1, (int)((i + 8)/8), fp)  ; i++ );


  


  
  if( feof( fp ) )      {
	printf( "End of file" );
  }
  
  if( ferror( fp ) )      {
	printf( "Read error" );
  }
  
  fclose(fp);
}
