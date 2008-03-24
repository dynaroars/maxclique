#include <iostream>
#include <cassert>
#include <cstdlib>
#include <ctime>
using namespace std;


const int arraySize = 20;
int mA[arraySize];


void sumUp(int mB[],const int mA[],const int arraySize){
  mB[0]=mA[0];
  printf("%d, ",mB[0]);
  for(int i = 1 ; i < arraySize ; ++i){
	mB[i]=mA[i]+mB[i-1];
	if (i!=0 && i%10==0) printf("\n");
	printf("%d, ",mB[i]);
  }
}


int main(){
  srand(time(0));
  for(int i = 0 ; i < arraySize ; ++i){
	mA[i]=rand()%100;

	if (i!=0 && i%10==0) printf("\n");
	printf("%d, ",mA[i]);
  }

    printf("\n\n");
  int mB[arraySize];
  sumUp(mB,mA,arraySize);
  
  printf("\n\n");
  double topPercent = 20.0;
  int topNum=(int)(20.0 * arraySize)/100;
  printf("Getting top %d numbers\n",topNum);

  for(int i = 0 ; i < arraySize ; ++i){
	if (topNum){
	  if (mB[i]>=rand()%mB[arraySize-1]){
		printf("%d, ",mA[i]);
		topNum--;
	  }
	}
  }

  printf("\n");
  return 0;
}



