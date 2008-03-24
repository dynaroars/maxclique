#include <cstdlib>
#include <iostream>
#include <cassert>
using namespace std;

int main(){

  int n=4;

  int sumi=0;
  long suml=0;
  double sumd=0;

  printf("before, sumi %d , suml %f, sumd %f \n",sumi,suml,sumd);

  for (int i = 0 ;i < n ; ++i){
	double randN = i/10;
	printf("%f, ",randN);
	sumi+=(int)randN;
	suml+=(long)randN;
	sumd+=(double)randN;
  }

  int mysum=0;
  for(int i = 0 ;i < n ; ++i){
	mysum+=i;
  }
  cout << "\n" << mysum << endl;

  printf("\nafter, sumi %d , suml %f, sumd %f, mysum %d\n",mysum,sumi,suml,sumd,mysum);

  assert((double)suml==sumd);

  cout <<suml-((long)sumi-2)<<endl;
  
  return 0;
}
