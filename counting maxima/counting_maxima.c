#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

int N = 10;
int K = 2;
int d = 0;
int Nsample = 1001;

int array_to_dec(int n, int v[n]);
int Global_Max(int d, double f[d]);
int Local_Maximum( int S, double f[d]);
int Count_Maxima( double f[d], double fmax);
void Landscape(int N, int K, int d, double fit[d]);
void localMaximaNumber();
void localMaximaNumber_graph();

int main(){
 localMaximaNumber_graph();
 return 0;
}

void localMaximaNumber_graph(){
  FILE *fptr;
  fptr = fopen("/home/larissa/local_maxima_avg_vs_k.txt", "w");
  fclose(fptr);
  N = 10;
  d = pow(2, N);
  localMaximaNumber();
  N = 12;
  d = pow(2, N);
  localMaximaNumber();
  N = 15;
  d = pow(2, N);
  localMaximaNumber();
  return;
}
int array_to_dec(int n, int v[n]){
  int i, dec;
  dec = 0;
  for (i = 0; i < n; i++){
    dec = dec + v[i]*pow(2,n-i-1);
    }
  return dec;
}
void Landscape(int N, int K, int d, double fit[d]){
  int dim = pow(2, K+1);
  int i, j, k, l; //contadores
  double phi[dim][N], avg[N];
  int num, v[N], B, neighboors[N+1], Matrix[d][N+1];

	//Landscape Matrix
	for (i = 0; i < d; i++){
	  B = i;
	  for (j = 0; j < N; j++){
		 v[N-j-1] = B%2;
		 B = B/2;
	  }
	  Matrix[i][0] = i ;
	  for (j = 0; j < N; j++){
		 for (l = 0; l< K+1; l++){
		   if(j+l < N){neighboors[l] = v[j+l];}
		   else {neighboors[l] = v[j+l-N];}
		 }
		 Matrix[i][j+1] = array_to_dec(K+1, neighboors);
	  }
	}
	//Fitness Components
	for (i = 0; i < dim; i++){
	  for (j = 0; j < N; j++){
		 phi[i][j] = rand();
		 phi[i][j] = phi[i][j]/RAND_MAX;

	  }
	}
	//Fitness
	for (i=0; i < d; i++){
	  fit[i] = 0;
	  for (j = 0; j < N; j++){
		 num = Matrix[i][j+1];
		 fit[i] = fit[i] + phi[num][j];
	  }

	  fit[i] = fit[i]/N;
    // printf("fit = %.2f \n", fit[i]);
	}
	//printf("%d \n", d);

//	for (i = 0; i < d; i++){
//	printf("%d, %.2f \n ",i, fit[i]);
//	}

	return;
}
int Global_Max(int d, double f[d]){
  int k = 0;
  double max = 0.;
  int i;
  for (i = 0; i < d; i++){
    if (f[i] > max){
      max = f[i];
      k = i;
    }
  }
  return k;
}
int Local_Maximum( int S, double f[d]){
  int i, state[N], dec;
  double fn, fs;
  int flag = 0;
  fs = f[S];
  for (i = 0; i < N; i++){
    state[N-i-1] = S%2;
    S = (S-state[N-i-1])/2;
  }
  for (i = 0; i < N; i++){
    if(i > 0){
      state[i-1] = 1 - state[i-1];
      state[i] = 1 - state[i];
      }
    else{state[i] = 1 - state[i];}
    dec = array_to_dec(N, state);
    fn = f[dec];
    if(fn > fs ){
      flag = 1;
      break;
		}
  }
  return flag;
}
int Count_Maxima(double f[d], double fmax){
  //This function Counts the number of local maxima and generates a file with the information about their relative fitness
  int i, count, flag;
  double fr;
  FILE *arq;
  count = 0;
  arq = fopen("/home/larissa/Landscape Information K = 4.txt", "w");
  fprintf(arq, "Landscape with N = 12 and K = %d and Local Maximus: \n", K);
  for(i = 0; i < d; i++){
    flag = Local_Maximum(i, f);
    if(flag == 0){
      count = count + 1;
      fr = f[i];
      fr = fr/fmax;
    //  fprintf(arq, "%d \t %f \n", count, fr);
    }
  }
  fclose(arq);
  return count;
}
