//Genetic Algorithm
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

int N = 10;
int K = 2;
int d = 0;
int Nsample = 501;

int array_to_dec(int n, int v[n]);
int Global_Max(int d, double f[d]);
int Local_Maximum( int S, double f[d]);
int Count_Maxima( double f[d], double fmax);
void Landscape(int N, int K, int d, double fit[d]);
void localMaximaNumber();
void complexityCatastrophy();
void averageFitnessMaxima(int N, int K);
double fitnessMaxima(int d, double f[d]);

int main(){
  complexityCatastrophy();
  return 0;
}

//Creating Data files
void complexityCatastrophy(){ //calcula a diferença entre o valor médio fitness dos máximos e valor médio do lanscape
  int  max, sample, i;
  double fitness_difference[N];
  int Nmax = 18;
  FILE *fptr;
  time_t t;
  srand((unsigned) time(&t));
  fptr = fopen("/home/larissa/Downloads/graph1.txt", "w");
  fclose(fptr);
  for (N = 4; N < Nmax; N++){
    K = N-1;
    averageFitnessMaxima(N, K);
    }
  return;
  }
//Calculates maxima average fitness:
void averageFitnessMaxima(int N, int K){
  d = pow(2, N);
  FILE *fptr;
  double fit[d], fit_maxima, fit_maximum, avg;
  int sample, i, max;
  Nsample = 50;
  fit_maxima = 0.;
  fit_maximum = 0.;
  avg = 0.;
  for( sample = 0;  sample < Nsample; sample++){
    Landscape(N, K, d, fit);
  //  avg_fit_local_max = averageFitnessMaxima(fit, max);
    fit_maxima += fitnessMaxima(d, fit);
    max = Global_Max(d, fit);
    fit_maximum += fit[max];
    for (i = 0; i < d; i++){
      avg += fit[i];
    }
  }
  //  printf("%.3f\n",delta_fit );
  fit_maxima = fit_maxima/Nsample;
  fit_maximum = fit_maximum/Nsample;
  avg = avg/(d*Nsample);

  fptr = fopen("/home/larissa/Downloads/graph1.txt", "a");
  printf("%d, %d, %.4f , %.4f, %.4f \n", K, N, fit_maxima, fit_maximum, avg);
  fprintf(fptr, "%d \t %d \t %f \t %f \n", K, N, fit_maxima, fit_maximum);
  fclose(fptr);
  return;
  }
double fitnessMaxima(int d, double f[d]){
 int flag, i;
 double avg_fit_max = 0.;
 int count = 0;
  for(i = 0; i < d; i++){
      flag = Local_Maximum( i, f);
      if(flag == 0){
        avg_fit_max += f[i];
        count += 1;
      }
    }
  //printf("%.3f, %d \n", avg_fit_max, count);
  avg_fit_max = avg_fit_max/count;
  return avg_fit_max;
}
//Generating the landscape
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
//Identifies Maxima
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
void localMaximaNumber(){
  int  i, max, Nmax, sample;
  double avg[N], fit[d];
  FILE *fptr;
  time_t t;
  srand((unsigned) time(&t));
  Nsample = 3000;
    for (K = 0; K < N; K++){
      fptr = fopen("/home/larissa/local_maxima_avg_vs_k.txt", "a");
      avg[K] = 0.;
      for( sample= 0;  sample < Nsample; sample++){
        Landscape(N, K, d, fit);
        max = Global_Max(d, fit); //Global Maximum
        Nmax = Count_Maxima(fit, fit[max]);
        avg[K] += Nmax;
        //printf("%d, %d, %.2f \n", N, K, avg[K]);
      }
    avg[K] = avg[K]/Nsample;
    fprintf(fptr, "%d \t %d \t %.3f \n", N, K, avg[K]);
    printf("%d, %d, %.2f \n", N, K, avg[K]);
    fclose(fptr);
    }
  return;
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
