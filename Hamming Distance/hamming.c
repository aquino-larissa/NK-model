//Genetic Algorithm
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
int HammingDistance(int state1, int state2);
double averageHammingDistance(double f[d], int max);
void localMaximaNumber();
void HammingDistanceMaxima(int N, int d);
void HammingDistanceMaximavsN(int N, int d);
void Graph1();
void Graph2();

int main(){
 //Graph1();
 Graph2();
 return 0;
}
void Graph1(){
  FILE *fptr;
 //complexityCatastrophy();
  fptr = fopen("/home/larissa/maxima_hamming_distance_avg_vs_k.txt", "w");
  fclose(fptr);
  N = 10;
  d = pow(2,N);
  HammingDistanceMaxima(N, d);
  N = 12;
  d = pow(2,N);
  HammingDistanceMaxima(N, d);
  N = 15;
  d = pow(2,N);
  HammingDistanceMaxima(N, d);
  return;
}
void Graph2(){
  N = 15;
  int d = pow(2,N);
  double fit[d];
  int Nsample = 500;
  int Dh, sample, i, k, max, flag, count[N];
  double avg_fit[N];
  double maxima_number[N];
  FILE *fptr;
  fptr = fopen("/home/larissa/Documentos/modelo NK/Hamming Distance/graph2-N15.txt", "w");
  fclose(fptr);
  for (i = 0; i < N; i++){
  		 avg_fit[i] = 0.;
    	 maxima_number[i] = 0.;
    	 count[i] = 0;
  }
  for (k = 0; k < 5; k++){
    K = 2*k;
    for (sample = 0; sample < Nsample; sample++){
      Landscape(N, K, d, fit);
      max = Global_Max(d, fit);
      for (i = 0; i < d; i++){
          Dh = HammingDistance(i, max);
          avg_fit[Dh] += fit[i];
          count[Dh] +=1;
          flag = Local_Maximum(i, fit);
          if(flag == 0){maxima_number[Dh] = maxima_number[Dh] + 1.;}
      }
    }
    fptr = fopen("/home/larissa/Documentos/modelo NK/Hamming Distance/graph2-N15.txt", "a");
    for(i = 0; i < N; i++){
		 avg_fit[i] = avg_fit[i]/(count[i]);
		 maxima_number[i] = maxima_number[i]/Nsample;
		 fprintf(fptr, "%d\t%d\t%d\t%f\t%f\n", N, K, i, avg_fit[i], maxima_number[i]);
		 printf("%d\t%d\t%d\t%f\t%f\n", N, K, i, avg_fit[i], maxima_number[i]);
		 avg_fit[i] = 0.;
    	 maxima_number[i] = 0.;
    	 count[i] = 0;
    }
  fclose(fptr);
  }
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
double averageHammingDistance(double f[d], int max){
  int i, count, flag;
  double avg_distance;
  count = 0;
  avg_distance = 0.;
  for(i = 0; i < d; i++){
    if(i!=max){
      flag = Local_Maximum( i, f);
      if(flag == 0){
        avg_distance += HammingDistance(i, max);
        count += 1;
      }
      //printf("%d\n", count);
    }
  }
  if(count >0){
    avg_distance = avg_distance/count;}

  return avg_distance;
}
int HammingDistance(int state1, int state2){
  int count, i;
  int v1[N], v2[N];
  count = 0;
  for (i = 0; i < N; i++){
    v1[N-i-1] = state1%2;
    v2[N-i-1] = state2%2;
    state1 = state1/2;
    state2 = state2/2;
    if(v1[N-i-1] != v2[N-i-1]){
      count += 1;
    }
  }
  return count;
}
double averageFitnessMaxima(double f[d], int max){
 int flag, i;
 double avg_fit = 0.;
 int count = 0;
  for(i = 0; i < d; i++){
    if(i!=max){
      flag = Local_Maximum( i, f);
      if(flag == 0){
        avg_fit += f[max]-f[i];
        count += 1;
      }
    }
  }
  if(count > 0){
    avg_fit = avg_fit/count;}
  return avg_fit;
}
double deltaFitnessMaxima(int d, double f[d], int max){
 int flag, i;
 double avg_fit = 0.;
 double avg_fit_max = 0.;
 double delta_fit = 0.;
 int count = 0;
  for(i = 0; i < d; i++){
      flag = Local_Maximum( i, f);
      avg_fit += f[i];
      if(flag == 0){
        avg_fit_max += f[i];
        count += 1;
      }
    }
  //printf("%.3f, %d \n", avg_fit_max, count);
  if(count > 0){
    delta_fit = (avg_fit_max/count) - (avg_fit/d);
  }
  else{delta_fit = avg_fit/d;}
  return delta_fit;
}
double fitnessMaxima(int N, int K){
  d = pow(2, N);
  double fit[d], avg_fit_local_max, delta_fit;
  int sample, i, max;
  Nsample = 1000;
  delta_fit = 0.;
  for( sample= 0;  sample < Nsample; sample++){
    Landscape(N, K, d, fit);
    max = Global_Max(d, fit); //Global Maximum
  //  avg_fit_local_max = averageFitnessMaxima(fit, max);
    avg_fit_local_max = deltaFitnessMaxima(d, fit, max);
    delta_fit += avg_fit_local_max;
  }
  printf("%.3f\n",delta_fit );
  delta_fit = delta_fit/Nsample;
  return delta_fit;
}
void complexityCatastrophy(){
  int  max, sample, i;
  double fitness_difference[N];
  int Nmax = 20;
  FILE *fptr;
  time_t t;
  srand((unsigned) time(&t));
  fptr = fopen("/home/larissa/complexity_catastrophy_graph3.txt", "w");
  fclose(fptr);
  for (i =0; i < 3; i++){
    if(i==0){K = 4;}
    if(i==1){K = 4;}
    if(i==2){K = 8;}
    for (N = 3; N < Nmax; N++){
      K = N-1;
      fptr = fopen("/home/larissa/complexity_catastrophy_graph3.txt", "a");
      fitness_difference[N] = fitnessMaxima(N, K);
      printf("%d, %d, %.2f \n", K, N, fitness_difference[N]);
      fprintf(fptr, "%d \t %d \t %.5f \n", K, N, fitness_difference[N]);
      fclose(fptr);
      }
  }
return;
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
void HammingDistanceMaxima(int N, int d){
  int  max, Nmax, sample;
  double avg_distance[N], fit[d];
  FILE *fptr;
  time_t t;
  srand((unsigned) time(&t));
  Nsample = 1000;
  for (K = 1; K < N; K++){
    fptr = fopen("/home/larissa/maxima_hamming_distance_avg_vs_k.txt", "a");
    avg_distance[K] = 0.;
      for( sample= 0;  sample < Nsample; sample++){
        Landscape(N, K, d, fit);
        max = Global_Max(d, fit); //Global Maximum
        avg_distance[K] += averageHammingDistance(fit, max); //Count_Maxima(N, K, d, fit, fit[max]);
    //    printf("%d, %d, %.2f\n", sample, K, avg_distance[K]);
      }
    avg_distance[K] = avg_distance[K]/Nsample;
    printf("%d, %d, %.2f\n",N, K, avg_distance[K]);
    fprintf(fptr, "%d \t %d \t %.3f \n", N, K, avg_distance[K]);
    fclose(fptr);
    }

  return;
}
void HammingDistanceMaximavsN(int N, int d){
  int  max, Nmax, sample;
  double avg_distance, fit[d];
  FILE *fptr;
  time_t t;
  K = 2;
  Nsample = 5000;
  srand((unsigned) time(&t));
  fptr = fopen("/home/larissa/maxima_hamming_distance_avg_vs_n.txt", "a");
  avg_distance = 0.;
  for( sample= 0;  sample < Nsample; sample++){
    Landscape(N, K, d, fit);
    max = Global_Max(d, fit); //Global Maximum
    avg_distance += averageHammingDistance(fit, max); //Count_Maxima(N, K, d, fit, fit[max]);
    //    printf("%d, %d, %.2f\n", sample, K, avg_distance[K]);
    }
    avg_distance = avg_distance/Nsample;
    printf("%d, %.2f\n", N, avg_distance);
    fprintf(fptr, "%d \t %.3f \n", N, avg_distance);
  fclose(fptr);
  return;
}
