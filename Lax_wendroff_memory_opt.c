
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define _CRT_SECURE_NO_WARNINGS
#define T 500
#define M_PIl 3.141592653589793238462643383279502884L /* pi */
#define grid_size_x 0.01
#define grid_size_y 0.01
#define step_ratio 0.34
#define c0 340
#include <stdint.h>
#include<time.h>
#include<pthread.h>
//#include<iostream>
 #include <sys/time.h>
//#include <sysinfoapi.h>

//pthread_mutex_t test_mutex;
#define N   3
#define AVE (17*12*2)/N
#define X AVE*N+1
#define Y AVE*N+1
#define Z 50

#define VARM(vi,vj) vi+(X+1)*vj

int t;
double time_step = step_ratio * grid_size_x / c0;

double Coeff_x = step_ratio;

double t_fun[T+1][N];


char* names[12];


double* u_n[N] ;
double* u[N] ;
double* u_p[N] ;
double* r_mat[N] ;
double* r_mat_n[N] ;
double* s_mat[N] ;
double* s_mat_n[N] ;
double* s_mat_p[N] ;
double* l_mat[N] ;
double* l_mat_n[N] ;


#pragma comment(lib,"pthreadVC2.lib")

void main()
{

    struct  timeval  start;
	struct  timeval  maclloc_start;
	struct  timeval  initialisation_start;
	struct  timeval  leap_frog_start;
	struct  timeval  p_thread_start;
	struct  timeval  reassign_start;
	struct  timeval  reassign_end;
    struct  timeval  end;
	long main_total_t=0, malloc_total_t=0,ini_total_t=0,leap_frog_total_t=0,p_thread_total_t=0,reassign_total_t=0;




	gettimeofday(&start,NULL);
	//Declare Function

	void* lax_wendroff_advection_all_in_1(void* arg);
	void* lax_wendroff_reassign(void* arg);

	
	float  viscosity;

	double    t0_g, tp_g, square_of_Coeff_x, Q_x, y;
	int var_c,var_ip1,var_jp1,var_im1,var_jm1,mid_k,mid_var,mid_var_ip1,mid_var_jp1,mid_var_im1,mid_var_jm1;
	int i, j, k, P_coeff, id,k_mid,pth;
	double k1,k2,k3,k1_mid,k2_mid,k3_mid;

	P_coeff = 1024 * 1024;

	pthread_t* pthread_id = NULL; //
	pthread_id = (pthread_t*)malloc(N * sizeof(pthread_t));

	square_of_Coeff_x = Coeff_x * Coeff_x;

	Q_x = (grid_size_x - c0 * step_ratio) / (grid_size_x + c0 * step_ratio);

	double mid_Y = 0.5 * Y;
	double mid_X = 0.5 * X;

	k2_mid = ((int)mid_Y - (double)AVE - 1);
	k3_mid = (double)AVE;
	k1_mid = k2_mid / k3_mid + 1;
	k_mid = (int)k1_mid;

	//mid_k=(((int)mid_Y-AVE-1) /AVE +1);
	mid_var=VARM((int)mid_X,(int)mid_Y)-(AVE+1)*(X+1)-(k_mid-1)*(AVE)*(X + 1);
	mid_var_ip1=VARM(((int)mid_X+1),(int)mid_Y)-(AVE+1)* (X + 1) -(k_mid-1)*(AVE)* (X + 1);
	mid_var_jp1=VARM((int)mid_X,((int)mid_Y+1))-(AVE+1)* (X + 1) -(k_mid-1)*(AVE)* (X + 1);
	mid_var_im1=VARM(((int)mid_X-1),(int)mid_Y)-(AVE+1)* (X + 1) -(k_mid-1)*(AVE)* (X + 1);
	mid_var_jm1=VARM((int)mid_X,((int)mid_Y-1))-(AVE+1)* (X + 1) -(k_mid-1)*(AVE)* (X + 1);

	
	char names1[] = { "wave_data_lax_wendroff_tm1.txt" };
	names[0] = names1;
	char names2[] = { "wave_data_lax_wendroff_tm2.txt" };
	names[1] = names2;
	char names3[] = { "wave_data_lax_wendroff_tm3.txt" };
	names[2] = names3;
	char names4[] = { "wave_data_lax_wendroff_tm4.txt" };
	names[3] = names4;
	char names5[] = { "wave_data_lax_wendroff_tm5.txt" };
	names[4] = names5;
	char names6[] = { "wave_data_lax_wendroff_tm6.txt" };
	names[5] = names6;
	char names7[] = { "wave_data_lax_wendroff_tm7.txt" };
	names[6] = names7;
	char names8[] = { "wave_data_lax_wendroff_tm8.txt" };
	names[7] = names8;
	char names9[] = { "wave_data_lax_wendroff_tm9.txt" };
    names[8] = names9;
  	char names10[] = { "wave_data_lax_wendroff_tm10.txt" };
  	names[9] = names10;
  	char names11[] = { "wave_data_lax_wendroff_tm11.txt" };
  	names[10] = names11;
  	char names12[] = { "wave_data_lax_wendroff_tm12.txt" };
  	names[11] = names12;
	gettimeofday(&maclloc_start,NULL);	

	for (i = 0; i < N; i++){
		if(i == 0 || i == N-1){
		u_n[i] = malloc(sizeof(double) * (X + 1)* (AVE + 1));
		u [i] = malloc(sizeof(double) * (X + 1)* (AVE + 1));
		u_p [i] = malloc(sizeof(double) * (X + 1)* (AVE + 1));
		r_mat [i] = malloc(sizeof(double) * (X + 1)* (AVE + 1));
		s_mat [i] = malloc(sizeof(double) * (X + 1)* (AVE + 1));
		l_mat [i] = malloc(sizeof(double) * (X + 1)* (AVE + 1));
		s_mat_n [i] = malloc(sizeof(double) * (X + 1)* (AVE + 1));
		r_mat_n [i] = malloc(sizeof(double) * (X + 1)* (AVE + 1));
		l_mat_n [i] = malloc(sizeof(double) * (X + 1)* (AVE + 1));
		s_mat_p [i] = malloc(sizeof(double) * (X + 1)* (AVE + 1));
		}
		else{
		u_n[i] = malloc(sizeof(double) * (X + 1)* (AVE));
		u [i] = malloc(sizeof(double) * (X + 1)* (AVE));
		u_p [i] = malloc(sizeof(double) * (X + 1)* (AVE));
		r_mat [i] = malloc(sizeof(double) * (X + 1)* (AVE));
		s_mat [i] = malloc(sizeof(double) * (X + 1)* (AVE));
		l_mat [i] = malloc(sizeof(double) * (X + 1)* (AVE));
		s_mat_n [i] = malloc(sizeof(double) * (X + 1)* (AVE));
		r_mat_n [i] = malloc(sizeof(double) * (X + 1)* (AVE));
		l_mat_n [i] = malloc(sizeof(double) * (X + 1)* (AVE));
		s_mat_p [i] = malloc(sizeof(double) * (X + 1)* (AVE));
		}
	}









	
	tp_g = 15;
	t0_g = 5;
	/*Matrix set up*/

	gettimeofday(&initialisation_start,NULL);

		for (j = 0; j <= Y; j++) {
			for (i = 0; i <= X; i++) {
				//
				if(j<=AVE){
					pth = 0;
				}
				else if (j >= (N-1)*AVE+1){
					pth = N-1;
				}
				else{
					k2 = (j - (double)AVE - 1);
					k3 = (double)AVE;
					k1 = k2 / k3 + 1;
					pth =(int)k1;
				}
				if (pth==0){
					var_c=VARM(i,j);
				}
				else{
					var_c=VARM(i,j)-(AVE+1)*(X+1)-(pth-1)*(AVE)* (X + 1);
				}
				//

			s_mat_p[pth][var_c] = r_mat_n[pth][var_c] = s_mat_n[pth][var_c] = l_mat_n[pth][var_c] = r_mat[pth][var_c] = s_mat[pth][var_c] = l_mat[pth][var_c] = u_p[pth][var_c] = u[pth][var_c] = u_n[pth][var_c] = 0;
			//printf("%lf  ", u_n[pth][var_c]);

			}

		}
	
	/*Main computation*/


	for (t = 0; t <= T; t++) {
		gettimeofday(&leap_frog_start,NULL);
		u[k_mid][mid_var] = P_coeff * exp(-pow((t - t0_g + 1) / tp_g, 2));
		 if (t == 0) {
			for (j = 1; j <= Y-1 ; j++) {
				for (i = 1; i <= X-1 ; i++) {
					//
					if (j <= AVE) {
						pth = 0;
					}
					else if (j >= (N - 1) * AVE + 1) {
						pth = N - 1;
					}
					else {
						k2 = (j - (double)AVE - 1);
						k3 = (double)AVE;
						k1 = k2 / k3 + 1;
						pth = (int)k1;
					}
				//

				switch (pth){
				case 0:
					var_c=VARM(i,j);
					var_ip1=VARM((i+1),j);
					var_im1=VARM((i-1),j);
					var_jm1=VARM(i,(j-1));
				if (j==AVE){
					var_jp1=VARM(i,(j+1))-(AVE+1)* (X + 1) -((pth+1)-1)*(AVE)* (X + 1);
					u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth+1][var_jp1] - 4 * u[pth][var_c]);
					//printf("%lf",u_n[pth][var_c]);
				}else{
					var_jp1=VARM(i,(j+1));
					u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1] - 4 * u[pth][var_c]);
					//printf("%lf",u_n[pth][var_c]);
				}
				//printf("[i=%d, j=%d] ", i, j);
					break;
				case N-1:
					var_c=VARM(i,j)-(AVE+1)* (X + 1) -(pth-1)*(AVE)* (X + 1);
					var_ip1=VARM((i+1),j)-(AVE+1)* (X + 1) -(pth-1)*(AVE)* (X + 1);
					var_im1=VARM((i-1),j)-(AVE+1)* (X + 1) -(pth-1)*(AVE)* (X + 1);
					var_jp1=VARM(i,(j+1))-(AVE+1)* (X + 1) -(pth-1)*(AVE)* (X + 1);
					if(j==pth*AVE+1){
						if (pth == 1) {
							var_jm1 = VARM(i, (j - 1));
						}
						else {
							var_jm1 = VARM(i, (j - 1)) - (AVE + 1) * (X+1) - ((pth - 1) - 1) * (AVE)*(X+1);
						}
						u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth-1][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1] - 4 * u[pth][var_c]);
						//printf("%lf",u_n[pth][var_c]);
					}else{
						var_jm1=VARM(i,(j-1))-(AVE+1)*(X+1)-(pth-1)*(AVE)*(X+1);
						u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1] - 4 * u[pth][var_c]);
						//printf("%lf",u_n[pth][var_c]);
					}
					//printf("[i=%d, j=%d] ", i, j);
					break;
				default:
					var_c=VARM(i,j)-(AVE+1)*(X+1)-(pth-1)*(AVE)*(X+1);
					var_ip1=VARM((i+1),j)-(AVE+1)*(X+1)-(pth-1)*(AVE)*(X+1);
					var_im1=VARM((i-1),j)-(AVE+1)*(X+1)-(pth-1)*(AVE)*(X+1);

				if (j==AVE*(pth+1)){
					var_jm1=VARM(i,(j-1))-(AVE+1)*(X+1)-(pth-1)*(AVE)*(X+1);
					var_jp1=VARM(i,(j+1))-(AVE+1)*(X+1)-((pth+1)-1)*(AVE)*(X+1);
					u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth+1][var_jp1] - 4 * u[pth][var_c]);
					//printf("%lf",u_n[pth][var_c]);
				}else if(j==AVE*pth+1){
					var_jp1=VARM(i,(j+1))-(AVE+1)*(X+1)-(pth-1)*(AVE)*(X+1);
					if (pth == 1) {
						var_jm1 = VARM(i, (j - 1));
					}
					else{
						var_jm1 = VARM(i, (j - 1)) - (AVE + 1) * (X+1) - ((pth - 1) - 1) * (AVE)*(X+1);
					}
					
					u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth-1][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1] - 4 * u[pth][var_c]);
					//printf("%lf",u_n[pth][var_c]);
				}else{
					var_jp1=VARM(i,(j+1))-(AVE+1)*(X+1)-(pth-1)*(AVE)*(X+1);
					var_jm1=VARM(i,(j-1))-(AVE+1)*(X+1)-(pth-1)*(AVE)*(X+1);
					u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1] - 4 * u[pth][var_c]);
	/*				printf("%lf",u_n[pth][var_c]);*/
				} 
				//printf("[i=%d, j=%d] ", i,j);
					break;
				}


				}

			}
		 }

		r_mat[k_mid][mid_var] = c0 * (u[k_mid][mid_var_ip1] - u[k_mid][mid_var_im1]) / (2.0 * grid_size_x);
		l_mat[k_mid][mid_var] = c0 * (u[k_mid][mid_var_jp1] - u[k_mid][mid_var_jm1]) / (2.0 * grid_size_x);
		s_mat[k_mid][mid_var] = 2 * (u[k_mid][mid_var] - u_p[k_mid][mid_var]) / time_step + s_mat_p[k_mid][mid_var];


		//--------------------ALL in one-------------------------------------

		gettimeofday(&p_thread_start,NULL);
		 for (id = 0; id < N; id++)
		 {
		 	pthread_create(pthread_id + id, NULL, lax_wendroff_advection_all_in_1, (void *) (intptr_t) id);
		 }

		 for (id = 0; id < N; id++)
		 {

		 	pthread_join(pthread_id[id], NULL);

		 }


		gettimeofday(&reassign_start,NULL);


		 for (id = 0; id < N; id++)
		 {
		 	pthread_create(pthread_id + id, NULL, lax_wendroff_reassign, (void *) (intptr_t) id);
		 }

		 for (id = 0; id < N; id++)
		 {

		 	pthread_join(pthread_id[id], NULL);

		 }

		gettimeofday(&reassign_end,NULL);
		leap_frog_total_t +=  ((long long)p_thread_start.tv_sec-(long long)leap_frog_start.tv_sec)*1000+ (p_thread_start.tv_usec-leap_frog_start.tv_usec)/1000;
		p_thread_total_t +=  ((long long)reassign_start.tv_sec-(long long)p_thread_start.tv_sec)*1000+ (reassign_start.tv_usec-p_thread_start.tv_usec)/1000;
		reassign_total_t +=  ((long long)reassign_end.tv_sec-(long long)reassign_start.tv_sec)*1000+ (reassign_end.tv_usec-reassign_start.tv_usec)/1000;
		//printf("%f s  and   %f s and %f s and %f s\n", (t2 - t1) / CLOCKS_PER_SEC, (t3 - t2) / CLOCKS_PER_SEC,(t3-t0)/CLOCKS_PER_SEC, (t3 - t1) / CLOCKS_PER_SEC);

	}
	gettimeofday(&end,NULL);
	main_total_t =  ((long long)end.tv_sec-(long long)start.tv_sec)*1000+ (end.tv_usec-start.tv_usec)/1000;
	malloc_total_t =  ((long long)initialisation_start.tv_sec-(long long)maclloc_start.tv_sec)*1000+ (initialisation_start.tv_usec-maclloc_start.tv_usec)/1000;
	//ini_total_t =  (end.tv_sec-initialisation_start.tv_sec)+ (end.tv_usec-start.tv_usec)*10e-6;

	printf("main  %ld s\n malloc  %ld\n leap frog  %ld\n pthread  %ld \n reassign  %ld \n",main_total_t,malloc_total_t,leap_frog_total_t,p_thread_total_t, reassign_total_t);
//	for (i=0,i<T+1,i++){
//			for (j=0,j<N,j++){
//					t_store=t_fun[VARM(i,j)];
//			printf("%p s \r", t_fun[VARM(i,j)]);
//
//			}

//	}
	// free(u_n);
	// free(u_p);
	// free(u);
	// free(r_mat);
	// free(r_mat_n);
	// free(s_mat);
	// free(s_mat_n);
	// free(l_mat);
	// free(l_mat_n);
	// free(s_mat_p);
	//pthread_mutex_destroy(&test_mutex);

	//for (j = 0; j <= Y; j++) {
	//	for (i = 0; i <= X; i++) {
	//		fprintf(fpWrite,"%lf ", u[VARM(i,j)]);
	//	}
	//	fprintf(fpWrite,"\n");
	//}




	//return 0;
}



void* lax_wendroff_advection_all_in_1(void* arg) {

	int            n =  (int) (intptr_t) arg;  //Nth thread
	//struct  timeval t_fun_start,t_fun_end;
	long t_fun_th;
	double      start = n * AVE+1;
	double      end = start + AVE - 1;

	long long      i, j, k;
	int var_c,var_ip1,var_jp1,var_im1,var_jm1;

	// gettimeofday(&t_fun_start,NULL);
	switch (n) {
	case 0:
	 	for (i = 1; i <= X - 1; i++){

	 	j=end;
	 	var_c=VARM(i,j);
	 	var_ip1=VARM((i+1),j);
	 	var_im1=VARM((i-1),j);
	 	var_jp1=VARM(i,(j+1))-(AVE+1)* (X + 1) -((n+1)-1)*(AVE)* (X + 1);
	 	var_jm1=VARM(i,(j-1));

	 	r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * (1 - Coeff_x) * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
	 	l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n+1][var_jp1] - s_mat[n][var_jm1])) + 0.5 * (1 - Coeff_x) * (l_mat[n+1][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
	 	s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_im1])))\
	 		+ (Coeff_x * (0.5 * (l_mat[n+1][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n+1][var_jp1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])));

	 } 
	 for (j = start; j <= end-1; j++) {
	 	for (i = 1; i <= X - 1; i++) {

	 	var_c=VARM(i,j);
	 	var_ip1=VARM((i+1),j);
	 	var_im1=VARM((i-1),j);
	 	var_jp1=VARM(i,(j+1));
	 	var_jm1=VARM(i,(j-1));

	 	r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * (1 - Coeff_x) * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
	 	l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n][var_jm1])) + 0.5 * (1 - Coeff_x) * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
	 	s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_im1])))\
	 		+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])));
				
	 	}
	 }
	 	for (j = 0; j <= end; j++) {
	 		for (i = 0; i <= X; i++) {
	 			var_c=VARM(i,j);
	 			u_n[n][var_c] = u[n][var_c]  + 0.5 * time_step * (s_mat_n[n][var_c]  - s_mat[n][var_c] );
	 		}
	 	}


		// if (t == 0) {
		// 	FILE* fpWrite = fopen(names[n], "w");
		// 	if (fpWrite == NULL)
		// 	{
		// 		printf("Error");
		// 		return 0;
		// 	}
		// 	for (j = 0; j <= end; j++) {
		// 		for (i = 0; i <= X; i++) {
		// 			var_c=VARM(i,j);
		// 			fprintf(fpWrite, "%lf ", u_n[n][var_c] );
		// 			//printf( "u_n[%d][%d]=%lf ",i,j, u_n[n][var_c]);
		// 		}
		// 		fprintf(fpWrite, "\n");
		// 	}
		// 	fclose(fpWrite);

		// }
		// else if (t % 5 == 0) {
		// 	FILE* fpWrite = fopen(names[n], "a");
		// 	if (fpWrite == NULL)
		// 	{
		// 		return 0;
		// 	}
		// 	for (j = 0; j <= end; j++) {
		// 		for (i = 0; i <= X; i++) {
		// 			var_c=VARM(i,j);
		// 			fprintf(fpWrite, "%lf ", u_n[n][var_c] );
		// 			//printf("u_n[%d][%d]=%lf ", i, j, u_n[n][var_c]);
		// 		}
		// 		fprintf(fpWrite, "\n");
		// 	}
		// 	fclose(fpWrite);
		// }

		break;

	case N-1:
	 for (i = 1; i <= X - 1; i++){

	 	j=start;
	 	var_c=VARM(i,j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 	var_ip1=VARM((i+1),j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 	var_im1=VARM((i-1),j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 	var_jp1=VARM(i,(j+1))-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 	var_jm1=VARM(i,(j-1))-(AVE+1)*(X + 1)-((n-1)-1)*(AVE)*(X + 1);

	 	r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * (1 - Coeff_x) * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
	 	l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n-1][var_jm1])) + 0.5 * (1 - Coeff_x) * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n-1][var_jm1]));
	 	s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_im1])))\
	 		+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n-1][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] - 2.0 * s_mat[n][var_c] + s_mat[n-1][var_jm1])));

	 } 
	 for (j = start+1; j <= end; j++) {
	 	for (i = 1; i <= X - 1; i++) {

	 			var_c=VARM(i,j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 			var_ip1=VARM((i+1),j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 			var_im1=VARM((i-1),j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 			var_jp1=VARM(i,(j+1))-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 			var_jm1=VARM(i,(j-1))-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);

	 			r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * (1 - Coeff_x) * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
	 			l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n][var_jm1])) + 0.5 * (1 - Coeff_x) * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
	 			s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_im1])))\
	 				+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])));
				
	 	}
	 }
	 	for (j = start; j <= Y; j++) {
	 		for (i = 0; i <= X; i++) {
	 			var_c=VARM(i,j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 			u_n[n][var_c] = u[n][var_c]  + 0.5 * time_step * (s_mat_n[n][var_c]  - s_mat[n][var_c] );
	 		}
		 }
		// if (t == 0) {
		// 	FILE* fpWrite = fopen(names[n], "w");
		// 	if (fpWrite == NULL)
		// 	{
		// 		printf("Error");
		// 		return 0;
		// 	}
		// 	for (j = start; j <= Y; j++) {
		// 		for (i = 0; i <= X; i++) {
		// 			var_c=VARM(i,j)-(AVE+1)*(X+1)-(n-1)*(AVE)* (X + 1);
		// 			fprintf(fpWrite, "%lf ", u_n[n][var_c] );
		// 			//printf("u_n[%d][%d]=%lf ", i, j, u_n[n][var_c]);
		// 		}
		// 		fprintf(fpWrite, "\n");
		// 	}
		// 	fclose(fpWrite);
		// }
		// else if (t % 5 == 0) {
		// 	FILE* fpWrite = fopen(names[n], "a");
		// 	if (fpWrite == NULL)
		// 	{
		// 		return 0;
		// 	}
		// 	for (j = start; j <= Y; j++) {
		// 		for (i = 0; i <= X; i++) {
		// 			var_c=VARM(i,j)-(AVE+1)* (X + 1) -(n-1)*(AVE)* (X + 1);
		// 			fprintf(fpWrite, "%lf ", u_n[n][var_c] );
		// 			//printf("u_n[%d][%d]=%lf ", i, j, u_n[n][var_c]);
		// 		}
		// 		fprintf(fpWrite, "\n");
		// 	}
		// 	fclose(fpWrite);
		// }
		break;

	default:
	 for (i = 1; i <= X - 1; i++){

	 	j=start;
	 	var_c=VARM(i,j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 	var_ip1=VARM((i+1),j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 	var_im1=VARM((i-1),j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 	var_jp1=VARM(i,(j+1))-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 	if (n == 1) {
	 		var_jm1 = VARM(i, (j - 1));
	 	}
	 	else {
	 		var_jm1=VARM(i,(j-1))-(AVE+1)*(X + 1)-((n-1)-1)*(AVE)*(X + 1);
	 	}
		
	 	r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * (1 - Coeff_x) * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
	 	l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n-1][var_jm1])) + 0.5 * (1 - Coeff_x) * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n-1][var_jm1]));
	 	s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_im1])))\
	 		+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n-1][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] - 2.0 * s_mat[n][var_c] + s_mat[n-1][var_jm1])));

	 } 
	 for (i = 1; i <= X - 1; i++){

	 	j=end;
	 	var_c=VARM(i,j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 	var_ip1=VARM((i+1),j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 	var_im1=VARM((i-1),j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 	var_jp1=VARM(i,(j+1))-(AVE+1)*(X + 1)-((n+1)-1)*(AVE)*(X + 1);
	 	var_jm1=VARM(i,(j-1))-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);

	 	r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * (1 - Coeff_x) * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
	 	l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n+1][var_jp1] - s_mat[n][var_jm1])) + 0.5 * (1 - Coeff_x) * (l_mat[n+1][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
	 	s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_im1])))\
	 		+ (Coeff_x * (0.5 * (l_mat[n+1][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n+1][var_jp1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])));

	 } 
	 for (j = start+1; j <= end-1; j++) {
	 	for (i = 1; i <= X - 1; i++) {

	 			var_c=VARM(i,j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 			var_ip1=VARM((i+1),j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 			var_im1=VARM((i-1),j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 			var_jp1=VARM(i,(j+1))-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 			var_jm1=VARM(i,(j-1))-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);

	 			r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * (1 - Coeff_x) * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
	 			l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n][var_jm1])) + 0.5 * (1 - Coeff_x) * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
	 			s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_im1])))\
	 				+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] - 2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])));
				
	 	}
	 }

	 	for (j = start; j <= end; j++) {
	 		for (i = 0; i <= X; i++) {
	 			var_c=VARM(i,j)-(AVE+1)*(X + 1)-(n-1)*(AVE)*(X + 1);
	 			u_n[n][var_c] = u[n][var_c]  + 0.5 * time_step * (s_mat_n[n][var_c]  - s_mat[n][var_c] );
	 		}
	 	}

		// if (t == 0) {
		// 	FILE* fpWrite = fopen(names[n], "w");
		// 	if (fpWrite == NULL)
		// 	{
		// 		printf("Error");
		// 		return 0;
		// 	}
		// 	for (j = start; j <= end; j++) {
		// 		for (i = 0; i <= X; i++) {
		// 			var_c=VARM(i,j)-(AVE+1)* (X + 1) -(n-1)*(AVE)* (X + 1);
		// 			fprintf(fpWrite, "%lf ", u_n[n][var_c] );
		// 			//printf("u_n[%d][%d]=%lf ", i, j, u_n[n][var_c]);
		// 		}
		// 		fprintf(fpWrite, "\n");
		// 	}
		// 	fclose(fpWrite);
		// }
		// else if (t % 5 == 0) {
		// 	FILE* fpWrite = fopen(names[n], "a");
		// 	if (fpWrite == NULL)
		// 	{
		// 		return 0;
		// 	}
		// 	for (j = start; j <= end; j++) {
		// 		for (i = 0; i <= X; i++) {
		// 			var_c=VARM(i,j)-(AVE+1)* (X + 1) -(n-1)*(AVE)* (X + 1);
		// 			fprintf(fpWrite, "%lf ", u_n[n][var_c] );
		// 			//printf("u_n[%d][%d]=%lf ", i, j, u_n[n][var_c]);
		// 		}
		// 		fprintf(fpWrite, "\n");
		// 	}
		// 	fclose(fpWrite);
		// }
		break;
	
	}

	// gettimeofday(&t_fun_end,NULL);
	// t_fun_th =  ((long long)t_fun_end.tv_sec-(long long)t_fun_start.tv_sec)*1000+ (t_fun_end.tv_usec-t_fun_start.tv_usec)/1000;
	//printf("t_fun[%d][%d]= %ld ms\n",t,n,t_fun_th );
	return 0;
}

void* lax_wendroff_reassign(void* arg){
	int		n =  (int) (intptr_t) arg;  //Nth thread
	//struct  timeval t_fun_start,t_fun_end;
	double	t_fun_th;
	double	start = n * AVE+1;
	double	end = start + AVE - 1;
	long long      i, j;
	int		var_c;
	 //u_p[n][0] = u[n][0];
	 switch (n){
	 case 0:
	 		for (j = 0; j <= end; j++) {
	 			for (i = 0; i <= X; i++) {
	 				var_c=VARM(i,j);
	 				u_p[n][var_c] = u[n][var_c];
	 				u[n][var_c]= u_n[n][var_c];
	 				r_mat[n][var_c] = r_mat_n[n][var_c];
	 				l_mat[n][var_c] = l_mat_n[n][var_c];
	 				s_mat[n][var_c] = s_mat_n[n][var_c];
	 				s_mat_p[n][var_c] = s_mat[n][var_c];
	 			}
	 		}
		
	 	break;

	 case N-1:
	 		for (j = start; j <= Y; j++) {
	 			for (i = 0; i <= X; i++) {
	 				var_c=VARM(i,j)-(AVE+1)* (X + 1) -(n-1)*(AVE)*(X + 1);
	 				u_p[n][var_c] = u[n][var_c];
	 				u[n][var_c]= u_n[n][var_c];
	 				r_mat[n][var_c] = r_mat_n[n][var_c];
	 				l_mat[n][var_c] = l_mat_n[n][var_c];
	 				s_mat[n][var_c] = s_mat_n[n][var_c];
	 				s_mat_p[n][var_c] = s_mat[n][var_c];
	 			}
	 		}
	 	break;
	 default:
	 		for (j = start; j <= end; j++) {
	 			for (i = 0; i <= X; i++) {
	 				var_c=VARM(i,j)-(AVE+1)* (X + 1) -(n-1)*(AVE)* (X + 1);
	 				u_p[n][var_c] = u[n][var_c];
	 				u[n][var_c]= u_n[n][var_c];
	 				r_mat[n][var_c] = r_mat_n[n][var_c];
	 				l_mat[n][var_c] = l_mat_n[n][var_c];
	 				s_mat[n][var_c] = s_mat_n[n][var_c];
	 				s_mat_p[n][var_c] = s_mat[n][var_c];
	 			}
	 		}
	 	break;
	 }
	return 0;
}
