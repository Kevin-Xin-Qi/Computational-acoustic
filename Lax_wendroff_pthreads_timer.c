
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define _CRT_SECURE_NO_WARNINGS
#define T 10
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
#define N   2
#define AVE (17*12*15)/N
#define X AVE*N+1
#define Y AVE*N+1
#define Z 50

int t;
double time_step = step_ratio * grid_size_x / c0;

double Coeff_x = step_ratio;

double t_fun[T+1][N];


char* names[8];


double** u_n = NULL;
double** u = NULL;
double** u_p = NULL;
double** r_mat = NULL;
double** r_mat_n = NULL;
double** s_mat = NULL;
double** s_mat_n = NULL;
double** s_mat_p = NULL;
double** l_mat = NULL;
double** l_mat_n = NULL;


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

	double   square_of_c, t0_g, tp_g, square_of_Coeff_x, Q_x, y;
	int i, j, k, P_coeff, id;

	P_coeff = 1024 * 1024;

	pthread_t* pthread_id = NULL; //
	pthread_id = (pthread_t*)malloc(N * sizeof(pthread_t));

	square_of_Coeff_x = Coeff_x * Coeff_x;

	Q_x = (grid_size_x - c0 * step_ratio) / (grid_size_x + c0 * step_ratio);
	
	double mid_Y = 0.5 * Y;
	double mid_X = 0.5 * X;

	
	char names1[] = { "wave_data_lax_wendroff_t1.txt" };
	names[0] = names1;
	char names2[] = { "wave_data_lax_wendroff_t2.txt" };
	names[1] = names2;
	char names3[] = { "wave_data_lax_wendroff_t3.txt" };
	names[2] = names3;
	char names4[] = { "wave_data_lax_wendroff_t4.txt" };
	names[3] = names4;
	char names5[] = { "wave_data_lax_wendroff_t5.txt" };
	names[4] = names5;
	char names6[] = { "wave_data_lax_wendroff_t6.txt" };
	names[5] = names6;
	char names7[] = { "wave_data_lax_wendroff_t7.txt" };
	names[6] = names7;
	char names8[] = { "wave_data_lax_wendroff_t8.txt" };
	names[7] = names8;


	





	gettimeofday(&maclloc_start,NULL);	
	u_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u_n[i] = malloc(sizeof(double) * (Y + 1));
	}
	u = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u[i] = malloc(sizeof(double) * (Y + 1));
	}
	u_p = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u_p[i] = malloc(sizeof(double) * (Y + 1));
	}
	r_mat = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		r_mat[i] = malloc(sizeof(double) * (Y + 1));
	}
	s_mat = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat[i] = malloc(sizeof(double) * (Y + 1));
	}
	l_mat = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		l_mat[i] = malloc(sizeof(double) * (Y + 1));
	}
	s_mat_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat_n[i] = malloc(sizeof(double) * (Y + 1));
	}
	r_mat_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		r_mat_n[i] = malloc(sizeof(double) * (Y + 1));
	}
	l_mat_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		l_mat_n[i] = malloc(sizeof(double) * (Y + 1));
	}
	s_mat_p = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat_p[i] = malloc(sizeof(double) * (Y + 1));
	}


	
	tp_g = 15;
	t0_g = 5;
	/*Matrix set up*/

	gettimeofday(&initialisation_start,NULL);
	for (j = 0; j <= Y; j++) {
		for (i = 0; i <= X; i++) {
			s_mat_p[i][j] = r_mat_n[i][j] = s_mat_n[i][j] = l_mat_n[i][j] = r_mat[i][j] = s_mat[i][j] = l_mat[i][j] = u_p[i][j] = u[i][j] = u_n[i][j] = 0;


		}
	}





	//pthread_mutex_init(&test_mutex, NULL);

	//time_t time_1 =  time(NULL);
	/*Main computation*/


	for (t = 0; t <= T; t++) {
		gettimeofday(&leap_frog_start,NULL);
		u[(int)mid_X][(int)mid_Y] = P_coeff * exp(-pow((t - t0_g + 1) / tp_g, 2));
		if (t == 0) {
			for (j = 1; j <= Y - 1; j++) {
				for (i = 1; i <= X - 1; i++) {
					u_n[i][j] = 2 * u[i][j] - u_p[i][j] + square_of_Coeff_x * (u[i - 1][j] + u[i][j - 1] + u[i + 1][j] + u[i][j + 1] - 4 * u[i][j]);
					//u_n[i][j] = leap_frog(u[i + 1][j], u[i - 1][j], u[i][j + 1], u[i][j - 1], u[i][j], u_p[i][j], square_of_Coeff_x);

				}
			}
		}
		r_mat[(int)mid_X][(int)mid_Y] = c0 * (u[(int)mid_X + 1][(int)mid_Y] - u[(int)mid_X - 1][(int)mid_Y]) / (2.0 * grid_size_x);
		l_mat[(int)mid_X][(int)mid_Y] = c0 * (u[(int)mid_X][(int)mid_Y + 1] - u[(int)mid_X][(int)mid_Y - 1]) / (2.0 * grid_size_x);
		s_mat[(int)mid_X][(int)mid_Y] = 2 * (u[(int)mid_X][(int)mid_Y] - u_p[(int)mid_X][(int)mid_Y]) / time_step + s_mat_p[(int)mid_X][(int)mid_Y];

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

		for (j = 0; j <= Y; j++) {
			for (i = 0; i <= X; i++) {
				u_p[i][j] = u[i][j];
				u[i][j] = u_n[i][j];
				r_mat[i][j] = r_mat_n[i][j];
				l_mat[i][j] = l_mat_n[i][j];
				s_mat[i][j] = s_mat_n[i][j];
				s_mat_p[i][j] = s_mat[i][j];
			}
		}
		// for (id = 0; id < N; id++)
		// {
		// 	pthread_create(pthread_id + id, NULL, lax_wendroff_reassign, (void *) (intptr_t) id);
		// }

		// for (id = 0; id < N; id++)
		// {

		// 	pthread_join(pthread_id[id], NULL);

		// }
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
//					t_store=t_fun[i][j];
//			printf("%p s \r", t_fun[i][j]);
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
	//		fprintf(fpWrite,"%lf ", u[i][j]);
	//	}
	//	fprintf(fpWrite,"\n");
	//}





}



void* lax_wendroff_advection_all_in_1(void* arg) {

	int            n =  (int) (intptr_t) arg;  //Nth thread
	struct  timeval t_fun_start,t_fun_end;
	long t_fun_th;
	double      start = n * AVE+1;
	double      end = start + AVE - 1;

	long long      i, j;
	gettimeofday(&t_fun_start,NULL);
	for (j = start; j <= end; j++) {
		for (i = 1; i <= X - 1; i++) {
				r_mat_n[i][j] = r_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i + 1][j] - s_mat[i - 1][j])) + 0.5 * (1 - Coeff_x) * (r_mat[i + 1][j] - 2.0 * r_mat[i][j] + r_mat[i - 1][j]));
				l_mat_n[i][j] = l_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i][j + 1] - s_mat[i][j - 1])) + 0.5 * (1 - Coeff_x) * (l_mat[i][j + 1] - 2.0 * l_mat[i][j] + l_mat[i][j - 1]));
				s_mat_n[i][j] = s_mat[i][j] + (Coeff_x * (0.5 * (r_mat[i + 1][j] - r_mat[i - 1][j]) + 0.5 * (1 - Coeff_x) * (s_mat[i + 1][j] - 2 * s_mat[i][j] + s_mat[i - 1][j])))\
					+ (Coeff_x * (0.5 * (l_mat[i][j + 1] - l_mat[i][j - 1]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j + 1] - 2 * s_mat[i][j] + s_mat[i][j - 1])));
				
		}
	}

	switch (n) {
	case 0:
		// for (j = 1; j <= end; j++) {
		// 	for (i = 1; i <= X - 1; i++) {
		// 		r_mat_n[i][j] = r_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i + 1][j] - s_mat[i - 1][j])) + 0.5 * (1 - Coeff_x) * (r_mat[i + 1][j] - 2.0 * r_mat[i][j] + r_mat[i - 1][j]));
		// 		l_mat_n[i][j] = l_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i][j + 1] - s_mat[i][j - 1])) + 0.5 * (1 - Coeff_x) * (l_mat[i][j + 1] - 2.0 * l_mat[i][j] + l_mat[i][j - 1]));
		// 		s_mat_n[i][j] = s_mat[i][j] + (Coeff_x * (0.5 * (r_mat[i + 1][j] - r_mat[i - 1][j]) + 0.5 * (1 - Coeff_x) * (s_mat[i + 1][j] - 2 * s_mat[i][j] + s_mat[i - 1][j])))\
		// 			+ (Coeff_x * (0.5 * (l_mat[i][j + 1] - l_mat[i][j - 1]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j + 1] - 2 * s_mat[i][j] + s_mat[i][j - 1])));
		// 		process_count++;
		// 	}
		// }
		for (j = 0; j <= end; j++) {
			for (i = 0; i <= X; i++) {
				u_n[i][j] = u[i][j] + 0.5 * time_step * (s_mat_n[i][j] - s_mat[i][j]);
			}
		}
		//for (j = 0; j <= end; j++) {
		//	for (i = 0; i <= X; i++) {
		//		u_p[i][j] = u[i][j];
		//		u[i][j] = u_n[i][j];
		//		r_mat[i][j] = r_mat_n[i][j];
		//		l_mat[i][j] = l_mat_n[i][j];
		//		s_mat[i][j] = s_mat_n[i][j];
		//		s_mat_p[i][j] = s_mat[i][j];

		//	}
		//}

		if (t == 0) {
			FILE* fpWrite = fopen(names[n], "w");
			if (fpWrite == NULL)
			{
				printf("Error");
				return 0;
			}
			for (j = 0; j <= end; j++) {
				for (i = 0; i <= X; i++) {
					fprintf(fpWrite, "%lf ", u[i][j]);
				}
				fprintf(fpWrite, "\n");
			}
			fclose(fpWrite);

		}
		else if (t % 5 == 0) {
			FILE* fpWrite = fopen(names[n], "a");
			if (fpWrite == NULL)
			{
				return 0;
			}
			for (j = 0; j <= end; j++) {
				for (i = 0; i <= X; i++) {
					fprintf(fpWrite, "%lf ", u[i][j]);
				}
				fprintf(fpWrite, "\n");
			}
			fclose(fpWrite);
		}

		break;

	case N-1:
		// for (j = start; j <= Y - 1; j++) {
		// 	for (i = 1; i <= X - 1; i++) {
		// 		r_mat_n[i][j] = r_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i + 1][j] - s_mat[i - 1][j])) + 0.5 * (1 - Coeff_x) * (r_mat[i + 1][j] - 2.0 * r_mat[i][j] + r_mat[i - 1][j]));
		// 		l_mat_n[i][j] = l_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i][j + 1] - s_mat[i][j - 1])) + 0.5 * (1 - Coeff_x) * (l_mat[i][j + 1] - 2.0 * l_mat[i][j] + l_mat[i][j - 1]));
		// 		s_mat_n[i][j] = s_mat[i][j] + (Coeff_x * (0.5 * (r_mat[i + 1][j] - r_mat[i - 1][j]) + 0.5 * (1 - Coeff_x) * (s_mat[i + 1][j] - 2 * s_mat[i][j] + s_mat[i - 1][j])))\
		// 			+ (Coeff_x * (0.5 * (l_mat[i][j + 1] - l_mat[i][j - 1]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j + 1] - 2 * s_mat[i][j] + s_mat[i][j - 1])));
		// 		process_count++;
		// 	}
		// }
		for (j = start; j <= Y; j++) {
			for (i = 0; i <= X; i++) {
				u_n[i][j] = u[i][j] + 0.5 * time_step * (s_mat_n[i][j] - s_mat[i][j]);
			}
		}
		//for (j = start; j <= Y; j++) {
		//	for (i = 0; i <= X; i++) {
		//		u_p[i][j] = u[i][j];
		//		u[i][j] = u_n[i][j];
		//		r_mat[i][j] = r_mat_n[i][j];
		//		l_mat[i][j] = l_mat_n[i][j];
		//		s_mat[i][j] = s_mat_n[i][j];
		//		s_mat_p[i][j] = s_mat[i][j];
		//	}
		//}
		if (t == 0) {
			FILE* fpWrite = fopen(names[n], "w");
			if (fpWrite == NULL)
			{
				printf("Error");
				return 0;
			}
			for (j = start; j <= Y; j++) {
				for (i = 0; i <= X; i++) {
					fprintf(fpWrite, "%lf ", u[i][j]);
				}
				fprintf(fpWrite, "\n");
			}
			fclose(fpWrite);
		}
		else if (t % 5 == 0) {
			FILE* fpWrite = fopen(names[n], "a");
			if (fpWrite == NULL)
			{
				return 0;
			}
			for (j = start; j <= Y; j++) {
				for (i = 0; i <= X; i++) {
					fprintf(fpWrite, "%lf ", u[i][j]);
				}
				fprintf(fpWrite, "\n");
			}
			fclose(fpWrite);
		}
		break;

	default:
		// for (j = start; j <= end; j++) {
		// 	for (i = 1; i <= X - 1; i++) {
		// 		r_mat_n[i][j] = r_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i + 1][j] - s_mat[i - 1][j])) + 0.5 * (1 - Coeff_x) * (r_mat[i + 1][j] - 2.0 * r_mat[i][j] + r_mat[i - 1][j]));
		// 		l_mat_n[i][j] = l_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i][j + 1] - s_mat[i][j - 1])) + 0.5 * (1 - Coeff_x) * (l_mat[i][j + 1] - 2.0 * l_mat[i][j] + l_mat[i][j - 1]));
		// 		s_mat_n[i][j] = s_mat[i][j] + (Coeff_x * (0.5 * (r_mat[i + 1][j] - r_mat[i - 1][j]) + 0.5 * (1 - Coeff_x) * (s_mat[i + 1][j] - 2 * s_mat[i][j] + s_mat[i - 1][j])))\
		// 			+ (Coeff_x * (0.5 * (l_mat[i][j + 1] - l_mat[i][j - 1]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j + 1] - 2 * s_mat[i][j] + s_mat[i][j - 1])));
		// 		process_count++;
		// 	}
		// }
		for (j = start; j <= end; j++) {
			for (i = 0; i <= X; i++) {
				u_n[i][j] = u[i][j] + 0.5 * time_step * (s_mat_n[i][j] - s_mat[i][j]);
			}
		}
		//for (j = start; j <= end; j++) {
		//	for (i = 0; i <= X; i++) {
		//		u_p[i][j] = u[i][j];
		//		u[i][j] = u_n[i][j];
		//		r_mat[i][j] = r_mat_n[i][j];
		//		l_mat[i][j] = l_mat_n[i][j];
		//		s_mat[i][j] = s_mat_n[i][j];
		//		s_mat_p[i][j] = s_mat[i][j];
		//	}
		//}
		if (t == 0) {
			FILE* fpWrite = fopen(names[n], "w");
			if (fpWrite == NULL)
			{
				printf("Error");
				return 0;
			}
			for (j = start; j <= end; j++) {
				for (i = 0; i <= X; i++) {
					fprintf(fpWrite, "%lf ", u[i][j]);
				}
				fprintf(fpWrite, "\n");
			}
			fclose(fpWrite);
		}
		else if (t % 5 == 0) {
			FILE* fpWrite = fopen(names[n], "a");
			if (fpWrite == NULL)
			{
				return 0;
			}
			for (j = start; j <= end; j++) {
				for (i = 0; i <= X; i++) {
					fprintf(fpWrite, "%lf ", u[i][j]);
				}
				fprintf(fpWrite, "\n");
			}
			fclose(fpWrite);
		}
		break;
	
	}
	//pthread_mutex_unlock(&test_mutex);
	gettimeofday(&t_fun_end,NULL);
	t_fun_th =  ((long long)t_fun_end.tv_sec-(long long)t_fun_start.tv_sec)*1000+ (t_fun_end.tv_usec-t_fun_start.tv_usec)/1000;
	printf("t_fun[%d][%d]= %ld ms\n",t,n,t_fun_th );
	
}

void* lax_wendroff_reassign(void* arg){
	int            n =  (int) (intptr_t) arg;  //Nth thread
	struct  timeval t_fun_start,t_fun_end;
	double t_fun_th;
	double      start = n * AVE+1;
	double      end = start + AVE - 1;
	long long      i, j;
	for (j = start; j <= end; j++) {
		for (i = 0; i <= X; i++) {
			u_p[i][j] = u[i][j];
			u[i][j] = u_n[i][j];
			r_mat[i][j] = r_mat_n[i][j];
			l_mat[i][j] = l_mat_n[i][j];
			s_mat[i][j] = s_mat_n[i][j];
			s_mat_p[i][j] = s_mat[i][j];
		}
	}
}
