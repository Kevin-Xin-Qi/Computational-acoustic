
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"Constants_pthread.h"

#include<time.h>
#include<pthread.h>
//#include<iostream>
//#include <sysinfoapi.h>

//pthread_mutex_t test_mutex;
#define N   3
#define AVE 17*4*3/N
#define X AVE*N+1
#define Y AVE*N+1
#define Z 50

int t;
double time_step = step_ratio * grid_size_x / c0;
int Number_of_thread = X * Y;
double Coeff_x = step_ratio;






char* names1[] = { "wave_data_lax_wendroff_t1.txt" };
char* names2[] = { "wave_data_lax_wendroff_t2.txt" };
char* names3[] = { "wave_data_lax_wendroff_t3.txt" };
char* names4[] = { "wave_data_lax_wendroff_t4.txt" };


char* names[4];

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

	//	clock_t start_time = time(NULL);
	double t0, t1, t2, t3, t_cal_total, t_file_total;
	t_cal_total = 0;
	t_file_total = 0;
	//Declare Function
	double leap_frog(double u_left, double u_right, double u_upper, double u_lower, double u_centre, double u_previous, double square_of_Coeff_x);

	double s_to_u(double u, double s_mat_n, double s_mat, double time_step);
	void* lax_wendroff_advection_all_in_1(void* arg);

	double Time_column[T + 1][2];
	float  viscosity;

	double   square_of_c, t0_g, tp_g, square_of_Coeff_x, Q_x, y;
	int i, j, k, P_coeff, id;

	P_coeff = 1024 * 1024;

	pthread_t* pthread_id = NULL; //
	pthread_id = (pthread_t*)malloc(N * sizeof(pthread_t));

	square_of_Coeff_x = Coeff_x * Coeff_x;

	Q_x = (grid_size_x - c0 * step_ratio) / (grid_size_x + c0 * step_ratio);
	y = Y;
	double mid_Y = 0.5 * Y;


	double mid_X = 0.5 * X;

	
	char names1[] = { "wave_data_lax_wendroff_t1.txt" };
	names[0] = &names1;
	
	char names2[] = { "wave_data_lax_wendroff_t2.txt" };
	names[1] = &names2;
	char names3[] = { "wave_data_lax_wendroff_t3.txt" };
	names[2] = &names3;
	char names4[] = { "wave_data_lax_wendroff_t4.txt" };
	names[3] = &names4;


	





	
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


	void* lax_wendroff_advection(void* arg);
	tp_g = 15;
	t0_g = 5;
	/*Matrix set up*/


	for (j = 0; j <= Y; j++) {
		for (i = 0; i <= X; i++) {
			s_mat_p[i][j] = r_mat_n[i][j] = s_mat_n[i][j] = l_mat_n[i][j] = r_mat[i][j] = s_mat[i][j] = l_mat[i][j] = u_p[i][j] = u[i][j] = u_n[i][j] = 0;


		}
	}





	//pthread_mutex_init(&test_mutex, NULL);

	//time_t time_1 =  time(NULL);
	/*Main computation*/
	t0 = clock();

	for (t = 0; t <= T; t++) {
		t1 = clock();
		u[(int)mid_X][(int)mid_Y] = P_coeff * exp(-pow((t - t0_g + 1) / tp_g, 2));
		if (t == 1) {
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
		for (id = 0; id < N; id++)
		{
			pthread_create(pthread_id + id, NULL, lax_wendroff_advection_all_in_1, id);
		}

		for (id = 0; id < N; id++)
		{

			pthread_join(pthread_id[id], NULL);

		}






		t2 = clock();
		// ----------------------FILE writing-------------------------------
		//if (t%10 == 1) {
			//if (t == 1) {
			//	FILE* fpWrite = fopen("wave_data_lax_wendroff_t.txt", "w");
			//	if (fpWrite == NULL)
			//	{
			//		return 0;
			//	}
			//	for (j = 0; j <= Y; j++) {
			//		for (i = 0; i <= X; i++) {
			//			fprintf(fpWrite, "%lf ", u[i][j]);
			//		}
			//		fprintf(fpWrite, "\n");
			//	}
			//	fclose(fpWrite);
			//}
			//else {
			//	FILE* fpWrite = fopen("wave_data_lax_wendroff_t.txt", "a");
			//	if (fpWrite == NULL)
			//	{
			//		return 0;
			//	}
			//	for (j = 0; j <= Y; j++) {
			//		for (i = 0; i <= X; i++) {
			//			fprintf(fpWrite, "%lf ", u[i][j]);
			//		}
			//		fprintf(fpWrite, "\n");
			//	}
			//	fclose(fpWrite);
			//}
		//}
		// ----------------------FILE writing END-------------------------------
		t3 = clock();

		t_cal_total += (t2 - t1) / CLOCKS_PER_SEC;
		t_file_total += (t3 - t2) / CLOCKS_PER_SEC;
		//printf("%f s  and   %f s and %f s and %f s\n", (t2 - t1) / CLOCKS_PER_SEC, (t3 - t2) / CLOCKS_PER_SEC,(t3-t0)/CLOCKS_PER_SEC, (t3 - t1) / CLOCKS_PER_SEC);

	}

	printf("%f s  and   %f s and %f s and %f s\n", (t2 - t1) / CLOCKS_PER_SEC, (t3 - t2) / CLOCKS_PER_SEC, (t3 - t0) / CLOCKS_PER_SEC, (t3 - t1) / CLOCKS_PER_SEC);
	free(u_n);
	free(u_p);
	free(u);
	free(r_mat);
	free(r_mat_n);
	free(s_mat);
	free(s_mat_n);
	free(l_mat);
	free(l_mat_n);
	free(s_mat_p);
	//pthread_mutex_destroy(&test_mutex);

	//for (j = 0; j <= Y; j++) {
	//	for (i = 0; i <= X; i++) {
	//		fprintf(fpWrite,"%lf ", u[i][j]);
	//	}
	//	fprintf(fpWrite,"\n");
	//}





}



void* lax_wendroff_advection_all_in_1(void* arg) {
	int            n = (int)arg;  //Nth thread

	double      start = n * AVE+1;
	double      end = start + AVE - 1;
	long long      i, j;
	switch (n) {
	case 0:
		for (j = 1; j <= end; j++) {
			for (i = 1; i <= X - 1; i++) {
				r_mat_n[i][j] = r_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i + 1][j] - s_mat[i - 1][j])) + 0.5 * (1 - Coeff_x) * (r_mat[i + 1][j] - 2.0 * r_mat[i][j] + r_mat[i - 1][j]));
				l_mat_n[i][j] = l_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i][j + 1] - s_mat[i][j - 1])) + 0.5 * (1 - Coeff_x) * (l_mat[i][j + 1] - 2.0 * l_mat[i][j] + l_mat[i][j - 1]));
				s_mat_n[i][j] = s_mat[i][j] + (Coeff_x * (0.5 * (r_mat[i + 1][j] - r_mat[i - 1][j]) + 0.5 * (1 - Coeff_x) * (s_mat[i + 1][j] - 2 * s_mat[i][j] + s_mat[i - 1][j])))\
					+ (Coeff_x * (0.5 * (l_mat[i][j + 1] - l_mat[i][j - 1]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j + 1] - 2 * s_mat[i][j] + s_mat[i][j - 1])));
			}
		}
		for (j = 0; j <= end; j++) {
			for (i = 0; i <= X; i++) {
				u_n[i][j] = u[i][j] + 0.5 * time_step * (s_mat_n[i][j] - s_mat[i][j]);
			}
		}
		for (j = 0; j <= end; j++) {
			for (i = 0; i <= X; i++) {
				u_p[i][j] = u[i][j];
				u[i][j] = u_n[i][j];
				r_mat[i][j] = r_mat_n[i][j];
				l_mat[i][j] = l_mat_n[i][j];
				s_mat[i][j] = s_mat_n[i][j];
				s_mat_p[i][j] = s_mat[i][j];

			}
		}

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
		else if (t % 25 == 0) {
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

	case 2:
		for (j = start; j <= Y - 1; j++) {
			for (i = 1; i <= X - 1; i++) {
				r_mat_n[i][j] = r_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i + 1][j] - s_mat[i - 1][j])) + 0.5 * (1 - Coeff_x) * (r_mat[i + 1][j] - 2.0 * r_mat[i][j] + r_mat[i - 1][j]));
				l_mat_n[i][j] = l_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i][j + 1] - s_mat[i][j - 1])) + 0.5 * (1 - Coeff_x) * (l_mat[i][j + 1] - 2.0 * l_mat[i][j] + l_mat[i][j - 1]));
				s_mat_n[i][j] = s_mat[i][j] + (Coeff_x * (0.5 * (r_mat[i + 1][j] - r_mat[i - 1][j]) + 0.5 * (1 - Coeff_x) * (s_mat[i + 1][j] - 2 * s_mat[i][j] + s_mat[i - 1][j])))\
					+ (Coeff_x * (0.5 * (l_mat[i][j + 1] - l_mat[i][j - 1]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j + 1] - 2 * s_mat[i][j] + s_mat[i][j - 1])));
			}
		}
		for (j = start; j <= Y; j++) {
			for (i = 0; i <= X; i++) {
				u_n[i][j] = u[i][j] + 0.5 * time_step * (s_mat_n[i][j] - s_mat[i][j]);
			}
		}
		for (j = start; j <= Y; j++) {
			for (i = 0; i <= X; i++) {
				u_p[i][j] = u[i][j];
				u[i][j] = u_n[i][j];
				r_mat[i][j] = r_mat_n[i][j];
				l_mat[i][j] = l_mat_n[i][j];
				s_mat[i][j] = s_mat_n[i][j];
				s_mat_p[i][j] = s_mat[i][j];
			}
		}
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
		else if (t % 25 == 0) {
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
		for (j = start; j <= end; j++) {
			for (i = 1; i <= X - 1; i++) {
				r_mat_n[i][j] = r_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i + 1][j] - s_mat[i - 1][j])) + 0.5 * (1 - Coeff_x) * (r_mat[i + 1][j] - 2.0 * r_mat[i][j] + r_mat[i - 1][j]));
				l_mat_n[i][j] = l_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i][j + 1] - s_mat[i][j - 1])) + 0.5 * (1 - Coeff_x) * (l_mat[i][j + 1] - 2.0 * l_mat[i][j] + l_mat[i][j - 1]));
				s_mat_n[i][j] = s_mat[i][j] + (Coeff_x * (0.5 * (r_mat[i + 1][j] - r_mat[i - 1][j]) + 0.5 * (1 - Coeff_x) * (s_mat[i + 1][j] - 2 * s_mat[i][j] + s_mat[i - 1][j])))\
					+ (Coeff_x * (0.5 * (l_mat[i][j + 1] - l_mat[i][j - 1]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j + 1] - 2 * s_mat[i][j] + s_mat[i][j - 1])));
			}
		}
		for (j = start; j <= end; j++) {
			for (i = 0; i <= X; i++) {
				u_n[i][j] = u[i][j] + 0.5 * time_step * (s_mat_n[i][j] - s_mat[i][j]);
			}
		}
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
		else if (t % 25 == 0) {
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





}
