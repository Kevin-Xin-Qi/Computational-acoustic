
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>
#define _CRT_SECURE_NO_WARNINGS
#define T 50
#define M_PIl 3.141592653589793238462643383279502884L /* pi */
#define grid_size_x 0.01
#define grid_size_y 0.01
#define step_ratio 0.26
#define c0 340
#include<pthread.h>
#include <stdint.h>

#define N   3
#define AVE 12*17/N
#define X AVE*N+1
#define Y AVE*N+1
#define Z AVE*N+1


double time_step = step_ratio * grid_size_x / c0;
double Coeff_x ;
double*** u_n = NULL;
double*** u = NULL;
double*** u_p = NULL;
double*** r_mat = NULL;
double*** r_mat_n = NULL;
double*** s_mat = NULL;
double*** s_mat_n = NULL;
double*** s_mat_p = NULL;
double*** l_mat = NULL;
double*** l_mat_n = NULL;
double*** q_mat = NULL;
double*** q_mat_n = NULL;

char* names[8] ;
int  t;

#pragma comment(lib,"pthreadVC2.lib")
void main()
{
	double t0, t1, t2, t3, t_cal_total, t_file_total;



	 Coeff_x = c0 * time_step / grid_size_x;

	void* lax_wendroff_advection_all_in_1(void* arg);
	void* lax_wendroff_reassign(void* arg);

	double   square_of_c,  time, t0_g, tp_g, square_of_Coeff_x, Q_x;
	int i, j, k, id;
	double P_coeff = 1024 * 1024;

	pthread_t* pthread_id = NULL; //
	pthread_id = (pthread_t*)malloc(N * sizeof(pthread_t));

	char names1[] = { "wave_data_lax_wendroff_3D_t1.txt" };
	names[0] = names1;
	char names2[] = { "wave_data_lax_wendroff_3D_t2.txt" };
	names[1] = names2;
	char names3[] = { "wave_data_lax_wendroff_3D_t3.txt" };
	names[2] = names3;
	char names4[] = { "wave_data_lax_wendroff_3D_t4.txt" };
	names[3] = names4;
	char names5[] = { "wave_data_lax_wendroff_3D_t5.txt" };
	names[4] = names5;
	char names6[] = { "wave_data_lax_wendroff_3D_t6.txt" };
	names[5] = names6;
	char names7[] = { "wave_data_lax_wendroff_3D_t7.txt" };
	names[6] = names7;
	char names8[] = { "wave_data_lax_wendroff_3D_t8.txt" };
	names[7] = names8;


	
	u_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u_n[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			u_n[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	
	u = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			u[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	
	u_p = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u_p[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			u_p[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	
	r_mat = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		r_mat[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			r_mat[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	
	r_mat_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		r_mat_n[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			r_mat_n[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}

	
	s_mat = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			s_mat[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	
	s_mat_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat_n[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			s_mat_n[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	
	s_mat_p = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat_p[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			s_mat_p[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	
	l_mat = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		l_mat[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			l_mat[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	
	l_mat_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		l_mat_n[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			l_mat_n[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	
	q_mat = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		q_mat[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			q_mat[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	
	q_mat_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		q_mat_n[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			q_mat_n[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}








	square_of_Coeff_x = Coeff_x * Coeff_x;

	Q_x = (grid_size_x - c0 * step_ratio) / (grid_size_x + c0 * step_ratio);

	double mid_Y = 0.5 * Y;
	double mid_X = 0.5 * X;
	double mid_Z = 0.5 * Z;





	tp_g = 15;
	t0_g = 5;
	/*Matrix initialisation*/

	for (k = 0; k <= Z; k++) {
		for (j = 0; j <= Y; j++) {
			for (i = 0; i <= X; i++) {
				q_mat[i][j][k]= q_mat_n[i][j][k]=s_mat_p[i][j][k] = r_mat_n[i][j][k] = s_mat_n[i][j][k] = l_mat_n[i][j][k] = r_mat[i][j][k] = s_mat[i][j][k] = l_mat[i][j][k] = u_p[i][j][k] = u[i][j][k] = u_n[i][j][k] = 0;
			}
		}
	}


	for (t = 0; t <= T; t++) {
		u[(int)mid_X][(int)mid_Y][(int)mid_Z] = P_coeff * exp(-pow((t - t0_g + 1) / tp_g, 2));
		if (t == 0) {
			for (k = 1; k <= Z - 1; k++) {
				for (j = 1; j <= Y - 1; j++) {
					for (i = 1; i <= X - 1; i++) {
						u_n[i][j][k] = 2 * u[i][j][k] - u_p[i][j][k] + square_of_Coeff_x * (u[i - 1][j][k] + u[i][j - 1][k] + u[i + 1][j][k] + u[i][j + 1][k] + u[i][j][k + 1] + u[i][j][k - 1] - 6 * u[i][j][k]);
						//printf("%f",u_n[i][j][k]);
					}
				}
			}
		}

		r_mat[(int)mid_X][(int)mid_Y][(int)mid_Z] = c0 * (u_n[(int)mid_X + 1][(int)mid_Y][(int)mid_Z] - u_n[(int)mid_X - 1][(int)mid_Y][(int)mid_Z]) / (2.0 * grid_size_x);
		l_mat[(int)mid_X][(int)mid_Y][(int)mid_Z] = c0 * (u_n[(int)mid_X][(int)mid_Y + 1][(int)mid_Z] - u_n[(int)mid_X][(int)mid_Y - 1][(int)mid_Z]) / (2.0 * grid_size_x);
		q_mat[(int)mid_X][(int)mid_Y][(int)mid_Z] = c0 * (u_n[(int)mid_X][(int)mid_Y][(int)mid_Z + 1] - u_n[(int)mid_X][(int)mid_Y][(int)mid_Z - 1]) / (2.0 * grid_size_x);
		s_mat[(int)mid_X][(int)mid_Y][(int)mid_Z] = 2 * (u_n[(int)mid_X][(int)mid_Y][(int)mid_Z] - u_p[(int)mid_X][(int)mid_Y][(int)mid_Z]) / time_step + s_mat_p[(int)mid_X][(int)mid_Y][(int)mid_Z];

		for (id = 0; id < N; id++)
		{
			pthread_create(pthread_id + id, NULL, lax_wendroff_advection_all_in_1, (void*)(intptr_t)id);
		}

		for (id = 0; id < N; id++)
		{
			pthread_join(pthread_id[id], NULL);
		}


		 for (id = 0; id < N; id++)
		 {
		 	pthread_create(pthread_id + id, NULL, lax_wendroff_reassign, (void *) (intptr_t) id);
		 }

		 for (id = 0; id < N; id++)
		 {

		 	pthread_join(pthread_id[id], NULL);

		 }
	}


}






void* lax_wendroff_advection_all_in_1(void* arg) {
	int            n = (int)(intptr_t)arg;  //Nth thread

	double      start = n * AVE + 1;
	double      end = start + AVE - 1;
	long long      i, j,k;
	double va, vb, vc, vd, ve,s_mat_c, s_mat_l, s_mat_r;
	switch (n) {
	case 0:
		for (k = 1; k <= end; k++) {
			for (j = 1; j <= Y - 1; j++) {
				for (i = 1; i <= X - 1; i++) {
					r_mat_n[i][j][k] = r_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i + 1][j][k] - s_mat[i - 1][j][k])) + 0.5 * Coeff_x * (r_mat[i + 1][j][k] - 2.0 * r_mat[i][j][k] + r_mat[i - 1][j][k]));
					l_mat_n[i][j][k] = l_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i][j + 1][k] - s_mat[i][j - 1][k])) + 0.5 * Coeff_x * (l_mat[i][j + 1][k] - 2.0 * l_mat[i][j][k] + l_mat[i][j - 1][k]));
					q_mat_n[i][j][k] = q_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i][j][k + 1] - s_mat[i][j][k - 1])) + 0.5 * Coeff_x * (q_mat[i][j][k + 1] - 2.0 * q_mat[i][j][k] + q_mat[i][j][k - 1]));
					s_mat_n[i][j][k] = s_mat[i][j][k] + (Coeff_x * (0.5 * (r_mat[i + 1][j][k] - r_mat[i - 1][j][k]) + 0.5 * (1 - Coeff_x) * (s_mat[i + 1][j][k] - 2 * s_mat[i][j][k] + s_mat[i - 1][j][k])))\
						+ (Coeff_x * (0.5 * (l_mat[i][j + 1][k] - l_mat[i][j - 1][k]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j + 1][k] - 2 * s_mat[i][j][k] + s_mat[i][j - 1][k])))\
						+ (Coeff_x * (0.5 * (q_mat[i][j][k + 1] - q_mat[i][j][k - 1]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j][k + 1] - 2 * s_mat[i][j][k] + s_mat[i][j][k - 1])));
				}
			}
		}
		for (k = 0; k <= end; k++) {
			for (j = 0; j <= Y; j++) {
				for (i = 0; i <= X; i++) {
					u_n[i][j][k] = u[i][j][k] + 0.5 * time_step * (s_mat_n[i][j][k] - s_mat[i][j][k]);
				}
			}
		}

		// if (t == 0) {
		// 	FILE* fpWrite = fopen(names[n], "w");
		// 	if (fpWrite == NULL)
		// 	{
		// 		printf("Error");
		// 		return 0;
		// 	}
		// 	for (k = 0; k <= end; k++) {
		// 		for (j = 0; j <= Y; j++) {
		// 			for (i = 0; i <= X; i++) {
		// 				fprintf(fpWrite, "%lf ", u_n[i][j][k]);
		// 			}
		// 			fprintf(fpWrite, "\n");
		// 		}
		// 	}
		// 	fclose(fpWrite);
		// }
		// else if (t % 5 == 0) {
		// 	FILE* fpWrite = fopen(names[n], "a");
		// 	if (fpWrite == NULL)
		// 	{
		// 		return 0;
		// 	}
		// 	for (k = 0; k <= end; k++) {
		// 		for (j = 0; j <= Y; j++) {
		// 			for (i = 0; i <= X; i++) {
		// 				fprintf(fpWrite, "%lf ", u_n[i][j][k]);
		// 			}
		// 			fprintf(fpWrite, "\n");
		// 		}
		// 	}
		// 	fclose(fpWrite);
		// }
		break;

	case N-1:
		for (k = start; k <= Z-1; k++) {
			for (j = 1; j <= Y - 1; j++) {
				for (i = 1; i <= X - 1; i++) {

					r_mat_n[i][j][k] = r_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i + 1][j][k] - s_mat[i - 1][j][k])) + 0.5 * Coeff_x * (r_mat[i + 1][j][k] - 2.0 * r_mat[i][j][k] + r_mat[i - 1][j][k]));
					l_mat_n[i][j][k] = l_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i][j + 1][k] - s_mat[i][j - 1][k])) + 0.5 * Coeff_x * (l_mat[i][j + 1][k] - 2.0 * l_mat[i][j][k] + l_mat[i][j - 1][k]));
					q_mat_n[i][j][k] = q_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i][j][k + 1] - s_mat[i][j][k - 1])) + 0.5 * Coeff_x * (q_mat[i][j][k + 1] - 2.0 * q_mat[i][j][k] + q_mat[i][j][k - 1]));
					s_mat_n[i][j][k] = s_mat[i][j][k] + (Coeff_x * (0.5 * (r_mat[i + 1][j][k] - r_mat[i - 1][j][k]) + 0.5 * (1 - Coeff_x) * (s_mat[i + 1][j][k] - 2 * s_mat[i][j][k] + s_mat[i - 1][j][k])))\
						+ (Coeff_x * (0.5 * (l_mat[i][j + 1][k] - l_mat[i][j - 1][k]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j + 1][k] - 2 * s_mat[i][j][k] + s_mat[i][j - 1][k])))\
						+ (Coeff_x * (0.5 * (q_mat[i][j][k + 1] - q_mat[i][j][k - 1]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j][k + 1] - 2 * s_mat[i][j][k] + s_mat[i][j][k - 1])));
				}
			}
		}
		for (k = start; k <= Z; k++) {
			for (j = 0; j <= Y; j++) {
				for (i = 0; i <= X; i++) {
					u_n[i][j][k] = u[i][j][k] + 0.5 * time_step * (s_mat_n[i][j][k] - s_mat[i][j][k]);
				}
			}
		}

		// if (t == 0) {
		// 	FILE* fpWrite = fopen(names[n], "w");
		// 	if (fpWrite == NULL)
		// 	{
		// 		printf("Error");
		// 		return 0;
		// 	}
		// 	for (k = start; k <= Z; k++) {
		// 		for (j = 0; j <= Y; j++) {
		// 			for (i = 0; i <= X; i++) {
		// 				fprintf(fpWrite, "%lf ", u_n[i][j][k]);
		// 			}
		// 			fprintf(fpWrite, "\n");
		// 		}
		// 	}
		// 	fclose(fpWrite);
		// }
		// else if (t % 5 == 0) {
		// 	FILE* fpWrite = fopen(names[n], "a");
		// 	if (fpWrite == NULL)
		// 	{
		// 		return 0;
		// 	}
		// 	for (k = start; k <= Z; k++) {
		// 		for (j = 0; j <= Y; j++) {
		// 			for (i = 0; i <= X; i++) {
		// 				fprintf(fpWrite, "%lf ", u_n[i][j][k]);
		// 			}
		// 			fprintf(fpWrite, "\n");
		// 		}
		// 	}
		// 	fclose(fpWrite);
		// }
		break;

	default:
		for (k = start; k <= end ; k++) {
			for (j = 1; j <= Y - 1; j++) {
				for (i = 1; i <= X - 1; i++) {

					r_mat_n[i][j][k] = r_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i + 1][j][k] - s_mat[i - 1][j][k])) + 0.5 * Coeff_x * (r_mat[i + 1][j][k] - 2.0 * r_mat[i][j][k] + r_mat[i - 1][j][k]));
                    l_mat_n[i][j][k] = l_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i][j + 1][k] - s_mat[i][j - 1][k])) + 0.5 * Coeff_x * (l_mat[i][j + 1][k] - 2.0 * l_mat[i][j][k] + l_mat[i][j - 1][k]));
					q_mat_n[i][j][k] = q_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i][j][k + 1] - s_mat[i][j][k - 1])) + 0.5 * Coeff_x * (q_mat[i][j][k + 1] - 2.0 * q_mat[i][j][k] + q_mat[i][j][k - 1]));
					s_mat_n[i][j][k] = s_mat[i][j][k] + (Coeff_x * (0.5 * (r_mat[i + 1][j][k] - r_mat[i - 1][j][k]) + 0.5 * (1 - Coeff_x) * (s_mat[i + 1][j][k] - 2 * s_mat[i][j][k] + s_mat[i - 1][j][k])))\
						+ (Coeff_x * (0.5 * (l_mat[i][j + 1][k] - l_mat[i][j - 1][k]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j + 1][k] - 2 * s_mat[i][j][k] + s_mat[i][j - 1][k])))\
						+ (Coeff_x * (0.5 * (q_mat[i][j][k + 1] - q_mat[i][j][k - 1]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j][k + 1] - 2 * s_mat[i][j][k] + s_mat[i][j][k - 1])));

				}
				
			}
		}
		for (k = start; k <= end; k++) {
			for (j = 0; j <= Y; j++) {
				for (i = 0; i <= X; i++) {
					u_n[i][j][k] = u[i][j][k] + 0.5 * time_step * (s_mat_n[i][j][k] - s_mat[i][j][k]);
					//if (i== 14 && j==14 && k==14) {
					//	printf("%f ", u_n[i][j][k]);
					//	printf("\n");
					//}
					
				}
			}
		}

		// if (t == 0) {
		// 	FILE* fpWrite = fopen(names[n], "w");
		// 	if (fpWrite == NULL)
		// 	{
		// 		printf("Error");
		// 		return 0;
		// 	}
		// 	for (k = start; k <= end; k++) {
		// 		for (j = 0; j <= Y; j++) {
		// 			for (i = 0; i <= X; i++) {
		// 				fprintf(fpWrite, "%lf ", u_n[i][j][k]);
		// 			}
		// 			fprintf(fpWrite, "\n");
		// 		}
		// 	}
		// 	fclose(fpWrite);
		// }
		// else if (t % 5 == 0) {
		// 	FILE* fpWrite = fopen(names[n], "a");
		// 	if (fpWrite == NULL)
		// 	{
		// 		return 0;
		// 	}
		// 	for (k = start; k <= end; k++) {
		// 		for (j = 0; j <= Y; j++) {
		// 			for (i = 0; i <= X; i++) {
		// 				fprintf(fpWrite, "%lf ", u_n[i][j][k]);
		// 			}
		// 			fprintf(fpWrite, "\n");
		// 		}
		// 	}
		// 	fclose(fpWrite);
		// }
		// break;
	}
	//pthread_mutex_unlock(&test_mutex);
}


void* lax_wendroff_reassign(void* arg){
	int		n =  (int) (intptr_t) arg;  //Nth thread
	//struct  timeval t_fun_start,t_fun_end;
	double	t_fun_th;
	double	start = n * AVE+1;
	double	end = start + AVE - 1;
	long long      i, j,k;
	int		var_c;
	 //u_p[n][0] = u[n][0];
	 switch (n){
	 case 0:
	 		for (k = 0; k <= end; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
					u_p[i][j][k] = u[i][j][k];
					u[i][j][k] = u_n[i][j][k];
					r_mat[i][j][k] = r_mat_n[i][j][k];
					l_mat[i][j][k] = l_mat_n[i][j][k];
					q_mat[i][j][k] = q_mat_n[i][j][k];
					s_mat[i][j][k] = s_mat_n[i][j][k];
					s_mat_p[i][j][k] = s_mat[i][j][k];
					}
				}
			}
		
	 	break;

	 case N-1:
	 for (k = start; k <= Z; k++) {
	 		for (j = 0; j <= Y; j++) {
	 			for (i = 0; i <= X; i++) {
					u_p[i][j][k] = u[i][j][k];
					u[i][j][k] = u_n[i][j][k];
					r_mat[i][j][k] = r_mat_n[i][j][k];
					l_mat[i][j][k] = l_mat_n[i][j][k];
					q_mat[i][j][k] = q_mat_n[i][j][k];
					s_mat[i][j][k] = s_mat_n[i][j][k];
					s_mat_p[i][j][k] = s_mat[i][j][k];
	 			}
	 		}
	 	}
	 	break;
	 default:
	 for (k = start; k <= end; k++) {
	 		for (j = 0; j <= Y; j++) {
	 			for (i = 0; i <= X; i++) {
					u_p[i][j][k] = u[i][j][k];
					u[i][j][k] = u_n[i][j][k];
					r_mat[i][j][k] = r_mat_n[i][j][k];
					l_mat[i][j][k] = l_mat_n[i][j][k];
					q_mat[i][j][k] = q_mat_n[i][j][k];
					s_mat[i][j][k] = s_mat_n[i][j][k];
					s_mat_p[i][j][k] = s_mat[i][j][k];
	 			}
	 		}
	 	}
	 	break;
	 }
	return 0;
}