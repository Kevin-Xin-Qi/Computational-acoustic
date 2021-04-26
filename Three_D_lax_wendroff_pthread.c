
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>
#include"Lax3D_constants.h"
#include<pthread.h>
#include<time.h>

#define N   3
#define AVE 17*4*3/N
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

char* names[4] ;
int i, j, k, t, P_coeff, id;

#pragma comment(lib,"pthreadVC2.lib")
void main()
{
	double t0, t1, t2, t3, t_cal_total, t_file_total;
	t0 = clock();



	double Coeff_x = c0 * time_step / grid_size_x;

	//float  viscosity;
	void* lax_wendroff_advection_all_in_1(void* arg);

	double   square_of_c,  time, t0_g, tp_g, square_of_Coeff_x, Q_x;

	P_coeff = 1024 * 1024;

	pthread_t* pthread_id = NULL; //
	pthread_id = (pthread_t*)malloc(N * sizeof(pthread_t));

	char names1[] = { "wave_data_lax_wendroff_3D_t1.txt" };
	names[0] = &names1;
	char names2[] = { "wave_data_lax_wendroff_3D_t2.txt" };
	names[1] = &names2;
	char names3[] = { "wave_data_lax_wendroff_3D_t3.txt" };
	names[2] = &names3;
	char names4[] = { "wave_data_lax_wendroff_3D_t4.txt" };
	names[3] = &names4;


	
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
	/*Matrix set up*/

	for (k = 0; k <= Z; k++) {
		for (j = 0; j <= Y; j++) {
			for (i = 0; i <= X; i++) {
				s_mat_p[i][j][k] = r_mat_n[i][j][k] = s_mat_n[i][j][k] = l_mat_n[i][j][k] = r_mat[i][j][k] = s_mat[i][j][k] = l_mat[i][j][k] = u_p[i][j][k] = u[i][j][k] = u_n[i][j][k] = 0;
			}
		}
	}


	/*Initial condition*/
	//u[(int)mid_X][(int)mid_Y] = exp(-pow((1  -t0_g) / tp_g, 2));
	//u[(int)mid_X][(int)mid_Y] = sin(2 * M_PIl * c0 * time_step * t * 100);

	//for (j = 1; j <= Y-1; j++) {
	//	for (i = 1; i <= X-1; i++) {

			//r_mat[i][j] = c0 * (u[i + 1][j] - u[i - 1][j]) / (2.0 * grid_size_x);
			//l_mat[i][j] = c0 * (u[i][j + 1] - u[i][j - 1]) / (2.0 * grid_size_y);
			//s_mat[i][j] = 2 * (u[i][j] - u_p[i][j]) / time_step + s_mat_p[i][j];
	//	}
	//}




	/*Main computation*/
	//FILE* fpWrite = fopen("wave_data_lax_wendroff_3d.txt", "w");
	//if (fpWrite == NULL)
	//{
	//	return 0;
	//}
	/*Initial condition*/
	for (t = 0; t <= T; t++) {
		u[(int)mid_X][(int)mid_Y][(int)mid_Z] = P_coeff* exp(-pow((t - t0_g + 1) / tp_g, 2));
		if (t == 0) {
			for (k = 1; k <= Z - 1; k++) {
				for (j = 1; j <= Y - 1; j++) {
					for (i = 1; i <= X - 1; i++) {
						u_n[i][j][k] = 2 * u[i][j][k] - u_p[i][j][k] + square_of_Coeff_x * (u[i - 1][j][k] + u[i][j - 1][k] + u[i + 1][j][k] + u[i][j + 1][k] + u[i][j][k + 1] + u[i][j][k - 1] - 6 * u[i][j][k]);
					}
				}
			}
		}

		r_mat[(int)mid_X][(int)mid_Y][(int)mid_Z] = c0 * (u[(int)mid_X + 1][(int)mid_Y][(int)mid_Z] - u[(int)mid_X - 1][(int)mid_Y][(int)mid_Z]) / (2.0 * grid_size_x);
		l_mat[(int)mid_X][(int)mid_Y][(int)mid_Z] = c0 * (u[(int)mid_X][(int)mid_Y + 1][(int)mid_Z] - u[(int)mid_X][(int)mid_Y - 1][(int)mid_Z]) / (2.0 * grid_size_x);
		q_mat[(int)mid_X][(int)mid_Y][(int)mid_Z] = c0 * (u[(int)mid_X][(int)mid_Y][(int)mid_Z + 1] - u[(int)mid_X][(int)mid_Y][(int)mid_Z - 1]) / (2.0 * grid_size_x);
		s_mat[(int)mid_X][(int)mid_Y][(int)mid_Z] = 2 * (u[(int)mid_X][(int)mid_Y][(int)mid_Z] - u_p[(int)mid_X][(int)mid_Y][(int)mid_Z]) / time_step + s_mat_p[(int)mid_X][(int)mid_Y][(int)mid_Z];

		for (id = 0; id < N; id++)
		{
			pthread_create(pthread_id + id, NULL, lax_wendroff_advection_all_in_1, id);
		}

		for (id = 0; id < N; id++)
		{
			pthread_join(pthread_id[id], NULL);
		}
	}

	t3 = clock();
	printf("%f s  \n",  (t3 - t0) / CLOCKS_PER_SEC);
	free(u_n);
	free(u_p);
	free(u);
	free(r_mat);
	free(r_mat_n);
	free(q_mat);
	free(q_mat_n);
	free(s_mat);
	free(s_mat_n);
	free(l_mat);
	free(l_mat_n);
	free(s_mat_p);


}

void* lax_wendroff_advection_all_in_1(void* arg) {
	int            n = (int)arg;  //Nth thread

	double      start = n * AVE + 1;
	double      end = start + AVE - 1;
	long long      i, j,k;
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
		if (t == 0) {
			FILE* fpWrite = fopen(names[n], "w");
			if (fpWrite == NULL)
			{
				printf("Error");
				return 0;
			}
			for (k = 0; k <= end; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
						fprintf(fpWrite, "%lf ", u[i][j][k]);
					}
					fprintf(fpWrite, "\n");
				}
			}
			fclose(fpWrite);
		}
		else if (t % 25 == 0) {
			FILE* fpWrite = fopen(names[n], "a");
			if (fpWrite == NULL)
			{
				return 0;
			}
			for (k = 0; k <= end; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
						fprintf(fpWrite, "%lf ", u[i][j][k]);
					}
					fprintf(fpWrite, "\n");
				}
			}
			fclose(fpWrite);
		}
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
		if (t == 0) {
			FILE* fpWrite = fopen(names[n], "w");
			if (fpWrite == NULL)
			{
				printf("Error");
				return 0;
			}
			for (k = start; k <= Z; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
						fprintf(fpWrite, "%lf ", u[i][j][k]);
					}
					fprintf(fpWrite, "\n");
				}
			}
			fclose(fpWrite);
		}
		else if (t % 25 == 0) {
			FILE* fpWrite = fopen(names[n], "a");
			if (fpWrite == NULL)
			{
				return 0;
			}
			for (k = start; k <= Z; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
						fprintf(fpWrite, "%lf ", u[i][j][k]);
					}
					fprintf(fpWrite, "\n");
				}
			}
			fclose(fpWrite);
		}
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
				}
			}
		}
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
		if (t == 0) {
			FILE* fpWrite = fopen(names[n], "w");
			if (fpWrite == NULL)
			{
				printf("Error");
				return 0;
			}
			for (k = start; k <= end; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
						fprintf(fpWrite, "%lf ", u[i][j][k]);
					}
					fprintf(fpWrite, "\n");
				}
			}
			fclose(fpWrite);
		}
		else if (t % 25 == 0) {
			FILE* fpWrite = fopen(names[n], "a");
			if (fpWrite == NULL)
			{
				return 0;
			}
			for (k = start; k <= end; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
						fprintf(fpWrite, "%lf ", u[i][j][k]);
					}
					fprintf(fpWrite, "\n");
				}
			}
			fclose(fpWrite);
		}
		break;
	}
	//pthread_mutex_unlock(&test_mutex);
}

