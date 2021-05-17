#define _CRT_SECURE_NO_WARNINGS
#define X 8*9*5*33
#define Y 8*9*5*33
#define Z 50
#define T 50
#define M_PIl 3.141592653589793238462643383279502884L /* pi */
#define grid_size_x 0.01
#define grid_size_y 0.01
#define step_ratio 0.34
#define c0 340
#include<stdio.h>
#include<math.h>
#include<stdlib.h>



long double time_step = step_ratio * grid_size_x / c0;




void main()
{

	double Coeff_x = c0 * time_step / grid_size_x;

	float  viscosity;

	long double   square_of_c, ut, ue, vt, ve, time,  t0_g, tp_g, square_of_Coeff_x, Q_x;
	int i, j, k, t;


	long double ** u_n;
		u_n = malloc(sizeof(long double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u_n[i] = malloc(sizeof(long double) * (Y + 1));
	}
	long double** u = malloc(sizeof(long double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u[i] = malloc(sizeof(long double) * (Y + 1));
	}
	long double** u_p = malloc(sizeof(long double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u_p[i] = malloc(sizeof(long double) * (Y + 1));
	}
	long double** r_mat = malloc(sizeof(long double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		r_mat[i] = malloc(sizeof(long double) * (Y + 1));
	}
	long double** s_mat = malloc(sizeof(long double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat[i] = malloc(sizeof(long double) * (Y + 1));
	}
	long double** l_mat = malloc(sizeof(long double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		l_mat[i] = malloc(sizeof(long double) * (Y + 1));
	}
	long double** s_mat_n = malloc(sizeof(long double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat_n[i] = malloc(sizeof(long double) * (Y + 1));
	}
	long double** r_mat_n = malloc(sizeof(long double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		r_mat_n[i] = malloc(sizeof(long double) * (Y + 1));
	}
	long double** l_mat_n = malloc(sizeof(long double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		l_mat_n[i] = malloc(sizeof(long double) * (Y + 1));
	}
	long double** s_mat_p = malloc(sizeof(long double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat_p[i] = malloc(sizeof(long double) * (Y + 1));
	}


	
	
	
	square_of_Coeff_x = Coeff_x * Coeff_x;

	Q_x = (grid_size_x-c0* step_ratio) / (grid_size_x + c0 * step_ratio);
	double mid_Y = round(Y / 2.0);

	double mid_X = round(X / 2.0);


	



	tp_g = 15;
	t0_g = 5;
	/*Matrix set up*/

	
	for (j = 0; j <= Y; j++) {
		for (i = 0; i <= X; i++) {
			s_mat_p[i][j]=r_mat_n[i][j] = s_mat_n[i][j] = l_mat_n[i][j] = r_mat[i][j] = s_mat[i][j] = l_mat[i][j] = u_p[i][j] = u[i][j] = u_n[i][j] = 0;


		}
	}




	



	/*Main computation*/
	// FILE* fpWrite = fopen("wave_data_lax_wendroff.txt", "w");
	// if (fpWrite == NULL)
	// {
	// 	return 0;
	// }
	for (t = 0; t <= T; t++) {
		u[(int)mid_X][(int)mid_Y] = exp(-pow((t - t0_g+1) / tp_g, 2));
		if (t == 1){
			for (j = 1; j <= Y - 1; j++) {
				for (i = 1; i <= X - 1; i++) {

					u_n[i][j] = 2 * u[i][j] - u_p[i][j] + square_of_Coeff_x * (u[i - 1][j] + u[i][j - 1] + u[i + 1][j] + u[i][j + 1] - 4 * u[i][j]);


				}
			}



		}
			r_mat[(int)mid_X][(int)mid_Y] = c0 * (u[(int)mid_X +1][(int)mid_Y] - u[(int)mid_X - 1][(int)mid_Y]) / (2.0 * grid_size_x);
			l_mat[(int)mid_X][(int)mid_Y] = c0 * (u[(int)mid_X][(int)mid_Y+1 ] - u[(int)mid_X][(int)mid_Y - 1]) / (2.0 * grid_size_x);
			s_mat[(int)mid_X][(int)mid_Y] = 2 * (u[(int)mid_X][(int)mid_Y] - u_p[(int)mid_X][(int)mid_Y]) / time_step + s_mat_p[(int)mid_X][(int)mid_Y];







		for (j = 1; j <= Y - 1; j++) {
			for (i = 1; i <= X - 1; i++) {

				r_mat_n[i][j] = r_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i + 1][j] - s_mat[i - 1][j])) + 0.5 * Coeff_x * (r_mat[i + 1][j] - 2.0 * r_mat[i][j] + r_mat[i - 1][j]));
				l_mat_n[i][j] = l_mat[i][j] + (Coeff_x * (0.5 * (s_mat[i][j + 1] - s_mat[i][j - 1])) + 0.5 * Coeff_x * (l_mat[i][j + 1] - 2.0 * l_mat[i][j] + l_mat[i][j - 1]));
				s_mat_n[i][j] = s_mat[i][j] + (Coeff_x * (0.5 * (r_mat[i + 1][j] - r_mat[i - 1][j]) + 0.5 * (1-Coeff_x) * (s_mat[i + 1][j] - 2 * s_mat[i][j] + s_mat[i - 1][j])))\
					+ (Coeff_x * (0.5 * (l_mat[i][j + 1] - l_mat[i][j - 1]) + 0.5 * (1-Coeff_x) * (s_mat[i][j + 1] - 2 * s_mat[i][j] + s_mat[i][j - 1])));
			}
		}

		for (j = 0; j <= Y ; j++) {
			for (i = 0; i <= X ; i++) {
				u_n[i][j] = u[i][j]+0.5* time_step *(s_mat_n[i][j]- s_mat[i][j]);
			

			}

		}
		/*Boundary condition*/
//u[(int)mid_X][(int)mid_Y] = exp(-pow((t* - t0_g) / tp_g, 2));

		// for (i = 0; i <= X; i++) {
		// 	u_n[0][i] = u[1][i] + Q_x * (u_n[1][i] - u[0][i]);
		// 	u_n[X][i] = u[X - 1][i] + Q_x * (u_n[X - 1][i] - u[X][i]);
		// 	u_n[i][0] = u[i][1] + Q_x * (u_n[i][1] - u[i][0]);
		// 	u_n[i][Y] = u[i][Y - 1] + Q_x * (u_n[i][Y - 1] - u[i][Y]);
		// }


		for (j = 0; j <= Y; j++) {
			for (i = 0; i <= X ; i++) {
				u_p[i][j] = u[i][j];
				u[i][j] = u_n[i][j];
				r_mat[i][j] = r_mat_n[i][j];
				l_mat[i][j] = l_mat_n[i][j];
				s_mat[i][j] = s_mat_n[i][j];
				s_mat_p[i][j] = s_mat[i][j];

			}

		}
		
	

		// for (j = 0; j <= Y; j++) {
		// 	for (i = 0; i <= X; i++) {
		// 		fprintf(fpWrite, "%lf ", u[i][j]);
		// 	}
		// 	fprintf(fpWrite, "\n");
		// }


	}



	//for (j = 0; j <= Y; j++) {
	//	for (i = 0; i <= X; i++) {
	//		fprintf(fpWrite,"%lf ", u[i][j]);
	//	}
	//	fprintf(fpWrite,"\n");
	//}



	// fclose(fpWrite);

	/*for (j = 0; j <= Y; j++) {
		for (i = 0; i <= X; i++) {
			printf("%lf ", u[i][j]);
		}
		printf("\n");
	}*/
}