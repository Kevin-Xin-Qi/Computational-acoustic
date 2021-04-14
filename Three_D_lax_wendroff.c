
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>
#include"Lax3D_constants.h"


double time_step = step_ratio * grid_size_x / c0;




void main()
{

	double Coeff_x = c0 * time_step / grid_size_x;

	//float  viscosity;

	double   square_of_c,  time, t0_g, tp_g, square_of_Coeff_x, Q_x;
	int i, j, k, t;


	double*** u_n;
	u_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u_n[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			u_n[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	double*** u;
	u = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			u[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	double*** u_p;
	u_p = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		u_p[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			u_p[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	double*** r_mat;
	r_mat = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		r_mat[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			r_mat[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	double*** r_mat_n;
	r_mat_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		r_mat_n[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			r_mat_n[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}

	double*** s_mat;
	s_mat = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			s_mat[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	double*** s_mat_n;
	s_mat_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat_n[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			s_mat_n[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	double*** s_mat_p;
	s_mat_p = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		s_mat_p[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			s_mat_p[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	double*** l_mat;
	l_mat = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		l_mat[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			l_mat[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	double*** l_mat_n;
	l_mat_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		l_mat_n[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			l_mat_n[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	double*** q_mat;
	q_mat = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		q_mat[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			q_mat[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}
	double*** q_mat_n;
	q_mat_n = malloc(sizeof(double*) * (X + 1));
	for (i = 0; i < X + 1; i++) {
		q_mat_n[i] = malloc(sizeof(double) * (Y + 1));
		for (j = 0; j < Y + 1; j++) {
			q_mat_n[i][j] = malloc(sizeof(double) * (Z + 1));
		}
	}


	//long double** u = malloc(sizeof(long double*) * (X + 1));
	//for (i = 0; i < X + 1; i++) {
	//	u[i] = malloc(sizeof(long double) * (Y + 1));
	//}
	//long double** u_p = malloc(sizeof(long double*) * (X + 1));
	//for (i = 0; i < X + 1; i++) {
	//	u_p[i] = malloc(sizeof(long double) * (Y + 1));
	//}
	//long double** r_mat = malloc(sizeof(long double*) * (X + 1));
	//for (i = 0; i < X + 1; i++) {
	//	r_mat[i] = malloc(sizeof(long double) * (Y + 1));
	//}
	//long double** s_mat = malloc(sizeof(long double*) * (X + 1));
	//for (i = 0; i < X + 1; i++) {
	//	s_mat[i] = malloc(sizeof(long double) * (Y + 1));
	//}
	//long double** l_mat = malloc(sizeof(long double*) * (X + 1));
	//for (i = 0; i < X + 1; i++) {
	//	l_mat[i] = malloc(sizeof(long double) * (Y + 1));
	//}
	//long double** s_mat_n = malloc(sizeof(long double*) * (X + 1));
	//for (i = 0; i < X + 1; i++) {
	//	s_mat_n[i] = malloc(sizeof(long double) * (Y + 1));
	//}
	//long double** r_mat_n = malloc(sizeof(long double*) * (X + 1));
	//for (i = 0; i < X + 1; i++) {
	//	r_mat_n[i] = malloc(sizeof(long double) * (Y + 1));
	//}
	//long double** l_mat_n = malloc(sizeof(long double*) * (X + 1));
	//for (i = 0; i < X + 1; i++) {
	//	l_mat_n[i] = malloc(sizeof(long double) * (Y + 1));
	//}
	//long double** s_mat_p = malloc(sizeof(long double*) * (X + 1));
	//for (i = 0; i < X + 1; i++) {
	//	s_mat_p[i] = malloc(sizeof(long double) * (Y + 1));
	//}





	square_of_Coeff_x = Coeff_x * Coeff_x;

	Q_x = (grid_size_x - c0 * step_ratio) / (grid_size_x + c0 * step_ratio);
	double mid_Y = round(Y / 2.0);

	double mid_X = round(X / 2.0);
	double mid_Z = round(Z / 2.0);





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
	FILE* fpWrite = fopen("wave_data_lax_wendroff_3d.txt", "w");
	if (fpWrite == NULL)
	{
		return 0;
	}
	/*Initial condition*/
	for (t = 0; t <= T; t++) {
		u[(int)mid_X][(int)mid_Y][(int)mid_Z] = exp(-pow((t - t0_g + 1) / tp_g, 2));
		if (t == 1) {
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





#pragma omp parallel 
		{
			for (k = 1; k <= Z - 1; k++) {
				for (j = 1; j <= Y - 1; j++) {
					for (i = 1; i <= X - 1; i++) {

						r_mat_n[i][j][k] = r_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i + 1][j][k] - s_mat[i - 1][j][k])) + 0.5 * Coeff_x * (r_mat[i + 1][j][k] - 2.0 * r_mat[i][j][k] + r_mat[i - 1][j][k]));
						l_mat_n[i][j][k] = l_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i][j + 1][k] - s_mat[i][j - 1][k])) + 0.5 * Coeff_x * (l_mat[i][j + 1][k] - 2.0 * l_mat[i][j][k] + l_mat[i][j - 1][k]));
						q_mat_n[i][j][k] = q_mat[i][j][k] + (Coeff_x * (0.5 * (s_mat[i][j][k + 1] - s_mat[i][j][k - 1])) + 0.5 * Coeff_x * (q_mat[i][j][k + 1] - 2.0 * q_mat[i][j][k] + q_mat[i][j][k - 1]));
					}
				}
			}









			for (k = 1; k <= Z - 1; k++) {
				for (j = 1; j <= Y - 1; j++) {
					for (i = 1; i <= X - 1; i++) {


						s_mat_n[i][j][k] = s_mat[i][j][k] + (Coeff_x * (0.5 * (r_mat_n[i + 1][j][k] - r_mat_n[i - 1][j][k]) + 0.5 * (1 - Coeff_x) * (s_mat[i + 1][j][k] - 2 * s_mat[i][j][k] + s_mat[i - 1][j][k])))\
							+ (Coeff_x * (0.5 * (l_mat_n[i][j + 1][k] - l_mat_n[i][j - 1][k]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j + 1][k] - 2 * s_mat[i][j][k] + s_mat[i][j - 1][k])))\
							+ (Coeff_x * (0.5 * (q_mat_n[i][j][k + 1] - q_mat_n[i][j][k - 1]) + 0.5 * (1 - Coeff_x) * (s_mat[i][j][k + 1] - 2 * s_mat[i][j][k] + s_mat[i][j][k - 1])));

					}
				}
			}
			for (k = 0; k <= Z; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
						u_n[i][j][k] = u[i][j][k] + 0.5 * time_step * (s_mat_n[i][j][k] - s_mat[i][j][k]);


					}

				}
			}
			/*Boundary condition*/
	//u[(int)mid_X][(int)mid_Y] = exp(-pow((t* - t0_g) / tp_g, 2));

			//for (i = 0; i <= X; i++) {
			//	u_n[0][i] = u[1][i] + Q_x * (u_n[1][i] - u[0][i]);
			//	u_n[X][i] = u[X - 1][i] + Q_x * (u_n[X - 1][i] - u[X][i]);
			//	u_n[i][0] = u[i][1] + Q_x * (u_n[i][1] - u[i][0]);
			//	u_n[i][Y] = u[i][Y - 1] + Q_x * (u_n[i][Y - 1] - u[i][Y]);
			//}

			for (k = 0; k <= Z; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
						u_p[i][j][k] = u[i][j][k];
						u[i][j][k] = u_n[i][j][k];
						r_mat[i][j][k] = r_mat_n[i][j][k];
						l_mat[i][j][k] = l_mat_n[i][j][k];
						s_mat[i][j][k] = s_mat_n[i][j][k];
						s_mat_p[i][j][k] = s_mat[i][j][k];

					}

				}
			}
		}

		for (k = 0; k <= Z; k++) {
			for (j = 0; j <= Y; j++) {
				for (i = 0; i <= X; i++) {
					fprintf(fpWrite, "%lf ", u[i][j][k]);
				}
				fprintf(fpWrite, "\n");
			}
	

		}
	}



	//for (j = 0; j <= Y; j++) {
	//	for (i = 0; i <= X; i++) {
	//		fprintf(fpWrite,"%lf ", u[i][j]);
	//	}
	//	fprintf(fpWrite,"\n");
	//}



	fclose(fpWrite);

	/*for (j = 0; j <= Y; j++) {
		for (i = 0; i <= X; i++) {
			printf("%lf ", u[i][j]);
		}
		printf("\n");
	}*/
}

