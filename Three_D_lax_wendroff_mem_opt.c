
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<omp.h>
#define _CRT_SECURE_NO_WARNINGS
#define T 500
#define M_PIl 3.141592653589793238462643383279502884L /* pi */
#define grid_size_x 0.01
#define grid_size_y 0.01
#define step_ratio 0.26
#define c0 340
#include<pthread.h>
#include<time.h>
#include <stdint.h>

#define N   3
#define AVE 12*7/N
#define X AVE*N+1
#define Y AVE*N+1
#define Z AVE*N+1
#define VARM(vi,vj,vk) vi+(X+1)*vj+vk*(X+1)*(X+1)

double time_step = step_ratio * grid_size_x / c0;
double Coeff_x ;
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
double* q_mat[N] ;
double* q_mat_n[N] ;

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
	int var_c,var_ip1,var_jp1,var_im1,var_jm1,var_qm1,var_qp1,mid_k,mid_var,mid_var_ip1,mid_var_jp1,mid_var_im1,mid_var_jm1,mid_var_qm1,mid_var_qp1;
	int i, j, k, id,  pth,k_mid;
	double k1,k2,k3,k1_mid,k2_mid,k3_mid;
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


	
	for (i = 0; i < N; i++){
		if(i == 0 || i == N-1){
		u_n[i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		u [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		u_p [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		r_mat [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		s_mat [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		l_mat [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		q_mat [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		s_mat_n [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		r_mat_n [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		l_mat_n [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		q_mat_n [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		s_mat_p [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE + 1));
		}
		else{
		u_n[i] = malloc(sizeof(double) * (X + 1)* (X + 1)* (AVE));
		u [i] = malloc(sizeof(double) * (X + 1)* (X + 1)* (AVE));
		u_p [i] = malloc(sizeof(double) * (X + 1)* (X + 1)* (AVE));
		r_mat [i] = malloc(sizeof(double) * (X + 1)* (X + 1)* (AVE));
		s_mat [i] = malloc(sizeof(double) * (X + 1)* (X + 1)* (AVE));
		l_mat [i] = malloc(sizeof(double) * (X + 1)* (X + 1)* (AVE));
		q_mat [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE));
		s_mat_n [i] = malloc(sizeof(double) * (X + 1)* (X + 1)* (AVE));
		r_mat_n [i] = malloc(sizeof(double) * (X + 1)* (X + 1)* (AVE));
		l_mat_n [i] = malloc(sizeof(double) * (X + 1)* (X + 1)* (AVE));
		q_mat_n [i] = malloc(sizeof(double) * (X + 1) * ( X + 1 ) *(AVE));
		s_mat_p [i] = malloc(sizeof(double) * (X + 1)* (X + 1)* (AVE));
		}
	}








	square_of_Coeff_x = Coeff_x * Coeff_x;

	Q_x = (grid_size_x - c0 * step_ratio) / (grid_size_x + c0 * step_ratio);

	double mid_Y = 0.5 * Y;
	double mid_X = 0.5 * X;
	double mid_Z = 0.5 * Z;

	k2_mid = ((int)mid_Z - (double)AVE - 1);
	k3_mid = (double)AVE;
	k1_mid = k2_mid / k3_mid + 1;
	k_mid = (int)k1_mid;

	
	mid_var=VARM((int)mid_X,(int)mid_Y,(int)mid_Z)-(AVE+1)*(X + 1)* (X + 1)-(k_mid-1)*(AVE)*(X + 1)* (X + 1);
	mid_var_ip1=VARM(((int)mid_X+1),(int)mid_Y,(int)mid_Z)-(AVE+1)* (X + 1)* (X + 1) -(k_mid-1)*(AVE)* (X + 1)* (X + 1);
	mid_var_jp1=VARM((int)mid_X,((int)mid_Y+1),(int)mid_Z)-(AVE+1)* (X + 1)* (X + 1) -(k_mid-1)*(AVE)* (X + 1)* (X + 1);
	mid_var_im1=VARM(((int)mid_X-1),(int)mid_Y,(int)mid_Z)-(AVE+1)* (X + 1)* (X + 1) -(k_mid-1)*(AVE)* (X + 1)* (X + 1);
	mid_var_jm1=VARM((int)mid_X,((int)mid_Y-1),(int)mid_Z)-(AVE+1)* (X + 1)* (X + 1) -(k_mid-1)*(AVE)* (X + 1)* (X + 1);
	mid_var_qp1=VARM((int)mid_X,(int)mid_Y,((int)mid_Z+1))-(AVE+1)*(X + 1)* (X + 1)-(k_mid-1)*(AVE)*(X + 1)* (X + 1);
	mid_var_qm1=VARM((int)mid_X,(int)mid_Y,((int)mid_Z-1))-(AVE+1)*(X + 1)* (X + 1)-(k_mid-1)*(AVE)*(X + 1)* (X + 1);





	tp_g = 15;
	t0_g = 5;
	/*Matrix initialisation*/

	for (k = 0; k <= Z; k++) {
		for (j = 0; j <= Y; j++) {
			for (i = 0; i <= X; i++) {
								//
				if(k<=AVE){
					pth = 0;
				}
				else if (k >= (N-1)*AVE+1){
					pth = N-1;
				}
				else{
					k2 = (k - (double)AVE - 1);
					k3 = (double)AVE;
					k1 = k2 / k3 + 1;
					pth =(int)k1;
				}
				if (pth==0){
					var_c=VARM(i,j,k);
				}
				else{
					var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
				}


				q_mat[pth][var_c]= q_mat_n[pth][var_c]=s_mat_p[pth][var_c] = r_mat_n[pth][var_c] = s_mat_n[pth][var_c] = l_mat_n[pth][var_c] = r_mat[pth][var_c] = s_mat[pth][var_c] = l_mat[pth][var_c] = u_p[pth][var_c] = u[pth][var_c] = u_n[pth][var_c] = 0;
					}
		}
	}


	for (t = 0; t <= T; t++) {
		u[k_mid][mid_var] = P_coeff * exp(-pow((t - t0_g + 1) / tp_g, 2));
		if (t == 0) {
			for (k = 1; k <= Z - 1; k++) {
				for (j = 1; j <= Y - 1; j++) {
					for (i = 1; i <= X - 1; i++) {
						if(k<=AVE){
							pth = 0;
						}
						else if (k >= (N-1)*AVE+1){
							pth = N-1;
						}
						else{
							k2 = (k - (double)AVE - 1);
							k3 = (double)AVE;
							k1 = k2 / k3 + 1;
							pth =(int)k1;
						}
						if (pth==0){
							var_c=VARM(i,j,k);
						}
						else{
							var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
						}

				switch (pth){
				case 0:
					var_c=VARM(i,j,k);
					var_ip1=VARM((i+1),j,k);
					var_im1=VARM((i-1),j,k);
					var_jm1=VARM(i,(j-1),k);
					var_jp1=VARM(i,(j+1),k);
					var_qm1=VARM(i,j,(k-1));
				if (k==AVE){
					var_qp1=VARM(i,j,(k+1))-(AVE+1)* (X + 1)* (X + 1) -((pth+1)-1)*(AVE)* (X + 1)* (X + 1);
					u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1]+u[pth+1][var_qp1]+u[pth][var_qm1] - 6 * u[pth][var_c]);
				}else{
					var_qp1=VARM(i,j,(k+1));
					u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1]+u[pth][var_qp1]+u[pth][var_qm1] - 6 * u[pth][var_c]);
				}
					break;

				case N-1:
					var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					var_ip1=VARM((i+1),j,k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					var_im1=VARM((i-1),j,k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					var_jm1=VARM(i,(j-1),k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					var_jp1=VARM(i,(j+1),k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					var_qp1=VARM(i,j,(k+1))-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					if(k==pth*AVE+1){var_qp1=VARM(i,j,(k+1))-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
						if (pth == 1) {
							var_qm1 = VARM(i,j, (k - 1));
						}
						else {
							var_qm1 = VARM(i,j, (k - 1)) -(AVE+1)*(X+1)* (X + 1)-((pth-1)-1)*(AVE)* (X + 1)* (X + 1);
						}
						u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1]+u[pth][var_qp1]+u[pth-1][var_qm1] - 6 * u[pth][var_c]);
						//printf("%lf",u_n[pth][var_c]);
					}else{
						var_qm1=VARM(i,j,(k-1)) -(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
						u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1]+u[pth][var_qp1]+u[pth][var_qm1] - 6 * u[pth][var_c]);
						//printf("%lf",u_n[pth][var_c]);
					}
					//printf("[i=%d, j=%d] ", i, j);
					break;
				default:
					var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					var_ip1=VARM((i+1),j,k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					var_im1=VARM((i-1),j,k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					var_jm1=VARM(i,(j-1),k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					var_jp1=VARM(i,(j+1),k)-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);

				if (k==AVE*(pth+1)){
					var_qm1=VARM(i,j,(k-1))-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					var_qp1=VARM(i,j,(k+1))-(AVE+1)*(X+1)* (X + 1)-((pth+1)-1)*(AVE)* (X + 1)* (X + 1);
					u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1]+u[pth+1][var_qp1]+u[pth][var_qm1] - 6 * u[pth][var_c]);
					//printf("%lf",u_n[pth][var_c]);
				}else if(k==AVE*pth+1){
					var_qp1=VARM(i,j,(k+1))-(AVE+1)* (X + 1)* (X + 1) -(pth-1)*(AVE)* (X + 1)* (X + 1);
					if (pth == 1) {
						var_qm1 = VARM(i,j, (k - 1));
					}
					else{
						var_qm1 = VARM(i,j, (k - 1)) -(AVE+1)*(X+1)* (X + 1)-((pth-1)-1)*(AVE)* (X + 1)* (X + 1);
					}
					
					u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1]+u[pth][var_qp1]+u[pth-1][var_qm1] - 6 * u[pth][var_c]);
					//printf("%lf",u_n[pth][var_c]);
				}else{
					var_qp1=VARM(i,j,(k+1))-(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					var_qm1=VARM(i,j,(k-1)) -(AVE+1)*(X+1)* (X + 1)-(pth-1)*(AVE)* (X + 1)* (X + 1);
					u_n[pth][var_c] = 2 * u[pth][var_c] - u_p[pth][var_c] + square_of_Coeff_x * (u[pth][var_im1] + u[pth][var_jm1] + u[pth][var_ip1] + u[pth][var_jp1]+u[pth][var_qp1]+u[pth][var_qm1] - 6 * u[pth][var_c]);
				} 
					break;
				}
					}
				}
			}
		}
		

		r_mat[k_mid][mid_var] = c0 * (u[k_mid][mid_var_ip1] - u[k_mid][mid_var_im1]) / (2.0 * grid_size_x);
		l_mat[k_mid][mid_var] = c0 * (u[k_mid][mid_var_jp1] - u[k_mid][mid_var_jm1]) / (2.0 * grid_size_x);
		q_mat[k_mid][mid_var] = c0 * (u[k_mid][mid_var_qp1] - u[k_mid][mid_var_qm1]) / (2.0 * grid_size_x);
		s_mat[k_mid][mid_var] = 2 * (u[k_mid][mid_var] - u_p[k_mid][mid_var]) / time_step + s_mat_p[k_mid][mid_var];

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
	int var_c,var_ip1,var_jp1,var_im1,var_jm1,var_qp1,var_qm1;
	double va, vb, vc, vd, ve,s_mat_c, s_mat_l, s_mat_r;
	switch (n) {
	case 0:
	for (j = 1; j <= Y - 1; j++) {
		for (i = 1; i <= X - 1; i++) {
			k=end;
			var_c=VARM(i,j,k);
			var_ip1=VARM((i+1),j,k);
			var_im1=VARM((i-1),j,k);
			var_jm1=VARM(i,(j-1),k);
			var_jp1=VARM(i,(j+1),k);
			var_qm1=VARM(i,j,(k-1));
			var_qp1=VARM(i,j,(k+1))-(AVE+1)* (X + 1)* (X + 1) -((n+1)-1)*(AVE)* (X + 1)* (X + 1);

			r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * Coeff_x * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
			l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n][var_jm1])) + 0.5 * Coeff_x * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
			q_mat_n[n][var_c] = q_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n+1][var_qp1] - s_mat[n][var_qm1])) + 0.5 * Coeff_x * (q_mat[n+1][var_qp1] - 2.0 * q_mat[n][var_c] + q_mat[n][var_qm1]));
			s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2 * s_mat[n][var_c] + s_mat[n][var_im1])))\
			+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] -  2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])))\
			+ (Coeff_x * (0.5 * (q_mat[n+1][var_qp1] - q_mat[n][var_qm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n+1][var_qp1] - 2.0* s_mat[n][var_c] + s_mat[n][var_qm1])));


		}
	}



	for (k = 1; k <= end-1; k++) {
		for (j = 1; j <= Y - 1; j++) {
			for (i = 1; i <= X - 1; i++) {
			var_c=VARM(i,j,k);
			var_ip1=VARM((i+1),j,k);
			var_im1=VARM((i-1),j,k);
			var_jm1=VARM(i,(j-1),k);
			var_jp1=VARM(i,(j+1),k);
			var_qm1=VARM(i,j,(k-1));
			var_qp1=VARM(i,j,(k+1));

			r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * Coeff_x * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
			l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n][var_jm1])) + 0.5 * Coeff_x * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
			q_mat_n[n][var_c] = q_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_qp1] - s_mat[n][var_qm1])) + 0.5 * Coeff_x * (q_mat[n][var_qp1] - 2.0 * q_mat[n][var_c] + q_mat[n][var_qm1]));
			s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2 * s_mat[n][var_c] + s_mat[n][var_im1])))\
			+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] -  2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])))\
			+ (Coeff_x * (0.5 * (q_mat[n][var_qp1] - q_mat[n][var_qm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_qp1] - 2.0* s_mat[n][var_c] + s_mat[n][var_qm1])));
			}
		}
	}
		for (k = 0; k <= end; k++) {
			for (j = 0; j <= Y; j++) {
				for (i = 0; i <= X; i++) {
					var_c=VARM(i,j,k);
					u_n[n][var_c] = u[n][var_c] + 0.5 * time_step * (s_mat_n[n][var_c] - s_mat[n][var_c]);
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
						var_c=VARM(i,j,k);
						fprintf(fpWrite, "%lf ", u_n[n][var_c]);
					}
					fprintf(fpWrite, "\n");
				}
			}
			fclose(fpWrite);
		}
		else if (t % 5 == 0) {
			FILE* fpWrite = fopen(names[n], "a");
			if (fpWrite == NULL)
			{
				return 0;
			}
			for (k = 0; k <= end; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
						var_c=VARM(i,j,k);
						fprintf(fpWrite, "%lf ", u_n[n][var_c]);
					}
					fprintf(fpWrite, "\n");
				}
			}
			fclose(fpWrite);
		}
		break;

	case N-1:
	for (j = 1; j <= Y - 1; j++) {
		for (i = 1; i <= X - 1; i++) {
			k=start;
			var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_ip1=VARM((i+1),j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_im1=VARM((i-1),j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_jm1=VARM(i,(j-1),k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_jp1=VARM(i,(j+1),k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_qp1=VARM(i,j,(k+1))-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			if (n == 1) {
				var_qm1 = VARM(i,j, (k - 1));
			}
			else {
				var_qm1 = VARM(i,j, (k - 1)) -(AVE+1)*(X+1)* (X + 1)-((n-1)-1)*(AVE)* (X + 1)* (X + 1);
			}
			r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * Coeff_x * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
			l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n][var_jm1])) + 0.5 * Coeff_x * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
			q_mat_n[n][var_c] = q_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_qp1] - s_mat[n-1][var_qm1])) + 0.5 * Coeff_x * (q_mat[n][var_qp1] - 2.0 * q_mat[n][var_c] + q_mat[n-1][var_qm1]));
			s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2 * s_mat[n][var_c] + s_mat[n][var_im1])))\
			+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] -  2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])))\
			+ (Coeff_x * (0.5 * (q_mat[n][var_qp1] - q_mat[n-1][var_qm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_qp1] - 2.0* s_mat[n][var_c] + s_mat[n-1][var_qm1])));
			
		}
	}
		for (k = start+1; k <= Z-1; k++) {
			for (j = 1; j <= Y - 1; j++) {
				for (i = 1; i <= X - 1; i++) {
					var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_ip1=VARM((i+1),j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_im1=VARM((i-1),j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_jm1=VARM(i,(j-1),k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_jp1=VARM(i,(j+1),k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_qp1=VARM(i,j,(k+1))-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_qm1=VARM(i,j,(k-1))-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);

					r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * Coeff_x * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
					l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n][var_jm1])) + 0.5 * Coeff_x * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
					q_mat_n[n][var_c] = q_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_qp1] - s_mat[n][var_qm1])) + 0.5 * Coeff_x * (q_mat[n][var_qp1] - 2.0 * q_mat[n][var_c] + q_mat[n][var_qm1]));
					s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2 * s_mat[n][var_c] + s_mat[n][var_im1])))\
					+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] -  2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])))\
					+ (Coeff_x * (0.5 * (q_mat[n][var_qp1] - q_mat[n][var_qm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_qp1] - 2.0* s_mat[n][var_c] + s_mat[n][var_qm1])));
				}
			}
		}
		for (k = start; k <= Z; k++) {
			for (j = 0; j <= Y; j++) {
				for (i = 0; i <= X; i++) {
					var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					u_n[n][var_c] = u[n][var_c] + 0.5 * time_step * (s_mat_n[n][var_c] - s_mat[n][var_c]);
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
						var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
						fprintf(fpWrite, "%lf ", u_n[n][var_c]);
					}
					fprintf(fpWrite, "\n");
				}
			}
			fclose(fpWrite);
		}
		else if (t % 5 == 0) {
			FILE* fpWrite = fopen(names[n], "a");
			if (fpWrite == NULL)
			{
				return 0;
			}
			for (k = start; k <= Z; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
						var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
						fprintf(fpWrite, "%lf ", u_n[n][var_c]);
					}
					fprintf(fpWrite, "\n");
				}
			}
			fclose(fpWrite);
		}
		break;

	default:
	for (j = 1; j <= Y - 1; j++) {
		for (i = 1; i <= X - 1; i++) {
			k = start;
		
			var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_ip1=VARM((i+1),j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_im1=VARM((i-1),j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_jm1=VARM(i,(j-1),k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_jp1=VARM(i,(j+1),k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_qp1=VARM(i,j,(k+1))-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			if (n == 1) {
				var_qm1 = VARM(i,j, (k - 1));
			}
			else {
				var_qm1 = VARM(i,j, (k - 1)) -(AVE+1)*(X+1)* (X + 1)-((n-1)-1)*(AVE)* (X + 1)* (X + 1);
			}
			r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * Coeff_x * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
			l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n][var_jm1])) + 0.5 * Coeff_x * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
			q_mat_n[n][var_c] = q_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_qp1] - s_mat[n-1][var_qm1])) + 0.5 * Coeff_x * (q_mat[n][var_qp1] - 2.0 * q_mat[n][var_c] + q_mat[n-1][var_qm1]));
			s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2 * s_mat[n][var_c] + s_mat[n][var_im1])))\
			+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] -  2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])))\
			+ (Coeff_x * (0.5 * (q_mat[n][var_qp1] - q_mat[n-1][var_qm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_qp1] - 2.0* s_mat[n][var_c] + s_mat[n-1][var_qm1])));
		}
	}
	for (j = 1; j <= Y - 1; j++) {
		for (i = 1; i <= X - 1; i++) {
			k=end;
			var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_ip1=VARM((i+1),j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_im1=VARM((i-1),j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_jm1=VARM(i,(j-1),k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_jp1=VARM(i,(j+1),k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_qm1=VARM(i,j,(k-1))-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
			var_qp1=VARM(i,j,(k+1))-(AVE+1)* (X + 1)* (X + 1) -((n+1)-1)*(AVE)* (X + 1)* (X + 1);

			r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * Coeff_x * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
			l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n][var_jm1])) + 0.5 * Coeff_x * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
			q_mat_n[n][var_c] = q_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n+1][var_qp1] - s_mat[n][var_qm1])) + 0.5 * Coeff_x * (q_mat[n+1][var_qp1] - 2.0 * q_mat[n][var_c] + q_mat[n][var_qm1]));
			s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2 * s_mat[n][var_c] + s_mat[n][var_im1])))\
			+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] -  2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])))\
			+ (Coeff_x * (0.5 * (q_mat[n+1][var_qp1] - q_mat[n][var_qm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n+1][var_qp1] - 2.0* s_mat[n][var_c] + s_mat[n][var_qm1])));
		}
	}
		for (k = start+1; k <= end-1 ; k++) {
			for (j = 1; j <= Y - 1; j++) {
				for (i = 1; i <= X - 1; i++) {
					var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_ip1=VARM((i+1),j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_im1=VARM((i-1),j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_jm1=VARM(i,(j-1),k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_jp1=VARM(i,(j+1),k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_qp1=VARM(i,j,(k+1))-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					var_qm1=VARM(i,j,(k-1))-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);

					r_mat_n[n][var_c] = r_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_ip1] - s_mat[n][var_im1])) + 0.5 * Coeff_x * (r_mat[n][var_ip1] - 2.0 * r_mat[n][var_c] + r_mat[n][var_im1]));
					l_mat_n[n][var_c] = l_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_jp1] - s_mat[n][var_jm1])) + 0.5 * Coeff_x * (l_mat[n][var_jp1] - 2.0 * l_mat[n][var_c] + l_mat[n][var_jm1]));
					q_mat_n[n][var_c] = q_mat[n][var_c] + (Coeff_x * (0.5 * (s_mat[n][var_qp1] - s_mat[n][var_qm1])) + 0.5 * Coeff_x * (q_mat[n][var_qp1] - 2.0 * q_mat[n][var_c] + q_mat[n][var_qm1]));
					s_mat_n[n][var_c] = s_mat[n][var_c] + (Coeff_x * (0.5 * (r_mat[n][var_ip1] - r_mat[n][var_im1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_ip1] - 2 * s_mat[n][var_c] + s_mat[n][var_im1])))\
					+ (Coeff_x * (0.5 * (l_mat[n][var_jp1] - l_mat[n][var_jm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_jp1] -  2.0 * s_mat[n][var_c] + s_mat[n][var_jm1])))\
					+ (Coeff_x * (0.5 * (q_mat[n][var_qp1] - q_mat[n][var_qm1]) + 0.5 * (1 - Coeff_x) * (s_mat[n][var_qp1] - 2.0* s_mat[n][var_c] + s_mat[n][var_qm1])));}
				
			}
		}
		for (k = start; k <= end; k++) {
			for (j = 0; j <= Y; j++) {
				for (i = 0; i <= X; i++) {
					var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
					u_n[n][var_c] = u[n][var_c] + 0.5 * time_step * (s_mat_n[n][var_c] - s_mat[n][var_c]);	
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
						var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
						fprintf(fpWrite, "%lf ", u_n[n][var_c]);
					}
					fprintf(fpWrite, "\n");
				}
			}
			fclose(fpWrite);
		}
		else if (t % 5 == 0) {
			FILE* fpWrite = fopen(names[n], "a");
			if (fpWrite == NULL)
			{
				return 0;
			}
			for (k = start; k <= end; k++) {
				for (j = 0; j <= Y; j++) {
					for (i = 0; i <= X; i++) {
						var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
						fprintf(fpWrite, "%lf ", u_n[n][var_c]);
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
						var_c=VARM(i,j,k);
						u_p[n][var_c] = u[n][var_c];
						u[n][var_c]= u_n[n][var_c];
						r_mat[n][var_c] = r_mat_n[n][var_c];
						l_mat[n][var_c] = l_mat_n[n][var_c];
						s_mat[n][var_c] = s_mat_n[n][var_c];
						s_mat_p[n][var_c] = s_mat[n][var_c];
						q_mat[n][var_c] = q_mat_n[n][var_c];
					}
				}
			}
		
	 	break;

	 case N-1:
	 for (k = start; k <= Z; k++) {
	 		for (j = 0; j <= Y; j++) {
	 			for (i = 0; i <= X; i++) {
	 				var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
	 				u_p[n][var_c] = u[n][var_c];
	 				u[n][var_c]= u_n[n][var_c];
	 				r_mat[n][var_c] = r_mat_n[n][var_c];
	 				l_mat[n][var_c] = l_mat_n[n][var_c];
	 				s_mat[n][var_c] = s_mat_n[n][var_c];
	 				s_mat_p[n][var_c] = s_mat[n][var_c];
					 q_mat[n][var_c] = q_mat_n[n][var_c];
	 			}
	 		}
	 	}
	 	break;
	 default:
	 for (k = start; k <= end; k++) {
	 		for (j = 0; j <= Y; j++) {
	 			for (i = 0; i <= X; i++) {
	 				var_c=VARM(i,j,k)-(AVE+1)*(X+1)* (X + 1)-(n-1)*(AVE)* (X + 1)* (X + 1);
	 				u_p[n][var_c] = u[n][var_c];
	 				u[n][var_c]= u_n[n][var_c];
	 				r_mat[n][var_c] = r_mat_n[n][var_c];
	 				l_mat[n][var_c] = l_mat_n[n][var_c];
	 				s_mat[n][var_c] = s_mat_n[n][var_c];
	 				s_mat_p[n][var_c] = s_mat[n][var_c];
					q_mat[n][var_c] = q_mat_n[n][var_c];
	 			}
	 		}
	 	}
	 	break;
	 }
	return 0;
}

