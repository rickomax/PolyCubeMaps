#include "adjust_orientation_LBFGS.h"

void evalfunc_ad_LBFGS(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer)
{
	double alpha = x[0]; double beta = x[1]; double gamma = x[2];
	double cos_alpha = std::cos(alpha); double sin_alpha = std::sin(alpha);
	double cos_beta = std::cos(beta); double sin_beta = std::sin(beta);
	double cos_gamma = std::cos(gamma); double sin_gamma = std::sin(gamma);
	double r00 = cos_alpha*cos_gamma - cos_beta*sin_alpha*sin_gamma;
	double r01 = sin_alpha*cos_gamma + cos_beta*cos_alpha*sin_gamma;
	double r02 = sin_beta*sin_gamma;

	double r10 = -cos_alpha*sin_gamma - cos_beta*sin_alpha*cos_gamma;
	double r11 = -sin_alpha*sin_gamma + cos_beta*cos_alpha*cos_gamma;
	double r12 = sin_beta*cos_gamma;

	double r20 = sin_beta*sin_alpha;
	double r21 = -sin_beta*cos_alpha;
	double r22 = cos_beta;

	user_pointer_ad* up = (user_pointer_ad*)(user_pointer);
	int nbf = up->bfn.size(); *f = 0.0;
	g[0] = 0; g[1] = 0; g[2] = 0;
	for (int i = 0; i < nbf; ++i)
	{
		OpenVolumeMesh::Geometry::Vec3d& n = up->bfn[i];
		double nx = r00 * n[0] + r01 * n[1] + r02 * n[2];
		double ny = r10 * n[0] + r11 * n[1] + r12 * n[2];
		double nz = r20 * n[0] + r21 * n[1] + r22 * n[2];

		*f += (nx*ny)*(nx*ny) + (ny*nz)*(ny*nz) + (nz*nx)*(nz*nx);

		double d_nx_alpha = (-sin_alpha*cos_gamma - cos_beta*cos_alpha*sin_gamma) * n[0] + (cos_alpha*cos_gamma - cos_beta*sin_alpha*sin_gamma) * n[1] + 0 * n[2];
		double d_ny_alpha = (sin_alpha*sin_gamma - cos_beta*cos_alpha*cos_gamma) * n[0] + (-cos_alpha*sin_gamma - cos_beta*sin_alpha*cos_gamma) * n[1] + 0 * n[2];
		double d_nz_alpha = sin_beta*cos_alpha * n[0] + sin_beta*sin_alpha * n[1] + 0 * n[2];

		double d_nx_beta = (sin_beta*sin_alpha*sin_gamma) * n[0] + (-sin_beta*cos_alpha*sin_gamma) * n[1] + cos_beta*sin_gamma * n[2];
		double d_ny_beta = (sin_beta*sin_alpha*cos_gamma) * n[0] + (-sin_beta*cos_alpha*cos_gamma) * n[1] + cos_beta*cos_gamma * n[2];
		double d_nz_beta = cos_beta*sin_alpha *n[0] - cos_beta*cos_alpha * n[1] - sin_beta * n[2];

		double d_nx_gamma = (-cos_alpha*sin_gamma - cos_beta*sin_alpha*cos_gamma) * n[0] + (-sin_alpha*sin_gamma + cos_beta*cos_alpha*cos_gamma) * n[1] + sin_beta*cos_gamma * n[2];
		double d_ny_gamma = (-cos_alpha*cos_gamma + cos_beta*sin_alpha*sin_gamma) * n[0] + (-sin_alpha*cos_gamma - cos_beta*cos_alpha*sin_gamma) * n[1] + (-sin_beta*sin_gamma) * n[2];
		double d_nz_gamma = 0.0;

		g[0] += 2 * (nx*ny)*(d_nx_alpha*ny + nx*d_ny_alpha) + 2 * (ny*nz)*(d_ny_alpha*nz + ny*d_nz_alpha) + 2 * (nz*nx)*(d_nz_alpha*nx + nz*d_nx_alpha);
		g[1] += 2 * (nx*ny)*(d_nx_beta*ny + nx*d_ny_beta) + 2 * (ny*nz)*(d_ny_beta*nz + ny*d_nz_beta) + 2 * (nz*nx)*(d_nz_beta*nx + nz*d_nx_beta);
		g[2] += 2 * (nx*ny)*(d_nx_gamma*ny + nx*d_ny_gamma) + 2 * (ny*nz)*(d_ny_gamma*nz + ny*d_nz_gamma) + 2 * (nz*nx)*(d_nz_gamma*nx + nz*d_nx_gamma);
	}
}

void newiteration_ad_LBFGS(int iter, int call_iter, double *x, double* f, double *g, double* gnorm, void* user_pointer)
{
	printf("%d %d; %4.3e, %4.3e\n", iter, call_iter, f[0], gnorm[0]);
}

bool ad_LBFGS(std::vector<OpenVolumeMesh::Geometry::Vec3d>& bfn, std::vector<double>& X)
{
	user_pointer_ad up;
	up.bfn = bfn;

	//LBFGS
	double parameter[20];
	int info[20];
	//initialize
	INIT_HLBFGS(parameter, info);
	parameter[2] = 0.9;
	info[4] = 100; //number of iteration
	info[6] = 0;
	info[7] = 0; //if with hessian 1, without 0
	info[10] = 0;

	int N = 3; int M = 7;

	printf("-------------------------------\n");
	printf("start LBFGS\n");

	HLBFGS(N, M, &X[0], evalfunc_ad_LBFGS, 0, HLBFGS_UPDATE_Hessian, newiteration_ad_LBFGS, parameter, info, &up);

	return true;
}
