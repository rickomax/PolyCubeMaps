#ifndef ADJUST_ORIENTATION_LBFGS_H
#define ADJUST_ORIENTATION_LBFGS_H

#include "..\VolumeMeshViewer\MeshDefinition.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

struct user_pointer_ad
{
	std::vector<OpenVolumeMesh::Geometry::Vec3d> bfn;
};

#include "BasicMath/MathSuite/LBFGS_Old_Web/HLBFGS.h"
void evalfunc_ad_LBFGS(int N, double* x, double *prev_x, double* f, double* g, void* user_pointer);
void newiteration_ad_LBFGS(int iter, int call_iter, double *x, double* f, double *g, double* gnorm, void* user_pointer);
bool ad_LBFGS(std::vector<OpenVolumeMesh::Geometry::Vec3d>& bfn, std::vector<double>& X);

#endif

