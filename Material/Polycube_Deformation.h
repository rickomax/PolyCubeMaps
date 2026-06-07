#ifndef POLYCUBE_DEFORMATION_H
#define POLYCUBE_DEFORMATION_H

#include "..\..\VolumeMeshViewer\MeshDefinition.h"
#include <Eigen\Dense>
#include <hash_map>

class polycube_deformation_interface
{
public:
	polycube_deformation_interface();
	~polycube_deformation_interface();

	bool get_prepare_ok(){ return prepare_ok; };
	void set_prepare_ok(bool ok){ prepare_ok = ok; };
	void load_deformation_result(const char* filename, VolumeMesh* mesh_);
	void save_deformation_result(const char* filename, VolumeMesh* mesh_);

	void prepare_for_deformation(VolumeMesh* mesh_);
	void compute_source_S(VolumeMesh* mesh_);
	void compute_distortion(VolumeMesh* mesh_);

	//2015.11.25 refine the mesh to polycube
	void exp_mips_deformation_refine_polycube_omp(int max_iter, int iter2, double energy_power, double angle_area_ratio, VolumeMesh * mesh_, bool use_xyz = true, bool use_RGNF = true, bool use_half_sigma_s = false);
	void exp_mips_deformation_refine_polycube_normal_omp(int max_iter, double energy_power, double angle_area_ratio, VolumeMesh * mesh_);
	void set_sigma_r(double r) { sigma_r = r; }
	void set_sigma_s(double s) { sigma_s = s; }
	void set_RGNF_iter_count(int c) { RGNF_iter_count = c; }
	void load_boundary_face_label(const char* filename, VolumeMesh* mesh_);
	void load_up_down_chart(const char* filename);
	void flatten_boundary_face_omp(int max_iter, double energy_power, double angle_area_ratio, VolumeMesh * mesh_);
	void flatten_boundary_face_omp_normal_position(int max_iter, double energy_power, double angle_area_ratio, VolumeMesh * mesh_);
	void find_all_corner(VolumeMesh * mesh_);
	void find_all_chart_value(VolumeMesh * mesh_);
	void save_polycube_constraint(const char* filename, VolumeMesh* mesh_);
	void save_polycube_int(const char* filename, VolumeMesh* mesh_);
	void deform_ARAP_polycube(VolumeMesh* mesh_);
	void boundary_mapping_polycube(VolumeMesh* mesh_); //using IVF
	void adjust_orientation(VolumeMesh* mesh_);
	
	void assign_pos_mesh(VolumeMesh* mesh_, bool r_order = false);
	void set_vertex_color(int n_color, const std::vector<int>& v_color, VolumeMesh * mesh_);

	std::vector<int>& get_change_big_diff(){ return change_big_flag; };
	void set_change_big_diff(std::vector<int>& change_big_flag_){ change_big_flag = change_big_flag_; };

	std::vector<OpenVolumeMesh::CellHandle>& get_flipped_cell(){ return flipped_cell; };

private:
	//2015.11.25 refine polycube
	bool compute_exp_misp_energy_refine_polycube(const int& vc_size, const double& old_e, double& new_e,
		const double* posx, const double* posy, const double* posz,
		const std::vector<std::vector<double> >& vc_S, double* exp_vec,
		const double& npx, const double& npy, const double& npz,
		double alpha, double beta,
		const double& ga, const int& vf_size, 
		const double* vf_px, const double* vf_py, const double* vf_pz, const double* mu );
	bool compute_exp_misp_energy_refine_polycube_normal(const int& vc_size, const double& old_e, double& new_e,
		const double* posx, const double* posy, const double* posz,
		const std::vector<std::vector<double> >& vc_S, double* exp_vec,
		const double& npx, const double& npy, const double& npz,
		double alpha, double beta,
		const double& mu1, const double& mu2, const int& vv_size,
		const double* vv_px, const double* vv_py, const double* vv_pz);
	bool compute_exp_misp_energy_refine_polycube_position(const int& vc_size, const double& old_e, double& new_e,
		const double* posx, const double* posy, const double* posz,
		const std::vector<std::vector<double> >& vc_S, double* exp_vec,
		const double& npx, const double& npy, const double& npz,
		double alpha, double beta,
		const double& ga, const int& vf_size, OpenVolumeMesh::Geometry::Vec3i& v_type,
		const double* vf_px, const double* vf_py, const double* vf_pz);

	bool local_check_negative_volume4(const int& vc_size, const double* posx, const double* posy, const double* posz,
		const double* vc_n_cross_x, const double* vc_n_cross_y, const double* vc_n_cross_z,
		const double& npx, const double& npy, const double& npz);

	//2015.11.25 refine the mesh to polycube
	void exp_mips_deformation_refine_one_polycube(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, double ga = 1.0, bool update_all = true);
	void exp_mips_deformation_refine_one_polycube_normal(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, double mu1, double mu2, bool update_all = true);//use normal
	//void exp_mips_deformation_refine_one_polycube_position(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, double ga = 1.0);

	void check_big_distortion_cell(
		const int& vc_size, std::vector<int>& large_distortion_flag,
		const double* posx, const double* posy, const double* posz,
		const std::vector<std::vector<double> >& vc_S,
		const double& npx, const double& npy, const double& npz);

private:
	bool prepare_ok;
	//basic data
	std::vector< double > cell_volume;
	std::vector< std::vector<int> > cell_vertex;
	std::vector< Eigen::Matrix3d > cell_S;
	std::vector< std::vector<std::vector<int> > > cell_vertex_vertex;
	std::vector< std::vector<int>> vertex_cell;
	std::vector< std::vector<std::vector<int>>> vertex_cell_vertex;
	std::vector< std::vector<std::vector<double> > > vcv_S;

	std::vector<double> dpx; std::vector<double> dpy; std::vector<double> dpz;
	int number_of_color; int max_vc_size;
	std::vector< std::vector<int> > vertex_diff_color;
	std::vector<std::vector<std::vector<double> >> vc_pos_x_omp;
	std::vector<std::vector<std::vector<double> >> vc_pos_y_omp;
	std::vector<std::vector<std::vector<double> >> vc_pos_z_omp;
	std::vector<std::vector<double >> vc_pos_x2_omp;
	std::vector<std::vector<double >> vc_pos_y2_omp;
	std::vector<std::vector<double >> vc_pos_z2_omp;
	std::vector<std::vector<double >> vc_pos_x3_omp;
	std::vector<std::vector<double >> vc_pos_y3_omp;
	std::vector<std::vector<double >> vc_pos_z3_omp;
	std::vector<std::vector<std::vector<double> >> vc_S_omp;
	std::vector<std::vector<double >> vc_n_cross_x_omp;
	std::vector<std::vector<double >> vc_n_cross_y_omp;
	std::vector<std::vector<double >> vc_n_cross_z_omp;
	std::vector<std::vector<double >> exp_vec_omp;
	std::vector<std::vector<double >> gx_vec_omp;
	std::vector<std::vector<double >> gy_vec_omp;
	std::vector<std::vector<double >> gz_vec_omp;
	std::vector<std::vector<double >> mu_vec_omp;
	std::vector<std::vector<int >> large_dis_flag_omp;

	std::vector<int> is_polycube_handles;

	std::vector<OpenVolumeMesh::CellHandle> flipped_cell;

	//2015.11.25 boundary vertex information
	std::vector<OpenVolumeMesh::Geometry::Vec3i> bfv_id;
	std::vector< std::vector<int> > bvf_id;
	std::vector< std::vector<int> > bvv_id;
	std::vector< OpenVolumeMesh::Geometry::Vec3d > target_bfn;
	std::vector<int> bf_chart;
	std::vector<std::vector<int>> polycube_chart; //each chart which stores the faces
	std::vector<int> polycube_chart_label;//0 for x, 1 for y, 2 for z, for each chart
	std::vector<double> chart_mean_value; //the mean value of each chart
	std::vector<int> up_chart_id; std::vector<int> down_chart_id; std::vector<double> diff_up_down;
	std::vector<std::vector<int> > bef_id;
	std::vector<OpenVolumeMesh::Geometry::Vec3i> vertex_type;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> src_pos;

	double avg_boundary_edge_length; int boundary_face_number; 
	std::vector<int> is_bv;
	double sigma_s; double sigma_r; int RGNF_iter_count;
	void filter_boundary_normal_RGNF(VolumeMesh * mesh_,double sigma_s_, double sigma_r_, int iter_count_, bool xyz_flag = false);
	void filter_boundary_normal_laplace(VolumeMesh * mesh_, int iter_count_);
	void filter_boundary_normal_ring(VolumeMesh * mesh_, double sigma_s_, int iter_count_, bool xyz_flag = false);

	bool check_normal_difference(double th = 45);
	void compute_face_normal(std::vector< OpenVolumeMesh::Geometry::Vec3d >& bfn);
	double compute_energy_ratio_filter_deform();
	double compute_energy_ratio_flatten_normal();
	std::vector<double> local_energy_ratio;
	std::vector<int> change_big_flag;
	std::vector< OpenVolumeMesh::Geometry::Vec3d > last_bfn;
	std::vector<int> distoriton_big_count;
	std::vector<int> distoriton_big_cell;
	double bound_K;
};


#endif // !TETRAHEDRAL_DEFORMATION_H