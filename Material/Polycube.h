#ifndef MESH_POLYCUBE_PATCH_H
#define MESH_POLYCUBE_PATCH_H

#include <QObject>
#include "..\VolumeMeshViewer\MeshDefinition.h"
#include <Eigen/Dense>
#include "Polycube_Deformation.h"
#include "Labeling.h"

class Polycube_Interface : public QObject
{
	Q_OBJECT
public:
	Polycube_Interface();
	Polycube_Interface(VolumeMesh* mesh);
	~Polycube_Interface();

	void SetMesh(VolumeMesh* mesh)
	{
		mesh_ = mesh;
		if (pc_label) delete pc_label;
		pc_label = new Labeling();
		reset_aLL_state();
	}

	void rotate_deformation();

	void save_boundary_de(const char* filename);
	void save_result(const char* filename);
	
	void load_vertex_color(const char* filename);
	void load_init_mesh(const char* filename);
	void do_polycube_AMIPS(int max_iter, double sigma_s);
	void filter_deform(double sigma_s);
	void flatten_normal();
	void label_face();
	void save_label(const char* filename);

public slots:
	void draw_for_this_situation();

signals:
	void updateGL_Manual_signal();
	void uppdate_vertex_position_signal();
	void save_opengl_screen_PC_signal(QString);
private:

private:
	void reset_state_once();
	void reset_aLL_state();

	void initlize_pc();
	VolumeMesh* mesh_;
	std::vector<Eigen::Matrix3d> R;
	std::vector<Eigen::Matrix3d> PT;
	std::vector<OpenVolumeMesh::Geometry::Vec4i> cv_id;
	std::vector<double> cell_volume;
	std::vector<OpenVolumeMesh::Geometry::Vec3i> bfv_id;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> bf_n;
	std::vector<int> fix_cell_boundary;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> src_pos;

	void compute_rotation();
	void compute_all_bf_normal();
	void compute_all_quaterion();

	void ARAP_deform();

	//polycube construction
	polycube_deformation_interface pd_normal; //for filtering-deformation
	polycube_deformation_interface pd_flatten; //for flattening boundary
	Labeling* pc_label;
	bool label_OK;
};

#endif