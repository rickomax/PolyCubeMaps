#include "..\VolumeMeshViewer\OpenGLHeaders.h"
#include "Polycube.h"
#include "..\..\BasicMath\MathSuite\SPARSE\Sparse_Matrix.h"
#include "..\..\BasicMath\MathSuite\SPARSE\Sparse_Solver.h"

Polycube_Interface::Polycube_Interface()
{
	mesh_ = NULL; pc_label = NULL;
	reset_state_once();
	reset_aLL_state();
}

Polycube_Interface::Polycube_Interface(VolumeMesh* mesh)
	: mesh_(mesh)
{
	reset_state_once();
	reset_aLL_state();
}

Polycube_Interface::~Polycube_Interface()
{
	if (pc_label) delete pc_label;
}

void Polycube_Interface::reset_state_once()
{

}

void Polycube_Interface::reset_aLL_state()
{
	R.clear(); PT.clear();
	cell_volume.clear(); fix_cell_boundary.clear();
	cv_id.clear(); bfv_id.clear(); bf_n.clear();
	initlize_pc(); label_OK = false;
}

void Polycube_Interface::draw_for_this_situation()
{
	//if (fix_cell_boundary.size() == 0) return;

	/*glColor3f(1.0, 0.5,0.5);
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < fix_cell_boundary.size(); ++i)
	{
		if (fix_cell_boundary[i] == 1)
		{
			for (int j = 0; j < 3; ++j)
			{
				OpenVolumeMesh::VertexHandle vh = OpenVolumeMesh::VertexHandle(bfv_id[i][j]);
				glVertex3dv(mesh_->vertex(vh).data());
			}
		}
	}
	glEnd();*/

	/*std::vector<int>& change_bd = pd_flatten.get_change_big_diff();
	glColor3f(1.0, 0.5, 0.5);
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < change_bd.size(); ++i)
	{
		if (change_bd[i] == 1)
		{
			for (int j = 0; j < 3; ++j)
			{
				OpenVolumeMesh::VertexHandle vh = OpenVolumeMesh::VertexHandle(bfv_id[i][j]);
				glVertex3dv(mesh_->vertex(vh).data());
			}
		}
	}
	glEnd();*/

	if (label_OK)
	{
		/*GLfloat mat_a_x[] = { 0.921, 0.1, 0.113, 1.0 };
		GLfloat mat_d_x[] = { 0.921, 0.1, 0.113, 1.0 };
		GLfloat mat_s_x[] = { 0.921, 0.1, 0.113, 1.0 };
		GLfloat mat_a_y[] = { 0.188, 0.667, 0.24, 1.0 };
		GLfloat mat_d_y[] = { 0.188, 0.667, 0.24, 1.0 };
		GLfloat mat_s_y[] = { 0.188, 0.667, 0.24, 1.0 };
		GLfloat mat_a_z[] = { 0, 0.35, 0.67, 1.0 };
		GLfloat mat_d_z[] = { 0, 0.35, 0.67, 1.0 };
		GLfloat mat_s_z[] = { 0, 0.35, 0.67, 1.0 };*/

		GLfloat mat_a_x[] = { 0.2, 0.2, 1.0, 1.0 };
		GLfloat mat_d_x[] = { 0.2, 0.2, 1.0, 1.0 };
		GLfloat mat_s_x[] = { 0.2, 0.2, 1.0, 1.0 };
		GLfloat mat_a_y[] = { 0.7, 0.7, 0.7, 1.0 };
		GLfloat mat_d_y[] = { 0.7, 0.7, 0.7, 1.0 };
		GLfloat mat_s_y[] = { 0.7, 0.7, 0.7, 1.0 };
		GLfloat mat_a_z[] = { 1.0, 0.2, 0.2, 1.0 };
		GLfloat mat_d_z[] = { 1.0, 0.2, 0.2, 1.0 };
		GLfloat mat_s_z[] = { 1.0, 0.2, 0.2, 1.0 };
		GLfloat shine[] = { 120.0f };
		GLfloat mat_a_flipped[] = { 0.95, 0.95, 0.0 };
		GLfloat mat_d_flipped[] = { 0.95, 0.95, 0.0 };
		GLfloat mat_s_flipped[] = { 0.95, 0.95, 0.0 };

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.5f, 2.0f);
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);

		/*glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_a_flipped);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_d_flipped);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_s_flipped);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
		std::vector<OpenVolumeMesh::Geometry::Vec3d> hfv(3); std::vector<OpenVolumeMesh::Geometry::Vec3d> shfv(3);
		std::vector<OpenVolumeMesh::CellHandle>& flipped_cell = pd_normal.get_flipped_cell();
		for (int i = 0; i < flipped_cell.size(); ++i)
		{
			OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(flipped_cell[i]);
			std::vector<OpenVolumeMesh::HalfFaceHandle> hfh_Vec = cell.halffaces();

			for (unsigned int j = 0; j < hfh_Vec.size(); ++j)
			{
				OpenVolumeMesh::HalfFaceVertexIter fhv_it = mesh_->hfv_iter(hfh_Vec[j]);
				hfv[0] = mesh_->vertex(*fhv_it); shfv[0] = src_pos[fhv_it->idx()];
				++fhv_it; hfv[1] = mesh_->vertex(*fhv_it); shfv[1] = src_pos[fhv_it->idx()];
				++fhv_it; hfv[2] = mesh_->vertex(*fhv_it); shfv[2] = src_pos[fhv_it->idx()];
				OpenVolumeMesh::Geometry::Vec3d hfn = OpenVolumeMesh::Geometry::cross(shfv[2] - shfv[0], shfv[1] - shfv[0]).normalize();
				glBegin(GL_TRIANGLES);
				glNormal3dv(hfn.data());
				glVertex3dv(hfv[2].data());
				glVertex3dv(hfv[1].data());
				glVertex3dv(hfv[0].data());
				glEnd();
			}
		}*/

		std::vector<int>& tet_label = pc_label->tetLabel; //+X,-X,+Y,-Y,+Z,-Z
		for (int i = 0; i < tet_label.size(); ++i)
		{
			if (tet_label[i] < 0) continue;

			OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
			/*OpenVolumeMesh::Geometry::Vec3d p0 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[0]));
			OpenVolumeMesh::Geometry::Vec3d p1 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[1]));
			OpenVolumeMesh::Geometry::Vec3d p2 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[2]));*/
			OpenVolumeMesh::Geometry::Vec3d p0 = src_pos[one_bfv_id[0]];
			OpenVolumeMesh::Geometry::Vec3d p1 = src_pos[one_bfv_id[1]];
			OpenVolumeMesh::Geometry::Vec3d p2 = src_pos[one_bfv_id[2]];
			OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0).normalize();

			p0 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[0]));
			p1 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[1]));
			p2 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[2]));

			int tl = tet_label[i];
			if (tl == 0 || tl == 3)//X
			{
				//glColor4f(0.921, 0.1, 0.113, 1.0);
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_a_x);
				glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_d_x);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_s_x);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
			}
			else if (tl == 1 || tl == 4)//Y
			{
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_a_y);
				glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_d_y);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_s_y);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
				//glColor4f(0.188, 0.667, 0.24, 1.0);
			}
			else if (tl == 2 || tl == 5)//Z
			{
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_a_z);
				glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_d_z);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_s_z);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
				//glColor4f(0, 0.35, 0.67, 1.0);
			}
			glBegin(GL_TRIANGLES);
			glNormal3dv(n.data());
			glVertex3dv(p0.data());
			glVertex3dv(p1.data());
			glVertex3dv(p2.data());
			glEnd();
		}

		glDisable(GL_POLYGON_OFFSET_FILL);
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		glLineWidth(2.0); glColor3f(0.0, 0.0, 0.0);
		std::vector<int> chart_edge = pc_label->EdgeBoundaryTag;
		for (int i = 0; i < chart_edge.size(); ++i)
		{
			glBegin(GL_LINES);
			if (chart_edge[i] == 1)
			{
				OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(OpenVolumeMesh::EdgeHandle(i));
				glVertex3dv(mesh_->vertex(edge.from_vertex()).data());
				glVertex3dv(mesh_->vertex(edge.to_vertex()).data());
			}
			glEnd();
		}
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
}

void Polycube_Interface::initlize_pc()
{
	//initialize the affine matrix
	if (!mesh_) return;
	int nc = mesh_->n_cells();
	R.resize(nc); PT.resize(nc); cell_volume.resize(nc); cv_id.resize(nc);
	OpenVolumeMesh::Geometry::Vec4i one_cv_id; Eigen::Matrix3d p;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> cv_p(4);
	for (OpenVolumeMesh::CellIter c_it = mesh_->cells_begin(); c_it != mesh_->cells_end(); ++c_it)
	{
		int c_id = c_it->idx(); int cv_count = 0;
		for (OpenVolumeMesh::CellVertexIter cv_it = mesh_->cv_iter(c_it.cur_handle()); cv_it; ++cv_it)
		{
			one_cv_id[cv_count] = cv_it->idx();
			cv_p[cv_count] = mesh_->vertex(*cv_it);
			++cv_count;
		}
		cv_id[c_id] = one_cv_id;
		p(0, 0) = cv_p[1][0] - cv_p[0][0]; p(1, 0) = cv_p[1][1] - cv_p[0][1]; p(2, 0) = cv_p[1][2] - cv_p[0][2];
		p(0, 1) = cv_p[2][0] - cv_p[0][0]; p(1, 1) = cv_p[2][1] - cv_p[0][1]; p(2, 1) = cv_p[2][2] - cv_p[0][2];
		p(0, 2) = cv_p[3][0] - cv_p[0][0]; p(1, 2) = cv_p[3][1] - cv_p[0][1]; p(2, 2) = cv_p[3][2] - cv_p[0][2];
		cell_volume[c_id] = std::abs(p.determinant());
		PT[c_id] = p.inverse();
	}
	bfv_id.clear(); bfv_id.resize(nc, OpenVolumeMesh::Geometry::Vec3i(-1, -1, -1));
	for (OpenVolumeMesh::FaceIter f_it = mesh_->faces_begin(); f_it != mesh_->faces_end(); ++f_it)
	{
		OpenVolumeMesh::HalfFaceHandle hfh0 = mesh_->halfface_handle(*f_it, 0);
		OpenVolumeMesh::HalfFaceHandle hfh1 = mesh_->halfface_handle(*f_it, 1);
		OpenVolumeMesh::CellHandle ch0 = mesh_->incident_cell(hfh0);
		OpenVolumeMesh::CellHandle ch1 = mesh_->incident_cell(hfh1);
		if (ch0 == VolumeMesh::InvalidCellHandle && ch1 != VolumeMesh::InvalidCellHandle)
		{
			int hfv_count = 0; int c_id = ch1.idx();
			for (OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh0); hfv_it; ++hfv_it)
			{
				int hfv_id = hfv_it->idx();
				bfv_id[c_id][hfv_count] = hfv_id;
				++hfv_count;
			}
		}
		else if (ch1 == VolumeMesh::InvalidCellHandle && ch0 != VolumeMesh::InvalidCellHandle)
		{
			int hfv_count = 0; int c_id = ch0.idx();
			for (OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh1); hfv_it; ++hfv_it)
			{
				int hfv_id = hfv_it->idx();
				bfv_id[c_id][hfv_count] = hfv_id;
				++hfv_count;
			}
		}
		else if (ch1 == VolumeMesh::InvalidCellHandle && ch0 == VolumeMesh::InvalidCellHandle)
		{
			printf("Error : Both two halffaces have no cells!\n");
		}
	}

	src_pos.resize(mesh_->n_vertices());
	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		src_pos[v_it->idx()] = mesh_->vertex(*v_it);
	}
}

void Polycube_Interface::rotate_deformation()
{
	initlize_pc();
	for (int i = 0; i < 3; ++i)
	{
		//
		compute_rotation();
		//
		ARAP_deform();

		emit uppdate_vertex_position_signal();
		emit updateGL_Manual_signal();
	}
}

void Polycube_Interface::compute_rotation()
{
	compute_all_bf_normal();
	compute_all_quaterion();
}

void Polycube_Interface::compute_all_bf_normal()
{
	int nc = mesh_->n_cells();
	if (bf_n.size() != nc) bf_n.resize(nc);
	if (fix_cell_boundary.size() != nc) fix_cell_boundary.resize(nc);
	
	int nv = mesh_->n_vertices();
	std::vector< std::vector<int> > vc_id(nv);
	std::vector< OpenVolumeMesh::Geometry::Vec3d > vn(nv, OpenVolumeMesh::Geometry::Vec3d(0,0,0));
	//nc = bfv_id.size()
	for (int i = 0; i < bfv_id.size(); ++i) //i == cell id
	{
		if (bfv_id[i][0] == -1)
		{
			fix_cell_boundary[i] = 0;
			continue;
		}
		
		fix_cell_boundary[i] = 1;

		OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
		OpenVolumeMesh::Geometry::Vec3d p0 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[0]));
		OpenVolumeMesh::Geometry::Vec3d p1 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[1]));
		OpenVolumeMesh::Geometry::Vec3d p2 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[2]));

		bf_n[i] = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0); //with area weight
		
		vc_id[one_bfv_id[0]].push_back(i); vn[one_bfv_id[0]] += bf_n[i];
		vc_id[one_bfv_id[1]].push_back(i); vn[one_bfv_id[1]] += bf_n[i];
		vc_id[one_bfv_id[2]].push_back(i); vn[one_bfv_id[2]] += bf_n[i];

		bf_n[i].normalize();
	}

	//sharp vertex
	for (int i = 0; i < nv; ++i) //i == vertex id
	{
		if (vc_id[i].size() != 0)
		{
			double ni_len = vn[i].norm();
			if (ni_len < 1e-10)
			{
				for (int j = 0; j < vc_id[i].size(); ++j)
				{
					fix_cell_boundary[vc_id[i][j]] = 0;
				}
			}
			else
			{
				vn[i] /= ni_len; //normalize
				double min_cos = 2; //maximal angle
				double max_cos = 2; //minimal angle
				for (int j = 0; j < vc_id[i].size(); ++j)
				{
					double tmp_cos = OpenVolumeMesh::Geometry::dot( vn[i], bf_n[ vc_id[i][j] ] );
					if (tmp_cos < min_cos)
					{
						min_cos = tmp_cos;
					}
					if (tmp_cos > max_cos)
					{
						max_cos = tmp_cos;
					}
				}
				if (min_cos < 0.95) //sharp vertex
				{
					for (int j = 0; j < vc_id[i].size(); ++j)
					{
						fix_cell_boundary[vc_id[i][j]] = 0;
					}
				}
			}
		}
	}
}

void Polycube_Interface::compute_all_quaterion()
{
	//find closest axis
	int nc = mesh_->n_cells(); 
	std::vector<OpenVolumeMesh::Geometry::Vec4d> all_quaterion(nc);
	std::vector<int> map_index(nc, -1); int var_count = 0;
	for (int i = 0; i < fix_cell_boundary.size(); ++i)//for all cells
	{
		if (fix_cell_boundary[i] == 1)
		{
			double max_component = std::abs(bf_n[i][0]); int max_id = 0;
			OpenVolumeMesh::Geometry::Vec3d axis(0.0, 0.0, 0.0);
			for (int j = 1; j<3; j++)
			{
				if (std::abs(bf_n[i][j])>max_component)
				{
					max_component = std::abs(bf_n[i][j]);
					max_id = j;
				}
			}
			axis[max_id] = bf_n[i][max_id] / max_component;
			OpenVolumeMesh::Geometry::Vec3d rot_n = OpenVolumeMesh::Geometry::cross(bf_n[i], axis);
			double sin_val = rot_n.norm();
			if (sin_val < 1e-10)
			{//angle = 0
				all_quaterion[i] = OpenVolumeMesh::Geometry::Vec4d(1,0,0,0);
			}
			else
			{
				rot_n /= sin_val;
				double cos_val = OpenVolumeMesh::Geometry::dot(bf_n[i], axis);
				double rot_angle = std::atan2(sin_val, cos_val);
				cos_val = std::cos(rot_angle *0.5); sin_val = std::sin(rot_angle *0.5);
				all_quaterion[i] = OpenVolumeMesh::Geometry::Vec4d(cos_val, sin_val*rot_n[0], sin_val*rot_n[1], sin_val*rot_n[2]);
			}
		}
		else
		{
			map_index[i] = var_count;
			++var_count;
		}
	}

	//smoothly propagate
	Sparse_Matrix A(var_count, var_count, NOSYM, CCS, 4);
	for (OpenVolumeMesh::CellIter c_it = mesh_->cells_begin(); c_it != mesh_->cells_end(); ++c_it)
	{
		int c_id = c_it->idx(); int var_id = map_index[c_id];
		if (var_id != -1)
		{
			OpenVolumeMesh::OpenVolumeMeshCell cell = mesh_->cell(*c_it);
			std::vector<OpenVolumeMesh::HalfFaceHandle>& hfh_vec = cell.get_halffaces();
			std::vector<int> cell_j; double cc_count = 0.0;
			for (int i = 0; i < hfh_vec.size(); ++i)
			{
				OpenVolumeMesh::HalfFaceHandle opp_hf = mesh_->opposite_halfface_handle(hfh_vec[i]);
				OpenVolumeMesh::CellHandle ch = mesh_->incident_cell(opp_hf);
				if (ch != VolumeMesh::InvalidCellHandle)
				{
					cc_count += 1.0;
					cell_j.push_back(ch.idx());
				}
			}
			double w = 1.0 / cc_count;
			A.fill_entry(var_id, var_id, 1.0);
			for (int i = 0; i < cell_j.size(); ++i)
			{
				int var_j = map_index[cell_j[i]];
				if (var_j == -1)
				{
					A.fill_rhs_entry(var_id, 0, w*all_quaterion[cell_j[i]][0]);
					A.fill_rhs_entry(var_id, 1, w*all_quaterion[cell_j[i]][1]);
					A.fill_rhs_entry(var_id, 2, w*all_quaterion[cell_j[i]][2]);
					A.fill_rhs_entry(var_id, 3, w*all_quaterion[cell_j[i]][3]);
				}
				else
				{
					A.fill_entry(var_id, var_j, -w);
				}
			}
		}
	}
	Sparse_Matrix* B = TransposeTimesSelf(&A, CCS, SYM_LOWER, true);
	solve_by_CHOLMOD(B);
	const std::vector<double>& tq = B->get_solution();

	for (int i = 0; i < nc; ++i)
	{
		int var_id = map_index[i];
		if (var_id != -1)
		{
			all_quaterion[i] = OpenVolumeMesh::Geometry::Vec4d(tq[var_id], tq[var_id + var_count], tq[var_id + var_count * 2], tq[var_id + var_count * 3]);
			all_quaterion[i].normalize();
		}
	}

	delete B;

	//convert quaternion to rotation matrix
	if (R.size() != nc) R.resize(nc);
	for (int i = 0; i < nc; ++i)
	{
		OpenVolumeMesh::Geometry::Vec4d& q = all_quaterion[i];
		Eigen::Matrix3d& Ri = R[i];
		double s = q[0]; double x = q[1]; double y = q[2]; double z = q[3];
		Ri(0, 0) = 1 - 2 * (y*y + z*z); Ri(0, 1) = 2 * x*y - 2 * s*z;   Ri(0, 2) = 2 * s*y + 2 * x*z;
		Ri(1, 0) = 2 * x*y + 2 * s*z;   Ri(1, 1) = 1 - 2 * (x*x + z*z); Ri(1, 2) = -2 * s*x + 2 * y*z;
		Ri(2, 0) = -2 * s*y + 2 * x*z;  Ri(2, 1) = 2 * s*x + 2 * y*z;   Ri(2, 2) = 1 - 2 * (y*y + x*x);
	}
}

//fix first point
void Polycube_Interface::ARAP_deform()
{
	OpenVolumeMesh::Geometry::Vec3d p0 = mesh_->vertex(OpenVolumeMesh::VertexHandle(0));//first vertex
	int nc = mesh_->n_cells(); int nv = mesh_->n_vertices();
	Sparse_Matrix C(nc*3, nv-1, NOSYM, CCS, 3);
	OpenVolumeMesh::Geometry::Vec4d cs;double px, py, pz;
	for (int c_id = 0; c_id < nc; ++c_id)
	{
		OpenVolumeMesh::Geometry::Vec4i& one_cv_id = cv_id[c_id];
		Eigen::Matrix3d& CP = PT[c_id]; double cv = std::sqrt(cell_volume[c_id]); cv = 1.0;
		Eigen::Matrix3d& CR = R[c_id];
		for (int j = 0; j < 3; ++j)
		{
			cs[1] = CP(0, j)*cv; cs[2] = CP(1, j)*cv; cs[3] = CP(2, j)*cv;
			cs[0] = -(cs[1] + cs[2] + cs[3]);
			for (int k = 0; k < 4; ++k)
			{
				if (one_cv_id[k] == 0) //0 means the first vertex in mesh, k means the k+1 vertex in the cell
				{
					C.fill_rhs_entry(3 * c_id + j, 0, -cs[k] * p0[0]);
					C.fill_rhs_entry(3 * c_id + j, 1, -cs[k] * p0[1]);
					C.fill_rhs_entry(3 * c_id + j, 2, -cs[k] * p0[2]);
				}
				else
				{
					C.fill_entry(3 * c_id + j, one_cv_id[k] - 1, cs[k]);
				}
			}
			C.fill_rhs_entry(3 * c_id + j, 0, cv*CR(0, j));
			C.fill_rhs_entry(3 * c_id + j, 1, cv*CR(1, j));
			C.fill_rhs_entry(3 * c_id + j, 2, cv*CR(2, j));
		}
	}
	Sparse_Matrix* D = TransposeTimesSelf(&C, CCS, SYM_LOWER, true);
	solve_by_CHOLMOD(D);
	const std::vector< double >& tv = D->get_solution();
	for (int i = 0; i < nv-1; ++i)
	{
		OpenVolumeMesh::VertexHandle vh( i + 1 );
		mesh_->set_vertex( vh, OpenVolumeMesh::Geometry::Vec3d( tv[i], tv[i + nv - 1], tv[i + 2 * nv - 2] ) );
	}

	//update PT
	Eigen::Matrix3d p; std::vector<OpenVolumeMesh::Geometry::Vec3d> cv_p(4);
	for (int c_id = 0; c_id < nc; ++c_id)
	{
		OpenVolumeMesh::Geometry::Vec4i& one_cv_id = cv_id[c_id];
		for (int j = 0; j < 4; ++j)
		{
			cv_p[j] = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_cv_id[j]));
		}
		p(0, 0) = cv_p[1][0] - cv_p[0][0]; p(1, 0) = cv_p[1][1] - cv_p[0][1]; p(2, 0) = cv_p[1][2] - cv_p[0][2];
		p(0, 1) = cv_p[2][0] - cv_p[0][0]; p(1, 1) = cv_p[2][1] - cv_p[0][1]; p(2, 1) = cv_p[2][2] - cv_p[0][2];
		p(0, 2) = cv_p[3][0] - cv_p[0][0]; p(1, 2) = cv_p[3][1] - cv_p[0][1]; p(2, 2) = cv_p[3][2] - cv_p[0][2];
		PT[c_id] = p.inverse();
	}
}

void Polycube_Interface::save_boundary_de(const char* filename)
{
	FILE* f_b_de = fopen(filename, "w");

	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		if (mesh_->is_boundary(*v_it))
		{
			int v_id = v_it->idx();
			OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
			fprintf(f_b_de, "%d %15.14f %15.14f %15.14f\n", v_id, p[0],p[1],p[2]);
		}
	}

	fclose(f_b_de);
}
void Polycube_Interface::save_result(const char* filename)
{
	FILE* f_r = fopen(filename, "w");

	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
		fprintf(f_r, "%15.14f %15.14f %15.14f\n", p[0], p[1], p[2]);
	}

	fclose(f_r);
}

void Polycube_Interface::load_vertex_color(const char* filename)
{
	FILE* f_vc = fopen(filename, "r");
	int nv = mesh_->n_vertices(); char buf[4096];
	int v_id_, c_; int number_of_color = 0;
	std::vector<int> vertex_color(nv, 0);
	char v_id[128]; char c[128]; int v_count = 0;
	while (!feof(f_vc))
	{
		fgets(buf, 4096, f_vc);
		if (v_count < nv)
		{
			sscanf(buf, "%s %s", v_id, c);
			v_id_ = atoi(v_id); c_ = atoi(c);
			vertex_color[v_id_] = c_;
			if (c_ > number_of_color) number_of_color = c_;
			++v_count;
		}
	}
	fclose(f_vc);
	printf("Number Of Colors: %d\n", number_of_color);

	pd_normal.set_vertex_color(number_of_color, vertex_color, mesh_);
	pd_flatten.set_vertex_color(number_of_color, vertex_color, mesh_);

	//prepare for deformation
	printf("Prepare Deformation.......\n");
	pd_normal.prepare_for_deformation(mesh_);
	pd_flatten.prepare_for_deformation(mesh_);
	pd_flatten.adjust_orientation(mesh_);
	pd_normal.assign_pos_mesh(mesh_, true);
	printf("Prepare Labeling.......\n");
	pc_label->prepare_data(mesh_);
	printf("Prepare Finishing.......\n");
	label_OK = false;

	emit uppdate_vertex_position_signal();
	emit updateGL_Manual_signal();
}

void Polycube_Interface::load_init_mesh(const char* filename)
{
	pd_normal.load_deformation_result(filename, mesh_);
	pd_flatten.load_deformation_result(filename, mesh_);

	emit uppdate_vertex_position_signal();
	emit updateGL_Manual_signal();
}

void Polycube_Interface::do_polycube_AMIPS(int max_iter, double sigma_s)
{
	pd_normal.set_sigma_s(sigma_s);

	//E:\\data\\Volume Mesh\\Polycube_IBCD\\Final Data\\elephant\\png
	//E:\\data\\Volume Mesh\\Polycube_IBCD\\Final Data\\rocker\\png
	/*QString mesh_name1("E:\\data\\Volume Mesh\\Polycube_IBCD\\Final Data\\rocker\\png\\png_00");
	QString mesh_name2("E:\\data\\Volume Mesh\\Polycube_IBCD\\Final Data\\rocker\\png\\png_0");
	QString mesh_name3("E:\\data\\Volume Mesh\\Polycube_IBCD\\Final Data\\rocker\\png\\png_");
	QString num = ""; QString extension = ".png"; QString image_name; int image_num = 0;
	num = QString::number(image_num);
	if (image_num < 10)
	{
		image_name = mesh_name1 + num + extension;
	}
	else if (image_num < 100)
	{
		image_name = mesh_name2 + num + extension;
	}
	else
	{
		image_name = mesh_name3 + num + extension;
	}

	emit save_opengl_screen_PC_signal(image_name);
	++image_num;*/

	long start_t = clock();
	for (int i = 0; i < max_iter; ++i)
	{
		//filter deformation
		// use_xyz = true,  use_RGNF = true,  use_half_sigma_s = false
		/*if(i==0) pd_normal.exp_mips_deformation_refine_polycube_omp(13, 15, 1.0, 0.5, mesh_, false, true, false);
		else pd_normal.exp_mips_deformation_refine_polycube_omp(12, 5, 1.0, 0.5, mesh_, false, true, false);*/

		//for animation
		if (i == 0)
		{
			for (int j = 0; j < 13; ++j)
			{
				pd_normal.exp_mips_deformation_refine_polycube_omp(1, 15, 1.0, 0.5, mesh_, false, true, false);
				pd_normal.assign_pos_mesh(mesh_);
				emit uppdate_vertex_position_signal();
				/*num = QString::number(image_num);
				if (image_num < 10)
				{
					image_name = mesh_name1 + num + extension;
				}
				else if (image_num < 100)
				{
					image_name = mesh_name2 + num + extension;
				}
				else
				{
					image_name = mesh_name3 + num + extension;
				}

				emit save_opengl_screen_PC_signal(image_name);
				++image_num;*/
			}
			
		}
		else
		{
			for (int j = 0; j < 12; ++j)
			{
				pd_normal.exp_mips_deformation_refine_polycube_omp(1, 5, 1.0, 0.5, mesh_, false, true, false);
				pd_normal.assign_pos_mesh(mesh_);
				emit uppdate_vertex_position_signal();
				/*num = QString::number(image_num);
				if (image_num < 10)
				{
				image_name = mesh_name1 + num + extension;
				}
				else if (image_num < 100)
				{
				image_name = mesh_name2 + num + extension;
				}
				else
				{
				image_name = mesh_name3 + num + extension;
				}

				emit save_opengl_screen_PC_signal(image_name);
				++image_num;*/
			}
		}

		pd_normal.assign_pos_mesh(mesh_);

		//labeling
		label_OK = pc_label->LabelingMain();
		printf("---------------------------------------------\n");
		if (label_OK)
		{
			printf("Labeling Success! CC: %d\n", pc_label->count_corner_number());
			break;
		}
		else
		{
			pc_label->get_charts_bound();
			pc_label->GetTetLabel();
			pc_label->GetEdgeChartTag();
			printf("Labeling Fail!\n");
		}
		printf("---------------------------------------------\n");

		pd_flatten.adjust_orientation(mesh_);
		//try to flatten boundary
		pd_flatten.compute_source_S(mesh_);
		//pd_flatten.assign_pos_mesh(mesh_, true);
		pd_flatten.exp_mips_deformation_refine_polycube_normal_omp(15, 1, 0.5, mesh_);
		pd_flatten.assign_pos_mesh(mesh_);
		pd_normal.assign_pos_mesh(mesh_, true);
	}

	long end_t = clock();
	long diff = end_t - start_t;
	double t = (double)(diff) / CLOCKS_PER_SEC;
	printf("Polycube Map : %f s\n", t);

	//labeling
	/*label_OK = pc_label->LabelingMain();
	printf("---------------------------------------------\n");
	if (label_OK)
	{
		printf("Labeling Success! CC: %d\n", pc_label->count_corner_number());
	}
	else
	{
		printf("Labeling Fail!\n");
	}
	printf("---------------------------------------------\n");*/

	pd_normal.assign_pos_mesh(mesh_);
	pd_normal.compute_distortion(mesh_);
	emit uppdate_vertex_position_signal();
}

void Polycube_Interface::filter_deform(double sigma_s)
{
	long start_t = clock();

	pd_normal.set_change_big_diff(pd_flatten.get_change_big_diff());
	pd_normal.set_sigma_s(sigma_s);
	/*for (int j = 0; j < 10; ++j)
	{*/
		pd_normal.exp_mips_deformation_refine_polycube_omp(13, 15, 1.0, 0.5, mesh_, false, true, false);
		//pd_normal.exp_mips_deformation_refine_polycube_omp(15, 1.0, 0.5, mesh_, true, true, false);
	//}
	//pd_normal.exp_mips_deformation_refine_polycube_omp(15, 1.0, 0.5, mesh_, true, true, true); //max
	long end_t = clock();
	long diff = end_t - start_t;
	double t = (double)(diff) / CLOCKS_PER_SEC;
	printf("Filter Deform : %f s\n", t);

	pd_normal.assign_pos_mesh(mesh_);
	pd_normal.compute_distortion(mesh_);
	emit uppdate_vertex_position_signal();

	label_OK = pc_label->LabelingMain();
	printf("---------------------------------------------\n");
	if (label_OK)
	{
		printf("Labeling Success! CC: %d\n", pc_label->count_corner_number());
	}
	else
	{
		pc_label->get_charts_bound();
		pc_label->GetTetLabel();
		pc_label->GetEdgeChartTag();
		printf("Labeling Fail!\n");
	}
	printf("---------------------------------------------\n");
}

void Polycube_Interface::flatten_normal()
{
	long start_t = clock();

	pd_flatten.adjust_orientation(mesh_);
	pd_flatten.compute_source_S(mesh_);
	//pd_flatten.assign_pos_mesh(mesh_, true);

	pd_flatten.exp_mips_deformation_refine_polycube_normal_omp(15, 1, 0.5, mesh_);
	pd_flatten.assign_pos_mesh(mesh_);

	long end_t = clock();
	long diff = end_t - start_t;
	double t = (double)(diff) / CLOCKS_PER_SEC;
	printf("Flatten Normal : %f s\n", t);

	pd_normal.assign_pos_mesh(mesh_, true);
	pd_normal.compute_distortion(mesh_);
	emit uppdate_vertex_position_signal();

	/*label_OK = pc_label->LabelingMain();
	printf("---------------------------------------------\n");
	if (label_OK)
	{
		printf("Labeling Success! CC: %d\n", pc_label->count_corner_number());
	}
	else
	{
		pc_label->get_charts_bound();
		pc_label->GetTetLabel();
		pc_label->GetEdgeChartTag();
		printf("Labeling Fail!\n");
	}
	printf("---------------------------------------------\n");*/
}


void Polycube_Interface::label_face()
{
	label_OK = pc_label->LabelingMain();
	//label_OK = pc_label->InitLabelMain();
	printf("---------------------------------------------\n");
	if (label_OK)
	{
		printf("Labeling Success! CC: %d\n", pc_label->count_corner_number());
	}
	else
	{
		printf("Labeling Fail!\n");
	}
	printf("---------------------------------------------\n");
}

void Polycube_Interface::save_label(const char* filename)
{
	FILE* f_label = fopen(filename, "w");
	std::vector<int>& tet_label = pc_label->tetLabel; //+X,-X,+Y,-Y,+Z,-Z
	std::vector<int>& chart_id = pc_label->tet_face_chart_id;
	for (int i = 0; i < tet_label.size(); ++i)
	{
		if (tet_label[i] < 0)
		{
			fprintf(f_label, "%d -1 -1\n", i);
		}
		else
		{
			fprintf(f_label, "%d %d %d\n", i, tet_label[i], chart_id[i]);
		}
	}
	fclose(f_label);
}