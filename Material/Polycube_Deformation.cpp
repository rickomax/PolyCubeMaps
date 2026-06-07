#include "Polycube_Deformation.h"
#include "omp.h"
#include <Qstring>
#include "..\..\BasicMath\MathSuite\SPARSE\Sparse_Matrix.h"
#include "..\..\BasicMath\MathSuite\SPARSE\Sparse_Solver.h"
#include "..\..\BasicMath\ConvexQuadOptimization.h"
#include "cholmod.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "mosek.h"
#include <Queue>
#include "..\..\BasicMath\fmath.hpp"
#define Hex_Meshing

polycube_deformation_interface::polycube_deformation_interface()
{
	prepare_ok = false;
	sigma_s = 1.0; sigma_r = 0.3; RGNF_iter_count = 1;
}

polycube_deformation_interface::~polycube_deformation_interface()
{
}

void polycube_deformation_interface::load_boundary_face_label(const char* filename, VolumeMesh* mesh_)
{
	int nc = mesh_->n_cells();
	target_bfn.resize(nc, OpenVolumeMesh::Geometry::Vec3d(0,0,0));
	bf_chart.clear(); bf_chart.resize(nc, -1);
	std::vector<int> chart_xyz(nc, -1);
	FILE* f_bfn = fopen(filename, "r");
	char buf[4096] = "";
	char label[128]; char c_id[128]; char chart[128]; int max_chart_id = 0;
	while (!feof(f_bfn))
	{
		fgets(buf, 4096, f_bfn);
		if (buf[0] == '\n' || buf[0] == '\r' || buf[0] == '\0') continue;

		sscanf(buf, "%s %s %s", c_id, label, chart);
		int l = atoi(label); int cid = atoi(c_id);
		switch (l)
		{
		case 0:
			target_bfn[cid][0] = 1.0;//x
			chart_xyz[cid] = 0;
			break;
		case 1:
			target_bfn[cid][1] = 1.0;//y
			chart_xyz[cid] = 2;
			break;
		case 2:
			target_bfn[cid][2] = 1.0;//z
			chart_xyz[cid] = 4;
			break;
		case 3:
			target_bfn[cid][0] = -1.0;//-x
			chart_xyz[cid] = 1;
			break;
		case 4:
			target_bfn[cid][1] = -1.0;//-y
			chart_xyz[cid] = 3;
			break;
		case 5:
			target_bfn[cid][2] = -1.0;//-z
			chart_xyz[cid] = 5;
			break;
		default:
			break;
		}

		bf_chart[cid] = atoi(chart);
		if (bf_chart[cid] > max_chart_id)
		{
			max_chart_id = bf_chart[cid];
		}
	}
	fclose(f_bfn);
	printf("Chart Number : %d\n", max_chart_id + 1);

	polycube_chart.clear(); polycube_chart.resize(max_chart_id+1);
	polycube_chart_label.clear(); polycube_chart_label.resize(max_chart_id + 1);
	chart_mean_value.clear(); chart_mean_value.resize(max_chart_id + 1, 0.0);
	for (int i = 0; i < nc; ++i)
	{
		int chart_id = bf_chart[i];
		if (chart_id < 0) continue;
		polycube_chart[chart_id].push_back(i);
		polycube_chart_label[chart_id] = chart_xyz[i];
	}
	int nv = bvf_id.size(); 
	vertex_type.clear(); vertex_type.resize(nv, OpenVolumeMesh::Geometry::Vec3i(-1,-1,-1));
	for (int i = 0; i < nv; ++i)
	{
		std::vector<int>& one_vf_id = bvf_id[i];
		int bvf_size = one_vf_id.size();
		if (bvf_size > 0)
		{
			for (int j = 0; j < bvf_size; ++j)
			{
				vertex_type[i][chart_xyz[one_vf_id[j]]/2] = bf_chart[ one_vf_id[j] ];
			}
		}
	}

	int ne = bef_id.size(); 
	std::vector<OpenVolumeMesh::Geometry::Vec3i> edge_two_chart;
	for (int i = 0; i< ne; ++i)
	{
		std::vector<int>& one_ef_id = bef_id[i];
		int bef_size = one_ef_id.size();
		if (bef_size == 2)
		{
			if (bf_chart[one_ef_id[0]] != bf_chart[one_ef_id[1]] && bf_chart[one_ef_id[1]] >= 0 && bf_chart[one_ef_id[0]] >= 0)
			{
				int l0 = bf_chart[one_ef_id[0]]; int l1 = bf_chart[one_ef_id[1]];
				if (l0 < l1)
				{
					edge_two_chart.push_back(OpenVolumeMesh::Geometry::Vec3i(l0, l1, i));
				}
				else
				{
					edge_two_chart.push_back(OpenVolumeMesh::Geometry::Vec3i(l1, l0, i));
				}
			}
		}
		else if (bef_size != 0)
		{
			printf("Error boundary edge face %d\n", bef_size);
		}
	}

	std::vector<int> empty_(1);
	std::vector<OpenVolumeMesh::Geometry::Vec2i> two_label; std::vector<std::vector<int>> edge_with_same_label;
	for (int i = 0; i < edge_two_chart.size(); ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& two_chart = edge_two_chart[i]; 
		int inserted_id = -1;
		for (int j = 0; j < two_label.size(); ++j)
		{
			if (two_chart[0] == two_label[j][0] && two_chart[1] == two_label[j][1]) //find label
			{
				inserted_id = j;
				edge_with_same_label[j].push_back(two_chart[2]);
			}
		}
		if (inserted_id < 0)
		{
			two_label.push_back(OpenVolumeMesh::Geometry::Vec2i(two_chart[0], two_chart[1]));
			empty_[0] = two_chart[2];
			edge_with_same_label.push_back(empty_);
		}
	}

	up_chart_id.clear(); down_chart_id.clear(); diff_up_down.clear();
	for (int i = 0; i < edge_with_same_label.size(); ++i)
	{
		std::vector<int> two_corner;
		for (int j = 0; j < edge_with_same_label[i].size(); ++j)
		{
			int edge_id = edge_with_same_label[i][j];
			OpenVolumeMesh::OpenVolumeMeshEdge edge = mesh_->edge(OpenVolumeMesh::EdgeHandle(edge_id));
			int fv = edge.from_vertex().idx(); int tv = edge.to_vertex().idx();
			if (vertex_type[fv][0] >= 0 && vertex_type[fv][1] >= 0 && vertex_type[fv][2] >= 0)
			{
				two_corner.push_back(fv);
			}
			if (vertex_type[tv][0] >= 0 && vertex_type[tv][1] >= 0 && vertex_type[tv][2] >= 0)
			{
				two_corner.push_back(tv);
			}
		}

		std::vector<OpenVolumeMesh::Geometry::Vec2i> final_tc;
		if (two_corner.size() != 2)
		{
			printf("Error corner number %d\n", two_corner.size());
			std::vector<int> visited(two_corner.size(), -1);
			for (int ijk = 0; ijk < two_corner.size(); ++ijk)
			{
				if (visited[ijk] == 1) continue;

				visited[ijk] = 1; int seed_v = two_corner[ijk]; int last_v = -1;//std::queue<int> Q; Q.push(two_corner[0]);
				std::vector<int>& one_vv = bvv_id[seed_v]; int label_count = 0;
				for (int j = 0; j < one_vv.size(); ++j)
				{
					label_count = 0;
					for (int k = 0; k < 3; ++k)
					{
						if (vertex_type[one_vv[j]][k] == two_label[i][0] || vertex_type[one_vv[j]][k] == two_label[i][1])
						{
							++label_count;
						}
					}
					if (label_count == 2)
					{
						last_v = seed_v;
						seed_v = one_vv[j];
						break;
					}
				}
				if (vertex_type[seed_v][0] >= 0 && vertex_type[seed_v][1] >= 0 && vertex_type[seed_v][2] >= 0)
				{
					final_tc.push_back(OpenVolumeMesh::Geometry::Vec2i(two_corner[ijk], seed_v));
					for (int kj = 0; kj < two_corner.size(); ++kj)
					{
						if (two_corner[kj] == seed_v) visited[kj] = 1;
					}
				}
				else
				{
					while (1)
					{
						one_vv = bvv_id[seed_v]; 
						for (int j = 0; j < one_vv.size(); ++j)
						{
							label_count = 0;
							for (int k = 0; k < 3; ++k)
							{
								if (vertex_type[one_vv[j]][k] == two_label[i][0] || vertex_type[one_vv[j]][k] == two_label[i][1])
								{
									++label_count;
								}
							}
							if (label_count == 2 && one_vv[j] != last_v)
							{
								last_v = seed_v;
								seed_v = one_vv[j];
								break;
							}
						}
						if (vertex_type[seed_v][0] >= 0 && vertex_type[seed_v][1] >= 0 && vertex_type[seed_v][2] >= 0)
						{
							final_tc.push_back(OpenVolumeMesh::Geometry::Vec2i(two_corner[ijk], seed_v));
							for (int kj = 0; kj < two_corner.size(); ++kj)
							{
								if (two_corner[kj] == seed_v) visited[kj] = 1;
							}
							break;
						}
					}
				}
			}
		}
		else
		{
			final_tc.push_back(OpenVolumeMesh::Geometry::Vec2i(two_corner[0], two_corner[1]));
		}

		for (int j = 0; j < final_tc.size(); ++j)
		{
			int v0 = final_tc[j][0]; int v1 = final_tc[j][1];
			int l0 = -1; int l1 = -1; int k0 = -1; int k1 = -1;
			for (int j = 0; j < 3; ++j)
			{
				if (vertex_type[v0][j] != two_label[i][0] && vertex_type[v0][j] != two_label[i][1])
				{
					l0 = vertex_type[v0][j]; k0 = j;
				}
				if (vertex_type[v1][j] != two_label[i][0] && vertex_type[v1][j] != two_label[i][1])
				{
					l1 = vertex_type[v1][j]; k1 = j;
				}
			}
			if (k0 != k1)
			{
				printf("Error corner label %d %d %d %d\n", v0, k0, v1, k1);
			}
			OpenVolumeMesh::Geometry::Vec3d p0 = mesh_->vertex(OpenVolumeMesh::VertexHandle(v0));
			OpenVolumeMesh::Geometry::Vec3d p1 = mesh_->vertex(OpenVolumeMesh::VertexHandle(v1));
			if (p0[k0] > p1[k1])
			{
				up_chart_id.push_back(l0); down_chart_id.push_back(l1); diff_up_down.push_back(p0[k0] - p1[k1]);
			}
			else
			{
				up_chart_id.push_back(l1); down_chart_id.push_back(l0); diff_up_down.push_back(p1[k1] - p0[k0]);
			}
		}
	}
}

void polycube_deformation_interface::load_up_down_chart(const char* filename)
{
	FILE* f_up_down = fopen(filename, "r");
	char buf[4096] = "";
	char up_chart[128]; char down_chart[128]; char diff_chart[128];
	fgets(buf, 4096, f_up_down); int up_down_count = atoi(buf);
	up_chart_id.clear(); up_chart_id.resize(up_down_count);
	down_chart_id.clear(); down_chart_id.resize(up_down_count);
	diff_up_down.clear(); diff_up_down.resize(up_down_count);
	int count = 0;
	while (!feof(f_up_down))
	{
		fgets(buf, 4096, f_up_down);
		sscanf(buf, "%s %s %s", up_chart, down_chart, diff_chart);
		up_chart_id[count] = atoi(up_chart);
		down_chart_id[count] = atoi(down_chart);
		diff_up_down[count] = atof(diff_chart);
		++count;
		if (count == up_down_count) break;
	}
	fclose(f_up_down);
}

void polycube_deformation_interface::compute_source_S(VolumeMesh* mesh_)
{
	unsigned nv = mesh_->n_vertices();
	Eigen::Matrix3d VS, IS; 
	dpx.resize(nv); dpy.resize(nv); dpz.resize(nv); 
	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		int v_id = v_it->idx();
		OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
		dpx[v_id] = p[0]; dpy[v_id] = p[1]; dpz[v_id] = p[2];

		int vc_size = vertex_cell[v_id].size();
		for (int i = 0; i < vc_size; ++i)
		{
			const int* vv_id = vertex_cell_vertex[v_id][i].data();
			OpenVolumeMesh::Geometry::Vec3d r = mesh_->vertex(OpenVolumeMesh::VertexHandle(vv_id[0])) - p;
			OpenVolumeMesh::Geometry::Vec3d s = mesh_->vertex(OpenVolumeMesh::VertexHandle(vv_id[1])) - p;
			OpenVolumeMesh::Geometry::Vec3d t = mesh_->vertex(OpenVolumeMesh::VertexHandle(vv_id[2])) - p;
			VS(0, 0) = r[0]; VS(1, 0) = r[1]; VS(2, 0) = r[2];
			VS(0, 1) = s[0]; VS(1, 1) = s[1]; VS(2, 1) = s[2];
			VS(0, 2) = t[0]; VS(1, 2) = t[1]; VS(2, 2) = t[2];
			IS = VS.inverse();
			std::vector<double>& cv_s = vcv_S[v_id][i];
			cv_s[0] = IS(0, 0); cv_s[1] = IS(0, 1); cv_s[2] = IS(0, 2);
			cv_s[3] = IS(1, 0); cv_s[4] = IS(1, 1); cv_s[5] = IS(1, 2);
			cv_s[6] = IS(2, 0); cv_s[7] = IS(2, 1); cv_s[8] = IS(2, 2);
		}
	}
}

void polycube_deformation_interface::prepare_for_deformation(VolumeMesh* mesh_)
{
	std::vector<OpenVolumeMesh::VertexHandle> hfv_vec(3);
	std::vector<int> hfv_vec_id(3);
	std::vector<OpenVolumeMesh::HalfFaceHandle> cell_hfh_vec(4);

	unsigned nv = mesh_->n_vertices(); unsigned nc = mesh_->n_cells();
	vertex_cell.clear(); cell_vertex.clear();
	vertex_cell.resize(nv); cell_vertex.resize(nc);
	cell_vertex_vertex.clear(); cell_S.resize(nc);
	cell_vertex_vertex.resize(nc); 
	vertex_cell_vertex.clear(); vertex_cell_vertex.resize(nv);
	vcv_S.clear(); vcv_S.resize(nv);
	cell_volume.clear(); cell_volume.resize(nc, 0.0);
	Eigen::Matrix3d VS, IS; std::vector<double> S_(9);
	for (OpenVolumeMesh::CellIter c_it = mesh_->cells_begin(); c_it != mesh_->cells_end(); ++c_it)
	{
		double cv_count = 0.0; int c_id = c_it->idx();
		const std::vector<OpenVolumeMesh::VertexHandle>& vertices_ = mesh_->cell(*c_it).vertices();
		cell_vertex[c_id].clear();
		for (unsigned i = 0; i < vertices_.size(); ++i)
		{
			int v_id = vertices_[i];
			cell_vertex[c_id].push_back(v_id);
		}
		OpenVolumeMesh::Geometry::Vec3d cp = mesh_->vertex(vertices_[0]);
		OpenVolumeMesh::Geometry::Vec3d cr = mesh_->vertex(vertices_[1]) - cp;
		OpenVolumeMesh::Geometry::Vec3d cs = mesh_->vertex(vertices_[2]) - cp;
		OpenVolumeMesh::Geometry::Vec3d ct = mesh_->vertex(vertices_[3]) - cp;
		VS(0, 0) = cr[0]; VS(1, 0) = cr[1]; VS(2, 0) = cr[2];
		VS(0, 1) = cs[0]; VS(1, 1) = cs[1]; VS(2, 1) = cs[2];
		VS(0, 2) = ct[0]; VS(1, 2) = ct[1]; VS(2, 2) = ct[2];
		IS = VS.inverse();
		cell_S[c_id] = IS;

		mesh_->get_halffaces_from_cell(*c_it, cell_hfh_vec);
		for (unsigned i = 0; i < vertices_.size(); ++i)
		{
			int v_id = vertices_[i];
			if (v_id < 0 || v_id > nv)
			{
				printf("%d\n", v_id);
			}

			if (vertex_cell[v_id].size() == 0) vertex_cell[v_id].reserve(15);
			vertex_cell[v_id].push_back(c_it->idx());

			for (unsigned j = 0; j < cell_hfh_vec.size(); ++j)
			{
				bool find_hf = true;
				mesh_->get_vertices_from_halfface(cell_hfh_vec[j], hfv_vec);
				for (unsigned jj = 0; jj < hfv_vec.size(); ++jj)
				{
					if (hfv_vec[jj] == v_id)
					{
						find_hf = false;
						break;
					}
				}
	
				if (find_hf)
				{
					hfv_vec_id[0] = hfv_vec[0].idx(); hfv_vec_id[1] = hfv_vec[1].idx(); hfv_vec_id[2] = hfv_vec[2].idx();
					cell_vertex_vertex[c_id].push_back(hfv_vec_id);
					vertex_cell_vertex[v_id].push_back(hfv_vec_id);
					OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(vertices_[i]);
					OpenVolumeMesh::Geometry::Vec3d r = mesh_->vertex(hfv_vec[0]) - p;
					OpenVolumeMesh::Geometry::Vec3d s = mesh_->vertex(hfv_vec[1]) - p;
					OpenVolumeMesh::Geometry::Vec3d t = mesh_->vertex(hfv_vec[2]) - p;
					VS(0, 0) = r[0]; VS(1, 0) = r[1]; VS(2, 0) = r[2];
					VS(0, 1) = s[0]; VS(1, 1) = s[1]; VS(2, 1) = s[2];
					VS(0, 2) = t[0]; VS(1, 2) = t[1]; VS(2, 2) = t[2];
					
					cell_volume[c_id] = -VS.determinant();
					IS = VS.inverse();
					S_[0] = IS(0, 0); S_[1] = IS(0, 1); S_[2] = IS(0, 2); 
					S_[3] = IS(1, 0); S_[4] = IS(1, 1); S_[5] = IS(1, 2);
					S_[6] = IS(2, 0); S_[7] = IS(2, 1); S_[8] = IS(2, 2);
					vcv_S[v_id].push_back(S_);
					break;
				}
			}
		}
	}

	dpx.resize(nv); dpy.resize(nv); dpz.resize(nv); max_vc_size = 0; src_pos.resize(nv);
	is_bv.clear(); is_bv.resize(nv, -1);
	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		int v_id = v_it->idx();
		OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
		dpx[v_id] = p[0]; dpy[v_id] = p[1]; dpz[v_id] = p[2];
		src_pos[v_id] = p;
		int vc_size = vertex_cell[v_id].size();
		if (vc_size > max_vc_size)
		{
			max_vc_size = vc_size;
		}

		if (mesh_->is_boundary(*v_it))
		{
			is_bv[v_id] = 1;
		}
	}

	int max_num_t = omp_get_max_threads();
	vc_pos_x_omp.resize(max_num_t); vc_pos_y_omp.resize(max_num_t); vc_pos_z_omp.resize(max_num_t);
	vc_pos_x2_omp.resize(max_num_t); vc_pos_y2_omp.resize(max_num_t); vc_pos_z2_omp.resize(max_num_t);
	vc_pos_x3_omp.resize(max_num_t); vc_pos_y3_omp.resize(max_num_t); vc_pos_z3_omp.resize(max_num_t);
	vc_n_cross_x_omp.resize(max_num_t); vc_n_cross_y_omp.resize(max_num_t); vc_n_cross_z_omp.resize(max_num_t);
	vc_S_omp.resize(max_num_t); exp_vec_omp.resize(max_num_t); large_dis_flag_omp.resize(max_num_t);
	gx_vec_omp.resize(max_num_t); gy_vec_omp.resize(max_num_t); gz_vec_omp.resize(max_num_t);
	mu_vec_omp.resize(max_num_t);
	for (int i = 0; i < max_num_t; ++i)
	{
		vc_pos_x_omp[i].resize(max_vc_size); vc_pos_y_omp[i].resize(max_vc_size); vc_pos_z_omp[i].resize(max_vc_size);
		vc_pos_x2_omp[i].resize(4 * max_vc_size); vc_pos_y2_omp[i].resize(4 * max_vc_size); vc_pos_z2_omp[i].resize(4 * max_vc_size);
		vc_pos_x3_omp[i].resize(4 * max_vc_size); vc_pos_y3_omp[i].resize(4 * max_vc_size); vc_pos_z3_omp[i].resize(4 * max_vc_size);
		vc_n_cross_x_omp[i].resize(max_vc_size); vc_n_cross_y_omp[i].resize(max_vc_size); vc_n_cross_z_omp[i].resize(max_vc_size);
		vc_S_omp[i].resize(max_vc_size);  exp_vec_omp[i].resize(max_vc_size); large_dis_flag_omp[i].resize(max_vc_size);
		gx_vec_omp[i].resize(max_vc_size); gy_vec_omp[i].resize(max_vc_size); gz_vec_omp[i].resize(max_vc_size);
		mu_vec_omp[i].resize(max_vc_size);
		for (int j = 0; j < max_vc_size; ++j)
		{
			vc_S_omp[i][j].resize(9);
		}
	}

	//2015.11.25 boundary vertex information
	bfv_id.clear(); bfv_id.resize(nc, OpenVolumeMesh::Geometry::Vec3i(-1, -1, -1));
	bvf_id.clear(); bvf_id.resize(nv); boundary_face_number = 0; avg_boundary_edge_length = 0.0;
	bef_id.clear(); bef_id.resize(mesh_->n_edges());
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
				bvf_id[hfv_id].push_back(c_id);
				++hfv_count;
			}

			OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(*f_it);
			std::vector<OpenVolumeMesh::HalfEdgeHandle>& heh_vec = face.get_halfedges();
			for (int i = 0; i < heh_vec.size(); ++i)
			{
				int e_id = mesh_->edge_handle(heh_vec[i]);
				bef_id[e_id].push_back(c_id);
			}

			++boundary_face_number;
		}
		else if (ch1 == VolumeMesh::InvalidCellHandle && ch0 != VolumeMesh::InvalidCellHandle)
		{
			int hfv_count = 0; int c_id = ch0.idx();
			for (OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh_->hfv_iter(hfh1); hfv_it; ++hfv_it)
			{
				int hfv_id = hfv_it->idx();
				bfv_id[c_id][hfv_count] = hfv_id;
				bvf_id[hfv_id].push_back(c_id);
				++hfv_count;
			}

			OpenVolumeMesh::OpenVolumeMeshFace face = mesh_->face(*f_it);
			std::vector<OpenVolumeMesh::HalfEdgeHandle>& heh_vec = face.get_halfedges();
			for (int i = 0; i < heh_vec.size(); ++i)
			{
				int e_id = mesh_->edge_handle(heh_vec[i]);
				bef_id[e_id].push_back(c_id);
			}

			++boundary_face_number;
		}
		else if (ch1 == VolumeMesh::InvalidCellHandle && ch0 == VolumeMesh::InvalidCellHandle)
		{
			printf("Error : Both two halffaces have no cells!\n");
		}
	}
	bvv_id.clear(); bvv_id.resize(nv);
	for (int i = 0; i < nv; ++i)
	{
		if (bvf_id[i].size() > 0) //boundary vertex
		{
			std::vector<int>& one_ring_f_id = bvf_id[i]; int one_ring_size = one_ring_f_id.size();
			std::vector<int> one_ring_v_id(one_ring_size); std::vector<int> one_ring_v_f_id(one_ring_size);
			std::vector<int> one_ring_order(one_ring_size);
			for (int j = 0; j < one_ring_size; ++j)
			{
				OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[one_ring_f_id[j]];
				for (int k = 0; k < 3; ++k)
				{
					if (i == fv_id[k])
					{
						one_ring_order[j] = k;
						break;
					}
				}
			}
			OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[one_ring_f_id[0]];
			int start_v = fv_id[(one_ring_order[0] + 1) % 3];
			int v2 = start_v; one_ring_v_id[0] = v2; one_ring_v_f_id[0] = one_ring_f_id[0];
			int v3 = fv_id[(one_ring_order[0] + 2) % 3];
			for (int j = 1; j < one_ring_size; ++j)
			{
				v2 = v3;
				if (v2 == start_v) break;
				one_ring_v_id[j] = v2;
				for (int k = 1; k < one_ring_size; ++k)// face k
				{
					OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[one_ring_f_id[k]];
					if (fv_id[(one_ring_order[k] + 1) % 3] == v2)
					{
						v3 = fv_id[(one_ring_order[k] + 2) % 3];
						one_ring_v_f_id[j] = one_ring_f_id[k];
						break;
					}
				}
			}
			//printf("%d\n", i);
			bvv_id[i] = one_ring_v_id;
			bvf_id[i] = one_ring_v_f_id;
		}
	}

	avg_boundary_edge_length = 0.0; OpenVolumeMesh::Geometry::Vec3d p0, p1, p2; double count = 0.0;
	for (int i = 0; i < bfv_id.size(); ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
		if (one_bfv_id[0] < 0) continue;

		p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
		p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
		p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];
		avg_boundary_edge_length += (p0 - p1).norm();
		avg_boundary_edge_length += (p1 - p2).norm();
		avg_boundary_edge_length += (p2 - p0).norm();
		count += 3.0;
	}
	avg_boundary_edge_length /= count;
	printf("Average Boundary Edge length : %e\n", avg_boundary_edge_length);

	change_big_flag.clear(); last_bfn.clear();
	change_big_flag.resize(nc, -1); last_bfn.resize(nc);
	bound_K = 1.0;

	prepare_ok = true;
}

void polycube_deformation_interface::assign_pos_mesh(VolumeMesh* mesh_, bool r_order)
{
	if (r_order)
	{
		for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
		{
			int v_id = v_it->idx();
			OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
			dpx[v_id] = p[0]; dpy[v_id] = p[1]; dpz[v_id] = p[2];
		}
	}
	else
	{
		for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
		{
			int v_id = v_it->idx();
			OpenVolumeMesh::Geometry::Vec3d p(dpx[v_id], dpy[v_id], dpz[v_id]);
			mesh_->set_vertex(*v_it, p);
		}
	}
}

void polycube_deformation_interface::compute_distortion(VolumeMesh* mesh_)
{
	Eigen::Matrix3d VS, A; Eigen::Vector3d a;
	double max_cd = 0.0; double min_cd = 1e30; double avg_cd = 0.0; int max_cd_cell_id = -1;
	double max_iso = 0.0; double min_iso = 1e30; double avg_iso = 0.0; int max_iso_cell_id = -1;
	double max_vol = 0.0; double min_vol = 1e30; double avg_vol = 0.0; int max_vol_cell_id = -1;
	double cd_count = 0.0; int flip_count = 0; flipped_cell.clear();
	std::vector<double> all_iso_d; std::vector<double> all_con_d; std::vector<double> all_vol_d;
	int nc = mesh_->n_cells();
	for (int c_id = 0; c_id < nc; ++c_id)
	{
		std::vector<int>& cv_id = cell_vertex[c_id];
		double cpx = dpx[cv_id[0]]; double cpy = dpy[cv_id[0]]; double cpz = dpz[cv_id[0]];
		VS(0, 0) = dpx[cv_id[1]] - cpx; VS(1, 0) = dpy[cv_id[1]] - cpy; VS(2, 0) = dpz[cv_id[1]] - cpz;
		VS(0, 1) = dpx[cv_id[2]] - cpx; VS(1, 1) = dpy[cv_id[2]] - cpy; VS(2, 1) = dpz[cv_id[2]] - cpz;
		VS(0, 2) = dpx[cv_id[3]] - cpx; VS(1, 2) = dpy[cv_id[3]] - cpy; VS(2, 2) = dpz[cv_id[3]] - cpz;
		A = VS * cell_S[c_id];

		if (A.determinant() < 0)
		{
			++flip_count;
			flipped_cell.push_back(OpenVolumeMesh::CellHandle(c_id) );
			all_iso_d.push_back(-1);
			all_con_d.push_back(-1);
			all_vol_d.push_back(-1);
		}
		else
		{
			Eigen::JacobiSVD<Eigen::Matrix3d> svd(A);
			a = svd.singularValues();
			double cd = a(0) / a(2);
			if (cd > max_cd) { max_cd = cd; max_cd_cell_id = c_id; }
			if (cd < min_cd) min_cd = cd;
			avg_cd += cd; cd_count += 1.0;
			all_con_d.push_back(cd);

			double iso_d = a(0) > 1.0 / a(2) ? a(0) : 1.0 / a(2);
			if (iso_d > max_iso){ max_iso = iso_d; max_iso_cell_id = c_id; }
			if (iso_d < min_iso) min_iso = iso_d;
			avg_iso += iso_d;
			all_iso_d.push_back(iso_d);

			double dJ = a(0)*a(1)*a(2);

			double vol_d = dJ > 1.0 / dJ ? dJ : 1.0 / dJ;
			if (vol_d > max_vol){ max_vol = vol_d; max_vol_cell_id = c_id; }
			if (vol_d < min_vol) min_vol = vol_d;
			avg_vol += vol_d;
			all_vol_d.push_back(vol_d);
		}
	}

	avg_cd /= cd_count; avg_iso /= cd_count; avg_vol /= cd_count;
	double std_id = 0.0; double std_cd = 0.0; double std_vd = 0.0;
	for (int i = 0; i < all_iso_d.size(); ++i)
	{
		if (all_iso_d[i]>0)
		{
			std_cd += (avg_cd - all_con_d[i])*(avg_cd - all_con_d[i]);
			std_id += (avg_iso - all_iso_d[i])*(avg_iso - all_iso_d[i]);
			std_vd += (avg_vol - all_vol_d[i])*(avg_vol - all_vol_d[i]);
		}
	}
	std_id = std::sqrt(std_id / (cd_count - 1)); std_cd = std::sqrt(std_cd / (cd_count - 1)); std_vd = std::sqrt(std_vd / (cd_count - 1));
	printf("---------------------------------------------------------\n");
	printf("Flip Count : %d\n", flip_count);
	printf("Isometric Distortion: %f/%f/%f/%f\n", max_iso, min_iso, avg_iso, std_id);
	printf("conformal Distortion: %f/%f/%f/%f\n", max_cd, min_cd, avg_cd, std_cd);
	printf("volume Distortion: %f/%f/%f/%f\n", max_vol, min_vol, avg_vol, std_vd);
}

double polycube_deformation_interface::compute_energy_ratio_filter_deform()
{
	int nc = cell_S.size(); int nv = bvv_id.size();
	//distortion
	int max_num_t = omp_get_max_threads();
	std::vector<Eigen::Matrix3d> VS(max_num_t), A(max_num_t), A2(max_num_t);
	std::vector<double> distortion_v(nc, -1.0); 
#pragma omp parallel for schedule(dynamic, 1)
	for (int c_id = 0; c_id < nc; ++c_id)
	{
		int id = omp_get_thread_num();
		//int id = 0;

		const int* cv_id = cell_vertex[c_id].data();
		int v_id_0 = cv_id[0]; int v_id_1 = cv_id[1]; int v_id_2 = cv_id[2]; int v_id_3 = cv_id[3];
		double x0 = dpx[v_id_0]; double y0 = dpy[v_id_0]; double z0 = dpz[v_id_0];

		VS[id](0, 0) = dpx[v_id_1] - x0; VS[id](1, 0) = dpy[v_id_1] - y0; VS[id](2, 0) = dpz[v_id_1] - z0;
		VS[id](0, 1) = dpx[v_id_2] - x0; VS[id](1, 1) = dpy[v_id_2] - y0; VS[id](2, 1) = dpz[v_id_2] - z0;
		VS[id](0, 2) = dpx[v_id_3] - x0; VS[id](1, 2) = dpy[v_id_3] - y0; VS[id](2, 2) = dpz[v_id_3] - z0;
		A[id] = VS[id] * cell_S[c_id];

		A2[id] = A[id].transpose()*A[id];
		double AF = A2[id](0, 0) + A2[id](1, 1) + A2[id](2, 2);
		double AF_05 = std::sqrt(AF);
		double AF2 = A2[id].squaredNorm();
		double AF_I = (AF*AF - AF2)*0.5; double AF_I_05 = std::sqrt(AF_I);
		double det_A = A[id].determinant();
		double i_det_A = 1.0 / det_A;

		double g = AF_05 * AF_I_05;
		double e = g*i_det_A;

		double k = ((e*e - 1.0)*0.125 + 0.5*(det_A + i_det_A)) * 0.5;
		if (k > 60) k = 60;
		distortion_v[c_id] = std::exp(k);
	}

	if (distoriton_big_count.size() != nv)
	{
		distoriton_big_count.clear();
		distoriton_big_count.resize(nv, 0);
	}
	if (distoriton_big_cell.size() != nc)
	{
		distoriton_big_cell.clear();
		distoriton_big_cell.resize(nc, -1);
	}
	double large_dis_th = std::exp(bound_K);
	//boundary normal energy
	double e_sum = 0.0;
	double n_diff = 0.0; double count_diff = 0.0; local_energy_ratio.clear(); local_energy_ratio.resize(nc, -1);
	OpenVolumeMesh::Geometry::Vec3d p0, p1, p2; double max_ratio = 0.0; double ratio = 0.0;
	std::vector<double> n_diff_v(nc, 10.0);
	for (int i = 0; i < nc; ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[i];
		if (fv_id[0] >= 0)
		{
			OpenVolumeMesh::Geometry::Vec3d& ln = target_bfn[i];
			p0[0] = dpx[fv_id[0]]; p0[1] = dpy[fv_id[0]]; p0[2] = dpz[fv_id[0]];
			p1[0] = dpx[fv_id[1]]; p1[1] = dpy[fv_id[1]]; p1[2] = dpz[fv_id[1]];
			p2[0] = dpx[fv_id[2]]; p2[1] = dpy[fv_id[2]]; p2[2] = dpz[fv_id[2]];

			OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0).normalize();
			n_diff += 0.5*(n - ln).sqrnorm();
			count_diff += 1.0;

			n_diff_v[i] = 0.5*(n - ln).sqrnorm();
			e_sum += distortion_v[i];// / (1e-10 + (n - ln).sqrnorm());
			local_energy_ratio[i] = distortion_v[i] / (1e-10 + 0.5*(n - ln).sqrnorm());
			if (local_energy_ratio[i] > 1e16) local_energy_ratio[i] = 1e16;
		}
		
		if (distortion_v[i] > large_dis_th)
		{
			distoriton_big_cell[i] = 1;
			const int* cv_id = cell_vertex[i].data();
			for (int j = 0; j < 4; ++j)
			{
				distoriton_big_count[cv_id[j]] += 1;
			}
		}
		
	}

	n_diff /= count_diff;
	e_sum /= count_diff;
	double s = e_sum / (n_diff + 1e-8);
	printf("1 %f %f %f\n", e_sum, n_diff, s);

	//return s;

	//2.0 and 0.1
	double up_diff_th = 2.0*n_diff; double down_diff_th = 0.1*n_diff;
	double avg_n_diff = 0.0; count_diff = 0.0; e_sum = 0.0;
	for (int i = 0; i<nc; ++i)
	{
		if (n_diff_v[i] < up_diff_th && n_diff_v[i] > down_diff_th)
		{
			avg_n_diff += n_diff_v[i]; count_diff += 1.0;
			e_sum += distortion_v[i];
		}
	} 
	avg_n_diff /= count_diff;
	e_sum /= count_diff;

	s = e_sum / (avg_n_diff + 1e-8);
	printf("2 %f %f %f\n", e_sum, avg_n_diff, s);
	printf("--------------------------------------------\n");
	if (s > 1e16) s = 1e16;
	//if (s < 100) s = 100;

	for (int i = 0; i < nc; ++i)
	{
		if (local_energy_ratio[i] > 2*s)//from 1.5 to 10
		{
			local_energy_ratio[i] = 2*s;
		}
		else if (local_energy_ratio[i] < 0.1* s)
		{
			local_energy_ratio[i] = s*0.1;
		}
	}

	return s;
}

double polycube_deformation_interface::compute_energy_ratio_flatten_normal()
{
	//distortion
	int max_num_t = omp_get_max_threads();
	std::vector<Eigen::Matrix3d> VS(max_num_t), A(max_num_t), A2(max_num_t);
	int nc = cell_S.size();
	std::vector<double> e_sum(max_num_t, 0);
#pragma omp parallel for schedule(dynamic, 1)
	for (int c_id = 0; c_id < nc; ++c_id)
	{
		int id = omp_get_thread_num();

		const int* cv_id = cell_vertex[c_id].data();
		int v_id_0 = cv_id[0]; int v_id_1 = cv_id[1]; int v_id_2 = cv_id[2]; int v_id_3 = cv_id[3];
		double x0 = dpx[v_id_0]; double y0 = dpy[v_id_0]; double z0 = dpz[v_id_0];

		VS[id](0, 0) = dpx[v_id_1] - x0; VS[id](1, 0) = dpy[v_id_1] - y0; VS[id](2, 0) = dpz[v_id_1] - z0;
		VS[id](0, 1) = dpx[v_id_2] - x0; VS[id](1, 1) = dpy[v_id_2] - y0; VS[id](2, 1) = dpz[v_id_2] - z0;
		VS[id](0, 2) = dpx[v_id_3] - x0; VS[id](1, 2) = dpy[v_id_3] - y0; VS[id](2, 2) = dpz[v_id_3] - z0;
		A[id] = VS[id] * cell_S[c_id];

		A2[id] = A[id].transpose()*A[id];
		double AF = A2[id](0, 0) + A2[id](1, 1) + A2[id](2, 2);
		double AF_05 = std::sqrt(AF);
		double AF2 = A2[id].squaredNorm();
		double AF_I = (AF*AF - AF2)*0.5; double AF_I_05 = std::sqrt(AF_I);
		double det_A = A[id].determinant();
		double i_det_A = 1.0 / det_A;

		double g = AF_05 * AF_I_05;
		double e = g*i_det_A;

		double k = ((e*e - 1.0)*0.125 + 0.5*(det_A + i_det_A)) * 0.5;
		if (k > 60) k = 60;
		e_sum[id] += std::exp(k);
	}

	for (int i = 1; i < max_num_t; ++i)
	{
		e_sum[0] += e_sum[i];
	}
	e_sum[0] /= nc;

	//flatten energy
	double e_flatten = 0.0; double count_flatten = 0.0;
	OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
	for (int i = 0; i < nc; ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[i];
		if (fv_id[0] >= 0)
		{
			p0[0] = dpx[fv_id[0]]; p0[1] = dpy[fv_id[0]]; p0[2] = dpz[fv_id[0]];
			p1[0] = dpx[fv_id[1]]; p1[1] = dpy[fv_id[1]]; p1[2] = dpz[fv_id[1]];
			p2[0] = dpx[fv_id[2]]; p2[1] = dpy[fv_id[2]]; p2[2] = dpz[fv_id[2]];
			OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0);
			double nx2 = n[0] * n[0]; double ny2 = n[1] * n[1]; double nz2 = n[2] * n[2];
			double inv_len_2 = 1.0 / (nx2 + ny2 + nz2);
			double cx = std::sqrt(inv_len_2 *nx2 + 0.001); double cy = std::sqrt(inv_len_2 *ny2 + 0.001); double cz = std::sqrt(inv_len_2 *nz2 + 0.001);
			e_flatten += cx + cy + cz;
			count_flatten += 1.0;
		}
	}
	e_flatten /= count_flatten;

	double s = e_sum[0] / (e_flatten + 1e-10);
	if (s > 1e20) return 1e20;
	return s;
}

void polycube_deformation_interface::set_vertex_color(int n_color, const std::vector<int>& v_color, VolumeMesh * mesh_)
{
	number_of_color = n_color; int nv = mesh_->n_vertices();
	vertex_diff_color.clear();
	if (number_of_color == 0)
	{
		vertex_diff_color.resize(nv);
	}
	else
	{
		vertex_diff_color.resize(number_of_color);
	}
	for (unsigned i = 0; i < nv;++i)
	{
		if (number_of_color == 0)
		{
			vertex_diff_color[i].push_back(i);
		}
		else
		{
			vertex_diff_color[v_color[i] - 1].push_back(i);
		}
	}
}

#include "..\Tetrahedrization\permutohedral.h"
void polycube_deformation_interface::filter_boundary_normal_RGNF(VolumeMesh * mesh_, double sigma_s_, double sigma_r_, int iter_count_, bool xyz_flag)
{
	//avg_boundary_edge_length = 0.0;
	//std::vector<double> face_center(3 * boundary_face_number);
	std::vector<double> six_pos(6 * boundary_face_number, 0.0);
	std::vector<double> N0(3 * boundary_face_number); //face area weight normal
	int count = 0; OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
	for (int i = 0; i < bfv_id.size(); ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
		if (one_bfv_id[0] < 0) continue;
		
		p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
		p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
		p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];
		/*OpenVolumeMesh::Geometry::Vec3d p0 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[0]));
		OpenVolumeMesh::Geometry::Vec3d p1 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[1]));
		OpenVolumeMesh::Geometry::Vec3d p2 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[2]));*/

		//avg_boundary_edge_length += (p0 - p1).norm();
		//avg_boundary_edge_length += (p1 - p2).norm();
		//avg_boundary_edge_length += (p2 - p0).norm();

		OpenVolumeMesh::Geometry::Vec3d pt = (p0 + p1 + p2) / 3.0;
		OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0)*0.5; //face weight normal
		//n.normalize();
		if (xyz_flag) //change normal to align the axis
		{
			double fa = n.norm()*0.5;
			//fa = 1;
			double max_component = std::abs(n[0]); int max_id = 0;
			OpenVolumeMesh::Geometry::Vec3d axis(0.0, 0.0, 0.0);
			for (int j = 1; j<3; j++)
			{
				if (std::abs(n[j])>max_component)
				{
					max_component = std::abs(n[j]);
					max_id = j;
				}
			}
			axis[max_id] = n[max_id] / max_component;
			n = axis * fa;
		}

		six_pos[6 * count + 0] = pt[0]; six_pos[6 * count + 1] = pt[1]; six_pos[6 * count + 2] = pt[2];

		N0[3 * count + 0] = n[0]; N0[3 * count + 1] = n[1]; N0[3 * count + 2] = n[2];

		++count;
	}
	//avg_boundary_edge_length /= (3.0*count);
	
	double inv_pos = 1.0 / (sigma_s_ *avg_boundary_edge_length);
	double inv_nor = 1.0 / sigma_r_;

	for (int i = 0; i < count; ++i)
	{
		six_pos[6 * i + 0] *= inv_pos; six_pos[6 * i + 1] *= inv_pos; six_pos[6 * i + 2] *= inv_pos;
	}

	//
	std::vector<double> N1(3 * boundary_face_number, 0.0);
	for (int i = 0; i < iter_count_; i++)
	{
		/*if (xyz_flag)
		{
			for (int i = 0; i < count; ++i)
			{
				six_pos[6 * i + 3] = N0[3 * i + 0] * inv_nor;
				six_pos[6 * i + 4] = N0[3 * i + 1] * inv_nor;
				six_pos[6 * i + 5] = N0[3 * i + 2] * inv_nor;
			}
		}*/
		PermutohedralLattice::filter(six_pos.data(), 6, N0.data(), 3, boundary_face_number, N1.data());

		for (int j = 0; j < boundary_face_number; j++)
		{
			OpenVolumeMesh::Geometry::Vec3d tn(N1[3 * j + 0], N1[3 * j + 1], N1[3 * j + 2]);
			double d = tn.norm();
			if (d < 1.0e-16)
			{
				OpenVolumeMesh::Geometry::Vec3d tn2(N0[3 * j + 0], N0[3 * j + 1], N0[3 * j + 2]);
				tn2.normalize();
				N1[3 * j + 0] = tn2[0]; N1[3 * j + 1] = tn2[1]; N1[3 * j + 2] = tn2[2];
				
				//printf("Small Normal.\n");
			}
			else
			{
				N1[3 * j + 0] = tn[0] / d; N1[3 * j + 1] = tn[1] / d; N1[3 * j + 2] = tn[2] / d;
			}

			six_pos[6 * j + 3] = N1[3 * j + 0] * inv_nor;
			six_pos[6 * j + 4] = N1[3 * j + 1] * inv_nor;
			six_pos[6 * j + 5] = N1[3 * j + 2] * inv_nor;
		}
	}

	count = 0;
	if (target_bfn.size() != bfv_id.size())
	{
		target_bfn.resize(bfv_id.size());
	}
	for (int i = 0; i < bfv_id.size(); ++i)
	{
		if (bfv_id[i][0] < 0) continue;
		/*if (change_big_flag[i] == 1)
		{
			double max_component = std::abs(N0[count * 3 + 0]); int max_id = 0;
			OpenVolumeMesh::Geometry::Vec3d axis(0.0, 0.0, 0.0);
			for (int j = 1; j<3; j++)
			{
				if (std::abs(N0[count * 3 + j])>max_component)
				{
					max_component = std::abs(N0[count * 3 + j]);
					max_id = j;
				}
			}
			axis[max_id] = N0[count * 3 + max_id] / max_component;
			target_bfn[i] = axis;
		}
		else*/
		{
			double max_component = std::abs(N1[count * 3 + 0]); int max_id = 0;
			OpenVolumeMesh::Geometry::Vec3d axis(0.0, 0.0, 0.0);
			for (int j = 1; j<3; j++)
			{
				if (std::abs(N1[count * 3 + j])>max_component)
				{
					max_component = std::abs(N1[count * 3 + j]);
					max_id = j;
				}
			}
			axis[max_id] = N1[count * 3 + max_id] / max_component;
			target_bfn[i] = axis;
		}
		/*else
		{
			target_bfn[i][0] = N1[count * 3 + 0];
			target_bfn[i][1] = N1[count * 3 + 1];
			target_bfn[i][2] = N1[count * 3 + 2];
		}*/
		++count;
	}
}

void polycube_deformation_interface::filter_boundary_normal_laplace(VolumeMesh * mesh_, int iter_count_)
{
	if (iter_count_ == 0) return;
	int nv = mesh_->n_vertices();
	std::vector<OpenVolumeMesh::Geometry::Vec3d> l_p(nv);
	double slen = avg_boundary_edge_length*avg_boundary_edge_length*0.6*0.6;
	for (int i = 0; i < bvv_id.size(); ++i)
	{
		std::vector<int>& one_bvv_id = bvv_id[i];
		if (one_bvv_id.size() == 0) continue;

		OpenVolumeMesh::Geometry::Vec3d fp(dpx[i], dpy[i], dpz[i]);
		OpenVolumeMesh::Geometry::Vec3d& fp_0 = l_p[i];
		fp_0 = fp;
		double w_sum = 1.0;
		for (int j = 0; j < one_bvv_id.size(); ++j)
		{
			OpenVolumeMesh::Geometry::Vec3d fq(dpx[one_bvv_id[j]], dpy[one_bvv_id[j]], dpz[one_bvv_id[j]]);
			double w = std::exp(-(fp - fq).sqrnorm() / slen); w_sum += w;
			fp_0 += w * fq;
		}
		fp_0 /= w_sum;
	}

	for (int k = 0; k < iter_count_ - 1; ++k)
	{
		for (int i = 0; i < bvv_id.size(); ++i)
		{
			std::vector<int>& one_bvv_id = bvv_id[i];
			if (one_bvv_id.size() == 0) continue;

			OpenVolumeMesh::Geometry::Vec3d& fp_0 = l_p[i];
			OpenVolumeMesh::Geometry::Vec3d fp = l_p[i];
			double w_sum = 1.0;
			//fp[0] = dpx[one_bvv_id[0]]; fp[1] = dpy[one_bvv_id[0]]; fp[2] = dpz[one_bvv_id[0]];
			for (int j = 0; j < one_bvv_id.size(); ++j)
			{
				double w = std::exp(-(fp - l_p[one_bvv_id[j]]).sqrnorm() / slen); w_sum += w;
				//fp[0] += dpx[one_bvv_id[j]]; fp[1] += dpy[one_bvv_id[j]]; fp[2] += dpz[one_bvv_id[j]];
				fp_0 += w*l_p[one_bvv_id[j]];
			}
			fp_0 /= w_sum;
		}
	}

	//normal
	if (target_bfn.size() != bfv_id.size())
	{
		target_bfn.resize(bfv_id.size());
	}
	OpenVolumeMesh::Geometry::Vec3d p0, p1, p2, N0, N1;
	for (int i = 0; i < bfv_id.size(); ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
		if (one_bfv_id[0] < 0) continue;

		if (change_big_flag[i] == 1)
		{
			p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
			p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
			p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];

			N0 = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0).normalize();

			double max_component = std::abs(N0[0]); int max_id = 0;
			OpenVolumeMesh::Geometry::Vec3d axis(0.0, 0.0, 0.0);
			for (int j = 1; j<3; j++)
			{
				if (std::abs(N0[j]) > max_component)
				{
					max_component = std::abs(N0[j]);
					max_id = j;
				}
			}
			axis[max_id] = N0[max_id] / max_component;
			target_bfn[i] = axis;
		}
		else
		{
			p0 = l_p[one_bfv_id[0]]; p1 = l_p[one_bfv_id[1]]; p2 = l_p[one_bfv_id[2]];
			N1 = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0).normalize();

			double max_component = std::abs(N1[0]); int max_id = 0;
			OpenVolumeMesh::Geometry::Vec3d axis(0.0, 0.0, 0.0);
			for (int j = 1; j<3; j++)
			{
				if (std::abs(N1[j])>max_component)
				{
					max_component = std::abs(N1[j]);
					max_id = j;
				}
			}
			axis[max_id] = N1[max_id] / max_component;
			target_bfn[i] = axis;
		}
	}
}

void polycube_deformation_interface::filter_boundary_normal_ring(VolumeMesh * mesh_, double sigma_s_, int iter_count_, bool xyz_flag /* = false */)
{
	avg_boundary_edge_length = 0.0;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> N0(boundary_face_number);//face area weight normal
	std::vector<OpenVolumeMesh::Geometry::Vec3d> fc(boundary_face_number);
	int count = 0; int nc = bfv_id.size(); OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
	std::vector<int> map_face(nc, -1);
	for (int i = 0; i < nc; ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
		if (one_bfv_id[0] < 0) continue;

		p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
		p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
		p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];

		/*OpenVolumeMesh::Geometry::Vec3d p0 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[0]));
		OpenVolumeMesh::Geometry::Vec3d p1 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[1]));
		OpenVolumeMesh::Geometry::Vec3d p2 = mesh_->vertex(OpenVolumeMesh::VertexHandle(one_bfv_id[2]));*/

		avg_boundary_edge_length += (p0 - p1).norm();
		avg_boundary_edge_length += (p1 - p2).norm();
		avg_boundary_edge_length += (p2 - p0).norm();

		OpenVolumeMesh::Geometry::Vec3d pt = (p0 + p1 + p2) / 3.0;
		OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0)*0.5; //face weight normal
		//n.normalize();
		if (xyz_flag) //change normal to align the axis
		{
			double fa = n.norm()*0.5;
			//fa = 1;
			double max_component = std::abs(n[0]); int max_id = 0;
			OpenVolumeMesh::Geometry::Vec3d axis(0.0, 0.0, 0.0);
			for (int j = 1; j<3; j++)
			{
				if (std::abs(n[j])>max_component)
				{
					max_component = std::abs(n[j]);
					max_id = j;
				}
			}
			axis[max_id] = n[max_id] / max_component;
			n = axis * fa;
		}

		fc[count] = pt; N0[count] = n;

		map_face[i] = count;
		++count;
	}
	avg_boundary_edge_length /= (3.0*count);
	double inv_pos = 1.0 / (sigma_s_ *avg_boundary_edge_length);
	double inv_pos2 = inv_pos*inv_pos*0.5;

	//
	std::vector<int> visited_f(nc, -1); std::vector<int> ring_face; ring_face.reserve(25);
	std::vector<OpenVolumeMesh::Geometry::Vec3d> N1(boundary_face_number);
	for (int i = 0; i < iter_count_; i++)
	{
		for (int f_id = 0; f_id < nc; ++f_id) //for each boundary face
		{
			OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[f_id];
			if (one_bfv_id[0] < 0) continue;

#if 1
			for (int j = 0; j < 3; ++j)
			{
				std::vector<int>& one_bvf_id = bvf_id[one_bfv_id[j]];
				for (int k = 0; k < one_bvf_id.size(); ++k)
				{
					if (visited_f[one_bvf_id[k]] < 0)
					{
						visited_f[one_bvf_id[k]] = 1;
						ring_face.push_back(one_bvf_id[k]);
					}
				}
			}
#else
			for (int j = 0; j < 3; ++j)
			{
				std::vector<int>& one_bvf_id = bvf_id[one_bfv_id[j]];
				for (int k = 0; k < one_bvf_id.size(); ++k)
				{
					if (visited_f[one_bvf_id[k]] < 0)
					{
						visited_f[one_bvf_id[k]] = 1;
						ring_face.push_back(one_bvf_id[k]);
					}
				}
			}
			for (int ijk = 0; ijk < ring_face.size(); ++ijk)
			{
				OpenVolumeMesh::Geometry::Vec3i& one_bfv_id_ = bfv_id[ring_face[ijk]];
				for (int j = 0; j < 3; ++j)
				{
					std::vector<int>& one_bvf_id = bvf_id[one_bfv_id_[j]];
					for (int k = 0; k < one_bvf_id.size(); ++k)
					{
						if (visited_f[one_bvf_id[k]] < 0)
						{
							visited_f[one_bvf_id[k]] = 1;
							ring_face.push_back(one_bvf_id[k]);
						}
					}
				}
			}
#endif
			
			OpenVolumeMesh::Geometry::Vec3d& fc_0 = fc[map_face[f_id]];
			OpenVolumeMesh::Geometry::Vec3d n(0,0,0);
			for (int j = 0; j < ring_face.size(); ++j)
			{
				double len2 = (fc[map_face[ring_face[j]]] - fc_0).sqrnorm();
				double w = std::exp(-len2*inv_pos2);
				n += w * N0[map_face[ring_face[j]]];
			}

			double d = n.norm();
			if (d < 1.0e-16)
			{
				N1[map_face[f_id]] = N0[map_face[f_id]];
				N1[map_face[f_id]].normalize();
			}
			else
			{
				N1[map_face[f_id]] = n / d;
			}

			//printf("%d %d %d\n", f_id, map_face[f_id], ring_face.size());
			for (int j = 0; j < ring_face.size(); ++j)
			{
				visited_f[ring_face[j]] = -1;
			}
			ring_face.clear(); ring_face.reserve(25);
		}
	}

	count = 0;
	if (target_bfn.size() != bfv_id.size())
	{
		target_bfn.resize(bfv_id.size());
	}
	for (int i = 0; i < bfv_id.size(); ++i)
	{
		if (bfv_id[i][0] < 0) continue;

		/*if (xyz_flag)
		{
			double max_component = std::abs(N1[count][0]); int max_id = 0;
			OpenVolumeMesh::Geometry::Vec3d axis(0.0, 0.0, 0.0);
			for (int j = 1; j<3; j++)
			{
				if (std::abs(N1[count][j])>max_component)
				{
					max_component = std::abs(N1[count][j]);
					max_id = j;
				}
			}
			axis[max_id] = N1[count][max_id] / max_component;
			target_bfn[i] = axis;
			double cos_v = OpenVolumeMesh::Geometry::dot(axis, N1[count]);
			if ( cos_v < std::cos(5) )
			{
				printf("1 %f %f %f\2 %f %f %f\n", axis[0], axis[1], axis[2], N1[count][0], N1[count][1], N1[count][2]);
			}
		}
		else*/
		{
			target_bfn[i] = N1[count];
		}

		++count;
	}
}

void polycube_deformation_interface::exp_mips_deformation_refine_polycube_omp(int max_iter, int iter2, double energy_power, double angle_area_ratio, VolumeMesh * mesh_, bool use_xyz, bool use_RGNF, bool use_half_sigma_s)
{
	double ga = 1; double lamda = 0.1;
	//printf("%f %f %d\n", sigma_s, sigma_r, RGNF_iter_count);
	double temp_ss = sigma_s;
	if (use_half_sigma_s) sigma_s *= 0.5;

	for (int j = 0; j < max_iter; ++j)
	{
		if (use_RGNF) filter_boundary_normal_RGNF(mesh_, sigma_s, sigma_r, RGNF_iter_count, use_xyz);
		else filter_boundary_normal_ring(mesh_, 2 * sigma_s, RGNF_iter_count, true);

		//filter_boundary_normal_laplace(mesh_, 5);

		use_half_sigma_s = true;
		/*if (j == 0) use_half_sigma_s = false; 
		else use_half_sigma_s = true;*/

		ga = compute_energy_ratio_filter_deform();
		for (int k = 0; k < iter2; ++k)
		{
			//printf("%d ", k);
			int n_color = mesh_->n_vertices();
			if (number_of_color != 0) n_color = number_of_color;

			for (unsigned i = 0; i < n_color; ++i)
			{
				int one_color_size = vertex_diff_color[i].size();
#pragma omp parallel for schedule(dynamic, 1)
				for (int j = 0; j < one_color_size; ++j)
				{
					int id = omp_get_thread_num();
					exp_mips_deformation_refine_one_polycube(vertex_diff_color[i][j], angle_area_ratio, energy_power, id, ga, use_half_sigma_s);
					//exp_mips_deformation_refine_one_polycube_normal(vertex_diff_color[i][j], angle_area_ratio, energy_power, id, mu1, mu2);
				}
			}
		}
		//compute_energy_ratio_filter_deform();

		//lamda *= 2;
		if (bound_K < 4.0)
		{
			bound_K += 0.1;
		}
	}

	//change_big_flag.clear(); change_big_flag.resize(mesh_->n_cells(), -1);

	sigma_s = temp_ss;
}

void polycube_deformation_interface::compute_face_normal(std::vector< OpenVolumeMesh::Geometry::Vec3d >& bfn)
{
	int nc = bfv_id.size();
	OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
	for (int i = 0; i < nc; ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[i];
		if (fv_id[0] >= 0)
		{
			p0[0] = dpx[fv_id[0]]; p0[1] = dpy[fv_id[0]]; p0[2] = dpz[fv_id[0]];
			p1[0] = dpx[fv_id[1]]; p1[1] = dpy[fv_id[1]]; p1[2] = dpz[fv_id[1]];
			p2[0] = dpx[fv_id[2]]; p2[1] = dpy[fv_id[2]]; p2[2] = dpz[fv_id[2]];

			bfn[i] = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0).normalize();
		}
	}
}

bool polycube_deformation_interface::check_normal_difference(double th)
{
	double cos_th = std::cos(th);int nc = bfv_id.size();
	OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
	if( change_big_flag.size() != nc) change_big_flag.resize(nc, -1);
	for (int i = 0; i < nc; ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& fv_id = bfv_id[i];
		change_big_flag[i] = -1;
		if (fv_id[0] >= 0)
		{
			OpenVolumeMesh::Geometry::Vec3d& ln = last_bfn[i];
			p0[0] = dpx[fv_id[0]]; p0[1] = dpy[fv_id[0]]; p0[2] = dpz[fv_id[0]];
			p1[0] = dpx[fv_id[1]]; p1[1] = dpy[fv_id[1]]; p1[2] = dpz[fv_id[1]];
			p2[0] = dpx[fv_id[2]]; p2[1] = dpy[fv_id[2]]; p2[2] = dpz[fv_id[2]];

			OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0).normalize();
			double cos_ = OpenVolumeMesh::Geometry::dot(n, ln);
			if (cos_ < cos_th) change_big_flag[i] = 1;

			ln = n;
		}
	}

	return true;
}

void polycube_deformation_interface::exp_mips_deformation_refine_polycube_normal_omp(int max_iter, double energy_power, double angle_area_ratio, VolumeMesh * mesh_)
{
	compute_face_normal(last_bfn);
	double mu1 = 1e4; double mu2 = 0.0;
	//double ga = 1e6;
	for (int iter = 0; iter < 4; ++iter)
	{
		//filter_boundary_normal_RGNF(mesh_, 1.0, sigma_r, 4, true);
		for (int k = 0; k < max_iter ; ++k)
		{
			//printf("%d ", k);
			int n_color = mesh_->n_vertices();
			if (number_of_color != 0) n_color = number_of_color;

			for (unsigned i = 0; i < n_color; ++i)
			{
				int one_color_size = vertex_diff_color[i].size();
#pragma omp parallel for schedule(dynamic, 1)
				for (int j = 0; j < one_color_size; ++j)
				{
					int id = omp_get_thread_num();
					//exp_mips_deformation_refine_one_polycube(vertex_diff_color[i][j], angle_area_ratio, energy_power, id, ga);
					exp_mips_deformation_refine_one_polycube_normal(vertex_diff_color[i][j], angle_area_ratio, energy_power, id, mu1, mu2);
				}
			}
		}
		//printf("\n");
		mu1 *= 5;
		//ga *= 5;
	}
	check_normal_difference(45);
}

void polycube_deformation_interface::flatten_boundary_face_omp(int max_iter, double energy_power, double angle_area_ratio, VolumeMesh * mesh_)
{
	double ga = 1.0;
	for (int iter3 = 0; iter3 < 1; ++iter3)
	{
		for (int k = 0; k < max_iter; ++k)
		{
			printf("%d ", k);
			int n_color = mesh_->n_vertices();
			if (number_of_color != 0) n_color = number_of_color;

			for (unsigned i = 0; i < n_color; ++i)
			{
				int one_color_size = vertex_diff_color[i].size();
#pragma omp parallel for schedule(dynamic, 1)
				for (int j = 0; j < one_color_size; ++j)
				{
					int id = omp_get_thread_num();
					exp_mips_deformation_refine_one_polycube(vertex_diff_color[i][j], angle_area_ratio, energy_power, id, ga);
				}
			}
		}
		ga *= 2;
		printf("\n");
	}
}

void polycube_deformation_interface::flatten_boundary_face_omp_normal_position(int max_iter, double energy_power, double angle_area_ratio, VolumeMesh * mesh_)
{
	/*double ga = 1;
	find_all_corner(mesh_);
	for (int iter3 = 0; iter3 < 1; ++iter3)
	{
		for (int k = 0; k < max_iter; ++k)
		{
			//printf("%d ", k);
			int n_color = mesh_->n_vertices();
			if (number_of_color != 0) n_color = number_of_color;

			for (unsigned i = 0; i < n_color; ++i)
			{
				int one_color_size = vertex_diff_color[i].size();
#pragma omp parallel for schedule(dynamic, 1)
				for (int j = 0; j < one_color_size; ++j)
				{
					int id = omp_get_thread_num();
					exp_mips_deformation_refine_one_polycube_position(vertex_diff_color[i][j], angle_area_ratio, energy_power, id, ga);
				}
			}
		}
		ga *= 2;
	}*/
	//ga *= 2;
	//printf("\n");
}

void polycube_deformation_interface::find_all_corner(VolumeMesh * mesh_)
{
	//compute the mean value for each chart
	int chart_size = polycube_chart_label.size(); double cube_len = sigma_r * avg_boundary_edge_length;
	printf("%f\n", cube_len);
	std::vector<int> with_big_diff_chart; std::vector<double> max_chart; std::vector<int> min_chart;
	for (int i = 0; i < chart_size; ++i)
	{
		int chart_label = polycube_chart_label[i]/2; double sum = 0.0;
		std::vector<int>& one_chart = polycube_chart[i]; double max_value = -1e30; double min_value = 1e30;
		for (int j = 0; j < one_chart.size(); ++j)
		{
			int bf_id = one_chart[j]; double one_face_sum = 0.0;
			for (int k = 0; k < 3; ++k)
			{
				one_face_sum += mesh_->vertex(OpenVolumeMesh::VertexHandle(bfv_id[bf_id][k]))[chart_label];
			}
			double xyz_value = one_face_sum / 3.0;
			sum += xyz_value;
			if (max_value < xyz_value) max_value = xyz_value;
			if (min_value > xyz_value) min_value = xyz_value;
		}
		//chart_mean_value[i] = sum/one_chart.size();
		double max_diff = max_value - min_value; 
		chart_mean_value[i] = 0.5*(max_value+min_value);

#ifdef Hex_Meshing
		double int_v = std::floor(chart_mean_value[i] / cube_len + 0.5);
		chart_mean_value[i] = int_v * cube_len;
#endif
		//chart_mean_value[i] = min_value;
		if (max_diff > 0.001)
		{
			with_big_diff_chart.push_back(i); max_chart.push_back(max_value); min_chart.push_back(min_value);
		}
	}
	
	/*for (int i = 0; i < with_big_diff_chart.size(); ++i)
	{

	}*/
}


void polycube_deformation_interface::find_all_chart_value(VolumeMesh * mesh_)
{
	//get the initial chart value as the objectives
	find_all_corner(mesh_);

	/* Bounds on constraints. */
	std::vector<MSKboundkeye> bkc; std::vector<double> blc; std::vector<double> buc;
	/* Bounds on variables. */
	std::vector<MSKboundkeye> bkx; std::vector<double> blx; std::vector<double> bux;
	/* sparse representation of the A matrix stored by column. */
	std::vector<MSKlidxt> aptrb; std::vector<MSKidxt> asub; std::vector<double> aval;
	//objectives
	std::vector<MSKidxt> qsubi; std::vector<MSKidxt> qsubj; std::vector<double> qval;
	
	//constraints
	double ratio_diff = 0.8; double cube_len = sigma_r * avg_boundary_edge_length;
	std::vector< Eigen::Triplet<double> > coef;
	int up_down_chart_count = up_chart_id.size();
	int chart_number = chart_mean_value.size();
	double diff_u_d; int up_chart, down_chart;
	for (int i = 0; i < up_down_chart_count; ++i)
	{
		up_chart = up_chart_id[i];
		down_chart = down_chart_id[i];
		diff_u_d = diff_up_down[i] * ratio_diff;

#ifdef Hex_Meshing
		if (diff_u_d < cube_len)
		{
			diff_u_d = cube_len;
		}
#endif
		coef.push_back(Eigen::Triplet<double>(i, up_chart, +1.0));
		coef.push_back(Eigen::Triplet<double>(i, down_chart, -1.0));
		bkc.push_back(MSK_BK_LO); blc.push_back(diff_u_d); buc.push_back(+MSK_INFINITY);
	}
	Eigen::SparseMatrix<double> A(up_down_chart_count, chart_number);
	A.setFromTriplets(coef.begin(), coef.end());
	int size = A.nonZeros(); int cols = A.cols();
	/* Below is the sparse representation of the A matrix stored by column. */
	aptrb.resize(cols + 1); asub.resize(size); aval.resize(size);
	aptrb[cols] = size;
	int ind = 0;
	for (int k = 0; k < cols; ++k)
	{
		aptrb[k] = ind;
		for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
		{
			aval[ind] = it.value();
			asub[ind] = it.row();
			++ind;
		}
	}
	//objectives
	std::vector<double> c(chart_number, 0.0);
	for (int i = 0; i < chart_number; ++i)
	{
		qsubi.push_back(i); qsubj.push_back(i); qval.push_back(2.0);
		c[i] = -2.0 * chart_mean_value[i];
	}

	//variables
	bkx.resize(chart_number, MSK_BK_FR); blx.resize(chart_number, -MSK_INFINITY); bux.resize(chart_number, +MSK_INFINITY);

	//using Mosek to solve the problem
	std::vector<double> x(chart_number, 0.0); int mosek_count = 0;
	bool solve_success = solveConvexQuadPorgramming_mosek(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
	while ( !solve_success)
	{
		ratio_diff *= 0.8;
		for (int i = 0; i < up_down_chart_count; ++i)
		{
			diff_u_d = diff_up_down[i] * ratio_diff;
#ifdef Hex_Meshing
			if (diff_u_d < cube_len)
			{
				diff_u_d = cube_len;
			}
#endif
			blc[i] = diff_u_d;
		}
		solve_success = solveConvexQuadPorgramming_mosek(bkc, blc, buc, bkx, blx, bux, aptrb, asub, aval, qsubi, qsubj, qval, c, x);
		++mosek_count;
		if (mosek_count > 20) break;
	}
	if (mosek_count > 20)
	{
		printf("-------------------------------------------\nNo Solution!!!!!\n");
	}
	else
	{
		for (int i = 0; i < chart_number; ++i)
		{
#ifdef Hex_Meshing
			chart_mean_value[i] = cube_len*std::floor(x[i] / cube_len + 0.5);
#else
			chart_mean_value[i] = x[i];
#endif
			
		}
	}
}

#include "..\Tetrahedrization\adjust_orientation_LBFGS.h"
void polycube_deformation_interface::adjust_orientation(VolumeMesh* mesh_)
{
	//return;

	OpenVolumeMesh::Geometry::Vec3d p0, p1, p2;
	std::vector<OpenVolumeMesh::Geometry::Vec3d> all_bfn;
	for (int i = 0; i < bfv_id.size(); ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& one_bfv_id = bfv_id[i];
		if (one_bfv_id[0] < 0) continue;

		p0[0] = dpx[one_bfv_id[0]]; p0[1] = dpy[one_bfv_id[0]]; p0[2] = dpz[one_bfv_id[0]];
		p1[0] = dpx[one_bfv_id[1]]; p1[1] = dpy[one_bfv_id[1]]; p1[2] = dpz[one_bfv_id[1]];
		p2[0] = dpx[one_bfv_id[2]]; p2[1] = dpy[one_bfv_id[2]]; p2[2] = dpz[one_bfv_id[2]];

		OpenVolumeMesh::Geometry::Vec3d pt = (p0 + p1 + p2) / 3.0;
		OpenVolumeMesh::Geometry::Vec3d n = OpenVolumeMesh::Geometry::cross(p1 - p0, p2 - p0); //face weight normal
		all_bfn.push_back( n.normalize() );
	}

	std::vector<double> X(3, 0.0);
	for (int i = 0; i < 1; ++i)
	{
		ad_LBFGS(all_bfn, X);
	}

	double alpha = X[0]; double beta = X[1]; double gamma = X[2];
	printf("%f %f %f\n", alpha, beta, gamma);
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

	OpenVolumeMesh::Geometry::Vec3d p;
	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		OpenVolumeMesh::Geometry::Vec3d q = mesh_->vertex(*v_it);
		p[0] = r00 * q[0] + r01 * q[1] + r02 * q[2];
		p[1] = r10 * q[0] + r11 * q[1] + r12 * q[2];
		p[2] = r20 * q[0] + r21 * q[1] + r22 * q[2];
		mesh_->set_vertex(*v_it, p);
		int v_id = v_it->idx();
		dpx[v_id] = p[0] /*+ avg_p[0]*/; dpy[v_id] = p[1]/* + avg_p[1]*/; dpz[v_id] = p[2]/* + avg_p[2]*/;
	}
}

#include "..\Tetrahedrization\Polycube_Boundary_Map_IVF.h"
void polycube_deformation_interface::boundary_mapping_polycube(VolumeMesh* mesh_)
{
	int nv = mesh_->n_vertices(); int nc = mesh_->n_cells();
	std::vector<int> map_v(nv, -1); std::vector<int> map_c(nc, -1);
	int chart_size = polycube_chart_label.size();
	for (int i = 0; i < chart_size; ++i) //for each chart
	{
		int var_v_count = 0; int var_c_count = 0;
		std::vector<int>& one_chart = polycube_chart[i];
		for (int j = 0; j < one_chart.size(); ++j) //for each face on chart
		{
			int bf_id = one_chart[j]; //
			for (int k = 0; k < 3; ++k)
			{
				int v_id = bfv_id[bf_id][k];
				bool fixed = ((is_polycube_handles[3 * v_id + 0] == 1) && (is_polycube_handles[3 * v_id + 1] == 1) && (is_polycube_handles[3 * v_id + 2] == 1));
				if (!fixed)
				{
					if (map_v[v_id] < 0)
					{
						map_v[v_id] = var_v_count; ++var_v_count;
					}
					if (map_c[bf_id] < 0)
					{
						map_c[bf_id] = var_c_count; ++var_c_count;
					}
				}
			}
		}

		//mapping by 2D IVF mapping
		int chart_label = polycube_chart_label[i];
		polycube_boundary_mapping_IVF pbm;
		pbm.optimize_IVF(mesh_, chart_label, src_pos, bfv_id, bvf_id, map_v, var_v_count, map_c, var_c_count);

		//clear map
		for (int j = 0; j < nv; ++j)
		{
			if (map_v[j] >= 0)
			{
				is_polycube_handles[3 * j + 0] = 1; is_polycube_handles[3 * j + 1] = 1; is_polycube_handles[3 * j + 2] = 1;
			}
			map_v[j] = -1;
		}
		for (int j = 0; j < nc; ++j)
		{
			map_c[j] = -1;
		}
	}
}

void polycube_deformation_interface::deform_ARAP_polycube(VolumeMesh* mesh_)
{
	//find_all_corner(mesh_);
	find_all_chart_value(mesh_);
	int nv = vertex_type.size();
	is_polycube_handles.clear(); is_polycube_handles.resize(3 * nv, -1);
	//for vertex on polycube edge
	std::vector<int> polycube_edge_v; std::vector<int> polycube_edge_v_neighbor;
	std::vector<int> polycube_edge_v_type;
	for (int i = 0; i < nv; ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[i];
		OpenVolumeMesh::VertexHandle vh = OpenVolumeMesh::VertexHandle(i);
		OpenVolumeMesh::Geometry::Vec3d np = mesh_->vertex(vh);
		if (v_type[0] >= 0)
		{
			np[0] = chart_mean_value[v_type[0]];
			dpx[i] = np[0];
			is_polycube_handles[3 * i + 0] = 1;
		}
		if (v_type[1] >= 0)
		{
			np[1] = chart_mean_value[v_type[1]];
			dpy[i] = np[1];
			is_polycube_handles[3 * i + 1] = 1;
		}
		if (v_type[2] >= 0)
		{
			np[2] = chart_mean_value[v_type[2]];
			dpz[i] = np[2];
			is_polycube_handles[3 * i + 2] = 1;
		}
		mesh_->set_vertex(vh, np);

		if (v_type[0] >= 0 && v_type[1] >= 0 && v_type[2] < 0)// x-y
		{
			polycube_edge_v.push_back(i);
			polycube_edge_v_type.push_back(2);
			is_polycube_handles[3 * i + 2] = 1;
		}
		else if (v_type[0] < 0 && v_type[1] >= 0 && v_type[2] >= 0)// y-z
		{
			polycube_edge_v.push_back(i);
			polycube_edge_v_type.push_back(0);
			is_polycube_handles[3 * i + 0] = 1;
		}
		else if (v_type[0] >= 0 && v_type[1] < 0 && v_type[2] >= 0)// x-z
		{
			polycube_edge_v.push_back(i);
			polycube_edge_v_type.push_back(1);
			is_polycube_handles[3 * i + 1] = 1;
		}
	}

	for (int i = 0; i < polycube_edge_v.size(); ++i)
	{
		int v_id = polycube_edge_v[i]; OpenVolumeMesh::Geometry::Vec3d np = mesh_->vertex(OpenVolumeMesh::VertexHandle(v_id));
		std::vector<int>& one_vv = bvv_id[v_id];
		int vv_size = one_vv.size(); int add_nv = 0;
		for (int j = 0; j < vv_size; ++j)
		{
			int vj_id = one_vv[j];
			OpenVolumeMesh::Geometry::Vec3i& vj_type = vertex_type[vj_id]; int chart_count = 0; bool can_add = false;
			if (vj_type[0] >= 0) ++chart_count;
			if (vj_type[1] >= 0) ++chart_count;
			if (vj_type[2] >= 0) ++chart_count;
			if (chart_count == 3)// 3
			{
				//polycube_edge_v_neighbor.push_back(vj_id);
				//printf("%d %d %d\n", i, chart_count, vj_id);
				can_add = true;
			}
			else if (chart_count == 2)
			{
				if (vj_type[polycube_edge_v_type[i]] < 0)
				{
					can_add = true;
					//printf("%d %d %d\n", i, chart_count, vj_id);
				}
			}
			if (can_add)
			{
				OpenVolumeMesh::Geometry::Vec3d npj = mesh_->vertex(OpenVolumeMesh::VertexHandle(vj_id));
				bool nv_ok = true;
				for (int k = 0; k < 3; ++k)
				{
					if (k != polycube_edge_v_type[i])
					{
						if (std::abs(npj[k] - np[k]) > 1e-10) nv_ok = false;
					}
				}
				if (nv_ok)
				{
					polycube_edge_v_neighbor.push_back(vj_id);
					//printf("%d %d %d\n", i, chart_count, vj_id);
					++add_nv;
				}
			}
		}
		if (add_nv > 0 && add_nv != 2)
		{
			printf("%d %d\n", v_id, add_nv);
			if (add_nv > 2)
			{
				int pevn_size = polycube_edge_v_neighbor.size();
				std::vector<int> delete_v;
				std::vector<int> visited_f(mesh_->n_cells(), 0);
				std::vector<int>& one_vf = bvf_id[v_id];
				for (int j = 0; j < one_vf.size(); ++j)
				{
					visited_f[one_vf[j]] += 1;
				}
				for (int j = pevn_size - add_nv; j < pevn_size; ++j)
				{
					int vj_id = polycube_edge_v_neighbor[j];
					std::vector<int>& one_vj_f = bvf_id[vj_id];
					int two_f[2]; int two_count = 0;
					for (int k = 0; k < one_vj_f.size(); ++k)
					{
						if (visited_f[one_vj_f[k]] == 1)
						{
							two_f[two_count] = one_vj_f[k]; ++two_count;
						}
					}

					printf("%d %d\n", bf_chart[two_f[0]], bf_chart[two_f[1]]);
					if (bf_chart[two_f[0]] == bf_chart[two_f[1]])
					{
						delete_v.push_back(j);
					}
				}

				for (int j = 0; j < delete_v.size(); ++j)
				{
					polycube_edge_v_neighbor.erase(polycube_edge_v_neighbor.begin() + delete_v[delete_v.size() - 1 - j] );
				}
			}
			/*else if(add_nv == 1)
			{
				polycube_edge_v_neighbor.pop_back();
			}*/
		}
	}
	//return;
	//average deformation of the chart boundary
	printf("%d %d\n", polycube_edge_v.size(), polycube_edge_v_neighbor.size());
	if (polycube_edge_v_neighbor.size() == 2 * polycube_edge_v.size())
	{
		for (int i = 0; i < 25; ++i)
		{
			for (int j = 0; j < polycube_edge_v.size(); ++j)
			{
				int v_id = polycube_edge_v[j];
				OpenVolumeMesh::VertexHandle vh = OpenVolumeMesh::VertexHandle(v_id);
				OpenVolumeMesh::Geometry::Vec3d np = mesh_->vertex(vh);

				OpenVolumeMesh::VertexHandle vh0 = OpenVolumeMesh::VertexHandle(polycube_edge_v_neighbor[2 * j + 0]);
				OpenVolumeMesh::Geometry::Vec3d p0 = mesh_->vertex(vh0);
				OpenVolumeMesh::VertexHandle vh1 = OpenVolumeMesh::VertexHandle(polycube_edge_v_neighbor[2 * j + 1]);
				OpenVolumeMesh::Geometry::Vec3d p1 = mesh_->vertex(vh1);

				int edge_type = polycube_edge_v_type[j];
				np[edge_type] = 0.5*(p0[edge_type] + p1[edge_type]);
				mesh_->set_vertex(vh, np);
				dpx[v_id] = np[0]; dpy[v_id] = np[1]; dpz[v_id] = np[2];
			}
		}
	}

	//return;

	//map to boundary for each chart
	boundary_mapping_polycube(mesh_);

	//return;

	//integer for hex meshing

	//ARAP
	int nc = mesh_->n_cells(); Eigen::Matrix3d VS, Q, U, V, R;
	std::vector < std::vector < double >> Cell_R(nc);
	for (int c_id = 0; c_id < nc; ++c_id)
	{
		Cell_R[c_id].resize(9);
	}
	for (int iter = 0; iter < 3; ++iter)
	{
		//local step
		for (int c_id = 0; c_id < nc; ++c_id)
		{
			std::vector<int>& cv_id = cell_vertex[c_id];
			double cpx = dpx[cv_id[0]]; double cpy = dpy[cv_id[0]]; double cpz = dpz[cv_id[0]];
			VS(0, 0) = dpx[cv_id[1]] - cpx; VS(1, 0) = dpy[cv_id[1]] - cpy; VS(2, 0) = dpz[cv_id[1]] - cpz;
			VS(0, 1) = dpx[cv_id[2]] - cpx; VS(1, 1) = dpy[cv_id[2]] - cpy; VS(2, 1) = dpz[cv_id[2]] - cpz;
			VS(0, 2) = dpx[cv_id[3]] - cpx; VS(1, 2) = dpy[cv_id[3]] - cpy; VS(2, 2) = dpz[cv_id[3]] - cpz;
			Q = VS * cell_S[c_id];
			Eigen::JacobiSVD<Eigen::Matrix3d> svd(Q, Eigen::ComputeFullU | Eigen::ComputeFullV);
			U = svd.matrixU(); V = svd.matrixV();
			if (Q.determinant() > 0)
			{
				//R = U*V.transpose();
				R = V*U.transpose();
			}
			else
			{
				//V(0, 2) = -V(0, 2); V(1, 2) = -V(1, 2); V(2, 2) = -V(2, 2);
				//R = U*V.transpose();
				U(0, 2) = -U(0, 2); U(1, 2) = -U(1, 2); U(2, 2) = -U(2, 2);
				R = V*U.transpose();
			}
			Cell_R[c_id][0] = R(0, 0); Cell_R[c_id][1] = R(0, 1); Cell_R[c_id][2] = R(0, 2);
			Cell_R[c_id][3] = R(1, 0); Cell_R[c_id][4] = R(1, 1); Cell_R[c_id][5] = R(1, 2);
			Cell_R[c_id][6] = R(2, 0); Cell_R[c_id][7] = R(2, 1); Cell_R[c_id][8] = R(2, 2);
		}
		//global step
		for (int i = 0; i < 3; ++i)//for x-y-z
		{
			int var_count = 0; std::vector<int> map_v(nv, -1);
			for (int j = 0; j < nv; ++j)
			{
				if (is_polycube_handles[3 * j + i] != 1)
				{
					map_v[j] = var_count;
					++var_count;
				}
			}
			Sparse_Matrix A = new Sparse_Matrix(nc * 3, var_count, NOSYM, CCS, 1);
			std::vector<double> cs(4);
			for (int c_id = 0; c_id < nc; ++c_id)// for each cell
			{
				double cv = std::sqrt(std::abs(cell_volume[c_id])) * 100;
				std::vector<int>& cv_id = cell_vertex[c_id];
				Eigen::Matrix3d& CS = cell_S[c_id];
				std::vector<double>& CR = Cell_R[c_id];
				for (int j = 0; j < 3; ++j) //for three column
				{
					cs[1] = CS(0, j)*cv; cs[2] = CS(1, j)*cv; cs[3] = CS(2, j)*cv;
					cs[0] = -(cs[1] + cs[2] + cs[3]);
					for (int k = 0; k < 4; ++k)
					{
						int var_id = map_v[cv_id[k]];
						if (var_id < 0) //k
						{
							A.fill_rhs_entry(3 * c_id + j, -cs[k] * mesh_->vertex(OpenVolumeMesh::VertexHandle(cv_id[k]))[i]);
						}
						else
						{
							A.fill_entry(3 * c_id + j, var_id, cs[k]);
						}
					}
					A.fill_rhs_entry(3 * c_id + j, cv*CR[i * 3 + j]);
				}
			}

			Sparse_Matrix* D = TransposeTimesSelf(&A, CCS, SYM_LOWER, true);
			solve_by_CHOLMOD(D);
			const std::vector<double>& xyz = D->get_solution();

			for (int j = 0; j < nv; ++j)
			{
				if (is_polycube_handles[3 * j + i] != 1)
				{
					int var_id = map_v[j];
					OpenVolumeMesh::VertexHandle vh = OpenVolumeMesh::VertexHandle(j);
					OpenVolumeMesh::Geometry::Vec3d np = mesh_->vertex(vh);
					np[i] = xyz[var_id]; mesh_->set_vertex(vh, np);
				}
			}

			delete D;
		}

		for (int j = 0; j < nv; ++j)
		{
			OpenVolumeMesh::VertexHandle vh = OpenVolumeMesh::VertexHandle(j);
			OpenVolumeMesh::Geometry::Vec3d np = mesh_->vertex(vh);
			dpx[j] = np[0]; dpy[j] = np[1]; dpz[j] = np[2];
		}
	}

}

void polycube_deformation_interface::save_polycube_constraint(const char* filename, VolumeMesh* mesh_)
{
	FILE* f_pcde = fopen(filename, "w");
	int nv = vertex_type.size();
	for (int i = 0; i < nv; ++i)
	{
		OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[i];
		OpenVolumeMesh::VertexHandle vh = OpenVolumeMesh::VertexHandle(i);
		OpenVolumeMesh::Geometry::Vec3d np = mesh_->vertex(vh);
		if (v_type[0] >= 0)
		{
			fprintf(f_pcde, "%d %15.14f\n", 3 * i + 0, np[0]);
		}
		if (v_type[1] >= 0)
		{
			fprintf(f_pcde, "%d %15.14f\n", 3 * i + 1, np[1]);
		}
		if (v_type[2] >= 0)
		{
			fprintf(f_pcde, "%d %15.14f\n", 3 * i + 2, np[2]);
		}
	}
	fclose(f_pcde);
}

void polycube_deformation_interface::save_polycube_int(const char* filename, VolumeMesh* mesh_)
{
	double cube_len = sigma_r * avg_boundary_edge_length;
	printf("%f\n", cube_len);
	FILE* f_der = fopen(filename, "w");
	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
		OpenVolumeMesh::Geometry::Vec3d q = p / cube_len;
		int v_id = v_it->idx(); OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[v_id];
		if (v_type[0] >= 0)
		{
			q[0] = std::floor(q[0] + 0.5);
			fprintf(f_der, "%d ", (int)(q[0]));
		}
		else
		{
			fprintf(f_der, "%.19f ", (q[0]));
		}
		if (v_type[1] >= 0)
		{
			q[1] = std::floor(q[1] + 0.5);
			fprintf(f_der, "%d ", (int)(q[1]));
		}
		else
		{
			fprintf(f_der, "%.19f ", (q[1]));
		}
		if (v_type[2] >= 0)
		{
			q[2] = std::floor(q[2] + 0.5);
			fprintf(f_der, "%d", (int)(q[2]));
		}
		else
		{
			fprintf(f_der, "%.19f", (q[2]));
		}

		fprintf(f_der, "\n");
	}
	fclose(f_der);
}

void polycube_deformation_interface::exp_mips_deformation_refine_one_polycube(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, double ga, bool update_all)
{
	bool bv_flag = false;
	if (bvf_id[v_id].size() > 0)
	{
		bv_flag = true;
	}
	else
	{
		if (distoriton_big_count[v_id] == 0)
			return;
	}

	const double* p_dpx = dpx.data(); const double* p_dpy = dpy.data(); const double* p_dpz = dpz.data();
	double x0 = p_dpx[v_id]; double y0 = p_dpy[v_id]; double z0 = p_dpz[v_id];
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double local_energy = 0.0; double gx = 0.0; double gy = 0.0; double gz = 0.0;
	double min_radius = 1e30; double len;
	unsigned jj = 0; double alpha = (1.0 - angle_area_ratio)*energy_power; double beta = energy_power - alpha;
	int vc_size = vertex_cell[v_id].size(); const int* vc_id = vertex_cell[v_id].data();
	int v_id_1, v_id_2, v_id_3, k, c_id;
	double D00, D10, D20, D01, D11, D21, D02, D12, D22;
	double C00, C01, C02, C10, C11, C12, C20, C21, C22;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
	double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k, exp_e;
	double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
	double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
	double dvex, dvey, dvez, g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e, d_mips_e_x, d_mips_e_y, d_mips_e_z;
	double* p_vc_pos_x = vc_pos_x2_omp[omp_id].data();
	double* p_vc_pos_y = vc_pos_y2_omp[omp_id].data();
	double* p_vc_pos_z = vc_pos_z2_omp[omp_id].data();
	double* p_vc_n_cross_x = vc_n_cross_x_omp[omp_id].data();
	double* p_vc_n_cross_y = vc_n_cross_y_omp[omp_id].data();
	double* p_vc_n_cross_z = vc_n_cross_z_omp[omp_id].data();
	double* exp_vec = exp_vec_omp[omp_id].data();
	double* gx_vec = gx_vec_omp[omp_id].data();
	double* gy_vec = gy_vec_omp[omp_id].data();
	double* gz_vec = gz_vec_omp[omp_id].data();
	double* mu_vec = mu_vec_omp[omp_id].data();
	std::vector<std::vector<double>>& vc_S = vc_S_omp[omp_id];
	for (unsigned i = 0; i < vc_size; ++i)
	{
		const int* vv_id = vertex_cell_vertex[v_id][i].data();

		memcpy(&vc_S[i][0], &vcv_S[v_id][i][0], 9 * sizeof(double));
		double* s_data = vc_S[i].data();
		k = 3 * i;
		v_id_1 = vv_id[0];
		v_id_2 = vv_id[1];
		v_id_3 = vv_id[2];
		x1 = p_dpx[v_id_1]; x2 = p_dpx[v_id_2]; x3 = p_dpx[v_id_3];
		y1 = p_dpy[v_id_1]; y2 = p_dpy[v_id_2]; y3 = p_dpy[v_id_3];
		z1 = p_dpz[v_id_1]; z2 = p_dpz[v_id_2]; z3 = p_dpz[v_id_3];
		p_vc_pos_x[k + 0] = x1; p_vc_pos_x[k + 1] = x2; p_vc_pos_x[k + 2] = x3;
		p_vc_pos_y[k + 0] = y1; p_vc_pos_y[k + 1] = y2; p_vc_pos_y[k + 2] = y3;
		p_vc_pos_z[k + 0] = z1; p_vc_pos_z[k + 1] = z2; p_vc_pos_z[k + 2] = z3;
		p_vc_n_cross_x[i] = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1);
		p_vc_n_cross_y[i] = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
		p_vc_n_cross_z[i] = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);

		D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
		D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
		D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;

		len = D00*D00 + D10*D10 + D20*D20; if (len < min_radius) min_radius = len;
		len = D01*D01 + D11*D11 + D21*D21; if (len < min_radius) min_radius = len;
		len = D02*D02 + D12*D12 + D22*D22; if (len < min_radius) min_radius = len;

		C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
		C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
		C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];

		A00 = C00*D00 + C10*D01 + C20*D02; d_A00_x = -C00 - C10 - C20;
		A10 = C00*D10 + C10*D11 + C20*D12; d_A10_y = -C00 - C10 - C20;
		A20 = C00*D20 + C10*D21 + C20*D22; d_A20_z = -C00 - C10 - C20;
		A01 = C01*D00 + C11*D01 + C21*D02; d_A01_x = -C01 - C11 - C21;
		A11 = C01*D10 + C11*D11 + C21*D12; d_A11_y = -C01 - C11 - C21;
		A21 = C01*D20 + C11*D21 + C21*D22; d_A21_z = -C01 - C11 - C21;
		A02 = C02*D00 + C12*D01 + C22*D02; d_A02_x = -C02 - C12 - C22;
		A12 = C02*D10 + C12*D11 + C22*D12; d_A12_y = -C02 - C12 - C22;
		A22 = C02*D20 + C12*D21 + C22*D22; d_A22_z = -C02 - C12 - C22;

		A2_00 = A00 *A00 + A10*A10 + A20 *A20;
		A2_01 = A00*A01 + A10*A11 + A20*A21;
		A2_02 = A00*A02 + A10*A12 + A20*A22;
		A2_11 = A01 *A01 + A11*A11 + A21 *A21;
		A2_12 = A01*A02 + A11*A12 + A21*A22;
		A2_22 = A02 *A02 + A12*A12 + A22 *A22;

		AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
		AF2 = A2_00*A2_00 + A2_01*A2_01 * 2 + A2_02*A2_02 * 2 + A2_11*A2_11 + A2_12*A2_12 * 2 + A2_22*A2_22;
		AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
		det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
		i_det_A = 1.0 / det_A;

		D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
		D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
		D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);

		D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00*d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00*d_A02_x)
			+ 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02*d_A01_x));

		D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10*d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10*d_A12_y)
			+ 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12*d_A11_y));

		D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20*d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20*d_A22_z)
			+ 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22*d_A21_z));

		D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
		D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;

		D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
		D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
		D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;

		inv_det_A_2_05 = 0.5 *(1.0 - i_det_A*i_det_A);
		dvex = D_det_A_x*inv_det_A_2_05;
		dvey = D_det_A_y*inv_det_A_2_05;
		dvez = D_det_A_z*inv_det_A_2_05;

		g = AF_05 * AF_I_05;
		dgx = D_AF_x * AF_I_05 + AF_05*D_AF_I_x;
		dgy = D_AF_y * AF_I_05 + AF_05*D_AF_I_y;
		dgz = D_AF_z * AF_I_05 + AF_05*D_AF_I_z;

		dex = (dgx *i_det_A - g * D_det_A_x *i_det_A*i_det_A);
		dey = (dgy *i_det_A - g * D_det_A_y *i_det_A*i_det_A);
		dez = (dgz *i_det_A - g * D_det_A_z *i_det_A*i_det_A);

		e = g*i_det_A;
		mips_e = (e*e - 1.0)*0.125;
		//mips_e = e/3.0;
		volume_e = 0.5*(det_A + i_det_A);
		exp_k = (alpha*mips_e + beta *volume_e);
		if (exp_k > 60) exp_k = 60;
		exp_vec[i] = exp_k;
		/*exp_e = std::exp(exp_k);
		local_energy += exp_e;*/

		d_mips_e_x = e * 0.25 * dex;
		d_mips_e_y = e * 0.25 * dey;
		d_mips_e_z = e * 0.25 * dez;
		/*d_mips_e_x = dex/3.0;
		d_mips_e_y = dey/3.0;
		d_mips_e_z = dez/3.0;*/

		gx_vec[i] = alpha * d_mips_e_x + beta * dvex;
		gy_vec[i] = alpha * d_mips_e_y + beta * dvey;
		gz_vec[i] = alpha * d_mips_e_z + beta * dvez;

		/*dex = alpha * d_mips_e_x + beta * dvex;
		dey = alpha * d_mips_e_y + beta * dvey;
		dez = alpha * d_mips_e_z + beta * dvez;

		dex *= exp_e; dey *= exp_e; dez *= exp_e;*/

		//gx += dex; gy += dey; gz += dez;
	}

#if NDEBUG
	fmath::expd_v(exp_vec, vc_size);
	for (int i = 0; i < vc_size; ++i)
	{
		exp_e = exp_vec[i];
		local_energy += exp_e;
		gx += gx_vec[i] * exp_e;
		gy += gy_vec[i] * exp_e;
		gz += gz_vec[i] * exp_e;
		/*gx += gx_vec[i] ;
		gy += gy_vec[i] ;
		gz += gz_vec[i] ;*/
	}
#else
	for (int i = 0; i < vc_size; ++i)
	{
		exp_e = std::exp(exp_vec[i]);
		local_energy += exp_e;
		gx += gx_vec[i] * exp_e;
		gy += gy_vec[i] * exp_e;
		gz += gz_vec[i] * exp_e;
	}
#endif
	
	double mu = 1.0; int bvf_size = 0; int bvv_size = 0;
	double* p_vf_pos_x = NULL; double* p_vf_pos_y = NULL; double* p_vf_pos_z = NULL;
	if (bv_flag)
	{
		double normal_e = 0.0; double smooth_e = 0;
		double ne_gx = 0.0; double ne_gy = 0.0; double ne_gz = 0.0;
		double se_gx = 0.0; double se_gy = 0.0; double se_gz = 0.0;

		p_vf_pos_x = vc_pos_x3_omp[omp_id].data();
		p_vf_pos_y = vc_pos_y3_omp[omp_id].data();
		p_vf_pos_z = vc_pos_z3_omp[omp_id].data();

		std::vector<int>& one_vv_id = bvv_id[v_id];
		std::vector<int>& one_vf_id = bvf_id[v_id];
		bvv_size = one_vv_id.size(); bvf_size = bvv_size;
		for (int i = 0; i < bvv_size; ++i)
		{
			int j = (i + 1) % bvv_size; int k = (i + 2) % bvv_size;

			x1 = p_dpx[one_vv_id[i]]; y1 = p_dpy[one_vv_id[i]]; z1 = p_dpz[one_vv_id[i]];
			x2 = p_dpx[one_vv_id[j]]; y2 = p_dpy[one_vv_id[j]]; z2 = p_dpz[one_vv_id[j]];
			x3 = p_dpx[one_vv_id[k]]; y3 = p_dpy[one_vv_id[k]]; z3 = p_dpz[one_vv_id[k]];

			p_vf_pos_x[4 * i + 0] = x1; p_vf_pos_x[4 * i + 1] = x2; p_vf_pos_x[4 * i + 2] = x3;
			p_vf_pos_y[4 * i + 0] = y1; p_vf_pos_y[4 * i + 1] = y2; p_vf_pos_y[4 * i + 2] = y3;
			p_vf_pos_z[4 * i + 0] = z1; p_vf_pos_z[4 * i + 1] = z2; p_vf_pos_z[4 * i + 2] = z3;

			OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[one_vf_id[i]];
			p_vf_pos_x[4 * i + 3] = tn[0]; p_vf_pos_y[4 * i + 3] = tn[1]; p_vf_pos_z[4 * i + 3] = tn[2];

			double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
			double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
			double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);

			double len2 = nx*nx + ny*ny + nz*nz;  double len = std::sqrt(len2);

			//double E_ne = 1.0 + 1.0e-8 - (nx*tn[0] + ny*tn[1] + nz*tn[2])/len;//polycube energy
			double E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));//polycube energy
			if (update_all)
			{
				mu_vec[i] = local_energy_ratio[one_vf_id[i]];
			}
			else
			{
				mu_vec[i] = 1.0;
			}
			////local_energy_ratio[one_vf_id[i]];
			normal_e += mu_vec[i] * E_ne;

			double d_nx_x0 = 0.0;     double d_nx_y0 = z1 - z2; double d_nx_z0 = y2 - y1;
			double d_ny_x0 = z2 - z1; double d_ny_y0 = 0.0;     double d_ny_z0 = x1 - x2;
			double d_nz_x0 = y1 - y2; double d_nz_y0 = x2 - x1; double d_nz_z0 = 0.0;

			double d_len2_x0 = (nx*d_nx_x0 + ny*d_ny_x0 + nz*d_nz_x0); //miss 2
			double d_len2_y0 = (nx*d_nx_y0 + ny*d_ny_y0 + nz*d_nz_y0);
			double d_len2_z0 = (nx*d_nx_z0 + ny*d_ny_z0 + nz*d_nz_z0);

			double d_len_x0 = d_len2_x0 / len; double d_len_y0 = d_len2_y0 / len; double d_len_z0 = d_len2_z0 / len;

			double d_nx_len_x0 = d_nx_x0 / len - nx*d_len_x0 / (len*len);
			double d_nx_len_y0 = d_nx_y0 / len - nx*d_len_y0 / (len*len);
			double d_nx_len_z0 = d_nx_z0 / len - nx*d_len_z0 / (len*len);

			double d_ny_len_x0 = d_ny_x0 / len - ny*d_len_x0 / (len*len);
			double d_ny_len_y0 = d_ny_y0 / len - ny*d_len_y0 / (len*len);
			double d_ny_len_z0 = d_ny_z0 / len - ny*d_len_z0 / (len*len);

			double d_nz_len_x0 = d_nz_x0 / len - nz*d_len_x0 / (len*len);
			double d_nz_len_y0 = d_nz_y0 / len - nz*d_len_y0 / (len*len);
			double d_nz_len_z0 = d_nz_z0 / len - nz*d_len_z0 / (len*len);

			ne_gx += mu_vec[i] * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0*(ny / len - tn[1]) + d_nz_len_x0*(nz / len - tn[2]));
			ne_gy += mu_vec[i] * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0*(ny / len - tn[1]) + d_nz_len_y0*(nz / len - tn[2]));
			ne_gz += mu_vec[i] * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0*(ny / len - tn[1]) + d_nz_len_z0*(nz / len - tn[2]));

			//
			/*double nx_2 = (y0 - y2)*(z0 - z3) - (y0 - y3)*(z0 - z2);
			double ny_2 = (x0 - x3)*(z0 - z2) - (x0 - x2)*(z0 - z3);
			double nz_2 = (x0 - x2)*(y0 - y3) - (x0 - x3)*(y0 - y2);

			double len_2 = std::sqrt(nx_2*nx_2 + ny_2*ny_2 + nz_2*nz_2);
			double nx_d = nx / len - nx_2 / len_2; double ny_d = ny / len - ny_2 / len_2; double nz_d = nz / len - nz_2 / len_2;
			double E_smooth = nx_d*nx_d + ny_d*ny_d + nz_d*nz_d; //smooth energy
			smooth_e += E_smooth;

			double d_nx_2_x0 = 0.0;     double d_nx_2_y0 = z2 - z3; double d_nx_2_z0 = y3 - y2;
			double d_ny_2_x0 = z3 - z2; double d_ny_2_y0 = 0.0;     double d_ny_2_z0 = x2 - x3;
			double d_nz_2_x0 = y2 - y3; double d_nz_2_y0 = x3 - x2; double d_nz_2_z0 = 0.0;

			double d_len_2_x0 = (ny_2*d_ny_2_x0 + nz_2*d_nz_2_x0) / len_2;
			double d_len_2_y0 = (nx_2*d_nx_2_y0 + nz_2*d_nz_2_y0) / len_2;
			double d_len_2_z0 = (nx_2*d_nx_2_z0 + ny_2*d_ny_2_z0) / len_2;
			double d_nx_d_x0 = d_nx_x0 / len - nx*d_len_x0 / (len*len) - d_nx_2_x0 / len_2 + nx_2*d_len_2_x0 / (len_2*len_2);
			double d_nx_d_y0 = d_nx_y0 / len - nx*d_len_y0 / (len*len) - d_nx_2_y0 / len_2 + nx_2*d_len_2_y0 / (len_2*len_2);
			double d_nx_d_z0 = d_nx_z0 / len - nx*d_len_z0 / (len*len) - d_nx_2_z0 / len_2 + nx_2*d_len_2_z0 / (len_2*len_2);

			double d_ny_d_x0 = d_ny_x0 / len - ny*d_len_x0 / (len*len) - d_ny_2_x0 / len_2 + ny_2*d_len_2_x0 / (len_2*len_2);
			double d_ny_d_y0 = d_ny_y0 / len - ny*d_len_y0 / (len*len) - d_ny_2_y0 / len_2 + ny_2*d_len_2_y0 / (len_2*len_2);
			double d_ny_d_z0 = d_ny_z0 / len - ny*d_len_z0 / (len*len) - d_ny_2_z0 / len_2 + ny_2*d_len_2_z0 / (len_2*len_2);

			double d_nz_d_x0 = d_nz_x0 / len - nz*d_len_x0 / (len*len) - d_nz_2_x0 / len_2 + nz_2*d_len_2_x0 / (len_2*len_2);
			double d_nz_d_y0 = d_nz_y0 / len - nz*d_len_y0 / (len*len) - d_nz_2_y0 / len_2 + nz_2*d_len_2_y0 / (len_2*len_2);
			double d_nz_d_z0 = d_nz_z0 / len - nz*d_len_z0 / (len*len) - d_nz_2_z0 / len_2 + nz_2*d_len_2_z0 / (len_2*len_2);

			double d_E_smooth_x0 = 2.0*(nx_d*d_nx_d_x0 + ny_d*d_ny_d_x0 + nz_d*d_nz_d_x0);
			double d_E_smooth_y0 = 2.0*(nx_d*d_nx_d_y0 + ny_d*d_ny_d_y0 + nz_d*d_nz_d_y0);
			double d_E_smooth_z0 = 2.0*(nx_d*d_nx_d_z0 + ny_d*d_ny_d_z0 + nz_d*d_nz_d_z0);
			se_gx += d_E_smooth_x0; se_gy += d_E_smooth_y0; se_gz += d_E_smooth_z0;*/
		}

		/*std::vector<int>& one_vf_id = bvf_id[v_id];
		bvf_size = one_vf_id.size();
		for (int i = 0; i < bvf_size; ++i)
		{
			OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[one_vf_id[i]];
			OpenVolumeMesh::Geometry::Vec3i& one_fv_id = bfv_id[one_vf_id[i]];
			for (int j = 0; j < 3; ++j)
			{
				if (one_fv_id[j] == v_id)
				{
					x1 = p_dpx[one_fv_id[(j + 1) % 3]]; y1 = p_dpy[one_fv_id[(j + 1) % 3]]; z1 = p_dpz[one_fv_id[(j + 1) % 3]];
					x2 = p_dpx[one_fv_id[(j + 2) % 3]]; y2 = p_dpy[one_fv_id[(j + 2) % 3]]; z2 = p_dpz[one_fv_id[(j + 2) % 3]];
					break;
				}
			}

			p_vf_pos_x[3 * i + 0] = x1; p_vf_pos_x[3 * i + 1] = x2; p_vf_pos_x[3 * i + 2] = tn[0];
			p_vf_pos_y[3 * i + 0] = y1; p_vf_pos_y[3 * i + 1] = y2; p_vf_pos_y[3 * i + 2] = tn[1];
			p_vf_pos_z[3 * i + 0] = z1; p_vf_pos_z[3 * i + 1] = z2; p_vf_pos_z[3 * i + 2] = tn[2];

			double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
			double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
			double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);

			double len2 = nx*nx + ny*ny + nz*nz;  double len = std::sqrt(len2); 

			//double E_ne = 1.0 + 1.0e-8 - (nx*tn[0] + ny*tn[1] + nz*tn[2])/len;//polycube energy
			double E_ne = 0.5*((nx / len - tn[0])*(nx / len - tn[0]) + (ny / len - tn[1])*(ny / len - tn[1]) + (nz / len - tn[2])*(nz / len - tn[2]));//polycube energy
			if (update_all)
			{
				mu_vec[i] = local_energy_ratio[one_vf_id[i]];
			}
			else
			{
				mu_vec[i] = 1.0;
			}
			////local_energy_ratio[one_vf_id[i]];
			normal_e += mu_vec[i] * E_ne;

			double d_nx_x0 = 0.0;     double d_nx_y0 = z1 - z2; double d_nx_z0 = y2 - y1;
			double d_ny_x0 = z2 - z1; double d_ny_y0 = 0.0;     double d_ny_z0 = x1 - x2;
			double d_nz_x0 = y1 - y2; double d_nz_y0 = x2 - x1; double d_nz_z0 = 0.0;

			double d_len2_x0 = (nx*d_nx_x0 + ny*d_ny_x0 + nz*d_nz_x0); //miss 2
			double d_len2_y0 = (nx*d_nx_y0 + ny*d_ny_y0 + nz*d_nz_y0);
			double d_len2_z0 = (nx*d_nx_z0 + ny*d_ny_z0 + nz*d_nz_z0);

			double d_len_x0 = d_len2_x0 / len; double d_len_y0 = d_len2_y0 / len; double d_len_z0 = d_len2_z0 / len;

			double d_nx_len_x0 = d_nx_x0 / len - nx*d_len_x0 / (len*len);
			double d_nx_len_y0 = d_nx_y0 / len - nx*d_len_y0 / (len*len);
			double d_nx_len_z0 = d_nx_z0 / len - nx*d_len_z0 / (len*len);

			double d_ny_len_x0 = d_ny_x0 / len - ny*d_len_x0 / (len*len);
			double d_ny_len_y0 = d_ny_y0 / len - ny*d_len_y0 / (len*len);
			double d_ny_len_z0 = d_ny_z0 / len - ny*d_len_z0 / (len*len);

			double d_nz_len_x0 = d_nz_x0 / len - nz*d_len_x0 / (len*len);
			double d_nz_len_y0 = d_nz_y0 / len - nz*d_len_y0 / (len*len);
			double d_nz_len_z0 = d_nz_z0 / len - nz*d_len_z0 / (len*len);

			ne_gx += mu_vec[i] * (d_nx_len_x0*(nx / len - tn[0]) + d_ny_len_x0*(ny / len - tn[1]) + d_nz_len_x0*(nz / len - tn[2]));
			ne_gy += mu_vec[i] * (d_nx_len_y0*(nx / len - tn[0]) + d_ny_len_y0*(ny / len - tn[1]) + d_nz_len_y0*(nz / len - tn[2]));
			ne_gz += mu_vec[i] * (d_nx_len_z0*(nx / len - tn[0]) + d_ny_len_z0*(ny / len - tn[1]) + d_nz_len_z0*(nz / len - tn[2]));
		}*/

		if (update_all)
		{
			local_energy += normal_e ;
			gx += ne_gx; gy += ne_gy; gz += ne_gz;
			/*local_energy += 5 * smooth_e;
			gx += 5 * se_gx; gy += 5 * se_gy; gz += 5 * se_gz;*/
		}
		else
		{
			double s = local_energy / (normal_e + 1e-8);
			//if (s < ga)
			{
				ga = s;
			}
			if (s > 1e16) ga = 1e16;
			//else ga *= s;
			for (int i = 0; i < bvf_size; ++i)
			{
				mu_vec[i] = ga;
			}

			local_energy += ga * normal_e;
			gx += ga*ne_gx; gy += ga*ne_gy; gz += ga*ne_gz;
		}
		
	}

	min_radius = std::sqrt(min_radius);
	double dx = gx, dy = gy, dz = gz;
	double move_d = sqrt(dx*dx + dy*dy + dz*dz);
	if (move_d > min_radius)
	{
		double t = min_radius / move_d;
		dx *= t; dy *= t; dz *= t;
		move_d = min_radius;
	}

	double npx = x0 - dx; double npy = y0 - dy; double npz = z0 - dz;
	while (!local_check_negative_volume4(vc_size, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, p_vc_n_cross_x, p_vc_n_cross_y, p_vc_n_cross_z, npx, npy, npz))
	{
		dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
		if (move_d < 1e-8)
		{
			npx = x0; npy = y0; npz = z0;
			move_d = 1e-8; break;
		}
		else
		{
			npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
		}
	}
#if 1
	double new_e = 0;
	//if (check_flag) new_e = 2.0*move_d*move_d;
	bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
	if (!e_ok)
	{
		while (!e_ok)
		{
			dx *= 0.2; dy *= 0.2; dz *= 0.2; move_d *= 0.2;
			if (move_d < 1e-8)
			{
				npx = x0; npy = y0; npz = z0;
				move_d = 1e-8; break;
			}
			else
			{
				npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
			}
			new_e = 0;
			e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
		}

		if (e_ok)
		{
			dx *= 5; dy *= 5; dz *= 5; move_d *= 5; e_ok = false;
			while (!e_ok)
			{
				dx *= 0.5; dy *= 0.5; dz *= 0.5; move_d *= 0.5;
				if (move_d < 1e-8)
				{
					npx = x0; npy = y0; npz = z0;
					move_d = 1e-8; break;
				}
				else
				{
					npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
				}
				new_e = 0;
				e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
			}

			if (e_ok)
			{
				dx *= 2; dy *= 2; dz *= 2; move_d *= 2; e_ok = false;
				while (!e_ok)
				{
					dx *= 0.875; dy *= 0.875; dz *= 0.875; move_d *= 0.875;
					if (move_d < 1e-8)
					{
						npx = x0; npy = y0; npz = z0;
						move_d = 1e-8; break;
					}
					else
					{
						npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
					}
					new_e = 0;
					e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z, mu_vec);
				}
			}
		}
	}
#else
	double new_e = 0;
	bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
	while (!e_ok)
	{
		dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
		if (move_d < 1e-8)
		{
			npx = x0; npy = y0; npz = z0;
			move_d = 1e-8; break;
		}
		else
		{
			npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
		}
		new_e = 0;
		e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
	}
#endif
	dpx[v_id] = npx; dpy[v_id] = npy; dpz[v_id] = npz;

	//
	check_big_distortion_cell(vc_size, large_dis_flag_omp[omp_id], p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, npx, npy, npz);
	for (int i = 0; i < vc_size; ++i)
	{
		int c_id = vertex_cell[v_id][i];
		if (large_dis_flag_omp[omp_id][i] == 1) // large
		{
			if (distoriton_big_cell[c_id] == 1) //original is large
			{

			}
			else //original is small, from small to large
			{
				distoriton_big_cell[c_id] = 1;
				//add 1 for all vertices
				for (int j = 0; j < 4; ++j)
				{
					//con_dis_above_th_vertex[cell_vertex[c_id][j]] += 1;
					distoriton_big_count[cell_vertex[c_id][j]] += 1;
				}
			}
		}
		else
		{
			if (distoriton_big_cell[c_id] == 1) //original is large, from large to small
			{
				distoriton_big_cell[c_id] = -1;
				//decrease 1 for all vertices
				for (int j = 0; j < 4; ++j)
				{
					//con_dis_above_th_vertex[cell_vertex[c_id][j]] -= 1;
					distoriton_big_count[cell_vertex[c_id][j]] -= 1;
				}
			}
			else//original is small
			{

			}
		}
	}
}

void polycube_deformation_interface::check_big_distortion_cell(
	const int& vc_size, std::vector<int>& large_distortion_flag,
	const double* posx, const double* posy, const double* posz,
	const std::vector<std::vector<double> >& vc_S,
	const double& npx, const double& npy, const double& npz)
{
	double D00, D10, D20, D01, D11, D21, D02, D12, D22;
	double C00, C01, C02, C10, C11, C12, C20, C21, C22;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
	double AF, AF2, AF_I, e2, mips_e, volume_e, k;
	for (unsigned i = 0; i < vc_size; ++i)
	{
		const double* s_data = vc_S[i].data();
		int j = 3 * i;
		D00 = posx[j + 0] - npx; D10 = posy[j + 0] - npy; D20 = posz[j + 0] - npz;
		D01 = posx[j + 1] - npx; D11 = posy[j + 1] - npy; D21 = posz[j + 1] - npz;
		D02 = posx[j + 2] - npx; D12 = posy[j + 2] - npy; D22 = posz[j + 2] - npz;

		C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
		C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
		C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];

		A00 = C00*D00 + C10*D01 + C20*D02;
		A10 = C00*D10 + C10*D11 + C20*D12;
		A20 = C00*D20 + C10*D21 + C20*D22;
		A01 = C01*D00 + C11*D01 + C21*D02;
		A11 = C01*D10 + C11*D11 + C21*D12;
		A21 = C01*D20 + C11*D21 + C21*D22;
		A02 = C02*D00 + C12*D01 + C22*D02;
		A12 = C02*D10 + C12*D11 + C22*D12;
		A22 = C02*D20 + C12*D21 + C22*D22;

		A2_00 = A00 *A00 + A10*A10 + A20 *A20;
		A2_01 = A00*A01 + A10*A11 + A20*A21;
		A2_02 = A00*A02 + A10*A12 + A20*A22;
		A2_11 = A01 *A01 + A11*A11 + A21 *A21;
		A2_12 = A01*A02 + A11*A12 + A21*A22;
		A2_22 = A02 *A02 + A12*A12 + A22 *A22;

		AF = A2_00 + A2_11 + A2_22;
		AF2 = A2_00*A2_00 + A2_01*A2_01 * 2 + A2_02*A2_02 * 2 + A2_11*A2_11 + A2_12*A2_12 * 2 + A2_22*A2_22;
		AF_I = (AF*AF - AF2)*0.5;

		det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
		i_det_A = 1/det_A;
		e2 = AF*AF_I * (i_det_A*i_det_A);
		mips_e = (e2 - 1.0)*0.125;
		volume_e = 0.5*(det_A + i_det_A);
		k = (mips_e + volume_e)*0.5;
		if (k > bound_K) large_distortion_flag[i] = 1;
		else large_distortion_flag[i] = -1;
	}
}

void polycube_deformation_interface::exp_mips_deformation_refine_one_polycube_normal(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, double mu1, double mu2, bool update_all)
{
	bool bv_flag = false;
	if (bvv_id[v_id].size() > 0) bv_flag = true;

	const double* p_dpx = dpx.data(); const double* p_dpy = dpy.data(); const double* p_dpz = dpz.data();
	double x0 = p_dpx[v_id]; double y0 = p_dpy[v_id]; double z0 = p_dpz[v_id];
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double local_energy = 0.0; double gx = 0.0; double gy = 0.0; double gz = 0.0;
	double min_radius = 1e30; double len;
	unsigned jj = 0; double alpha = (1.0 - angle_area_ratio)*energy_power; double beta = energy_power - alpha;
	int vc_size = vertex_cell[v_id].size(); const int* vc_id = vertex_cell[v_id].data();
	int v_id_1, v_id_2, v_id_3, k, c_id;
	double D00, D10, D20, D01, D11, D21, D02, D12, D22;
	double C00, C01, C02, C10, C11, C12, C20, C21, C22;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
	double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k, exp_e;
	double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
	double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
	double dvex, dvey, dvez, g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e, d_mips_e_x, d_mips_e_y, d_mips_e_z;
	double* p_vc_pos_x = vc_pos_x2_omp[omp_id].data();
	double* p_vc_pos_y = vc_pos_y2_omp[omp_id].data();
	double* p_vc_pos_z = vc_pos_z2_omp[omp_id].data();
	double* p_vc_n_cross_x = vc_n_cross_x_omp[omp_id].data();
	double* p_vc_n_cross_y = vc_n_cross_y_omp[omp_id].data();
	double* p_vc_n_cross_z = vc_n_cross_z_omp[omp_id].data();
	double* exp_vec = exp_vec_omp[omp_id].data();
	double* gx_vec = gx_vec_omp[omp_id].data();
	double* gy_vec = gy_vec_omp[omp_id].data();
	double* gz_vec = gz_vec_omp[omp_id].data();
	std::vector<std::vector<double>>& vc_S = vc_S_omp[omp_id];
	for (unsigned i = 0; i < vc_size; ++i)
	{
		const int* vv_id = vertex_cell_vertex[v_id][i].data();

		memcpy(&vc_S[i][0], &vcv_S[v_id][i][0], 9 * sizeof(double));
		double* s_data = vc_S[i].data();
		k = 3 * i;
		v_id_1 = vv_id[0];
		v_id_2 = vv_id[1];
		v_id_3 = vv_id[2];
		x1 = p_dpx[v_id_1]; x2 = p_dpx[v_id_2]; x3 = p_dpx[v_id_3];
		y1 = p_dpy[v_id_1]; y2 = p_dpy[v_id_2]; y3 = p_dpy[v_id_3];
		z1 = p_dpz[v_id_1]; z2 = p_dpz[v_id_2]; z3 = p_dpz[v_id_3];
		p_vc_pos_x[k + 0] = x1; p_vc_pos_x[k + 1] = x2; p_vc_pos_x[k + 2] = x3;
		p_vc_pos_y[k + 0] = y1; p_vc_pos_y[k + 1] = y2; p_vc_pos_y[k + 2] = y3;
		p_vc_pos_z[k + 0] = z1; p_vc_pos_z[k + 1] = z2; p_vc_pos_z[k + 2] = z3;
		p_vc_n_cross_x[i] = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1);
		p_vc_n_cross_y[i] = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
		p_vc_n_cross_z[i] = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);

		D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
		D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
		D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;

		len = D00*D00 + D10*D10 + D20*D20; if (len < min_radius) min_radius = len;
		len = D01*D01 + D11*D11 + D21*D21; if (len < min_radius) min_radius = len;
		len = D02*D02 + D12*D12 + D22*D22; if (len < min_radius) min_radius = len;

		C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
		C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
		C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];

		A00 = C00*D00 + C10*D01 + C20*D02; d_A00_x = -C00 - C10 - C20;
		A10 = C00*D10 + C10*D11 + C20*D12; d_A10_y = -C00 - C10 - C20;
		A20 = C00*D20 + C10*D21 + C20*D22; d_A20_z = -C00 - C10 - C20;
		A01 = C01*D00 + C11*D01 + C21*D02; d_A01_x = -C01 - C11 - C21;
		A11 = C01*D10 + C11*D11 + C21*D12; d_A11_y = -C01 - C11 - C21;
		A21 = C01*D20 + C11*D21 + C21*D22; d_A21_z = -C01 - C11 - C21;
		A02 = C02*D00 + C12*D01 + C22*D02; d_A02_x = -C02 - C12 - C22;
		A12 = C02*D10 + C12*D11 + C22*D12; d_A12_y = -C02 - C12 - C22;
		A22 = C02*D20 + C12*D21 + C22*D22; d_A22_z = -C02 - C12 - C22;

		A2_00 = A00 *A00 + A10*A10 + A20 *A20;
		A2_01 = A00*A01 + A10*A11 + A20*A21;
		A2_02 = A00*A02 + A10*A12 + A20*A22;
		A2_11 = A01 *A01 + A11*A11 + A21 *A21;
		A2_12 = A01*A02 + A11*A12 + A21*A22;
		A2_22 = A02 *A02 + A12*A12 + A22 *A22;

		AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
		AF2 = A2_00*A2_00 + A2_01*A2_01 * 2 + A2_02*A2_02 * 2 + A2_11*A2_11 + A2_12*A2_12 * 2 + A2_22*A2_22;
		AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
		det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
		i_det_A = 1.0 / det_A;

		D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
		D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
		D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);

		D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00*d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00*d_A02_x)
			+ 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02*d_A01_x));

		D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10*d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10*d_A12_y)
			+ 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12*d_A11_y));

		D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20*d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20*d_A22_z)
			+ 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22*d_A21_z));

		D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
		D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;

		D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
		D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
		D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;

		inv_det_A_2_05 = 0.5 *(1.0 - i_det_A*i_det_A);
		dvex = D_det_A_x*inv_det_A_2_05;
		dvey = D_det_A_y*inv_det_A_2_05;
		dvez = D_det_A_z*inv_det_A_2_05;

		g = AF_05 * AF_I_05;
		dgx = D_AF_x * AF_I_05 + AF_05*D_AF_I_x;
		dgy = D_AF_y * AF_I_05 + AF_05*D_AF_I_y;
		dgz = D_AF_z * AF_I_05 + AF_05*D_AF_I_z;

		dex = (dgx *i_det_A - g * D_det_A_x *i_det_A*i_det_A);
		dey = (dgy *i_det_A - g * D_det_A_y *i_det_A*i_det_A);
		dez = (dgz *i_det_A - g * D_det_A_z *i_det_A*i_det_A);

		e = g*i_det_A;
		mips_e = (e*e - 1.0)*0.125;
		//mips_e = e/3.0;
		volume_e = 0.5*(det_A + i_det_A);
		exp_k = (alpha*mips_e + beta *volume_e);
		if (exp_k > 60) exp_k = 60;
		exp_vec[i] = exp_k;
		/*exp_e = std::exp(exp_k);
		local_energy += exp_e;*/

		d_mips_e_x = e * 0.25 * dex;
		d_mips_e_y = e * 0.25 * dey;
		d_mips_e_z = e * 0.25 * dez;
		/*d_mips_e_x = dex/3.0;
		d_mips_e_y = dey/3.0;
		d_mips_e_z = dez/3.0;*/

		gx_vec[i] = alpha * d_mips_e_x + beta * dvex;
		gy_vec[i] = alpha * d_mips_e_y + beta * dvey;
		gz_vec[i] = alpha * d_mips_e_z + beta * dvez;

		/*dex = alpha * d_mips_e_x + beta * dvex;
		dey = alpha * d_mips_e_y + beta * dvey;
		dez = alpha * d_mips_e_z + beta * dvez;

		dex *= exp_e; dey *= exp_e; dez *= exp_e;

		gx += dex; gy += dey; gz += dez;*/
	}

#if NDEBUG
	fmath::expd_v(exp_vec, vc_size);
	for (int i = 0; i < vc_size; ++i)
	{
		exp_e = exp_vec[i];
		local_energy += exp_e;
		gx += gx_vec[i] * exp_e;
		gy += gy_vec[i] * exp_e;
		gz += gz_vec[i] * exp_e;
		/*gx += gx_vec[i] ;
		gy += gy_vec[i] ;
		gz += gz_vec[i] ;*/
	}
#else
	for (int i = 0; i < vc_size; ++i)
	{
		exp_e = std::exp(exp_vec[i]);
		local_energy += exp_e;
		gx += gx_vec[i] * exp_e;
		gy += gy_vec[i] * exp_e;
		gz += gz_vec[i] * exp_e;
	}
#endif

	int bvv_size = 0; 
	double* p_vv_pos_x = NULL; double* p_vv_pos_y = NULL; double* p_vv_pos_z = NULL;
	if (bv_flag)
	{
		double polycube_e = 0.0; double ga = 0.001;
		double pe_gx = 0.0; double pe_gy = 0.0; double pe_gz = 0.0;
		double se_gx = 0.0; double se_gy = 0.0; double se_gz = 0.0;
		p_vv_pos_x = vc_pos_x3_omp[omp_id].data();
		p_vv_pos_y = vc_pos_y3_omp[omp_id].data();
		p_vv_pos_z = vc_pos_z3_omp[omp_id].data();

		std::vector<int>& one_vv_id = bvv_id[v_id];
		bvv_size = one_vv_id.size();
		for (int i = 0; i < bvv_size; ++i)
		{
			int j = (i + 1) % bvv_size; int k = (i + 2) % bvv_size;

			x1 = p_dpx[one_vv_id[i]]; y1 = p_dpy[one_vv_id[i]]; z1 = p_dpz[one_vv_id[i]];
			x2 = p_dpx[one_vv_id[j]]; y2 = p_dpy[one_vv_id[j]]; z2 = p_dpz[one_vv_id[j]];
			x3 = p_dpx[one_vv_id[k]]; y3 = p_dpy[one_vv_id[k]]; z3 = p_dpz[one_vv_id[k]];

			p_vv_pos_x[3 * i + 0] = x1; p_vv_pos_x[3 * i + 1] = x2; p_vv_pos_x[3 * i + 2] = x3;
			p_vv_pos_y[3 * i + 0] = y1; p_vv_pos_y[3 * i + 1] = y2; p_vv_pos_y[3 * i + 2] = y3;
			p_vv_pos_z[3 * i + 0] = z1; p_vv_pos_z[3 * i + 1] = z2; p_vv_pos_z[3 * i + 2] = z3;

			double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1); double nx2 = nx*nx;
			double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2); double ny2 = ny*ny;
			double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1); double nz2 = nz*nz;

			double len2 = nx*nx + ny*ny + nz*nz; double len4 = len2*len2;
			double len = std::sqrt(len2); double inv_len_4 = 1.0 / len4; double inv_len_2 = 1.0 / len2;
			double inv_len = 1.0 / len;

			/*double pe = (nx2*ny2 + ny2*nz2 + nz2*nx2);
			double E_pe = pe*inv_len_4;//polycube energy*/
			
			double cx = std::sqrt(inv_len_2 *nx2 + ga); double cy = std::sqrt(inv_len_2 *ny2 + ga); double cz = std::sqrt(inv_len_2 *nz2 + ga);

			double E_pe = cx + cy + cz;//polycube energy
			polycube_e += E_pe;

			double d_nx_x0 = 0.0;     double d_nx_y0 = z1 - z2; double d_nx_z0 = y2 - y1;
			double d_ny_x0 = z2 - z1; double d_ny_y0 = 0.0;     double d_ny_z0 = x1 - x2;
			double d_nz_x0 = y1 - y2; double d_nz_y0 = x2 - x1; double d_nz_z0 = 0.0;
			
			double d_len2_x0 = (nx*d_nx_x0 + ny*d_ny_x0 + nz*d_nz_x0); //miss 2
			double d_len2_y0 = (nx*d_nx_y0 + ny*d_ny_y0 + nz*d_nz_y0);
			double d_len2_z0 = (nx*d_nx_z0 + ny*d_ny_z0 + nz*d_nz_z0);

			/*double d_pe_x0 = (ny2 * 2 * nx*d_nx_x0 + nx2 * 2 * ny*d_ny_x0 
				+ nz2 * 2 * ny*d_ny_x0 + ny2 * 2 * nz*d_nz_x0 
				+ nz2 * 2 * nx*d_nx_x0 + nx2 * 2 * nz*d_nz_x0);
			double d_pe_y0 = (ny2 * 2 * nx*d_nx_y0 + nx2 * 2 * ny*d_ny_y0
				+ nz2 * 2 * ny*d_ny_y0 + ny2 * 2 * nz*d_nz_y0
				+ nz2 * 2 * nx*d_nx_y0 + nx2 * 2 * nz*d_nz_y0);
			double d_pe_z0 = (ny2 * 2 * nx*d_nx_z0 + nx2 * 2 * ny*d_ny_z0
				+ nz2 * 2 * ny*d_ny_z0 + ny2 * 2 * nz*d_nz_z0
				+ nz2 * 2 * nx*d_nx_z0 + nx2 * 2 * nz*d_nz_z0);

			double d_E_pe_x0 = d_pe_x0*inv_len_4 - 4.0*d_len2_x0*inv_len_2*inv_len_4*pe;
			double d_E_pe_y0 = d_pe_y0*inv_len_4 - 4.0*d_len2_y0*inv_len_2*inv_len_4*pe;
			double d_E_pe_z0 = d_pe_z0*inv_len_4 - 4.0*d_len2_z0*inv_len_2*inv_len_4*pe;
*/

			/*double nx_2 = (y0 - y2)*(z0 - z3) - (y0 - y3)*(z0 - z2);
			double ny_2 = (x0 - x3)*(z0 - z2) - (x0 - x2)*(z0 - z3);
			double nz_2 = (x0 - x2)*(y0 - y3) - (x0 - x3)*(y0 - y2);

			double len_2 = std::sqrt(nx_2*nx_2 + ny_2*ny_2 + nz_2*nz_2);
			double nx_d = nx / len - nx_2 / len_2; double ny_d = ny / len - ny_2 / len_2; double nz_d = nz / len - nz_2 / len_2;
			double E_smooth = nx_d*nx_d + ny_d*ny_d + nz_d*nz_d; //smooth energy
			smooth_e += E_smooth;*/

			/*double d_nx_2_x0 = 0.0;     double d_nx_2_y0 = z2 - z3; double d_nx_2_z0 = y3 - y2;
			double d_ny_2_x0 = z3 - z2; double d_ny_2_y0 = 0.0;     double d_ny_2_z0 = x2 - x3;
			double d_nz_2_x0 = y2 - y3; double d_nz_2_y0 = x3 - x2; double d_nz_2_z0 = 0.0;*/

			double d_len_x0 = d_len2_x0 *inv_len; double d_len_y0 = d_len2_y0 *inv_len; double d_len_z0 = d_len2_z0 *inv_len;
			/*double d_len_2_x0 = (ny_2*d_ny_2_x0 + nz_2*d_nz_2_x0) / len_2;
			double d_len_2_y0 = (nx_2*d_nx_2_y0 + nz_2*d_nz_2_y0) / len_2;
			double d_len_2_z0 = (nx_2*d_nx_2_z0 + ny_2*d_ny_2_z0) / len_2;*/

			double d_E_pe_x0 = (nx *inv_len)*(d_nx_x0 *inv_len - nx*d_len_x0 *inv_len_2) / cx + (ny *inv_len)*(d_ny_x0*inv_len - ny*d_len_x0 *inv_len_2) / cy + (nz*inv_len)*(d_nz_x0*inv_len - nz*d_len_x0*inv_len_2) / cz;
			double d_E_pe_y0 = (nx *inv_len)*(d_nx_y0 *inv_len - nx*d_len_y0 *inv_len_2) / cx + (ny *inv_len)*(d_ny_y0*inv_len - ny*d_len_y0 *inv_len_2) / cy + (nz*inv_len)*(d_nz_y0*inv_len - nz*d_len_y0*inv_len_2) / cz;
			double d_E_pe_z0 = (nx *inv_len)*(d_nx_z0 *inv_len - nx*d_len_z0 *inv_len_2) / cx + (ny *inv_len)*(d_ny_z0*inv_len - ny*d_len_z0 *inv_len_2) / cy + (nz*inv_len)*(d_nz_z0*inv_len - nz*d_len_z0*inv_len_2) / cz;

			/*double d_nx_d_x0 = d_nx_x0 / len - nx*d_len_x0 / (len*len) - d_nx_2_x0 / len_2 + nx_2*d_len_2_x0 / (len_2*len_2);
			double d_nx_d_y0 = d_nx_y0 / len - nx*d_len_y0 / (len*len) - d_nx_2_y0 / len_2 + nx_2*d_len_2_y0 / (len_2*len_2);
			double d_nx_d_z0 = d_nx_z0 / len - nx*d_len_z0 / (len*len) - d_nx_2_z0 / len_2 + nx_2*d_len_2_z0 / (len_2*len_2);

			double d_ny_d_x0 = d_ny_x0 / len - ny*d_len_x0 / (len*len) - d_ny_2_x0 / len_2 + ny_2*d_len_2_x0 / (len_2*len_2);
			double d_ny_d_y0 = d_ny_y0 / len - ny*d_len_y0 / (len*len) - d_ny_2_y0 / len_2 + ny_2*d_len_2_y0 / (len_2*len_2);
			double d_ny_d_z0 = d_ny_z0 / len - ny*d_len_z0 / (len*len) - d_ny_2_z0 / len_2 + ny_2*d_len_2_z0 / (len_2*len_2);

			double d_nz_d_x0 = d_nz_x0 / len - nz*d_len_x0 / (len*len) - d_nz_2_x0 / len_2 + nz_2*d_len_2_x0 / (len_2*len_2);
			double d_nz_d_y0 = d_nz_y0 / len - nz*d_len_y0 / (len*len) - d_nz_2_y0 / len_2 + nz_2*d_len_2_y0 / (len_2*len_2);
			double d_nz_d_z0 = d_nz_z0 / len - nz*d_len_z0 / (len*len) - d_nz_2_z0 / len_2 + nz_2*d_len_2_z0 / (len_2*len_2);

			double d_E_smooth_x0 = 2.0*(nx_d*d_nx_d_x0 + ny_d*d_ny_d_x0 + nz_d*d_nz_d_x0);
			double d_E_smooth_y0 = 2.0*(nx_d*d_nx_d_y0 + ny_d*d_ny_d_y0 + nz_d*d_nz_d_y0);
			double d_E_smooth_z0 = 2.0*(nx_d*d_nx_d_z0 + ny_d*d_ny_d_z0 + nz_d*d_nz_d_z0);*/

			pe_gx += d_E_pe_x0; pe_gy += d_E_pe_y0; pe_gz += d_E_pe_z0;
			//se_gx += d_E_smooth_x0; se_gy += d_E_smooth_y0; se_gz += d_E_smooth_z0;
		}

		double s = local_energy / (polycube_e + 1e-8);
		if (s > 1e16) mu1 *= 1e16;
		else mu1 *= s;

		local_energy += mu1 * polycube_e;
		gx += mu1*pe_gx ;
		gy += mu1*pe_gy;
		gz += mu1*pe_gz;
	}

	min_radius = std::sqrt(min_radius);
	double dx = gx, dy = gy, dz = gz;
	double move_d = sqrt(dx*dx + dy*dy + dz*dz);
	if (move_d > min_radius)
	{
		double t = min_radius / move_d;
		dx *= t; dy *= t; dz *= t;
		move_d = min_radius;
	}

	double npx = x0 - dx; double npy = y0 - dy; double npz = z0 - dz;
	while (!local_check_negative_volume4(vc_size, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, p_vc_n_cross_x, p_vc_n_cross_y, p_vc_n_cross_z, npx, npy, npz))
	{
		dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
		if (move_d < 1e-8)
		{
			npx = x0; npy = y0; npz = z0;
			move_d = 1e-8; break;
		}
		else
		{
			npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
		}
	}
#if 1
	double new_e = 0;
	//if (check_flag) new_e = 2.0*move_d*move_d;
	bool e_ok = compute_exp_misp_energy_refine_polycube_normal(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, mu1, mu2, bvv_size, p_vv_pos_x, p_vv_pos_y, p_vv_pos_z);
	if (!e_ok)
	{
		while (!e_ok)
		{
			dx *= 0.2; dy *= 0.2; dz *= 0.2; move_d *= 0.2;
			if (move_d < 1e-8)
			{
				npx = x0; npy = y0; npz = z0;
				move_d = 1e-8; break;
			}
			else
			{
				npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
			}
			new_e = 0;
			//if (check_flag) new_e = 2.0*move_d*move_d;
			e_ok = compute_exp_misp_energy_refine_polycube_normal(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, mu1, mu2, bvv_size, p_vv_pos_x, p_vv_pos_y, p_vv_pos_z);
		}

		if (e_ok)
		{
			dx *= 5; dy *= 5; dz *= 5; move_d *= 5; e_ok = false;
			while (!e_ok)
			{
				dx *= 0.5; dy *= 0.5; dz *= 0.5; move_d *= 0.5;
				if (move_d < 1e-8)
				{
					npx = x0; npy = y0; npz = z0;
					move_d = 1e-8; break;
				}
				else
				{
					npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
				}
				new_e = 0;
				//if (check_flag) new_e = 2.0*move_d*move_d;
				e_ok = compute_exp_misp_energy_refine_polycube_normal(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, mu1, mu2, bvv_size, p_vv_pos_x, p_vv_pos_y, p_vv_pos_z);
			}

			if (e_ok)
			{
				dx *= 2; dy *= 2; dz *= 2; move_d *= 2; e_ok = false;
				while (!e_ok)
				{
					dx *= 0.875; dy *= 0.875; dz *= 0.875; move_d *= 0.875;
					if (move_d < 1e-8)
					{
						npx = x0; npy = y0; npz = z0;
						move_d = 1e-8; break;
					}
					else
					{
						npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
					}
					new_e = 0;
					//if (check_flag) new_e = 2.0*move_d*move_d;
					e_ok = compute_exp_misp_energy_refine_polycube_normal(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, mu1, mu2, bvv_size, p_vv_pos_x, p_vv_pos_y, p_vv_pos_z);
				}
			}
		}
	}
#else
	double new_e = 0;
	//bool e_ok = compute_exp_misp_energy4(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta);
	bool e_ok = compute_exp_misp_energy_refine_polycube_normal(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, mu1, mu2, bvv_size, p_vv_pos_x, p_vv_pos_y, p_vv_pos_z);
	while (!e_ok)
	{
		dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
		if (move_d < 1e-8)
		{
			npx = x0; npy = y0; npz = z0;
			move_d = 1e-8; break;
		}
		else
		{
			npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
		}
		new_e = 0;
		e_ok = compute_exp_misp_energy_refine_polycube_normal(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, mu1, mu2, bvv_size, p_vv_pos_x, p_vv_pos_y, p_vv_pos_z);
	}
#endif
	dpx[v_id] = npx; dpy[v_id] = npy; dpz[v_id] = npz;
}

/*
void polycube_deformation_interface::exp_mips_deformation_refine_one_polycube_position(int v_id, double angle_area_ratio, double energy_power, const int& omp_id, double ga)
{
	bool bv_flag = false;
	if (bvf_id[v_id].size() > 0) bv_flag = true;

	const double* p_dpx = dpx.data(); const double* p_dpy = dpy.data(); const double* p_dpz = dpz.data();
	double x0 = p_dpx[v_id]; double y0 = p_dpy[v_id]; double z0 = p_dpz[v_id];
	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double local_energy = 0.0; double gx = 0.0; double gy = 0.0; double gz = 0.0;
	double min_radius = 1e30; double len;
	unsigned jj = 0; double alpha = (1.0 - angle_area_ratio)*energy_power; double beta = energy_power - alpha;
	int vc_size = vertex_cell[v_id].size(); const int* vc_id = vertex_cell[v_id].data();
	int v_id_1, v_id_2, v_id_3, k, c_id;
	double D00, D10, D20, D01, D11, D21, D02, D12, D22;
	double C00, C01, C02, C10, C11, C12, C20, C21, C22;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
	double AF, AF_05, I_AF_05, AF2, AF_I, AF_I_05, I_AF_I_05, exp_k, exp_e;
	double D_AF_x, D_AF_y, D_AF_z, D_AF2_x, D_AF2_y, D_AF2_z, D_AF_I_x, D_AF_I_y, D_AF_I_z, D_det_A_x, D_det_A_y, D_det_A_z, inv_det_A_2_05;
	double d_A00_x, d_A10_y, d_A20_z, d_A01_x, d_A11_y, d_A21_z, d_A02_x, d_A12_y, d_A22_z;
	double dvex, dvey, dvez, g, dgx, dgy, dgz, dex, dey, dez, e, mips_e, volume_e, d_mips_e_x, d_mips_e_y, d_mips_e_z;
	double* p_vc_pos_x = vc_pos_x2_omp[omp_id].data();
	double* p_vc_pos_y = vc_pos_y2_omp[omp_id].data();
	double* p_vc_pos_z = vc_pos_z2_omp[omp_id].data();
	double* p_vc_n_cross_x = vc_n_cross_x_omp[omp_id].data();
	double* p_vc_n_cross_y = vc_n_cross_y_omp[omp_id].data();
	double* p_vc_n_cross_z = vc_n_cross_z_omp[omp_id].data();
	double* exp_vec = exp_vec_omp[omp_id].data();
	double* gx_vec = gx_vec_omp[omp_id].data();
	double* gy_vec = gy_vec_omp[omp_id].data();
	double* gz_vec = gz_vec_omp[omp_id].data();
	std::vector<std::vector<double>>& vc_S = vc_S_omp[omp_id];
	for (unsigned i = 0; i < vc_size; ++i)
	{
		const int* vv_id = vertex_cell_vertex[v_id][i].data();

		memcpy(&vc_S[i][0], &vcv_S[v_id][i][0], 9 * sizeof(double));
		double* s_data = vc_S[i].data();
		k = 3 * i;
		v_id_1 = vv_id[0];
		v_id_2 = vv_id[1];
		v_id_3 = vv_id[2];
		x1 = p_dpx[v_id_1]; x2 = p_dpx[v_id_2]; x3 = p_dpx[v_id_3];
		y1 = p_dpy[v_id_1]; y2 = p_dpy[v_id_2]; y3 = p_dpy[v_id_3];
		z1 = p_dpz[v_id_1]; z2 = p_dpz[v_id_2]; z3 = p_dpz[v_id_3];
		p_vc_pos_x[k + 0] = x1; p_vc_pos_x[k + 1] = x2; p_vc_pos_x[k + 2] = x3;
		p_vc_pos_y[k + 0] = y1; p_vc_pos_y[k + 1] = y2; p_vc_pos_y[k + 2] = y3;
		p_vc_pos_z[k + 0] = z1; p_vc_pos_z[k + 1] = z2; p_vc_pos_z[k + 2] = z3;
		p_vc_n_cross_x[i] = (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1);
		p_vc_n_cross_y[i] = (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
		p_vc_n_cross_z[i] = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);

		D00 = x1 - x0; D10 = y1 - y0; D20 = z1 - z0;
		D01 = x2 - x0; D11 = y2 - y0; D21 = z2 - z0;
		D02 = x3 - x0; D12 = y3 - y0; D22 = z3 - z0;

		len = D00*D00 + D10*D10 + D20*D20; if (len < min_radius) min_radius = len;
		len = D01*D01 + D11*D11 + D21*D21; if (len < min_radius) min_radius = len;
		len = D02*D02 + D12*D12 + D22*D22; if (len < min_radius) min_radius = len;

		C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
		C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
		C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];

		A00 = C00*D00 + C10*D01 + C20*D02; d_A00_x = -C00 - C10 - C20;
		A10 = C00*D10 + C10*D11 + C20*D12; d_A10_y = -C00 - C10 - C20;
		A20 = C00*D20 + C10*D21 + C20*D22; d_A20_z = -C00 - C10 - C20;
		A01 = C01*D00 + C11*D01 + C21*D02; d_A01_x = -C01 - C11 - C21;
		A11 = C01*D10 + C11*D11 + C21*D12; d_A11_y = -C01 - C11 - C21;
		A21 = C01*D20 + C11*D21 + C21*D22; d_A21_z = -C01 - C11 - C21;
		A02 = C02*D00 + C12*D01 + C22*D02; d_A02_x = -C02 - C12 - C22;
		A12 = C02*D10 + C12*D11 + C22*D12; d_A12_y = -C02 - C12 - C22;
		A22 = C02*D20 + C12*D21 + C22*D22; d_A22_z = -C02 - C12 - C22;

		A2_00 = A00 *A00 + A10*A10 + A20 *A20;
		A2_01 = A00*A01 + A10*A11 + A20*A21;
		A2_02 = A00*A02 + A10*A12 + A20*A22;
		A2_11 = A01 *A01 + A11*A11 + A21 *A21;
		A2_12 = A01*A02 + A11*A12 + A21*A22;
		A2_22 = A02 *A02 + A12*A12 + A22 *A22;

		AF = A2_00 + A2_11 + A2_22; AF_05 = std::sqrt(AF); I_AF_05 = 0.5 / AF_05;
		AF2 = A2_00*A2_00 + A2_01*A2_01 * 2 + A2_02*A2_02 * 2 + A2_11*A2_11 + A2_12*A2_12 * 2 + A2_22*A2_22;
		AF_I = (AF*AF - AF2)*0.5; AF_I_05 = std::sqrt(AF_I); I_AF_I_05 = 0.5 / AF_I_05;
		det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
		i_det_A = 1.0 / det_A;

		D_AF_x = (2.0*A00*d_A00_x + 2.0*A01*d_A01_x + 2.0*A02*d_A02_x);
		D_AF_y = (2.0*A10*d_A10_y + 2.0*A11*d_A11_y + 2.0*A12*d_A12_y);
		D_AF_z = (2.0*A20*d_A20_z + 2.0*A21*d_A21_z + 2.0*A22*d_A22_z);

		D_AF2_x = (2.0*A2_00*2.0*A00*d_A00_x + 4.0*A2_01*(A01*d_A00_x + A00*d_A01_x) + 4.0*A2_02*(A02*d_A00_x + A00*d_A02_x)
			+ 2.0*A2_11*2.0*A01*d_A01_x + 2.0*A2_22*2.0*A02*d_A02_x + 4.0*A2_12*(A01*d_A02_x + A02*d_A01_x));

		D_AF2_y = (2.0*A2_00*2.0*A10*d_A10_y + 4.0*A2_01*(A11*d_A10_y + A10*d_A11_y) + 4.0*A2_02*(A12*d_A10_y + A10*d_A12_y)
			+ 2.0*A2_11*2.0*A11*d_A11_y + 2.0*A2_22*2.0*A12*d_A12_y + 4.0*A2_12*(A11*d_A12_y + A12*d_A11_y));

		D_AF2_z = (2.0*A2_00*2.0*A20*d_A20_z + 4.0*A2_01*(A21*d_A20_z + A20*d_A21_z) + 4.0*A2_02*(A22*d_A20_z + A20*d_A22_z)
			+ 2.0*A2_11*2.0*A21*d_A21_z + 2.0*A2_22*2.0*A22*d_A22_z + 4.0*A2_12*(A21*d_A22_z + A22*d_A21_z));

		D_AF_I_x = (AF*D_AF_x - 0.5*D_AF2_x)*I_AF_I_05; D_AF_I_y = (AF*D_AF_y - 0.5*D_AF2_y)*I_AF_I_05; D_AF_I_z = (AF*D_AF_z - 0.5*D_AF2_z)*I_AF_I_05;
		D_AF_x *= I_AF_05; D_AF_y *= I_AF_05; D_AF_z *= I_AF_05;

		D_det_A_x = d_A00_x * A11 * A22 + d_A01_x * A12 * A20 + d_A02_x * A10 * A21 - d_A02_x * A11 * A20 - d_A01_x * A10 * A22 - d_A00_x * A12 * A21;
		D_det_A_y = A00 * d_A11_y * A22 + A01 * d_A12_y * A20 + A02 * d_A10_y * A21 - A02 * d_A11_y * A20 - A01 * d_A10_y * A22 - A00 * d_A12_y * A21;
		D_det_A_z = A00 * A11 * d_A22_z + A01 * A12 * d_A20_z + A02 * A10 * d_A21_z - A02 * A11 * d_A20_z - A01 * A10 * d_A22_z - A00 * A12 * d_A21_z;

		inv_det_A_2_05 = 0.5 *(1.0 - i_det_A*i_det_A);
		dvex = D_det_A_x*inv_det_A_2_05;
		dvey = D_det_A_y*inv_det_A_2_05;
		dvez = D_det_A_z*inv_det_A_2_05;

		g = AF_05 * AF_I_05;
		dgx = D_AF_x * AF_I_05 + AF_05*D_AF_I_x;
		dgy = D_AF_y * AF_I_05 + AF_05*D_AF_I_y;
		dgz = D_AF_z * AF_I_05 + AF_05*D_AF_I_z;

		dex = (dgx *i_det_A - g * D_det_A_x *i_det_A*i_det_A);
		dey = (dgy *i_det_A - g * D_det_A_y *i_det_A*i_det_A);
		dez = (dgz *i_det_A - g * D_det_A_z *i_det_A*i_det_A);

		e = g*i_det_A;
		mips_e = (e*e - 1.0)*0.125;
		//mips_e = e/3.0;
		volume_e = 0.5*(det_A + i_det_A);
		exp_k = (alpha*mips_e + beta *volume_e);
		if (exp_k > 60) exp_k = 60;
		exp_vec[i] = exp_k;

		/ *exp_e = std::exp(exp_k);
		local_energy += exp_e;* /

		d_mips_e_x = e * 0.25 * dex;
		d_mips_e_y = e * 0.25 * dey;
		d_mips_e_z = e * 0.25 * dez;
		/ *d_mips_e_x = dex/3.0;
		d_mips_e_y = dey/3.0;
		d_mips_e_z = dez/3.0;* /

		gx_vec[i] = alpha * d_mips_e_x + beta * dvex;
		gy_vec[i] = alpha * d_mips_e_y + beta * dvey;
		gz_vec[i] = alpha * d_mips_e_z + beta * dvez;

		/ *dex = alpha * d_mips_e_x + beta * dvex;
		dey = alpha * d_mips_e_y + beta * dvey;
		dez = alpha * d_mips_e_z + beta * dvez;

		dex *= exp_e; dey *= exp_e; dez *= exp_e;

		gx += dex; gy += dey; gz += dez;* /
	}

#if NDEBUG
	fmath::expd_v(exp_vec, vc_size);
	for (int i = 0; i < vc_size; ++i)
	{
		exp_e = exp_vec[i];
		local_energy += exp_e;
		gx += gx_vec[i] * exp_e;
		gy += gy_vec[i] * exp_e;
		gz += gz_vec[i] * exp_e;
	}
#else
	for (int i = 0; i < vc_size; ++i)
	{
		exp_e = std::exp(exp_vec[i]);
		local_energy += exp_e;
		gx += gx_vec[i] * exp_e;
		gy += gy_vec[i] * exp_e;
		gz += gz_vec[i] * exp_e;
	}
#endif

	double mu = 1.0; int bvf_size = 0;
	double* p_vf_pos_x = NULL; double* p_vf_pos_y = NULL; double* p_vf_pos_z = NULL;
	OpenVolumeMesh::Geometry::Vec3i& v_type = vertex_type[v_id];
	if (bv_flag)
	{
		double normal_e = 0.0; double Ex = 1.0, Ey = 1.0, Ez = 1.0;
		double ne_gx = 0.0; double ne_gy = 0.0; double ne_gz = 0.0;
		p_vf_pos_x = vc_pos_x3_omp[omp_id].data();
		p_vf_pos_y = vc_pos_y3_omp[omp_id].data();
		p_vf_pos_z = vc_pos_z3_omp[omp_id].data();
		std::vector<int>& one_vf_id = bvf_id[v_id];
		bvf_size = one_vf_id.size();
		for (int i = 0; i < bvf_size; ++i)
		{
			OpenVolumeMesh::Geometry::Vec3d& tn = target_bfn[one_vf_id[i]];
			OpenVolumeMesh::Geometry::Vec3i& one_fv_id = bfv_id[one_vf_id[i]];
			for (int j = 0; j < 3; ++j)
			{
				if (one_fv_id[j] == v_id)
				{
					x1 = p_dpx[one_fv_id[(j + 1) % 3]]; y1 = p_dpy[one_fv_id[(j + 1) % 3]]; z1 = p_dpz[one_fv_id[(j + 1) % 3]];
					x2 = p_dpx[one_fv_id[(j + 2) % 3]]; y2 = p_dpy[one_fv_id[(j + 2) % 3]]; z2 = p_dpz[one_fv_id[(j + 2) % 3]];
					break;
				}
			}

			p_vf_pos_x[3 * i + 0] = x1; p_vf_pos_x[3 * i + 1] = x2; p_vf_pos_x[3 * i + 2] = tn[0];
			p_vf_pos_y[3 * i + 0] = y1; p_vf_pos_y[3 * i + 1] = y2; p_vf_pos_y[3 * i + 2] = tn[1];
			p_vf_pos_z[3 * i + 0] = z1; p_vf_pos_z[3 * i + 1] = z2; p_vf_pos_z[3 * i + 2] = tn[2];

			double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1);
			double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2);
			double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1);

			double len2 = nx*nx + ny*ny + nz*nz; double len = std::sqrt(len2); 
			double inv_len2 = 1.0 / len2; double inv_len = 1.0 / len;

			//double E_ne = 1.0 + 1.0e-8 - (nx*tn[0] + ny*tn[1] + nz*tn[2])/len;//polycube energy
			double E_ne = 0.5*((nx *inv_len - tn[0])*(nx *inv_len - tn[0]) + (ny *inv_len - tn[1])*(ny *inv_len - tn[1]) + (nz *inv_len - tn[2])*(nz *inv_len - tn[2]));//polycube energy
			normal_e += E_ne;

			double d_nx_x0 = 0.0;     double d_nx_y0 = z1 - z2; double d_nx_z0 = y2 - y1;
			double d_ny_x0 = z2 - z1; double d_ny_y0 = 0.0;     double d_ny_z0 = x1 - x2;
			double d_nz_x0 = y1 - y2; double d_nz_y0 = x2 - x1; double d_nz_z0 = 0.0;

			double d_len2_x0 = (nx*d_nx_x0 + ny*d_ny_x0 + nz*d_nz_x0); //miss 2
			double d_len2_y0 = (nx*d_nx_y0 + ny*d_ny_y0 + nz*d_nz_y0);
			double d_len2_z0 = (nx*d_nx_z0 + ny*d_ny_z0 + nz*d_nz_z0);

			double d_len_x0 = d_len2_x0 *inv_len; double d_len_y0 = d_len2_y0 *inv_len; double d_len_z0 = d_len2_z0*inv_len;

			double d_nx_len_x0 = d_nx_x0 *inv_len - nx*d_len_x0 *inv_len2;
			double d_nx_len_y0 = d_nx_y0 *inv_len - nx*d_len_y0 *inv_len2;
			double d_nx_len_z0 = d_nx_z0 *inv_len - nx*d_len_z0 *inv_len2;

			double d_ny_len_x0 = d_ny_x0 *inv_len - ny*d_len_x0 *inv_len2;
			double d_ny_len_y0 = d_ny_y0 *inv_len - ny*d_len_y0 *inv_len2;
			double d_ny_len_z0 = d_ny_z0 *inv_len - ny*d_len_z0 *inv_len2;

			double d_nz_len_x0 = d_nz_x0 *inv_len - nz*d_len_x0 *inv_len2;
			double d_nz_len_y0 = d_nz_y0 *inv_len - nz*d_len_y0 *inv_len2;
			double d_nz_len_z0 = d_nz_z0 *inv_len - nz*d_len_z0 *inv_len2;

			ne_gx += (d_nx_len_x0*(nx*inv_len - tn[0]) + d_ny_len_x0*(ny*inv_len - tn[1]) + d_nz_len_x0*(nz*inv_len - tn[2]));
			ne_gy += (d_nx_len_y0*(nx*inv_len - tn[0]) + d_ny_len_y0*(ny*inv_len - tn[1]) + d_nz_len_y0*(nz*inv_len - tn[2]));
			ne_gz += (d_nx_len_z0*(nx*inv_len - tn[0]) + d_ny_len_z0*(ny*inv_len - tn[1]) + d_nz_len_z0*(nz*inv_len - tn[2]));
		}

		double pos_e = 0.0;
		double pe_gx = 0.0; double pe_gy = 0.0; double pe_gz = 0.0;
		if (v_type[0] >= 0)
		{
			double diff_x = (x0 - chart_mean_value[v_type[0]]);
			pos_e += diff_x*diff_x;
			pe_gx = 2 * diff_x;
		}
		if (v_type[1] >= 0)
		{
			double diff_y = (y0 - chart_mean_value[v_type[1]]);
			pos_e += diff_y*diff_y;
			pe_gy = 2 * diff_y;
		}
		if (v_type[2] >= 0)
		{
			double diff_z = (z0 - chart_mean_value[v_type[2]]);
			pos_e += diff_z*diff_z;
			pe_gz = 2 * diff_z;
		}

		double s = local_energy / (normal_e + pos_e + 1e-8);
		if (s > 1e16) ga *= 1e16;
		else ga *= s;

		local_energy += ga * (normal_e + pos_e);
		//gx += ga*(ne_gx + pe_gx); gy += ga*(ne_gy + pe_gy); gz += ga*(ne_gz + pe_gz);
		gx += ga*(ne_gx); gy += ga*(ne_gy ); gz += ga*(ne_gz );
	}

	min_radius = std::sqrt(min_radius);
	double dx = gx, dy = gy, dz = gz;
	double move_d = sqrt(dx*dx + dy*dy + dz*dz);
	if (move_d > min_radius)
	{
		double t = min_radius / move_d;
		dx *= t; dy *= t; dz *= t;
		move_d = min_radius;
	}

	double npx = x0 - dx; double npy = y0 - dy; double npz = z0 - dz;
	while (!local_check_negative_volume4(vc_size, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, p_vc_n_cross_x, p_vc_n_cross_y, p_vc_n_cross_z, npx, npy, npz))
	{
		dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
		if (move_d < 1e-8)
		{
			npx = x0; npy = y0; npz = z0;
			move_d = 1e-8; break;
		}
		else
		{
			npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
		}
	}
#if 1
	double new_e = 0;
	//if (check_flag) new_e = 2.0*move_d*move_d;
	//bool e_ok = compute_exp_misp_energy_refine_polycube_position(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, v_type, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
	bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
	if (!e_ok)
	{
		while (!e_ok)
		{
			dx *= 0.2; dy *= 0.2; dz *= 0.2; move_d *= 0.2;
			if (move_d < 1e-8)
			{
				npx = x0; npy = y0; npz = z0;
				move_d = 1e-8; break;
			}
			else
			{
				npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
			}
			new_e = 0;
			//if (check_flag) new_e = 2.0*move_d*move_d;
			//e_ok = compute_exp_misp_energy_refine_polycube_position(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, v_type, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
			e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
		}

		if (e_ok)
		{
			dx *= 5; dy *= 5; dz *= 5; move_d *= 5; e_ok = false;
			while (!e_ok)
			{
				dx *= 0.5; dy *= 0.5; dz *= 0.5; move_d *= 0.5;
				if (move_d < 1e-8)
				{
					npx = x0; npy = y0; npz = z0;
					move_d = 1e-8; break;
				}
				else
				{
					npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
				}
				new_e = 0;
				//if (check_flag) new_e = 2.0*move_d*move_d;
				//e_ok = compute_exp_misp_energy_refine_polycube_position(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, v_type, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
				e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
			}

			if (e_ok)
			{
				dx *= 2; dy *= 2; dz *= 2; move_d *= 2; e_ok = false;
				while (!e_ok)
				{
					dx *= 0.875; dy *= 0.875; dz *= 0.875; move_d *= 0.875;
					if (move_d < 1e-8)
					{
						npx = x0; npy = y0; npz = z0;
						move_d = 1e-8; break;
					}
					else
					{
						npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
					}
					new_e = 0;
					//if (check_flag) new_e = 2.0*move_d*move_d;
					//e_ok = compute_exp_misp_energy_refine_polycube_position(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, v_type, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
					e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
				}
			}
		}
	}
#else
	double new_e = 0;
	bool e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
	while (!e_ok)
	{
		dx *= 0.8; dy *= 0.8; dz *= 0.8; move_d *= 0.8;
		if (move_d < 1e-8)
		{
			npx = x0; npy = y0; npz = z0;
			move_d = 1e-8; break;
		}
		else
		{
			npx = x0 - dx; npy = y0 - dy; npz = z0 - dz;
		}
		new_e = 0;
		e_ok = compute_exp_misp_energy_refine_polycube(vc_size, local_energy, new_e, p_vc_pos_x, p_vc_pos_y, p_vc_pos_z, vc_S, exp_vec, npx, npy, npz, alpha, beta, ga, bvf_size, p_vf_pos_x, p_vf_pos_y, p_vf_pos_z);
	}
#endif
	dpx[v_id] = npx; dpy[v_id] = npy; dpz[v_id] = npz;
}
*/

bool polycube_deformation_interface::local_check_negative_volume4(const int& vc_size, const double* posx, const double* posy, const double* posz,
	const double* vc_n_cross_x, const double* vc_n_cross_y, const double* vc_n_cross_z,
	const double& npx, const double& npy, const double& npz)
{
	int j = 0;
	for (unsigned i = 0; i < vc_size; ++i)
	{
		j = 3 * i;
		double c_v = vc_n_cross_x[i] * (npx - posx[j]) + vc_n_cross_y[i] * (npy - posy[j]) + vc_n_cross_z[i] * (npz - posz[j]);
		if (c_v < 1e-20) return false;
	}
	return true;
}

bool polycube_deformation_interface::compute_exp_misp_energy_refine_polycube(
	const int& vc_size, const double& old_e, double& new_e,
	const double* posx, const double* posy, const double* posz,
	const std::vector<std::vector<double> >& vc_S, double* exp_vec,
	const double& npx, const double& npy, const double& npz,
	double alpha, double beta,
	const double& ga, const int& vf_size,
	const double* vf_px, const double* vf_py, const double* vf_pz, const double* mu)
{
	//new_e = 0.0; 
	double D00, D10, D20, D01, D11, D21, D02, D12, D22;
	double C00, C01, C02, C10, C11, C12, C20, C21, C22;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
	double AF, AF2, AF_I, e2, mips_e, volume_e, k, exp_e;
	for (unsigned i = 0; i < vc_size; ++i)
	{
		const double* s_data = vc_S[i].data();
		int j = 3 * i;
		D00 = posx[j + 0] - npx; D10 = posy[j + 0] - npy; D20 = posz[j + 0] - npz;
		D01 = posx[j + 1] - npx; D11 = posy[j + 1] - npy; D21 = posz[j + 1] - npz;
		D02 = posx[j + 2] - npx; D12 = posy[j + 2] - npy; D22 = posz[j + 2] - npz;

		C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
		C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
		C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];

		A00 = C00*D00 + C10*D01 + C20*D02;
		A10 = C00*D10 + C10*D11 + C20*D12;
		A20 = C00*D20 + C10*D21 + C20*D22;
		A01 = C01*D00 + C11*D01 + C21*D02;
		A11 = C01*D10 + C11*D11 + C21*D12;
		A21 = C01*D20 + C11*D21 + C21*D22;
		A02 = C02*D00 + C12*D01 + C22*D02;
		A12 = C02*D10 + C12*D11 + C22*D12;
		A22 = C02*D20 + C12*D21 + C22*D22;

		det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;

		A2_00 = A00 *A00 + A10*A10 + A20 *A20;
		A2_01 = A00*A01 + A10*A11 + A20*A21;
		A2_02 = A00*A02 + A10*A12 + A20*A22;
		A2_11 = A01 *A01 + A11*A11 + A21 *A21;
		A2_12 = A01*A02 + A11*A12 + A21*A22;
		A2_22 = A02 *A02 + A12*A12 + A22 *A22;
		i_det_A = 1.0 / det_A;

		AF = A2_00 + A2_11 + A2_22;
		AF2 = A2_00*A2_00 + A2_01*A2_01 * 2 + A2_02*A2_02 * 2 + A2_11*A2_11 + A2_12*A2_12 * 2 + A2_22*A2_22;
		AF_I = (AF*AF - AF2)*0.5;

		e2 = AF*AF_I * (i_det_A*i_det_A);
		mips_e = (e2 - 1.0)*0.125;
		//mips_e = std::sqrt(e2) / 3.0;
		volume_e = 0.5*(det_A + i_det_A);

		k = alpha*mips_e + beta *volume_e;
		if (k > 60) k = 60;
		exp_vec[i] = k;
		/*exp_e = std::exp(k);
		new_e += exp_e;
		if (new_e > old_e) return false;*/
	}

#if NDEBUG
	fmath::expd_v(exp_vec, vc_size);
	for (int i = 0; i < vc_size; ++i)
	{
		new_e += exp_vec[i];
	}
#else
	for (int i = 0; i < vc_size; ++i)
	{
		exp_e = std::exp(exp_vec[i]);
		new_e += exp_e;
	}
#endif
	if (new_e > old_e) return false;

	double tnx, tny, tnz; double x1, x2, y1, y2, z1, z2,x3,y3,z3;
	for (int i = 0; i < vf_size; ++i)
	{
		int j = 4 * i;

		x1 = vf_px[j + 0]; x2 = vf_px[j + 1]; x3 = vf_px[j + 2]; tnx = vf_px[j + 3];
		y1 = vf_py[j + 0]; y2 = vf_py[j + 1]; y3 = vf_py[j + 2]; tny = vf_py[j + 3];
		z1 = vf_pz[j + 0]; z2 = vf_pz[j + 1]; z3 = vf_pz[j + 2]; tnz = vf_pz[j + 3];

		double nx = (npy - y1)*(npz - z2) - (npy - y2)*(npz - z1);
		double ny = (npx - x2)*(npz - z1) - (npx - x1)*(npz - z2);
		double nz = (npx - x1)*(npy - y2) - (npx - x2)*(npy - y1);
		double len = std::sqrt(nx*nx + ny*ny + nz*nz); double inv_len = 1 / len;
		double nx_d = nx*inv_len - tnx; double ny_d = ny*inv_len - tny; double nz_d = nz*inv_len - tnz;
		double E_ne = nx_d*nx_d + ny_d*ny_d + nz_d*nz_d;

		//new_e += 0.5*E_ne*ga;
		new_e += 0.5*E_ne*mu[i];
		if (new_e > old_e) return false;

		//no smooth energy
		/*double nx_2 = (npy - y2)*(npz - z3) - (npy - y3)*(npz - z2);
		double ny_2 = (npx - x3)*(npz - z2) - (npx - x2)*(npz - z3);
		double nz_2 = (npx - x2)*(npy - y3) - (npx - x3)*(npy - y2);

		double len_2 = std::sqrt(nx_2*nx_2 + ny_2*ny_2 + nz_2*nz_2);
		nx_d = nx *inv_len - nx_2 / len_2; ny_d = ny *inv_len - ny_2 / len_2; nz_d = nz *inv_len - nz_2 / len_2;
		double E_smooth = nx_d*nx_d + ny_d*ny_d + nz_d*nz_d; //smooth energy

		new_e += E_smooth * 5;
		if (new_e > old_e) return false;*/
	}

	return true;
}

bool polycube_deformation_interface::compute_exp_misp_energy_refine_polycube_normal(
	const int& vc_size, const double& old_e, double& new_e,
	const double* posx, const double* posy, const double* posz,
	const std::vector<std::vector<double> >& vc_S, double* exp_vec,
	const double& x0, const double& y0, const double& z0,
	double alpha, double beta,
	const double& mu1, const double& mu2, const int& vv_size,
	const double* vv_px, const double* vv_py, const double* vv_pz)
{
	//new_e = 0.0; 
	double D00, D10, D20, D01, D11, D21, D02, D12, D22;
	double C00, C01, C02, C10, C11, C12, C20, C21, C22;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
	double AF, AF2, AF_I, e2, mips_e, volume_e, k, exp_e;
	for (unsigned i = 0; i < vc_size; ++i)
	{
		const double* s_data = vc_S[i].data();
		int j = 3 * i;
		D00 = posx[j + 0] - x0; D10 = posy[j + 0] - y0; D20 = posz[j + 0] - z0;
		D01 = posx[j + 1] - x0; D11 = posy[j + 1] - y0; D21 = posz[j + 1] - z0;
		D02 = posx[j + 2] - x0; D12 = posy[j + 2] - y0; D22 = posz[j + 2] - z0;

		C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
		C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
		C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];

		A00 = C00*D00 + C10*D01 + C20*D02;
		A10 = C00*D10 + C10*D11 + C20*D12;
		A20 = C00*D20 + C10*D21 + C20*D22;
		A01 = C01*D00 + C11*D01 + C21*D02;
		A11 = C01*D10 + C11*D11 + C21*D12;
		A21 = C01*D20 + C11*D21 + C21*D22;
		A02 = C02*D00 + C12*D01 + C22*D02;
		A12 = C02*D10 + C12*D11 + C22*D12;
		A22 = C02*D20 + C12*D21 + C22*D22;

		det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
		i_det_A = 1.0 / det_A;

		A2_00 = A00 *A00 + A10*A10 + A20 *A20;
		A2_01 = A00*A01 + A10*A11 + A20*A21;
		A2_02 = A00*A02 + A10*A12 + A20*A22;
		A2_11 = A01 *A01 + A11*A11 + A21 *A21;
		A2_12 = A01*A02 + A11*A12 + A21*A22;
		A2_22 = A02 *A02 + A12*A12 + A22 *A22;

		AF = A2_00 + A2_11 + A2_22;
		AF2 = A2_00*A2_00 + A2_01*A2_01 * 2 + A2_02*A2_02 * 2 + A2_11*A2_11 + A2_12*A2_12 * 2 + A2_22*A2_22;
		AF_I = (AF*AF - AF2)*0.5;

		e2 = AF*AF_I * (i_det_A*i_det_A);
		mips_e = (e2 - 1.0)*0.125;
		//mips_e = std::sqrt(e2) / 3.0;
		volume_e = 0.5*(det_A + i_det_A);

		k = alpha*mips_e + beta *volume_e;
		if (k > 60) k = 60;
		exp_vec[i] = k;
		/*exp_e = std::exp(k);
		new_e += exp_e;
		if (new_e > old_e) return false;*/
	}

#if NDEBUG
	fmath::expd_v(exp_vec, vc_size);
	for (int i = 0; i < vc_size; ++i)
	{
		new_e += exp_vec[i];
	}
#else
	for (int i = 0; i < vc_size; ++i)
	{
		exp_e = std::exp(exp_vec[i]);
		new_e += exp_e;
	}
#endif
	if (new_e > old_e) return false;

	double x1, x2, x3, y1, y2, y3, z1, z2, z3; double ga = 0.001;
	for (int i = 0; i < vv_size; ++i)
	{
		int j = 3 * i;

		x1 = vv_px[j + 0]; x2 = vv_px[j + 1]; x3 = vv_px[j + 2]; 
		y1 = vv_py[j + 0]; y2 = vv_py[j + 1]; y3 = vv_py[j + 2];
		z1 = vv_pz[j + 0]; z2 = vv_pz[j + 1]; z3 = vv_pz[j + 2];

		double nx = (y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1); double nx2 = nx*nx;
		double ny = (x0 - x2)*(z0 - z1) - (x0 - x1)*(z0 - z2); double ny2 = ny*ny;
		double nz = (x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1); double nz2 = nz*nz;

		double len2 = nx2 + ny2 + nz2; 
		double len4 = len2*len2;
		double len = std::sqrt(len2);
		//double inv_len_4 = 1.0 / len4;
		double inv_len_2 = 1.0 / len2;

		/*double pe = (nx2*ny2 + ny2*nz2 + nz2*nx2);
		double E_pe = pe*inv_len_4;//polycube energy*/
		double cx = std::sqrt(inv_len_2 *nx2 + ga); double cy = std::sqrt(inv_len_2 *ny2 + ga); double cz = std::sqrt(inv_len_2 *nz2 + ga);
		double E_pe = cx + cy + cz;//polycube energy
		new_e += E_pe*mu1;
		if (new_e > old_e) return false;

		//no smooth energy
		/*double nx_2 = (y0 - y2)*(z0 - z3) - (y0 - y3)*(z0 - z2);
		double ny_2 = (x0 - x3)*(z0 - z2) - (x0 - x2)*(z0 - z3);
		double nz_2 = (x0 - x2)*(y0 - y3) - (x0 - x3)*(y0 - y2);

		double len_2 = std::sqrt(nx_2*nx_2 + ny_2*ny_2 + nz_2*nz_2);
		double nx_d = nx / len - nx_2 / len_2; double ny_d = ny / len - ny_2 / len_2; double nz_d = nz / len - nz_2 / len_2;
		double E_smooth = nx_d*nx_d + ny_d*ny_d + nz_d*nz_d; //smooth energy
		
		new_e += E_smooth*mu2;*/
		if (new_e > old_e) return false;
	}

	return true;
}

bool polycube_deformation_interface::compute_exp_misp_energy_refine_polycube_position(
	const int& vc_size, const double& old_e, double& new_e,
	const double* posx, const double* posy, const double* posz,
	const std::vector<std::vector<double> >& vc_S, double* exp_vec,
	const double& npx, const double& npy, const double& npz,
	double alpha, double beta,
	const double& ga, const int& vf_size, OpenVolumeMesh::Geometry::Vec3i& v_type,
	const double* vf_px, const double* vf_py, const double* vf_pz)
{
	//new_e = 0.0; 
	double D00, D10, D20, D01, D11, D21, D02, D12, D22;
	double C00, C01, C02, C10, C11, C12, C20, C21, C22;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	double det_A, i_det_A, A2_00, A2_01, A2_02, A2_11, A2_12, A2_22;
	double AF, AF2, AF_I, e2, mips_e, volume_e, k, exp_e;
	for (unsigned i = 0; i < vc_size; ++i)
	{
		const double* s_data = vc_S[i].data();
		int j = 3 * i;
		D00 = posx[j + 0] - npx; D10 = posy[j + 0] - npy; D20 = posz[j + 0] - npz;
		D01 = posx[j + 1] - npx; D11 = posy[j + 1] - npy; D21 = posz[j + 1] - npz;
		D02 = posx[j + 2] - npx; D12 = posy[j + 2] - npy; D22 = posz[j + 2] - npz;

		C00 = s_data[0]; C01 = s_data[1]; C02 = s_data[2];
		C10 = s_data[3]; C11 = s_data[4]; C12 = s_data[5];
		C20 = s_data[6]; C21 = s_data[7]; C22 = s_data[8];

		A00 = C00*D00 + C10*D01 + C20*D02;
		A10 = C00*D10 + C10*D11 + C20*D12;
		A20 = C00*D20 + C10*D21 + C20*D22;
		A01 = C01*D00 + C11*D01 + C21*D02;
		A11 = C01*D10 + C11*D11 + C21*D12;
		A21 = C01*D20 + C11*D21 + C21*D22;
		A02 = C02*D00 + C12*D01 + C22*D02;
		A12 = C02*D10 + C12*D11 + C22*D12;
		A22 = C02*D20 + C12*D21 + C22*D22;

		det_A = A00 * A11 * A22 + A01 * A12 * A20 + A02 * A10 * A21 - A02 * A11 * A20 - A01 * A10 * A22 - A00 * A12 * A21;
		i_det_A = 1.0 / det_A;

		A2_00 = A00 *A00 + A10*A10 + A20 *A20;
		A2_01 = A00*A01 + A10*A11 + A20*A21;
		A2_02 = A00*A02 + A10*A12 + A20*A22;
		A2_11 = A01 *A01 + A11*A11 + A21 *A21;
		A2_12 = A01*A02 + A11*A12 + A21*A22;
		A2_22 = A02 *A02 + A12*A12 + A22 *A22;

		AF = A2_00 + A2_11 + A2_22;
		AF2 = A2_00*A2_00 + A2_01*A2_01 * 2 + A2_02*A2_02 * 2 + A2_11*A2_11 + A2_12*A2_12 * 2 + A2_22*A2_22;
		AF_I = (AF*AF - AF2)*0.5;

		e2 = AF*AF_I * (i_det_A*i_det_A);
		mips_e = (e2 - 1.0)*0.125;
		//mips_e = std::sqrt(e2) / 3.0;
		volume_e = 0.5*(det_A + i_det_A);

		k = alpha*mips_e + beta *volume_e;
		if (k > 60) k = 60;
		exp_vec[i] = k;
		/*exp_e = std::exp(k);
		new_e += exp_e;
		if (new_e > old_e) return false;*/
	}

#if NDEBUG
	fmath::expd_v(exp_vec, vc_size);
	for (int i = 0; i < vc_size; ++i)
	{
		new_e += exp_vec[i];
	}
#else
	for (int i = 0; i < vc_size; ++i)
	{
		exp_e = std::exp(exp_vec[i]);
		new_e += exp_e;
	}
#endif
	if (new_e > old_e) return false;

	double tnx, tny, tnz, ne; double x1, x2, y1, y2, z1, z2;
	for (int i = 0; i < vf_size; ++i)
	{
		int j = 3 * i;

		x1 = vf_px[j + 0]; x2 = vf_px[j + 1]; tnx = vf_px[j + 2];
		y1 = vf_py[j + 0]; y2 = vf_py[j + 1]; tny = vf_py[j + 2];
		z1 = vf_pz[j + 0]; z2 = vf_pz[j + 1]; tnz = vf_pz[j + 2];

		double nx = (npy - y1)*(npz - z2) - (npy - y2)*(npz - z1);
		double ny = (npx - x2)*(npz - z1) - (npx - x1)*(npz - z2);
		double nz = (npx - x1)*(npy - y2) - (npx - x2)*(npy - y1);
		double len = std::sqrt(nx*nx + ny*ny + nz*nz);

		double pos_e = 0.0;
		if (v_type[0] >= 0)
		{
			double diff_x = (npx - chart_mean_value[v_type[0]]);
			pos_e += diff_x*diff_x;
		}
		if (v_type[1] >= 0)
		{
			double diff_y = (npy - chart_mean_value[v_type[1]]);
			pos_e += diff_y*diff_y;
		}
		if (v_type[2] >= 0)
		{
			double diff_z = (npz - chart_mean_value[v_type[2]]);
			pos_e += diff_z*diff_z;
		}

		double inv_len = 1.0 / len;
		double E_ne = 0.5*((nx *inv_len - tnx)*(nx*inv_len - tnx) + (ny *inv_len - tny)*(ny*inv_len - tny) + (nz *inv_len - tnz)*(nz *inv_len - tnz));

		new_e += (pos_e + E_ne)*ga;
		if (new_e > old_e) return false;
	}

	return true;
}

void polycube_deformation_interface::load_deformation_result(const char* filename, VolumeMesh* mesh_)
{
	prepare_for_deformation(mesh_);

	FILE* f_der = fopen(filename, "r");
	char buf[4096];
	char x[128]; char y[128]; char z[128];
	int v_count = 0; int nv = mesh_->n_vertices();
	while (!feof(f_der))
	{
		fgets(buf, 4096, f_der);
		sscanf(buf, "%s %s %s", x, y, z);
		
		OpenVolumeMesh::VertexHandle vh(v_count);
		OpenVolumeMesh::Geometry::Vec3d np = OpenVolumeMesh::Geometry::Vec3d(atof(x), atof(y), atof(z));

		mesh_->set_vertex(vh, np);
		dpx[v_count] = np[0]; dpy[v_count] = np[1]; dpz[v_count] = np[2]; 
		++v_count;
		if (v_count == nv) break;
	}

	fclose(f_der);

	compute_distortion(mesh_);
}

void polycube_deformation_interface::save_deformation_result(const char* filename, VolumeMesh* mesh_)
{
#if 1
	FILE* f_der = fopen(filename, "w");
	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
		fprintf(f_der, "%.19f %.19f %.19f\n", p[0], p[1], p[2]);
	}
	fclose(f_der);
#else
	FILE* f_der = fopen(filename, "w");
	for (OpenVolumeMesh::VertexIter v_it = mesh_->vertices_begin(); v_it != mesh_->vertices_end(); ++v_it)
	{
		OpenVolumeMesh::Geometry::Vec3d p = mesh_->vertex(*v_it);
		fprintf(f_der, "%f %f %f %d\n", p[0], p[1], p[2], is_deformation_constraint[v_it->idx()]);
	}
	fclose(f_der);
#endif
	
}
