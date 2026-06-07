#include "Labeling.h"
#include "..\..\BasicMath\ConvexQuadOptimization.h"
#include<vector>
#include<map>
#include<iostream>
//#define TEST_CHART
//#define TEST_EDGE
//#define TEST_ALLEXTREMAS
//#define TEST_SURFACE
//#define SOLVE_CONFLICT

using std::vector;
using std::map;
using std::pair;
using namespace std;
Labeling::Labeling()
{

}
void Labeling::prepare_data(VolumeMesh*mesh_)
{
	//VolumeMeshViewer = Viewer;
	mesh = mesh_;
	small_chart_threshold = 10;
	qpsigma = 0.8;
	smooth_iteration = 5;
	area_threshold = 5;
	CorrespondSurVol();
	GetIncidence();
}


Labeling::~Labeling()
{
}

bool Labeling::LabelingMain(double area_threshold_ , int small_chart_threshold_)
{
	//CorrespondSurVol();
	//GetIncidence();
	small_chart_threshold = small_chart_threshold_; area_threshold = area_threshold_;
	if(!GetPolyCubeStructure())return false;
	GetTetLabel();
	GetEdgeChartTag();
	//if(!GetChartCoor())return false;
	return true;
	
}
void Labeling::GetEdgeChartTag()
{
	EdgeBoundaryTag.clear(); EdgeBoundaryTag.resize(mesh->n_edges(),0);
	for (size_t i = 0; i < charts_boundary.size(); ++i)
	{
		for (size_t j = 0; j < charts_boundary[i].size(); ++j)
		{
			EdgeBoundaryTag[charts_boundary[i][j].idx()] = 1;
		}
	}
}
//bool Labeling::checkStatus()
//{
//	swallowed_chart = true;
//	get_charts_neighbor_();
//	for (size_t i = 0; i < chart_neighbor_.size(); ++i)
//	{
//		if (chart_neighbor_[i].size() < 4)
//		{
//			swallowed_chart = false;
//			printf("chart %d, neighbors %d \n", i, chart_neighbor_[i].size());
//		}
//	}
//	if (swallowed_chart && filtered_extrema)return true;
//	else
//		return false;
//}
int Labeling::s_vertex_id(VertexHandle v)
{
	return v2s_vertexid[v.idx()];
}
int Labeling::s_edge_id(EdgeHandle e)
{
	map<int, int>::iterator m_it = mapv2s_edge.find(e.idx());
	if (m_it == mapv2s_edge.end())printf("sedgeid error!\n");
	return m_it->second;
}
int Labeling::s_face_id(FaceHandle f)
{
	return v2s_face_id[f.idx()];
}
std::vector<FaceHandle> Labeling::f_3faces(FaceHandle f)
{
	HalfFaceHandle hf = mesh->halfface_handle(f, 0);
	std::vector<FaceHandle>faces;
	std::vector<HalfEdgeHandle> hes;
	mesh->get_halfedges_from_halfface(hf, hes);
	for (size_t i = 0; i < hes.size(); i++)
	{
		EdgeHandle edge_ = mesh->edge_handle(hes[i]);
		int s_edgeid = s_edge_id(edge_);
		assert(s_edgeid != -1);
		assert(s_face_id(b_edge2faces_[s_edgeid][0]) != -1);
		assert(s_face_id(b_edge2faces_[s_edgeid][1]) != -1);

		if (b_edge2faces_[s_edgeid][0] == f)
		{
			faces.push_back(b_edge2faces_[s_edgeid][1]);
		}
		else
		{
			faces.push_back(b_edge2faces_[s_edgeid][0]);
		}
	}
	return faces;
}
Geometry::Vec3d Labeling::calculateNormals(std::vector<Geometry::Vec3d>one_face)
{
	OpenVolumeMesh::Geometry::Vec3d normal(0.0, 0.0, 0.0);
	OpenVolumeMesh::Geometry::Vec3d v1; OpenVolumeMesh::Geometry::Vec3d v2;
	v1 = one_face[0] - one_face[1];
	v2 = one_face[2] - one_face[1];
	normal = (v2 % v1);
	normal /= normal.norm();
	return normal;
}
void Labeling:: updateFaceNormals()
{
	s_face_normals.clear(); s_face_normals.resize(n_boundary_faces);

	for (size_t i = 0; i < s2v_face_.size(); ++i)
	{
		FaceHandle f = s2v_face_[i];
		HalfFaceHandle hf1 = mesh->halfface_handle(f,0);
		HalfFaceHandle hf2 = mesh->halfface_handle(f, 1);
		vector<Geometry::Vec3d>one_face;
		if (mesh->incident_cell(hf1).idx() >= 0)hf1 = hf2;
		assert(mesh->incident_cell(hf1).idx() < 0);
		OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh->hfv_iter(hf1);
		for (hfv_it; hfv_it.valid(); ++hfv_it)
		{
			one_face.push_back(mesh->vertex(*hfv_it));
		}
		s_face_normals[i] = calculateNormals(one_face);
	}
}
void Labeling::CorrespondSurVol()
{
	s2v_face_.clear();
	boundaryMesh.clear();
	s2v_vertex.clear();
	v2s_vertexid.clear();
	v2s_vertexid.assign(mesh->n_vertices(), -1);
	v2s_face_id.clear();
	v2s_face_id.assign(mesh->n_faces(), -1);
	b_face2edges.clear();
	maps2v_edge.clear();
	mapv2s_edge.clear();
	b_face2cell.clear();
	int b_edge_count = 0;
	int b_vertex_count = 0;
	int b_face_count = 0;

	
	for (OpenVolumeMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{
		if (!mesh->is_boundary(*f_it)) continue;

		OpenVolumeMesh::HalfFaceHandle hfHandle = mesh->halfface_handle(*f_it, 0);
		OpenVolumeMesh::HalfFaceHandle hfHandle2 = mesh->halfface_handle(*f_it, 1);
		int cell_id = mesh->incident_cell(hfHandle2).idx();
		if (cell_id < 0)
		{
			hfHandle2 = hfHandle;
			cell_id = mesh->incident_cell(hfHandle2).idx();
			hfHandle = mesh->halfface_handle(*f_it, 1);
		}
		if (mesh->incident_cell(hfHandle).idx() < 0)
		{
			hfHandle = hfHandle2;
		}
		assert(mesh->incident_cell(hfHandle).idx()>=0);
		b_face2cell.push_back(mesh->incident_cell(hfHandle));
		
		HalfFaceHandle hf = mesh->halfface_handle(*f_it, 0);
		std::vector<HalfEdgeHandle> hes;
		std::vector<EdgeHandle>es;
		std::vector<VertexHandle>vs;
		mesh->get_halfedges_from_halfface(hf, hes);
		for (size_t i = 0; i < 3; ++i)
		{
			es.push_back(mesh->edge_handle(hes[i]));
			bool is_visited_to = false;
			bool is_visited_from = false;
			OpenVolumeMeshEdge e = mesh->edge(mesh->edge_handle(hes[i]));
			for (size_t j = 0; j < vs.size(); ++j)
			{
				if (vs[j] == e.to_vertex())is_visited_to = true;
				if (vs[j] == e.from_vertex())is_visited_from = true;
			}
			if (!is_visited_from)vs.push_back(e.from_vertex());
			if (!is_visited_to)vs.push_back(e.to_vertex());
		}
		if (vs[0] == vs[1] || vs[1] == vs[2] || vs[2] == vs[0])printf("he->vs error\n");
		b_face2edges.push_back(es);
		b_face2vertices.push_back(vs);
		//construct s2v_facehandle
		FaceHandle f = mesh->face_handle(hfHandle);
		s2v_face_.push_back(f);
		v2s_face_id[f.idx()] = b_face_count;
		b_face_count++;

		std::vector<HalfEdgeHandle>hf_edges;
		mesh->get_halfedges_from_halfface(hfHandle, hf_edges);
		for (size_t i = 0; i < hf_edges.size(); ++i)
		{
			EdgeHandle edge = mesh->edge_handle(hf_edges[i]);
			map<int, int>::iterator m_it = mapv2s_edge.find(edge.idx());
			if (m_it == mapv2s_edge.end())
			{
				//edge isn't be vitised
				mapv2s_edge.insert(pair<int, int>(edge.idx(), b_edge_count));
				maps2v_edge.insert(pair<int, EdgeHandle>(b_edge_count, edge));
				b_edge_count++;
			}
		}
		//end
		OpenVolumeMesh::HalfFaceVertexIter hfv_it = mesh->hfv_iter(hfHandle);
		for (hfv_it; hfv_it; ++hfv_it)
		{
			if (-1 == v2s_vertexid[hfv_it->idx()])
			{
				v2s_vertexid[hfv_it->idx()] = b_vertex_count;
				b_vertex_count++;
				s2v_vertex.push_back(*hfv_it);
			}
			//end of construct
		}
	}
	n_boundary_vertices = s2v_vertex.size();
	n_boundary_faces = s2v_face_.size();
	n_boundary_edges = maps2v_edge.size();
	//printf("construct boundary mesh! F %d E %d V %d\n", n_boundary_faces, n_boundary_edges, n_boundary_vertices);
}
void Labeling::GetIncidence()
{
	//v->faces
	b_vertex2facesid.resize(n_boundary_vertices);
	for (size_t i = 0; i < b_vertex2facesid.size(); ++i)b_vertex2facesid[i].clear();
	for (size_t i = 0; i < s2v_face_.size(); ++i)
	{
		VertexHandle v;
		for (size_t j = 0; j < 3; ++j)
		{
			v = b_face2vertices[i][j];
			int s_vid = s_vertex_id(v);
			b_vertex2facesid[s_vid].push_back(i);
		}
	}
	//edge->faces
	b_edge2faces_.resize(n_boundary_edges);
	for (size_t i = 0; i < b_edge2faces_.size(); ++i)b_edge2faces_[i].clear();
	for (size_t i = 0; i < s2v_face_.size(); ++i)
	{
		EdgeHandle e;
		for (size_t j = 0; j < 3; ++j)
		{
			e = b_face2edges[i][j];
			int s_eid = s_edge_id(e);
			b_edge2faces_[s_eid].push_back(s2v_face_[i]);
		}
	}
#ifdef TEST_SURFACE
	VolumeMeshViewer->specialFaces.clear();
	int bfcount = 0;
	for (size_t i = 0; i < v2s_face_id.size(); ++i)
	{
		if (v2s_face_id[i] != -1)
		{
			++bfcount;
			int s_id = v2s_face_id[i];
			VolumeMeshViewer->specialFaces.push_back(s2v_face_[s_id]);
		}
	}
	printf("%d\n", bfcount);
	VolumeMeshViewer->special_points.clear();
	for (size_t i = 0; i < s2v_vertex.size(); ++i)
	{
		VolumeMeshViewer->special_points.push_back(s2v_vertex[i]);
	}
	VolumeMeshViewer->updateGL();
	VolumeMeshViewer->specialEdges.clear();
	std::map<int, EdgeHandle>::iterator mit = maps2v_edge.begin();
	for (mit; mit != maps2v_edge.end();++mit)
	{
		VolumeMeshViewer->specialEdges.push_back(mit->second);
	}
#endif
}
bool Labeling::GetPolyCubeStructure()
{ 
	charts.clear();
	charts_boundary.clear();
	charts_boundary_start_edges.clear();
	chart_bound_faceid.clear();
	chart_neighbor_.clear();
	
	compute_average_bface_area();
	get_s_face_axis();
	/*GetAxisAlignTetLabel();
	GetInitBoundaryEdge();*/

	straighten_bound(false);
	//printf("straighten!\n");
	get_charts();
	get_charts_neighbor_();
	if (!solve_chart_neighbor())return false;;
#ifdef TEST_CHART
	printf("charts number :%d\n", charts.size());
	VolumeMeshViewer->FaceLabel.clear();
	VolumeMeshViewer->FaceLabel.resize(n_boundary_faces);
	for (size_t i = 0; i < charts.size(); ++i)
	{
		int label = chart_axis(i);
		for (size_t j = 0; j < charts[i].size(); ++j)
		{
			VolumeMeshViewer->FaceLabel[s_face_id(charts[i][j])] = label;
		}
	}
#endif
	get_charts_bound();
	GetVChartid();
	ChangeCEfaceLabel();
	//ChangeCornerFaceLabel();
	if(!extract_vbound())return false;
	TurningPoints(smooth_iteration);
	if (!filt_weak_extrema(area_threshold))return false;
#ifdef TEST_CHART
	VolumeMeshViewer->specialEdges.clear();
	for (size_t i = 0; i < charts_boundary.size(); ++i)
	{
		for (size_t j = 0; j < charts_boundary[i].size(); ++j)
		{
			VolumeMeshViewer->specialEdges.push_back(charts_boundary[i][j]);
		}
	}
	VolumeMeshViewer->updateGL();
#endif
	return true;
}
bool Labeling::GetChartCoor()
{
	GetVChartid();
	if(!extract_vbound())return false;
	if(!GetEdgeStructure())return false;
	//rule 1
	EdgeChartsPair();
	//rule2
	ParallelEdgeChartsPair();
	//get triplet
	ChartPairDist.clear();
#ifdef TEST_EDGE
	//VolumeMeshViewer->special_points.clear();
	for (size_t i = 0; i < ChartsEdges.size(); ++i)
	{
		if (!(is_corner(ChartsEdges[i][0]) && is_corner(ChartsEdges[i][ChartsEdges[i].size() - 1])))printf("EdgeError!\n");
	}
#endif
	for (size_t i = 0; i < chartspair_dist.size(); ++i)
	{
		for (size_t j = 0; j < chartspair_dist[i].size(); ++j)
		{
			if (chartspair_dist[i][j] != 10000)
			{
				//get triplets
				ChartPairDist.push_back(Eigen::Triplet<double>(i, j, chartspair_dist[i][j]));

				up_chart_id.push_back(i);
				down_chart_id.push_back(j);
				diff_up_down.push_back(chartspair_dist[i][j]);
			}
		}
	}
	InitChartMeanValue();
	compute_all_chart_value(qpsigma);

	//conflict average
	for (size_t i = 0; i < conflict_pair.size(); ++i)
	{
		int c1 = conflict_pair[i].first;
		int c2 = conflict_pair[i].second;
		double value = 0.5*(chart_mean_value[c1] + chart_mean_value[c2]);
		chart_mean_value[c1] = value;
		chart_mean_value[c2] = value;
	}
	return true;
}
bool Labeling::InitLabelMain()
{
	get_s_face_axis();
	GetAxisAlignTetLabel();
	GetInitBoundaryEdge();

	tetLabel = axis_align_tetLabel;
	EdgeBoundaryTag = Init_EdgeBoudnaryTag;

	return true;
}
//polycube structure funs
void Labeling::get_s_face_axis()
{
	s_face_axis.clear();
	updateFaceNormals();
	for (size_t i = 0; i < s_face_normals.size(); ++i)
	{
		s_face_axis.push_back(find_closest_axis(s_face_normals[i]));
	}
	//test
	//VolumeMeshViewer->FaceLabel.clear();
	
}
axis Labeling::find_closest_axis(Geometry::Vec3d normal)
{
	double max_ = 0;
	int max_id = 0;
	axis axe;
	for (int i = 0; i < 3; i++)
	{
		if (abs(normal[i])>max_)
		{
			max_ = abs(normal[i]);
			max_id = i;
		}
	}
	if (normal[max_id] > 0)
	{
		switch (max_id)
		{
		case 0:
			axe = px;
			break;
		case 1:
			axe = py;
			break;
		case 2:
			axe = pz;
			break;
		}
	}
	else
	{
		switch (max_id)
		{
		case 0:
			axe = nx;
			break;
		case 1:
			axe = ny;
			break;
		case 2:
			axe = nz;
			break;
		}
	}
	return axe;
}
//void Labeling::get_chartbound_faceid()
//{
//	chart_bound_faceid.clear();
//	for (size_t i = 0; i < s2v_face_.size(); ++i)
//	{
//		FaceHandle f = s2v_face_[i];
//		int s_fid = s_face_id(f);
//		assert(s_fid == i);
//		assert(s_fid >= 0 && s_fid< n_boundary_faces);
//		axis axe = s_face_axis[s_fid];
//		std::vector<FaceHandle> &faces = f_3faces(f);
//		assert(faces.size() == 3);
//		for (int k = 0; k < faces.size(); ++k)
//		{
//			int s_faceid_ = s_face_id(faces[k]);
//			assert(s_faceid_>= 0 && s_faceid_ < n_boundary_faces);
//			if (s_face_axis[s_faceid_] != axe)
//			{
//				chart_bound_faceid.push_back(s_fid);
//				break;
//			}
//		}
//	}
//	//printf("%d\n", chart_bound_faceid.size());
//}
void Labeling::delete_chart(int c_i)
{
	vector<vector<FaceHandle>>::iterator it = charts.begin() + c_i;
	charts.erase(it);
	//update face_chart_id
	for (size_t i = 0; i < face_chart_id.size(); ++i)
	{
		assert(face_chart_id[i] != c_i);
		if (face_chart_id[i]>c_i)
		{
			--face_chart_id[i];
		}
	}
}
void Labeling::straighten_bound(bool charts_flag)
{
	//find the faces which are on the bound of charts
	//get_chartbound_faceid();
	//straighten
	
	bool flag_straighten = true;
	int count = 0;
	while (1)
	{
		//straighten
		//get_chartbound_faceid();
		if (!flag_straighten)break;
		flag_straighten = false;
		for (size_t i = 0; i < s2v_face_.size(); ++i)
		//for (size_t i = 0; i < chart_bound_faceid.size(); ++i)
		{
			//FaceHandle f = s2v_face_[chart_bound_faceid[i]];
			FaceHandle f = s2v_face_[i];
			vector<int>axis_num;
			axis_num.resize(6);
			//int cs_faceid = chart_bound_faceid[i];
			axis axe0 = s_face_axis[i];
			axis axe1 = axe0;
			//axis_num[axe0]++;

			std::vector<FaceHandle> &cb_faces = f_3faces(f);
			for (size_t k = 0; k < 3; ++k)
			{
				axis axe = s_face_axis[s_face_id(cb_faces[k])];
				//int new_chartbound_faceid;
				//if (axe == axe0)
				//{
				//	new_chartbound_faceid = s_face_id(cb_faces[k]);
				//}
				axis_num[axe]++;
				if (axis_num[axe] == 2 && axe != axe0)
				{
					flag_straighten = true;
					axe1 = axe;
					s_face_axis[i] = axe1;
					if (charts_flag)
					{
						// we already have charts, add f to the chart and update sth.
						int c1_id = face_chart_id[s_face_id(cb_faces[k])];
						int c2_id = face_chart_id[i];
						charts[c1_id].push_back(f);
						face_chart_id[i] = c1_id;
						for (vector<FaceHandle>::iterator fit = charts[c2_id].begin(); fit != charts[c2_id].end(); ++fit)
						{
							if (*fit == f)
							{
								charts[c2_id].erase(fit);
								break;
							}
						}
						if (charts[c2_id].empty())
						{
							cout << "straighten_chart empty!" << endl;
							delete_chart(c2_id);
							get_charts_neighbor_();
						}

					}
					break;
				}
			}
		}
	}
}
void Labeling::get_charts()
{
	using namespace OpenVolumeMesh::Geometry;
	face_chart_id.clear();
	face_chart_id.resize(n_boundary_faces);
	std::vector<bool>visited;
	visited.assign(n_boundary_faces, false);
	//for (int i = 0; i < boundary_faces_num; ++i)visited[i] = false;
	vector<FaceHandle>::iterator f_it = s2v_face_.begin();
	int face_visited_num = 0;
	int chart_num = 0;
	for (f_it; f_it != s2v_face_.end(); ++f_it)
	{
		FaceHandle f = *f_it;
		int s_fid = s_face_id(f);
		if (visited[s_fid])continue;
		Vec3d start_normal = s_face_normals[s_fid];
		std::vector<FaceHandle> chart;
		//if (f_it->idx() < 0)cout << "error" << endl;
		//BFS
		std::list<FaceHandle>queue;
		queue.push_back(f);
		chart.push_back(f);
		face_chart_id[s_fid] = chart_num;
		face_visited_num++;
		visited[s_face_id(f)] = true;
		axis axe = s_face_axis[s_fid];
		while (!queue.empty())
		{
			f = queue.front();
			queue.pop_front();
			vector<FaceHandle>&faces = f_3faces(f);
			for (size_t k = 0; k < faces.size(); ++k)
			{
				int s_faceid = s_face_id(faces[k]);
				Vec3d normal = s_face_normals[s_faceid];
				double cos_angle = dot(start_normal, normal);
				if ((!visited[s_faceid]) && s_face_axis[s_faceid] == axe/* && cos_angle>cos(90.0/180.0*PI)*/)
				{
					chart.push_back(faces[k]);
					face_chart_id[s_faceid] = chart_num;
					queue.push_back(faces[k]);
					visited[s_faceid] = true;
					++face_visited_num;
				}
			}
			if (face_visited_num >= n_boundary_faces)break;
		}
		charts.push_back(chart);
		chart_num++;
	}
#ifdef TEST_CHART
	/*printf("charts number :%d\n", charts.size());
	VolumeMeshViewer->FaceLabel.clear();
	VolumeMeshViewer->FaceLabel.resize(n_boundary_faces);
	for (size_t i = 0; i < charts.size(); ++i)
	{
		int label = chart_axis(i) ;
		for (size_t j = 0; j < charts[i].size(); ++j)
		{
			VolumeMeshViewer->FaceLabel[s_face_id(charts[i][j])] = label;
		}
	}*/
	
#endif
}
axis Labeling::chart_axis(int c_i)
{
	return s_face_axis[s_face_id(charts[c_i][0])];
}
int Labeling::faceChartId(FaceHandle f)
{
	return face_chart_id[s_face_id(f)];
}
void Labeling::get_charts_neighbor_()
{
	vector<vector<bool>>visited; visited.resize(charts.size());
	chart_neighbor_.resize(charts.size());
	for (size_t i = 0; i < charts.size(); ++i)chart_neighbor_[i].clear();
	for (size_t i = 0; i < visited.size(); ++i)visited[i].resize(charts.size(), false);
	for (size_t i = 0; i < charts.size(); ++i)
	{
		for (size_t j = 0; j < charts[i].size(); ++j)
		{
			FaceHandle f = charts[i][j];
			vector<FaceHandle>&faces = f_3faces(f);
			for (size_t k = 0; k < 3; ++k)
			{
				int f_ci = face_chart_id[s_face_id(faces[k])];
				if (f_ci != i && (!visited[i][f_ci]))
				{
					visited[i][f_ci] = true;
					chart_neighbor_[i].push_back(f_ci);
				}
			}
		}
	}
	
}
bool Labeling::solve_chart_neighbor()
{
	solve_small_chart();//remove charts with 1 faces
	swallowed_chart = true;
	//1 nieghbor chart
	//printf("small chart thres %d\n", small_chart_threshold);
	while (1)
	{
		bool flag1 = false;
		for (size_t i = 0; i < chart_neighbor_.size(); ++i)
		{
			if (chart_neighbor_[i].size() == 1 && charts[i].size()<small_chart_threshold)
			{
				//for (size_t j = 0; j < charts[i].size(); ++j)SpecialFaces->push_back(charts[i][j]);
				del_1neighbor_chart(i);
				--i;
				flag1 = true;
			}
		}
		if ((!flag1))break;
	}
	////printf("%d %d\n", charts.size(), chart_neighbor_.size());
	//2-neighbor chart
	while (1)
	{
		bool flag2 = false;
		for (size_t i = 0; i < chart_neighbor_.size(); ++i)
		{
			if (chart_neighbor_[i].size() == 2 && charts[i].size()<small_chart_threshold)
			{
				//	cout << i << " 2nei " << endl;
				//for (size_t j = 0; j < charts[i].size(); ++j)SpecialFaces->push_back(charts[i][j]);
				del_2neighbor_chart(i, chart_neighbor_[i][0], chart_neighbor_[i][1]);
				--i;
				flag2 = true;
			}
		}
		//draw_charts();
		if ((!flag2))break;
	}
	//3 neighbor charts:swallowed by the biggest neighbor chart
	while (1)
	{
		bool flag3 = false;
		for (size_t i = 0; i < chart_neighbor_.size(); ++i)
		{
			if (chart_neighbor_[i].size() <4 && charts[i].size()<small_chart_threshold)
			{
				swallow_small_chart(i);
				--i;
				flag3 = true;
			}
		}
		if ((!flag3))break;
	}
	//after all swallow,some ajacent charts may have the same axe,merge them
	while (1)
	{
		bool flag4 = false;
		for (size_t i = 0; i < chart_neighbor_.size(); ++i)
		{
			axis c_axe = chart_axis(i);
			for (size_t j = 0; j < chart_neighbor_[i].size(); ++j)
			{
				int n_i = chart_neighbor_[i][j];
				if (c_axe == chart_axis(n_i))
				{
					merge_charts(i, n_i);
					--i;
					flag4 = true;
				}
			}
		}
		if ((!flag4))break;
	}
	for (size_t i = 0; i < chart_neighbor_.size(); ++i)
	{
		if (chart_neighbor_[i].size() < 4)
		{
			swallowed_chart = false;
			printf("chart %d, faces %d, neighbors %d \n", i,charts[i].size(), chart_neighbor_[i].size());
		}
	}
	if (!swallowed_chart)return false;
	straighten_bound(true);
	return true;
}
void Labeling::merge_charts(int c_i, int n_i)
{
	for (size_t i = 0; i < charts[c_i].size(); ++i)
	{
		update_axis_fcid_charts(charts[c_i][i], c_i, n_i);
	}
	delete_chart(c_i);
	//update neighbors
	get_charts_neighbor_();
}
void Labeling::update_axis_fcid_charts(FaceHandle f, int c_i, int new_ci)
//change f from chart c_i to chart new_ci
{
	axis axe = s_face_axis[s_face_id(charts[new_ci][0])];
	//update s_face_axis
	s_face_axis[s_face_id(f)] = axe;
	//update charts
	charts[new_ci].push_back(f);
	//update face_chart_id
	assert(new_ci != c_i);
	face_chart_id[s_face_id(f)] = new_ci;
}
void Labeling::swallow_small_chart(int c_i)
{
	int max_chart_num = 0;
	int max_ci = -1;
	for (size_t i = 0; i < chart_neighbor_[c_i].size(); ++i)
	{
		int n_ci = chart_neighbor_[c_i][i];
		if (charts[n_ci].size()>max_chart_num)
		{
			max_chart_num = charts[n_ci].size();
			max_ci = n_ci;
		}
	}
	//if (max_ci != -1)printf("find max neic\n");
	for (size_t i = 0; i < charts[c_i].size(); ++i)
	{
		update_axis_fcid_charts(charts[c_i][i], c_i, max_ci);
	}
	delete_chart(c_i);
	get_charts_neighbor_();
}
void Labeling::del_1neighbor_chart(int c_i)
{
	int neighbor_ci = chart_neighbor_[c_i][0];
	for (size_t j = 0; j < charts[c_i].size(); ++j)
	{
		update_axis_fcid_charts(charts[c_i][j], c_i, neighbor_ci);
	}
	delete_chart(c_i);
	//update neighbors
	get_charts_neighbor_();
}
void Labeling::del_2neighbor_chart(int c_i, int nei1, int nei2)
{
	if (charts[c_i].size() == 2)
	{
		FaceHandle f1 = charts[c_i][0];
		FaceHandle f2 = charts[c_i][1];
		int chart1, chart2;
		vector<FaceHandle>&faces_1 = f_3faces(f1);
		vector<FaceHandle>&faces_2 = f_3faces(f2);

		int e_num[2];
		int new_c;
		for (int k = 0; k < faces_1.size(); ++k)
		{
			int c1_id = face_chart_id[s_face_id(faces_1[k])];
			int c2_id = face_chart_id[s_face_id(faces_2[k])];
			if (c1_id == nei1)e_num[0]++;
			else
				e_num[1]++;
			if (c2_id == nei1)e_num[0]++;
			else
				e_num[1]++;
		}
		if (e_num[0] > e_num[1])
		{
			new_c = nei1;
		}
		else
		{
			new_c = nei2;
		}
		update_axis_fcid_charts(charts[c_i][0], c_i, new_c);
		update_axis_fcid_charts(charts[c_i][1], c_i, new_c);
		delete_chart(c_i);
		get_charts_neighbor_();
	
	}
	else
	{
		axis axe1 = chart_axis(nei1);
		axis axe2 = chart_axis(nei2);
	
		FaceHandle f1 = -1;
		FaceHandle f2 = -1;

		for (size_t i = 0; i < chart_bound_faceid.size(); ++i)
		{
			int s_faceid = chart_bound_faceid[i];
			//cout << endl << s_faceid;
			if (face_chart_id[s_faceid] == c_i)
				//find a chart bound of c_i
			{
				FaceHandle f = s2v_face_[s_faceid];
				vector<FaceHandle>&faces = f_3faces(f);
				assert(faces.size() == 3);
				for (int k = 0; k < faces.size(); ++k)
				{
					int chart_id = face_chart_id[s_face_id(faces[k])];
					if (chart_id == nei1 && f1 == -1)
					{
						
						f1 = f;
						break;
					}
					else if (chart_id == nei2 &&f2 == -1)
					{
						f2 = f;
						break;
					}
				}
			}
			if (f1 != -1 && f2 != -1)break;
		}
		if (f1 == -1)//there is only 1 face incident to chart c_i in chart nei1 
		{
			for (size_t i = 0; i < charts[c_i].size(); ++i)
			{
				FaceHandle f = charts[c_i][i];
				update_axis_fcid_charts(f, c_i, nei2);
			}
		}
		else if (f2 == -1)
		{
			for (size_t i = 0; i < charts[c_i].size(); ++i)
			{
				FaceHandle f = charts[c_i][i];
				update_axis_fcid_charts(f, c_i, nei1);
			}
		}
		else
		{
			//flood filling
			//test
			//cout << s_face_id(f1) << " " << s_face_id(f2) << endl;
			vector<int>fill_faces;// -1:hasn't been filled     nei1:add to neil  nei2:add to nei2
			fill_faces.assign(n_boundary_faces, -1);//initial hasn't been filled
			int count_1 = 0;
			int count_2 = 0;
			list<FaceHandle>queue1;
			queue1.push_back(f1);
			fill_faces[s_face_id(f1)] = nei1;
			count_1++;
			update_axis_fcid_charts(f1, c_i, nei1);
			list<FaceHandle>queue2;
			queue2.push_back(f2);
			fill_faces[s_face_id(f2)] = nei2;
			count_2++;
			update_axis_fcid_charts(f2, c_i, nei2);

			//	facefill(queue1, c_i, nei1, count_1, fill_faces);
			//facefill(queue2, c_i, nei2, count_2, fill_faces);
			assert(count_1 + count_2 <= charts[c_i].size());
			int count_while = 0;
			while (count_1 + count_2 < charts[c_i].size() && count_while<1000)
			{
				count_while++;
				assert(!(queue1.empty() && queue2.empty()));
				if (count_1>count_2)
				{
					if (!queue2.empty())
					{
						facefill(queue2, c_i, nei2, count_2, fill_faces);
						//update_axis_fcid_charts(f2, c_i, nei2);
					}
					if (!queue1.empty())
					{
						facefill(queue1, c_i, nei1, count_1, fill_faces);
						//update_axis_fcid_charts(f1, c_i, nei1);
					}
				}
				else
				{
					if (!queue1.empty())
					{
						facefill(queue1, c_i, nei1, count_1, fill_faces);
						//update_axis_fcid_charts(f1, c_i, nei1);
					}
					if (!queue2.empty())
					{
						facefill(queue2, c_i, nei2, count_2, fill_faces);
						//update_axis_fcid_charts(f2, c_i, nei2);
					}
				}
			}
		}
		//cout << count_while << endl;
		delete_chart(c_i);
		get_charts_neighbor_();
	}
}
void  Labeling::facefill(list<FaceHandle>&queue, int c_i, int nei, int &count, std::vector<int>&fill_faces)
//pop the front of the queue and fill his neighbors who haven't been visited
{
	FaceHandle f = queue.front();
	queue.pop_front();
	vector<FaceHandle>&faces = f_3faces(f);
	for (int k = 0; k < faces.size(); ++k)
	{
		int s_faceid = s_face_id(faces[k]);
		if (face_chart_id[s_faceid] == c_i &&fill_faces[s_faceid] == -1)
			//hasn't been filled
		{
			fill_faces[s_face_id(faces[k])] = nei;
			count++;
			update_axis_fcid_charts(faces[k], c_i, nei);
			queue.push_back(faces[k]);
		}
	}
}
void Labeling::solve_small_chart()
{
	bool flag = true;
	while (1)
	{
		if (!flag)break;
		flag = false;
		for (size_t i = 0; i < charts.size(); ++i)
		{
			vector<int>chart_num; chart_num.resize(charts.size());
			if (charts[i].size() == 1)
			{
				flag = true;
				FaceHandle f = charts[i][0];
				vector<FaceHandle >&faces = f_3faces(f);
				for (int k = 0; k < 3; ++k)
				{
					int c_id = face_chart_id[s_face_id(faces[k])];
					chart_num[c_id]++;
				}
				int max = 0;
				int cid = -1;
				for (int k = 0; k < charts.size(); ++k)
				{
					if (chart_num[k]>max)
					{
						max = chart_num[k];
						cid = k;
					}
				}
				update_axis_fcid_charts(f, i, cid);
				delete_chart(i);
				--i;
				get_charts_neighbor_();
			}
		}
	}
}
void Labeling::get_charts_bound()
{
	charts_boundary.clear();
	charts_boundary.resize(charts.size());
	charts_boundary_start_edges.clear();
	charts_boundary_start_edges.resize(charts.size());
	for (size_t i = 0; i < charts_boundary.size(); ++i)
	{
		get_chart_bound(i);
		//	cout << charts_boundary[i].size() << " ";
	}
	cout << endl;
}
void Labeling::get_chart_bound(int c_i)
{
	std::vector<FaceHandle>chart = charts[c_i];
	std::vector<EdgeHandle> &chart_bound = charts_boundary[c_i];
	std::vector<EdgeHandle>chart_bound_noorder;
	//vector<bool>visited_neighbor;
	//visited_neighbor.resize(charts.size());
	int c_ = 0;
	axis axe = s_face_axis[s_face_id(chart[0])];
	for (size_t j = 0; j < chart.size(); ++j)
	{
		FaceHandle f = chart[j];
		int s_faceid = s_face_id(f);
		std::vector<FaceHandle> &faces = f_3faces(f);
		assert(faces.size() == 3);
		for (int k = 0; k < 3; ++k)
		{
			if (s_face_axis[s_face_id(faces[k])] != axe)
			{
				c_++;
				EdgeHandle e = b_face2edges[s_faceid][k];
				chart_bound_noorder.push_back(e);
			}
		}
	}
	//reorder the boundary edges
	chart_bound.push_back(chart_bound_noorder[0]);
	charts_boundary_start_edges[c_i].push_back(0);
	vector<EdgeHandle>::iterator it = chart_bound_noorder.begin();
	chart_bound_noorder.erase(it);
	assert(chart_bound.size() == 1);
	//HalfEdgeHandle he0 = mesh->halfedge_handle(chart_bound[0], 0);
	//int start = mesh->from_vertex_of_halfedge(he0).idx();
	int start = mesh->edge(chart_bound[0]).from_vertex().idx();
	//int end = mesh->to_vertex_of_halfedge(he0).idx();
	int end = mesh->edge(chart_bound[0]).to_vertex().idx();

	assert(start != end);
	//count_ is for test
	int count_ = 1;
	//cout << start << " " <<end<< endl;
	while (chart_bound.size() < c_ && count_ < 100000)
	{
		count_++;
		for (size_t i = 0; i < chart_bound_noorder.size(); ++i)
		{
			//HalfEdgeHandle he0 = mesh->halfedge_handle(chart_bound_noorder[i], 0);
			OpenVolumeMeshEdge edge_ = mesh->edge(chart_bound_noorder[i]);
			if (edge_.from_vertex().idx() == end)
			{
				assert(start != end);
				chart_bound.push_back(chart_bound_noorder[i]);
				end = edge_.to_vertex().idx();
				vector<EdgeHandle>::iterator it = chart_bound_noorder.begin() + i;
				chart_bound_noorder.erase(it);
				break;
			}
			else if (edge_.to_vertex().idx() == end)
			{
				assert(start != end);
				chart_bound.push_back(chart_bound_noorder[i]);
				end =edge_.from_vertex().idx();
				vector<EdgeHandle>::iterator it = chart_bound_noorder.begin() + i;
				chart_bound_noorder.erase(it);
				break;
			}
		}

		if (start == end && chart_bound.size()<c_)
			//one bound of a chart is finded, then find other bounds of the chart
		{
			EdgeHandle e_start = chart_bound_noorder[0];
			OpenVolumeMeshEdge e_start_ = mesh->edge(e_start);
			//HalfEdgeHandle he_start = mesh->halfedge_handle(e_start, 0);

			charts_boundary_start_edges[c_i].push_back(chart_bound.size());
			chart_bound.push_back(e_start);

			vector<EdgeHandle>::iterator it = chart_bound_noorder.begin();
			chart_bound_noorder.erase(it);
			start = e_start_.from_vertex();
			end =e_start_.to_vertex();
			assert(start != end);
		}
	}

	//cout <<c_i<<" "<< count_ << endl;
}
//Bound , Edge structures
bool Labeling::is_visit_vchart(int c_i, vector<int>&v_chart)
{
	for (size_t i = 0; i < v_chart.size(); ++i)
	{
		if (c_i == v_chart[i])
			return true;
	}
	return false;
}
void Labeling::GetVChartid()
{
	v_chartsid.clear();
	v_chartsid.resize(n_boundary_vertices);
	for (size_t i = 0; i < charts.size(); ++i)
	{
		for (size_t j = 0; j < charts[i].size(); ++j)
		{
			FaceHandle f = charts[i][j];
			int s_faceid = s_face_id(f);
			assert(b_face2vertices[s_faceid].size() == 3);
			for (size_t k = 0; k < b_face2vertices[s_faceid].size(); ++k)
			{
				VertexHandle v = b_face2vertices[s_faceid][k];
				int s_vid = s_vertex_id(v);
				if (!is_visit_vchart(i, v_chartsid[s_vid]))
				{
					v_chartsid[s_vid].push_back(i);
				}
			}
		}
	}
	for (size_t i = 0; i < v_chartsid.size(); ++i)
	{
		if (v_chartsid[i].size() >3)
		{
			printf("corner %d valence >3 \n", s2v_vertex[i]);
			//VolumeMeshViewer->special_points.push_back(s2v_vertex[i]);
		}
	}
}

bool Labeling::extract_vbound()
{
	charts_vbound.clear();
	bound2chartid.clear();
	vbound_startid.resize(charts_boundary.size());
	charts_vbounds.resize(charts_boundary.size());
	chart2boundid.resize(charts_boundary.size());
	for (size_t i = 0; i < charts_boundary.size(); ++i)
	{
		vbound_startid[i].clear();
		charts_vbounds[i].clear();
		chart2boundid[i].clear();
	}

	for (size_t k = 0; k < charts_boundary.size(); ++k)
	{
		int c_i = k;
		vector<int>&chart_bound_startedgeid = charts_boundary_start_edges[c_i];
		vector<EdgeHandle>&chart_ebound = charts_boundary[c_i];
		//cout << chart_bound_startedgeid.size()<< "startsize chart_boundsize chart_size c_i"<<chart_bound.size()<<" "<<charts[c_i].size()<<" "<<c_i;
		for (size_t i = 0; i < chart_bound_startedgeid.size(); ++i)
		{
			int start_id = chart_bound_startedgeid[i];
			int end_id;
			if (i == chart_bound_startedgeid.size() - 1)
			{
				end_id = chart_ebound.size();
			}
			else
			{
				end_id = chart_bound_startedgeid[i + 1];
			}
			vbound_startid[c_i].push_back(start_id);
			if (end_id - start_id <= 2)
			{
				cout << "extract_chart_bound fun error,too few edges!" << " " << c_i << endl;
				return false;
			}
			vector<VertexHandle>bound;
			OpenVolumeMeshEdge e = mesh->edge(chart_ebound[start_id]);
			OpenVolumeMeshEdge nexte = mesh->edge(chart_ebound[start_id + 1]);
			VertexHandle v1 = e.from_vertex();
			VertexHandle v2 = e.to_vertex();
			VertexHandle nv1 = nexte.from_vertex();
			VertexHandle nv2 = nexte.to_vertex();

			//find the frist 3 vertices
			if (v1 == nv1)
			{
				charts_vbounds[c_i].push_back(v2);
				charts_vbounds[c_i].push_back(v1);
				charts_vbounds[c_i].push_back(nv2);
				bound.push_back(v2);
				bound.push_back(v1);
				bound.push_back(nv2);
			}
			else if (v1 == nv2)
			{
				charts_vbounds[c_i].push_back(v2);
				charts_vbounds[c_i].push_back(v1);
				charts_vbounds[c_i].push_back(nv1);
				bound.push_back(v2);
				bound.push_back(v1);
				bound.push_back(nv1);
			}
			else if (v2 == nv1)
			{
				charts_vbounds[c_i].push_back(v1);
				charts_vbounds[c_i].push_back(v2);
				charts_vbounds[c_i].push_back(nv2);
				bound.push_back(v1);
				bound.push_back(v2);
				bound.push_back(nv2);
			}
			else
			{
				charts_vbounds[c_i].push_back(v1);
				charts_vbounds[c_i].push_back(v2);
				charts_vbounds[c_i].push_back(nv1);
				bound.push_back(v1);
				bound.push_back(v2);
				bound.push_back(nv1);
			}
			for (int j = start_id + 2; j < end_id - 1; ++j)
			{
				OpenVolumeMeshEdge e = mesh->edge(chart_ebound[j]);
				vector<VertexHandle>::iterator vit = charts_vbounds[c_i].end() - 1;
				VertexHandle v1 = e.from_vertex();
				VertexHandle v2 = e.to_vertex();
				if (*vit == e.from_vertex())
				{
					charts_vbounds[c_i].push_back(v2);
					bound.push_back(v2);
				}
				else
				{
					assert(*vit == e.to_vertex());
					charts_vbounds[c_i].push_back(v1);
					bound.push_back(v1);
				}
			}
			charts_vbound.push_back(bound);
			chart2boundid[c_i].push_back(charts_vbound.size() - 1);
			bound2chartid.push_back(c_i);
		}
	}
	return true;
}
void Labeling::ChangeCEfaceLabel()
{
	printf("****************check corner_extrema faces' label*************\n");
	bool HasChangeLabel = false;
	int iter_change_count = 0;
	do
	{
		++iter_change_count;
		printf("itertaion %d \n", iter_change_count);
		HasChangeLabel = false;
		//turning point
		extract_vbound();
		TurningPoints(0);
		GetVChartid();
		changeface.clear();
		newchartid.clear();
		//GetEdgeStructure();
		//vector<bool>visit_corner; visit_corner.resize(n_boundary_vertices, false);
		for (size_t i = 0; i < extremas.size(); ++i)
		for (size_t j = 0; j < extremas[i].size(); ++j)
		{
			corner_n_extrema(i, j);
		}

		if (changeface.size() != 0)HasChangeLabel = true;
		//delete changed faces
		for (size_t i = 0; i < changeface.size(); ++i)
		{
			FaceHandle f = changeface[i];
			int s_fid = s_face_id(f);

			int oc_i = face_chart_id[s_fid];
			int nc_i = newchartid[i];
			s_face_axis[s_fid] = chart_axis(nc_i);
			face_chart_id[s_fid] = nc_i;
			charts[nc_i].push_back(f);
			//delete face
			vector<FaceHandle>::iterator it = charts[oc_i].begin();
			while (*it != f && it != charts[oc_i].end())++it;
			if (*it == f)
				charts[oc_i].erase(it);
			else
				cout << "delete face from chart error splitchart class!" << endl;
		}
		//chart_bound_faceid.clear();
		//chart_neighbor_.clear();
		//get_charts();
		get_charts_bound();
		straighten_bound(true);
	} while (HasChangeLabel && iter_change_count<10);

	if (iter_change_count == 10)printf("change label iteration error!\n");
}
void Labeling::corner_n_extrema(int c_i, int jj)
{
	int bounds_vid = extremas[c_i][jj];
	//vector<int>&bounds_startid = charts_boundary_start_edges[ii];
	VertexHandle extrema, corner, ano_v, ano_corner;
	extrema = charts_vbounds[c_i][bounds_vid];
	//test
	//VolumeMeshViewer->special_points.push_back(extrema);
	EdgeHandle extrema_corner_edge;
	int bid = -1; int eid_onb = -1; int cid_onb = -1;
	for (size_t i = 0; i < chart2boundid[c_i].size(); ++i)
	{
		int boundid = chart2boundid[c_i][i];
		for (size_t j = 0; j < charts_vbound[boundid].size(); ++j)
		{
			if (extrema == charts_vbound[boundid][j])
			{
				bid = boundid;
				eid_onb = j;
				break;
			}
		}
	}
	if (bid == -1)printf("cannot found extrema in bound\n");
	//find neighbor corner and c_e edge
	vector<VertexHandle>&bound = charts_vbound[bid];
	vector<int>e_n_cornersid_onb;
	int tmpid = (eid_onb + 1) % bound.size();
	if (is_corner(bound[tmpid]))e_n_cornersid_onb.push_back(tmpid);
	tmpid = (eid_onb - 1 + bound.size()) % bound.size();
	if (is_corner(bound[tmpid]))e_n_cornersid_onb.push_back(tmpid);
	if (e_n_cornersid_onb.empty())return;

	//printf("find corners %d neighboring extrema\n", bound[e_n_cornersid_onb[0]]);
	//test
	//VolumeMeshViewer->special_points.push_back(bound[e_n_cornersid_onb[0]]);
	//VolumeMeshViewer->special_points.push_back(extrema);
	//end test

	//find e_c_edge
	for (size_t ii = 0; ii < e_n_cornersid_onb.size(); ++ii)
	{
		cid_onb = e_n_cornersid_onb[ii];
		corner = bound[cid_onb];
		int s_corner_id = s_vertex_id(corner);
		bool found_edge = false;
		for (size_t j = 0; j < b_vertex2facesid[s_corner_id].size(); ++j)
		{
			int c_faceid = b_vertex2facesid[s_corner_id][j];
			extrema_corner_edge = -1;
			for (size_t k = 0; k < b_face2edges[c_faceid].size(); ++k)
			{
				OpenVolumeMeshEdge e = mesh->edge(b_face2edges[c_faceid][k]);
				if (e.from_vertex() == extrema && e.to_vertex() == corner)
				{
					extrema_corner_edge = b_face2edges[c_faceid][k];
					found_edge = true;
					break;
				}
				else if (e.to_vertex() == extrema &&e.from_vertex() == corner)
				{
					extrema_corner_edge = b_face2edges[c_faceid][k];
					found_edge = true;
					break;
				}
			}
			if (found_edge)break;
		}
		if (!found_edge)printf("error! cannot find edge\n");
		int s_edgeid = s_edge_id(extrema_corner_edge);
		FaceHandle f;
		FaceHandle anof;
		if (c_i == faceChartId(b_edge2faces_[s_edgeid][0]))
		{
			f = b_edge2faces_[s_edgeid][0]; anof = b_edge2faces_[s_edgeid][1];
		}
		else if (c_i == faceChartId(b_edge2faces_[s_edgeid][1]))
		{
			f = b_edge2faces_[s_edgeid][1]; anof = b_edge2faces_[s_edgeid][0];
		}
		//find goal_axe
		axis this_axe, ano_axe;
		int goal_axe = -1;
		this_axe = chart_axis(c_i);
		ano_axe = s_face_axis[s_face_id(anof)];
		vector<FaceHandle>&faces = f_3faces(f);

		for (size_t i = 0; i < faces.size(); ++i)
		{
			axis axe = s_face_axis[s_face_id(faces[i])];
			if ((axe != this_axe) && (axe != ano_axe))
			{
				goal_axe = axe; break;
			}
		}
		if (goal_axe == -1)return;// the face isn't on 2 bounds
		printf("extrema %d corner %d\n", extrema.idx(), corner.idx());
		//VolumeMeshViewer->specialFaces.push_back(f);
		//find another vertex of the triangle
		for (size_t k = 0; k < b_face2vertices[s_face_id(f)].size(); ++k)
		{
			VertexHandle tmp = b_face2vertices[s_face_id(f)][k];
			if ((tmp != extrema) && (tmp != corner))
			{
				ano_v = tmp; break;
			}
		}
		Geometry::Vec3d FutureNormal, CurrentNormal;
		//estimate future normal
		int c1id_onb = -1;
		int c2id_onb = -1;
		find_2neighbor_corners(bid, cid_onb, c1id_onb, c2id_onb);
		if (c1id_onb == -1 || c2id_onb == -1)
		{
			printf("find 2 neighbor corners error!\n");
		}
		int e1dir = 3 - (this_axe % 3) - (goal_axe % 3);//c->ano_v
		int e2dir = 3 - (this_axe % 3) - (ano_axe % 3);//c->extrema
		Geometry::Vec3d e1(0.0, 0.0, 0.0);
		Geometry::Vec3d e2(0.0, 0.0, 0.0);

		double tmp1 = (mesh->vertex(bound[c1id_onb]) - mesh->vertex(corner))[e1dir];//cid+
		double tmp2 = (mesh->vertex(bound[c2id_onb]) - mesh->vertex(corner))[e2dir];//cid-
		if (eid_onb == (cid_onb + 1) % bound.size())
		{
			//cid+ ; c->extrema
			e1[e1dir] = tmp2; e2[e2dir] = tmp1;
		}
		else if (eid_onb == (cid_onb - 1 + bound.size()) % bound.size())
		{
			//cid- ; c->extrema
			e1[e1dir] = tmp1; e2[e2dir] = tmp2;
		}
		FutureNormal = e2%e1; //c->extrema % c->ano_v
		if (e1[0] == 0 && e1[1] == 0 && e1[2] == 0)printf("error in find future normal\n");
		//find current normal
		Geometry::Vec3d ce1, ce2;
		ce1 = mesh->vertex(extrema) - mesh->vertex(corner);
		ce2 = mesh->vertex(ano_v) - mesh->vertex(corner);

		/*if (eid_onb == (cid_onb + 1) % bound.size()){ CurrentNormal = ce1 % ce2; }
		else if (eid_onb == (cid_onb - 1 + bound.size()) % bound.size()){ CurrentNormal= ce2 % ce1; }*/
		CurrentNormal = ce1 % ce2;
		if ((CurrentNormal | FutureNormal) < 0)
		{
			printf("change face  %d label! %d -> %d\n", f.idx(), this_axe, goal_axe);
			int new_ci = -1;
			int s_cid = s_vertex_id(corner);
			for (size_t i = 0; i < v_chartsid[s_cid].size(); ++i)
			{
				if ((v_chartsid[s_cid][i] != faceChartId(f)) && (v_chartsid[s_cid][i] != faceChartId(anof)))
				{
					new_ci = v_chartsid[s_cid][i];
				}
			}
			changeface.push_back(f);
			newchartid.push_back(new_ci);
		}
	}


}
void Labeling::find_2neighbor_corners(int bid, int cid_onb, int &c1id_onb, int &c2id_onb)
{
	int whilecount = 10000;
	int tmp = cid_onb;
	int iter = 0;
	do{
		iter++;
		tmp = (tmp + 1) % charts_vbound[bid].size();
		VertexHandle v = charts_vbound[bid][tmp];
		if (is_corner(v))
		{
			c1id_onb = tmp; break;
		}
	} while (iter < whilecount);
	tmp = cid_onb; iter = 0;
	do{
		tmp = (tmp - 1 + charts_vbound[bid].size()) % charts_vbound[bid].size();
		VertexHandle v = charts_vbound[bid][tmp];
		if (is_corner(v))
		{
			c2id_onb = tmp; break;
		}
	} while (iter < whilecount);
}
void Labeling::corner_tri(VertexHandle corner)
{
	
}
void Labeling::ChangeCornerFaceLabel()
{
	printf("****************check corner faces' label*************\n");
	
	bool HasChangeLabel = false;
	int iter_change_count = 0;
	do
	{
		++iter_change_count;
		printf("itertaion %d \n", iter_change_count);
		HasChangeLabel = false;
		//turning point
		//TurningPoints();
		GetVChartid();
		if (!extract_vbound())return;
		GetEdgeStructure();
		vector<bool>visit_corner; visit_corner.resize(n_boundary_vertices, false);
		for (size_t i = 0; i < NorepeatEdges.size(); ++i)
		{
			vector<VertexHandle>&Edge = ChartsEdges[NorepeatEdges[i]];
			VertexHandle c1, c2;
			c1 = Edge[0]; c2 = Edge[Edge.size() - 1];
			if (!visit_corner[s_vertex_id(c1)])
			{
				visit_corner[s_vertex_id(c1)] = true;
				corner_tri(c1);
			}
			if (!visit_corner[s_vertex_id(c2)])
			{
				visit_corner[s_vertex_id(c2)] = true;
				corner_tri(c2);
			}
			
		}
		if (changeface.size() == 0)HasChangeLabel = false;
		else
			HasChangeLabel = true;
		//delete changed faces
		for (size_t i = 0; i < changeface.size(); ++i)
		{
			FaceHandle f = changeface[i];
			int s_fid = s_face_id(f);
			int oc_i = face_chart_id[s_fid];
			int nc_i = newchartid[i];
			face_chart_id[s_fid] = nc_i;
			charts[nc_i].push_back(f);
			//delete face
			vector<FaceHandle>::iterator it = charts[oc_i].begin();
			while (*it != f && it != charts[oc_i].end())++it;
			if (*it == f)
				charts[oc_i].erase(it);
			else
				cout << "delete face from chart error splitchart class!" << endl;
		}
		chart_bound_faceid.clear();
		chart_neighbor_.clear();
		//get_charts();
		get_charts_bound();
		straighten_bound(true);
	} while (HasChangeLabel && iter_change_count<10);
	
	if (iter_change_count == 10)printf("change label iteration error!\n");
}
int Labeling::EstimateCoorDir(bool increase, int corner_bound_vid, int start_id, int end_id, int c_id, int dir)
{
	vector<VertexHandle>&bound = charts_vbound[c_id];
	VertexHandle v, corner;
	int size_ = end_id - start_id;
	corner = bound[corner_bound_vid];
	int id = corner_bound_vid;
	if (increase)
	{
		do
		{
			id = (id - start_id + 1) % size_ + start_id;
			v = bound[id];
		} while (v_chartsid[s_vertex_id(v)].size() != 3);
	}
	else
	{
		do
		{
			id = (id - start_id - 1 + size_) % size_ + start_id;
			v = bound[id];
		} while (v_chartsid[s_vertex_id(v)].size() != 3);
	}
	printf("nextcorner %d,corner %d\n", v, corner);
	if (mesh->vertex(v)[dir] > mesh->vertex(corner)[dir])return 1;
	else return -1;
}
void Labeling::EstimateEdgeDir(VertexHandle v1, VertexHandle v2, Geometry::Vec3d &e, int &c1id, int &c2id)
//get the edge direction unit vector e and 2 ajacent charts of the edge
{
	//get the Edge from v1 to v2
	vector<VertexHandle>*Edge;
	int signal=0;
	int length;
	bool test;
	for (size_t i = 0; i < ChartsEdges.size(); ++i)
	{
		Edge = &ChartsEdges[i];
		length = Edge->size();
		if (Edge->at(0) == v1 || Edge->at(length - 1) == v1)test = true;
		//Edge = &ChartsEdges[NorepeatEdges[i]];
		if ((Edge->at(0) == v1 )&& (Edge->at(1) == v2))
		{
			signal = 1; break;
		}
		
		if ((Edge->at(length - 1) == v1) && (Edge->at(length - 2) == v2))
		{
			signal = -1; break;
		}
	}
	if (!test)printf("ERROR_TEST\n");
	if (signal == 0)
	{
		//VolumeMeshViewer->special_points.push_back(v1);
		//VolumeMeshViewer->special_points.push_back(v2);
		printf("v1 or v2 isn't at bound\n"); return; 
	}
	VertexHandle start = Edge->at(0);
	VertexHandle end = Edge->at(length- 1);
	vector<int>&s_charts = v_chartsid[s_vertex_id(start)];
	vector<int>&e_charts = v_chartsid[s_vertex_id(end)];
	vector<int>sc;
	assert(s_charts.size() == 3 && e_charts.size() == 3);
	c1id = -1; c2id = -1;
	for (size_t i = 0; i<s_charts.size(); ++i)
	{
		for (size_t j = 0; j<e_charts.size(); ++j)
		{
			if (s_charts[i] == e_charts[j])
			{
				sc.push_back(s_charts[i]);
				break;
			}
		}
	}
	if (sc.size() != 2)printf("same chart find error %d in EstimateEdgeDir fun !\n", sc.size());
	c1id = sc[0]; c2id = sc[1];
	e = mesh->vertex(end) - mesh->vertex(start);
	int dir1, dir2, dir;
	dir1 = chart_axis(c1id) % 3; dir2 = chart_axis(c2id) % 3; dir = 3 - dir1 - dir2;
	e[dir1] = 0; e[dir2] = 0;
	if (signal == -1)e[dir] = 0 - e[dir];
}
void Labeling::TurningPoints(int smooth_iteration)
{
	//cout << "***********************is Editing split_nonplanar_chart!*********************" << endl;
	//cout << "chartssize" << pos_->charts.size() << endl;
	extremas.resize(charts.size());
	is_extrema.resize(charts.size());
	for (size_t i = 0; i < is_extrema.size(); ++i)is_extrema[i].resize(charts_boundary[i].size(), false);
	//isvisit_v.clear(); isvisit_v.assign(v_chartsid.size(), false);
	for (size_t i = 0; i < extremas.size(); ++i)
	{
		extremas[i].clear();
		for (size_t j = 0; j < vbound_startid[i].size(); ++j)
		{
			int start_id = vbound_startid[i][j];
			int end_id;
			if (j == vbound_startid[i].size() - 1)
			{
				end_id = charts_vbounds[i].size();
			}
			else
			{
				end_id = vbound_startid[i][j + 1];
			}
			get_extrema(start_id, end_id, i, smooth_iteration);
		}
	}

	
}
bool Labeling::filt_weak_extrema(double area_threshold)
{
	filtered_extrema = true;
	double  max = -1;
	VertexHandle max_area_extrema, va, vb;
	vector<bool>isvisit_extrema;
	isvisit_extrema.resize(n_boundary_vertices, false);
	printf("area percentage threshold:%f  \nextremas:\n", area_threshold);
	for (size_t i = 0; i < extremas.size(); ++i)
	{
		for (size_t j = 0; j < extremas[i].size(); ++j)
		{
			int bound_vid = extremas[i][j];
			int start_id = -1;
			int end_id = -1;
			find_se_id(bound_vid, i, start_id, end_id);
			//printf("s e %d %d\n", start_id, end_id);
			if (start_id == -1 || end_id == -1)cout << "error find se id --filt_weak_extrema fun" << endl;
			VertexHandle extrema = charts_vbounds[i][bound_vid];
			if (!isvisit_extrema[s_vertex_id(extrema)])
			{
				isvisit_extrema[s_vertex_id(extrema)] = true;
			}
			//find another two points
			VertexHandle v1 = -1; VertexHandle v2 = -1;
			int id = bound_vid;
			int size_ = end_id - start_id;
			do
			{
				id = (id - start_id + 1) % size_ + start_id;
				v1 = charts_vbounds[i][id];
			} while (v_chartsid[s_vertex_id(v1)].size() != 3 && !is_extrema[i][id]);
			id = bound_vid;
			do
			{
				//		printf("id %d\n", id);
				id = (id - start_id - 1 + size_) % size_ + start_id;
				v2 = charts_vbounds[i][id];
			} while (v_chartsid[s_vertex_id(v2)].size() != 3 && !is_extrema[i][id]);
			if (v1 == -1 || v2 == -1)cout << "v1,v2 find error --filt_weak_extrema" << endl;
			if (v1 == v2)cout << "v1==v2 --fild_weak_extrema" << endl;
			Geometry::Vec3d _e, _v1, _v2;
			_e = mesh->vertex(extrema);
			_v1 = mesh->vertex(v1);
			_v2 = mesh->vertex(v2);
			double area = 0.5*((_v1 - _e) % (_v2 - _e)).norm();
			if (area /avg_bface_area > area_threshold)
			{
				printf("extrema %d isn't filtered\n", charts_vbounds[i][bound_vid].idx());
				filtered_extrema = false;
				//VolumeMeshViewer->special_points.push_back(charts_vbounds[i][bound_vid].idx());
			}
			if (area > max)
			{
				va = v1;
				vb = v2;
				max = area;
				max_area_extrema = extrema;
			}
		}
	}
	if (max > 0)
	{
		printf("max extrema area %f\n", max);
		printf("max percentage %f\n", max / avg_bface_area);
	}
	if (!filtered_extrema)return false;
	else
		return true;
}
void Labeling::compute_average_bface_area()
{
	Geometry::Vec3d p_0, p_1, p_2;
	double max = -1;
	double sum = 0;
	for (size_t i = 0; i < b_face2vertices.size(); ++i)
	{
		p_0 = mesh->vertex(b_face2vertices[i][0]);
		p_1 = mesh->vertex(b_face2vertices[i][1]);
		p_2 = mesh->vertex(b_face2vertices[i][2]);
		double area = (OpenVolumeMesh::Geometry::cross(p_0 - p_1, p_2 - p_1)).norm();
		sum += area;
		if (max <area)max = area;
	}
	avg_bface_area = 0.5*sum / (double)b_face2vertices.size();
	printf("average boundary triange area %f\n", avg_bface_area);
}
void Labeling::find_se_id(int bound_vid, int c_i, int &start_id, int &end_id)
{
	vector<int>&_bound_startid = vbound_startid[c_i];
	bool found = false;
	for (size_t i = 0; i < _bound_startid.size(); ++i)
	{
		if (_bound_startid[i] >bound_vid)
		{
			start_id = _bound_startid[i - 1];
			end_id = _bound_startid[i];
			found = true;
			break;
		}
	}
	if (!found)
	{
		start_id = _bound_startid[_bound_startid.size() - 1];
		end_id = charts_boundary[c_i].size();
	}
}
void Labeling::get_extrema(int start_id, int end_id, int c_i, int smooth_iteration)
{
	MatrixXd smooth_vertices;
	//cout << "boundsize" << smoothed_bound.size() << endl;
	smooth_chart_loopbound(c_i, start_id, end_id, smooth_vertices, smooth_iteration);
	int row = smooth_vertices.rows();
	//first vertex
	int vid = s_vertex_id(charts_vbounds[c_i][start_id]);
	int previd = s_vertex_id(charts_vbounds[c_i][end_id - 1]);
	int nextvid = s_vertex_id(charts_vbounds[c_i][start_id + 1]);
	if (is_vertex_extrema(previd, vid, nextvid, 0, smooth_vertices))
	{
		//int b_vid = charts_vbounds[c_i][start_id];
		//if (!isvisit_v[b_vid])
		//{
		//isvisit_v[b_vid] = true;
		extremas[c_i].push_back(start_id);
		is_extrema[c_i][start_id] = true;
#ifdef TEST_ALLEXTREMAS
		VolumeMeshViewer->special_points.push_back(charts_vbounds[c_i][start_id]);
#endif
		//}
	}

	for (int i = 1; i < row - 1; ++i)
	{
		vid = s_vertex_id(charts_vbounds[c_i][start_id + i]);
		previd = s_vertex_id(charts_vbounds[c_i][start_id + i - 1]);
		nextvid = s_vertex_id(charts_vbounds[c_i][start_id + i + 1]);
		if (is_vertex_extrema(previd, vid, nextvid, i, smooth_vertices))
		{
			//int b_vid = charts_vbounds[c_i][start_id + i];
			//	if (!isvisit_v[b_vid])
			//{
			//isvisit_v[b_vid] = true;
			extremas[c_i].push_back(start_id + i);
			is_extrema[c_i][start_id + i] = true;
#ifdef TEST_ALLEXTREMAS
			VolumeMeshViewer->special_points.push_back(charts_vbounds[c_i][start_id + i]);
#endif
			//}
		}
	}

	//last vertex
	vid = s_vertex_id(charts_vbounds[c_i][end_id - 1]);
	previd = s_vertex_id(charts_vbounds[c_i][end_id - 2]);
	nextvid = s_vertex_id(charts_vbounds[c_i][start_id]);
	if (is_vertex_extrema(previd, vid, nextvid, row - 1, smooth_vertices))
	{
		//int b_vid = charts_vbounds[c_i][end_id - 1];
		//if (!isvisit_v[b_vid])
		//{
		//isvisit_v[b_vid] = true;
		extremas[c_i].push_back(end_id - 1);
		is_extrema[c_i][end_id - 1] = true;
		//pos_->specialPoints->push_back(bound_vertex[c_i][end_id-1]);
#ifdef TEST_ALLEXTREMAS
		VolumeMeshViewer->special_points.push_back(charts_vbounds[c_i][end_id - 1]);
#endif
		//}
	}
}
bool Labeling::is_vertex_extrema(int previd, int vid, int nextvid, int v_matrixid, MatrixXd &SmoothVertices)
{
	if (v_chartsid[vid].size() != 2)return false;
	int v_c1id = v_chartsid[vid][0];
	int v_c2id = v_chartsid[vid][1];
	bool has_same_chart1 = false;
	bool has_same_chart2 = false;
	for (int i = 0; i < v_chartsid[previd].size(); ++i)
	{
		if (v_chartsid[previd][i] == v_c1id)
		{
			has_same_chart1 = true;
		}
		else if (v_chartsid[previd][i] == v_c2id)
		{
			has_same_chart2 = true;
		}
	}
	if (!(has_same_chart1&&has_same_chart2))return false;
	has_same_chart1 = false;
	has_same_chart2 = false;
	for (int i = 0; i < v_chartsid[nextvid].size(); ++i)
	{
		if (v_chartsid[nextvid][i] == v_c1id)
		{
			has_same_chart1 = true;
		}
		else if (v_chartsid[nextvid][i] == v_c2id)
		{
			has_same_chart2 = true;
		}
	}
	if (!(has_same_chart1&&has_same_chart2))return false;

	int direction1 = (chart_axis(v_c1id)) % 3;
	int direction2 = (chart_axis(v_c2id)) % 3;
	int edge_direction = 3 - direction1 - direction2;
	int row = SmoothVertices.rows();
	double  prev_coor = SmoothVertices((v_matrixid + row - 1) % row, edge_direction);
	double  v_coor = SmoothVertices(v_matrixid, edge_direction);
	double nextv_coor = SmoothVertices((v_matrixid + 1) % row, edge_direction);
	if ((v_coor - prev_coor)*(nextv_coor - v_coor) < 0)return true;
	else return false;
}
void Labeling::smooth_chart_loopbound(int c_i, int start_id, int end_id, MatrixXd &smooth_vertices, int smooth_iteration)
{
	int n_vertices = end_id - start_id;
	smooth_vertices = MatrixXd::Zero(n_vertices, 3);
	for (int j = 0; j < n_vertices; ++j)
	{
		Geometry::Vec3d v = mesh->vertex(charts_vbounds[c_i][start_id + j]);
		smooth_vertices(j, 0) = v[0]; smooth_vertices(j, 1) = v[1]; smooth_vertices(j, 2) = v[2];
	}
	MatrixXd A_smooth = MatrixXd::Zero(n_vertices, n_vertices);
	for (int j = 0; j < n_vertices; ++j)
	{
		A_smooth(j, (j - 1 + n_vertices) % n_vertices) = 0.25;
		A_smooth(j, (j + 1) % n_vertices) = 0.25;
		A_smooth(j, j) = 0.5;
		vector<VertexHandle>::iterator v_it = charts_vbounds[c_i].begin() + j;
		Geometry::Vec3d v = mesh->vertex(*v_it);
	}
	for (int j = 0; j < smooth_iteration; ++j)
	{
		smooth_vertices = A_smooth*smooth_vertices;
	}
}
bool Labeling::GetEdgeStructure()
{
	//bound<->edge
	vector<VertexHandle>Edge;
	for (size_t i = 0; i < ChartsEdges.size(); ++i)ChartsEdges[i].clear();
	ChartsEdges.clear();
	NorepeatEdges.clear();
	bound2edgeid.clear(); bound2edgeid.resize(charts_vbound.size());
	edge2boundid.clear();
	Edges2Norepeatid.clear();
	for (size_t i = 0; i < charts_vbound.size(); ++i)
	{
		vector<VertexHandle>&bound = charts_vbound[i];
		vector<int>corners;
		//get all corners of the bound
		for (size_t j = 0; j < bound.size(); ++j)
		{
			if (is_corner(bound[j]))corners.push_back(j);
		}
		if (corners.size()<4)
		{
			printf("error! bound corners number <4!\n");
			return false;
			//VolumeMeshViewer->special_points.push_back(bound[corners[0]]);
		}
		for (size_t j = 0; j < corners.size() - 1; ++j)
		{
			int pid = corners[j];
			int nid = corners[j + 1];
			Edge.clear();
			for (int k = pid; k <= nid; ++k)
			{
				Edge.push_back(bound[k]);
			}
			ChartsEdges.push_back(Edge);
			edge2boundid.push_back(i);
			bound2edgeid[i].push_back(ChartsEdges.size() - 1);
			//check repeatation
			int status = has_thisEdge(Edge);
			if (status == -1)//no repeat
			{
				NorepeatEdges.push_back(ChartsEdges.size() - 1);
				Edges2Norepeatid.push_back(NorepeatEdges.size() - 1);
			}
			else
			{
				Edges2Norepeatid.push_back(status);
			}
		}
		Edge.clear();
		for (int k = corners[corners.size() - 1]; k < bound.size(); ++k)Edge.push_back(bound[k]);
		for (int k = 0; k <= corners[0]; ++k)Edge.push_back(bound[k]);
		ChartsEdges.push_back(Edge);
		edge2boundid.push_back(i);
		bound2edgeid[i].push_back(ChartsEdges.size() - 1);
		//check repeatation
		int status = has_thisEdge(Edge);
		if (status == -1)//no repeat
		{
			NorepeatEdges.push_back(ChartsEdges.size() - 1);
			Edges2Norepeatid.push_back(NorepeatEdges.size() - 1);
		}
		else
		{
			Edges2Norepeatid.push_back(status);
		}
	}
	return true;
}
bool Labeling::is_corner(VertexHandle v)
{
	int s_vid = s_vertex_id(v);
	if (v_chartsid[s_vid].size() == 3)return true;
	else return false;
}
int Labeling::has_thisEdge(vector<VertexHandle>&Edge)
{
	for (size_t i = 0; i < NorepeatEdges.size(); ++i)
	{
		int eid = NorepeatEdges[i];
		vector<VertexHandle>::iterator it1, it2, it3, it4;
		it1 = Edge.begin(); it2 = Edge.end() - 1;
		it3 = ChartsEdges[eid].begin(); it4 = ChartsEdges[eid].end() - 1;
		if ((*it1 == *it3 && *it2 == *it4) || (*it1 == *it4 && *it2 == *it3))
		{
			return i;
		}
	}
	return -1;
}
void Labeling::init_chartspair_dist()
{
	chartspair_dist.resize(charts.size());
	chartpair_edges.resize(charts.size());
	for (size_t i = 0; i < chartspair_dist.size(); ++i)
	{
		chartspair_dist[i].assign(charts.size(), 10000);
		chartpair_edges[i].assign(charts.size(), -1);
	}
}
void Labeling::GetTetLabel()
{
	tetLabel.resize(mesh->n_cells(), -1);
	tet_face_chart_id.resize(mesh->n_cells(), -1);
	
	for (size_t i = 0; i < b_face2cell.size(); ++i)
	{
		int c_id = b_face2cell[i].idx();

		if (c_id < 0)printf("face2cell error!\n");

		int label = s_face_axis[i];
		tetLabel[c_id] = label;
		tet_face_chart_id[c_id] = face_chart_id[i];
	}
}
void Labeling::GetAxisAlignTetLabel()
{
	axis_align_tetLabel.clear();
	axis_align_tetLabel.resize(mesh->n_cells(), -1);
	for (size_t i = 0; i < b_face2cell.size(); ++i)
	{
		int c_id = b_face2cell[i].idx();
		if (c_id < 0)printf("face2cell error!\n");
		int label = s_face_axis[i];
		axis_align_tetLabel[c_id] = label;
	}
}
void Labeling::GetInitBoundaryEdge()
{
	Init_EdgeBoudnaryTag.clear();
	Init_EdgeBoudnaryTag.resize(mesh->n_edges(), 0);
	for (size_t i = 0; i < b_edge2faces_.size(); ++i)
	{
		int first_axe = s_face_axis[ s_face_id(b_edge2faces_[i][0]) ];
		if (b_edge2faces_[i].size() == 1)
		{
			printf("edge with only 1 surface triangles!\n");
			continue;
		}
		for (size_t j = 1; j < b_edge2faces_[i].size(); ++j)
		{
			if (s_face_axis[s_face_id(b_edge2faces_[i][j])] != first_axe)
			{
				std::map<int, EdgeHandle>::iterator e_it = maps2v_edge.find(i);
				if (e_it == maps2v_edge.end())printf("edge find error!\n");
				EdgeHandle e = e_it->second;
				Init_EdgeBoudnaryTag[e.idx()] = 1;
				break;
			}
		}
	}
}
void Labeling::EdgeChartsPair()
{
	init_chartspair_dist();
	for (size_t i = 0; i <NorepeatEdges.size(); ++i)
	{
		vector<VertexHandle>Edge = ChartsEdges[NorepeatEdges[i]];
		VertexHandle s = Edge[0];
		VertexHandle e = Edge[Edge.size() - 1];
		vector<int>&s_charts = v_chartsid[s_vertex_id(s)];
		vector<int>&e_charts = v_chartsid[s_vertex_id(e)];
		vector<int>sc; sc.clear();
		int c1 = -1; int c2 = -1;
		for (size_t ii = 0; ii<3; ++ii)
		{
			for (size_t j = 0; j<3; ++j)
			{
				if (s_charts[ii] == e_charts[j])
				{
					sc.push_back(s_charts[ii]);
					break;
				}
			}
		}
		if (sc.size() != 2)printf("%d same chart of c-c Edge!\n", sc.size());
		for (size_t ii = 0; ii < 3; ++ii)
		{
			if (s_charts[ii] != sc[0] && s_charts[ii] != sc[1])c1 = s_charts[ii];
			if (e_charts[ii] != sc[0] && e_charts[ii] != sc[1])c2 = e_charts[ii];
		}
		if (c1 == -1 || c2 == -1)printf("diff chart find error!\n");
		int dir1 = chart_axis(sc[0]) % 3;
		int dir2 = chart_axis(sc[1]) % 3;
		int dir = 3 - dir1 - dir2;
		double l = mesh->vertex(s)[dir] - mesh->vertex(e)[dir];
		if (l>0)
		{
			chartspair_dist[c1][c2] = l;
			chartpair_edges[c1][c2] = NorepeatEdges[i];

		}
		else
		{
			chartspair_dist[c2][c1] = 0 - l;
			chartpair_edges[c2][c1] = NorepeatEdges[i];
		}
	}
	for (size_t i = 0; i < chartspair_dist.size(); ++i)
	{
		for (size_t j = i + 1; j < chartspair_dist[i].size(); ++j)
		{
			if (chartspair_dist[i][j] != 10000 && chartspair_dist[j][i] != 10000)
			{
				printf("Conflict in C2C rule!************* %d-%d****%lf  %lf\n", i, j, chartspair_dist[i][j], chartspair_dist[j][i]);
				chartspair_dist[i][j] = 0; chartspair_dist[j][i] = 0;
				conflict_pair.push_back(pair<int, int>(i, j));
#ifdef SOLVE_CONFLICT
				solve_conflict(i, j);
				solve_conflict(j, i);
#else
				if(chartspair_dist[i][j]>chartspair_dist[j][i])
				{
					chartspair_dist[j][i] = 10000;
				}
				else
				{
					chartspair_dist[i][j] = 10000;
				}
#endif
			}
		}
	}
}
void Labeling::solve_conflict(int c_i,int c_j)
{
	//int edgeid = chartpair_edges[c_i][c_j];
	//vector<VertexHandle> &Edge1 = ChartsEdges[edgeid];
	////vector<VertexHandle> Edge2 = ChartsEdges[chartpair_edges[c_j][c_i]];
	//VertexHandle cs, ce, remove_c, corner;
	//cs = Edge1[0];
	//ce = Edge1[Edge1.size() - 1];
	//VolumeMeshViewer->special_points.push_back(cs);	VolumeMeshViewer->special_points.push_back(ce);
	//int c1, c2;
	//int changing_c, orthogonal_c, added_c, deleted_c;
	//Get2chartsOfEdge(edgeid, c1, c2);//get charts sharing the edge

	////get changing chart , remove_c ,added_c
	//vector<int>&FacesIdOfcs = b_v2fid.at(s_vertex_id(cs));
	//vector<int>&FacesIdOfce = b_v2fid.at(s_vertex_id(ce));

	//vector<int>cs_c1_fnum; vector<int>cs_c2_fnum;
	//vector<int>ce_c1_fnum; vector<int>ce_c2_fnum;
	////vector<int>&FacesIdOfce = b_v2fid->at(s_vertex_id(cs));
	//for (size_t i = 0; i < FacesIdOfcs.size(); ++i)
	//{
	//	int f_id = FacesIdOfcs[i];
	//	int c_id = face_chart_id[f_id];
	//	if (c_id == c1)cs_c1_fnum.push_back(f_id);
	//	else if (c_id == c2)cs_c2_fnum.push_back(f_id);
	//}
	//if (cs_c1_fnum.size() == 1 && cs_c2_fnum.size() == 1)printf("two few faces!try another corner!\n");
	//for (size_t i = 0; i < FacesIdOfce.size(); ++i)
	//{
	//	int f_id = FacesIdOfce[i];
	//	int c_id = face_chart_id[f_id];
	//	if (c_id == c1)ce_c1_fnum.push_back(f_id);
	//	else if (c_id == c2)ce_c2_fnum.push_back(f_id);
	//}
	//if (ce_c1_fnum.size() == 1 && ce_c2_fnum.size() == 1)printf("two few faces!try another corner!\n");
	////check whether we can change c1 faces
	//bool find_changing_c = false;
	//if (cs_c1_fnum.size() >= ce_c1_fnum.size())
	//{
	//	if (cs_c1_fnum.size() >= 3 /*&& ce_c1_fnum.size() >= 2 */ && cs_c2_fnum.size() >= 2)
	//	{
	//		changing_c = c1; orthogonal_c = c2;
	//		find_changing_c = true;
	//		remove_c = ce; corner = cs;
	//	}
	//}
	//else
	//{
	//	if (ce_c1_fnum.size() >= 3 /*&& cs_c1_fnum.size() >= 2*/ && ce_c2_fnum.size() >= 2)
	//	{
	//		changing_c = c1; orthogonal_c = c2;
	//		find_changing_c = true;
	//		remove_c = cs; corner = ce;
	//	}
	//}
	//if (!find_changing_c)
	//{
	//	if (cs_c2_fnum.size() >= ce_c2_fnum.size())
	//	{
	//		if (cs_c2_fnum.size() >= 3 /*&& ce_c2_fnum.size() >= 2*/ && cs_c1_fnum.size() >= 2)
	//		{
	//			changing_c = c2; orthogonal_c = c1;
	//			find_changing_c = true;
	//			remove_c = ce; corner = cs;
	//		}
	//	}
	//	else
	//	{
	//		if (ce_c2_fnum.size() >= 3 /*&& cs_c2_fnum.size() >= 2*/ && ce_c1_fnum.size() >= 2)
	//		{
	//			changing_c = c2; orthogonal_c = c1;
	//			find_changing_c = true;
	//			remove_c = cs; corner = ce;
	//		}
	//	}
	//}
	//if (!find_changing_c)
	//{
	//	printf("cannot find changing chart directly!\n");
	//	return;
	//}
	//printf("remove_c %d corner %d\n", remove_c, corner);
	//// get chart which will add some face
	//for (size_t i = 0; i < v_chartsid[s_vertex_id(remove_c)].size(); ++i)
	//{
	//	if (v_chartsid[s_vertex_id(remove_c)][i] != c1 && v_chartsid[s_vertex_id(remove_c)][i] != c2)
	//	{
	//		added_c = v_chartsid[s_vertex_id(remove_c)][i];
	//		break;
	//	}
	//}

	//vector<int>facesIdOfVertex;
	//for (size_t i = 0; i < Edge1.size(); ++i)
	//{

	//	if (Edge1[i] != corner)
	//	{
	//		facesIdOfVertex = b_v2fid->at(s_vertex_id(Edge1[i]));
	//		for (size_t j = 0; j < facesIdOfVertex.size(); ++j)
	//		{
	//			int s_fid = facesIdOfVertex[j];
	//			FaceHandle f = s2v_face_->at(s_fid);
	//			if (face_chart_id[s_fid] == changing_c)
	//			{
	//				//change face label
	//				/*if (Edge1[i].idx() == 10050)
	//				{
	//				printf("added_c %d ");
	//				}*/
	//				s_face_axis[s_fid] = chart_axis(added_c);
	//				charts[added_c].push_back(f);
	//				delete_face_from_chart(f, changing_c);
	//				face_chart_id[s_fid] = added_c;
	//				//SpecialFaces->push_back(f);
	//			}
	//		}
	//	}
	//}


	//printf("solve conflict!\n");
	////delete another chart
	///*for (size_t i = 0; i < charts[deleted_c].size(); ++i)
	//{
	//FaceHandle f = charts[deleted_c][i];
	//charts[added_c].push_back(f);
	//face_chart_id[s_face_id(f)] = added_c;
	//}
	//vector<vector<FaceHandle>>::iterator it = charts.begin() + deleted_c;*/
	////charts.erase(it);
	////printf()

}
void Labeling::ParallelEdgeChartsPair()
{
	vector<vector<double>>mindist_edge;
	int c11, c12, c21, c22;
	int c_1, c_2;
	//for (size_t i = 0; i < mindist_edge.size(); ++i)mindist_edge[i].assign(2,-1);
	for (size_t i = 0; i < chart2boundid.size(); ++i)
	{
		int c_dir = (chart_axis(i)) % 3;
		//in each chart
		for (size_t j1 = 0; j1 < chart2boundid[i].size(); ++j1)
		{
			//find all c-c Edges of a chart
			int boundid1 = chart2boundid[i][j1];
			for (size_t k1 = 0; k1 < bound2edgeid[boundid1].size(); ++k1)
			{
				//find an edge in the chart
				int edgeid1 = bound2edgeid[boundid1][k1];
				//	int nore_eid1;
				double min_distant = 10000;
				int mindis_edgeid_of_edge1 = -1;
				int signal = 1;
				for (size_t j2 = 0; j2 < chart2boundid[i].size(); ++j2)
				{
					int boundid2 = chart2boundid[i][j2];
					for (size_t k2 = 0; k2 < bound2edgeid[boundid2].size(); ++k2)
					{
						//find another edge in the chart
						//		int nore_eid2;
						int edgeid2 = bound2edgeid[boundid2][k2];
						//find the minimum distance edge 
						//if (!is_same_Edge(edgeid1, edgeid2))
						if (edgeid1 != edgeid2)
						{
							//nore_eid1 = Edges2Norepeatid[edgeid1];
							//nore_eid2 = Edges2Norepeatid[edgeid2];
							//int dir1 = Nore_edgeDirection[nore_eid1];
							//int dir2 = Nore_edgeDirection[nore_eid2];
							int dir1 = EdgeDirection(edgeid1);
							int dir2 = EdgeDirection(edgeid2);
							if (dir1 == dir2)
							{
								int project_dir = 3 - dir1 - c_dir;
								//double dist = Nore_edgecoor[nore_eid1][project_dir] - Nore_edgecoor[nore_eid2][project_dir];
								double dist = EstimateEdgeCoor(edgeid1, project_dir) - EstimateEdgeCoor(edgeid2, project_dir);
								if (dist < 0 && (0 - dist)<min_distant)
								{
									signal = -1;
									min_distant = 0 - dist;
									//mindis_edgeid_of_edge1 = nore_eid2;
									mindis_edgeid_of_edge1 = edgeid2;
								}
								else if (dist>0 && dist < min_distant)
								{
									signal = 1;
									min_distant = dist;
									//mindis_edgeid_of_edge1 = nore_eid2;
									mindis_edgeid_of_edge1 = edgeid2;
								}
							}
						}
					}
				}
			
				if (mindis_edgeid_of_edge1 != -1)
				{
					//c11 = Nore_2chartsofEdge[nore_eid1][0]; c12 = Nore_2chartsofEdge[nore_eid1][1];
					//c21 = Nore_2chartsofEdge[mindis_edgeid_of_edge1][0]; c22 = Nore_2chartsofEdge[mindis_edgeid_of_edge1][1];
					Get2chartsOfEdge(edgeid1, c11, c12);
					Get2chartsOfEdge(mindis_edgeid_of_edge1, c21, c22);
					bool test = false;
					if (c11 == c21)
					{
						test = true;
						c_1 = c12; c_2 = c22;
					}
					else if (c11 == c22)
					{
						test = true;
						c_1 = c12; c_2 = c21;
					}
					else if (c12 == c21)
					{
						test = true;
						c_1 = c11; c_2 = c22;
					}
					else if (c12 == c22)
					{
						test = true;
						c_1 = c11; c_2 = c21;
					}
				
					assert(test == true);
					if (c_1 != c_2 && has_common_area(edgeid1, mindis_edgeid_of_edge1))
					{
						if (min_distant < chartspair_dist[c_1][c_2])
						{
							if (signal == 1)
							{
								chartspair_dist[c_1][c_2] = min_distant;
							}
							else
							{
								assert(signal == -1);
								chartspair_dist[c_2][c_1] = min_distant;
							}
						}
					}
				}
			}
		}
	}
	
}
int Labeling::EdgeDirection(int edgeid)
{
	vector<VertexHandle>&Edge = ChartsEdges[edgeid];
	VertexHandle s = Edge[0];
	VertexHandle e = Edge[Edge.size() - 1];
	vector<int>&s_charts = v_chartsid[s_vertex_id(s)];
	vector<int>&e_charts = v_chartsid[s_vertex_id(e)];
	vector<int>sc;
	for (size_t i = 0; i < 3; ++i)
	{
		for (size_t j = 0; j < 3; ++j)
		{
			if (s_charts[i] == e_charts[j])
			{
				sc.push_back(s_charts[i]);
				break;
			}
		}
	}
	if (sc.size() != 2)printf("%d same chart of c-c Edge!\n", sc.size());
	int dir1 = (chart_axis(sc[0])) % 3;
	int dir2 = (chart_axis(sc[1])) % 3;
	return 3 - dir1 - dir2;
}
double Labeling::EstimateEdgeCoor(int edgeid, int dir)
{
	vector<VertexHandle>&Edge = ChartsEdges[edgeid];
	double sum = 0;
	for (size_t i = 0; i < Edge.size(); ++i)
	{
		Geometry::Vec3d _v = mesh->vertex(Edge[i]);
		sum += _v[dir];
	}
	sum /= (double)(Edge.size());
	return sum;
}
void Labeling::Get2chartsOfEdge(int edgeid, int &c1, int &c2)
{
	vector<VertexHandle>Edge = ChartsEdges[edgeid];
	VertexHandle s = Edge[0];
	VertexHandle e = Edge[Edge.size() - 1];
	vector<int>&s_charts = v_chartsid[s_vertex_id(s)];
	vector<int>&e_charts = v_chartsid[s_vertex_id(e)];
	vector<int>sc;
	c1 = -1; c2 = -1;
	for (size_t i = 0; i<3; ++i)
	{
		for (size_t j = 0; j<3; ++j)
		{
			if (s_charts[i] == e_charts[j])
			{
				sc.push_back(s_charts[i]);
				break;
			}
		}
	}
	if (sc.size() != 2)printf("%d same chart of c-c Edge!\n", sc.size());
	c1 = sc[0]; c2 = sc[1];
}
bool Labeling::has_common_area(int edgeid1, int edgeid2)
{
	//int dir1 = Nore_edgeDirection[nore_eid1];
	//int dir2 = Nore_edgeDirection[nore_eid2];
	//int edgeid1, edgeid2;
	//edgeid1 = NorepeatEdges[nore_eid1]; edgeid2 = NorepeatEdges[nore_eid2];
	int dir1 = EdgeDirection(edgeid1);
	int dir2 = EdgeDirection(edgeid2);
	assert(dir1 == dir2);
	//judge if the 2 parallel edge has common areas
	double v11, v12, v21, v22;
	v11 = mesh->vertex(ChartsEdges[edgeid1][0])[dir1];
	v12 = mesh->vertex(ChartsEdges[edgeid1][ChartsEdges[edgeid1].size() - 1])[dir1];
	v21 = mesh->vertex(ChartsEdges[edgeid2][0])[dir2];
	v22 = mesh->vertex(ChartsEdges[edgeid2][ChartsEdges[edgeid2].size() - 1])[dir2];
	double max1, max2, min1, min2;
	if (v11 < v12)
	{
		max1 = v12; min1 = v11;
	}
	else
	{
		max1 = v11; min1 = v12;
	}
	if (v21 < v22)
	{
		max2 = v22; min2 = v21;
	}
	else
	{
		max2 = v21; min2 = v22;
	}
	if (max1 < min2 || max2 < min1)
	{
		
		return false;
	}
	else{
		return true;
	}
}
void Labeling::InitChartMeanValue()
{
	chart_mean_value.resize(charts.size());
	VertexHandle v;
	Geometry::Vec3d _v;
	vector<double>chartmins;
	vector<double>chartmaxs;
	chartmins.assign(chart_mean_value.size(), 10000);
	chartmaxs.assign(chart_mean_value.size(), -10000);
	for (size_t i = 0; i < v_chartsid.size(); ++i)
	{
		v = s2v_vertex.at(i);
		_v = mesh->vertex(v);
		for (size_t j = 0; j < v_chartsid[i].size(); ++j)
		{
			int c_i = v_chartsid[i][j];
			int dir = chart_axis(c_i) % 3;
			if (_v[dir]>chartmaxs[c_i])chartmaxs[c_i] = _v[dir];
			if (_v[dir] < chartmins[c_i])chartmins[c_i] = _v[dir];
		}
	}
	//printf("%f %f\n", chartmins[0], chartmaxs[0]);
	for (size_t i = 0; i < chart_mean_value.size(); ++i)
	{
		chart_mean_value[i] = 0.5*(chartmins[i] + chartmaxs[i]);
	}
	//	//if (chart_mean_value[i]<-7)printf("********** %d\n", i);
	//	printf("init %f\n", chart_mean_value[i]);
	//}
}
void Labeling::compute_all_chart_value(double sigma_r)
{
	using namespace Eigen;
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
	while (!solve_success)
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
			//printf("%f\n", x[i]);
#endif

		}
	}
}
int Labeling::count_corner_number()
{
	int cc_count = 0;
	for (int i = 0; i < v_chartsid.size(); ++i)
	{
		if (v_chartsid[i].size() >= 3)
		{
			++cc_count;
		}
	}
	return cc_count;
}