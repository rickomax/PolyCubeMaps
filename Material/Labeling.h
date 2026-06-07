#ifndef LABELING_H
#define LABELING_H
#include "VolumeMeshViewer\MeshDefinition.h"
#include "VolumeMeshViewer\InteractiveWidget.h"
#include<Eigen/Dense>
#include<Eigen/Sparse>
using namespace Eigen;
using namespace OpenVolumeMesh;
enum axis
{
	px,
	py,
	pz,
	nx,
	ny,
	nz
};
class Labeling
{
public:
	Labeling();
	//Labeling(VolumeMesh* mesh_);
	~Labeling();
public:
	void prepare_data(VolumeMesh* mesh_);
	bool InitLabelMain();//directly get labels and chart boundaries
	bool LabelingMain(double area_threshold_ = 5.0, int small_chart_threshold_ = 10);

	void CorrespondSurVol();
	//corresponding of vertex, edge, face between volume and surface
	Geometry::Vec3d calculateNormals(std::vector<Geometry::Vec3d>face);
	void updateFaceNormals();
	void GetIncidence();
	//incidence of surface mesh
	bool GetPolyCubeStructure();
	bool GetChartCoor();
	
public:
	//polycube structure
	axis find_closest_axis(Geometry::Vec3d normal);
	void get_s_face_axis();
	void get_charts();
	void draw_charts();
	void straighten_bound(bool charts_flag);
	void get_charts_neighbor_();
	bool solve_chart_neighbor();
	void swallow_small_chart(int c_i);
	void merge_charts(int c_i, int n_i);
	void get_charts_bound();
	void get_chart_bound(int c_i);
	void solve_small_chart();
	void get_chartbound_faceid();
	void delete_chart(int c_i);
	void del_1neighbor_chart(int c_i);
	void del_2neighbor_chart(int c_i, int nei1, int nei2);
	void update_axis_fcid_charts(FaceHandle f, int c_i, int new_ci);
	void  facefill(std::list<FaceHandle>&queue, int c_i, int nei, int &count, std::vector<int>&fill_faces);
	axis chart_axis(int c_i);
	int faceChartId(FaceHandle f);
	//tet label
	void GetAxisAlignTetLabel();
	void GetTetLabel();
	//edge chartboundary tag
	void GetEdgeChartTag();

	//get boundary edge before chart
	void GetInitBoundaryEdge();

	//PolyCube Edge Structure
	bool extract_vbound();
	bool is_visit_vchart(int c_i, std::vector<int>&v_chart);
	void GetVChartid();
	
	bool GetEdgeStructure();
	int has_thisEdge(std::vector<VertexHandle>&Edge);
	bool is_corner(VertexHandle v);
public:
	//change conrner face label
	void ChangeCornerFaceLabel();
	void corner_tri(VertexHandle corner);
	void EstimateEdgeDir(VertexHandle v1, VertexHandle v2, Geometry::Vec3d &e, int &c1id, int &c2id);
	//change corner-extrema face label
	void ChangeCEfaceLabel();
	void corner_n_extrema(int i, int j);//just correspond as extrema[i][j]=x, extrema=charts_vbounds[i][x]
	void find_2neighbor_corners(int bid, int cid_onb, int &c1id_onb, int &c2id_onb);
	//bool WillChangeLabel(int s_fid, VertexHandle corner, VertexHandle extrema, Geometry::Vec3d projection_vec);
	int EstimateCoorDir(bool increase, int corner_bound_vid, int start_id, int end_id, int c_id, int dir);
private:
	std::vector<FaceHandle>changeface;
	std::vector<int>newchartid;
public:
	//get turning points(extremas)
	void TurningPoints(int smooth_iteration);
	bool is_vertex_extrema(int previd, int vid, int nextvid, int v_matrixid, MatrixXd &SmoothVertices);
	void get_extrema(int start_id, int end_id, int c_id, int smooth_interation);
	void smooth_chart_loopbound(int c_i,int start_id,int end_id,MatrixXd &smooth_vertices,int smooth_iteration);
	bool filt_weak_extrema(double area_threshold);
	void find_se_id(int bound_vid, int c_i, int &start_id, int &end_id);
	int count_corner_number();
private:
	//polycube structure vectors
	std::vector< axis> s_face_axis;
	std::vector<std::vector<FaceHandle> >charts;
	std::vector<std::vector<EdgeHandle> >charts_boundary;
	std::vector<std::vector<int> >charts_boundary_start_edges;
	std::vector<int>chart_bound_faceid;
	std::vector<std::vector<int>>chart_neighbor_;
	std::vector<int>face_chart_id;

	//v->charts
	std::vector<std::vector<int>>v_chartsid;
	//bound<->chart
	std::vector<std::vector<int>>chart2boundid;
	std::vector <int>bound2chartid;
	std::vector<std::vector<VertexHandle>>charts_vbounds;
	std::vector<std::vector<VertexHandle>>charts_vbound;
	//bound<->edge
	std::vector<int>edge2boundid;
	std::vector<std::vector<int>>bound2edgeid;
	std::vector<int>Nore_edgeDirection;
	std::vector<Geometry::Vec3d>Nore_edgecoor;
	std::vector<std::vector<int>>Nore_2chartsofEdge;

	std::vector<std::vector<VertexHandle>>ChartsEdges;
	std::vector<int>Edges2Norepeatid;
	std::vector<int>NorepeatEdges;
	//extremas
	
	std::vector<std::vector<int>>vbound_startid;
	std::vector<std::vector<int>>extremas;
	std::vector<std::vector<bool>>is_extrema;
	std::vector<bool>isvisit_v;
	double extrema_threshold;
public:
	//find chart pair and solve QP to get chart coordinates
	void init_chartspair_dist();
	void EdgeChartsPair();
	void ParallelEdgeChartsPair();
	int EdgeDirection(int edgeid);
	double EstimateEdgeCoor(int edgeid, int dir);
	void Get2chartsOfEdge(int edgeid, int &c1, int &c2);
	bool has_common_area(int edgeid1, int edgeid2);
	void InitChartMeanValue();
	void compute_all_chart_value(double sigma_r);
	//solve conflict
	void solve_conflict(int c_i, int c_j);
	EdgeHandle surfaceEdge_of_v12(VertexHandle v1, VertexHandle v2);

	int pairs_num;
private:
	std::vector<std::vector<double>>chartspair_dist;
	std::vector<std::vector<int>>chartpair_edges;
	std::vector<std::pair<int, int>>conflict_pair;
	//QP
	double avg_boundary_edge_length;
	std::vector<int>up_chart_id;
	std::vector<int>down_chart_id;
	std::vector<double>diff_up_down;
public:
	int s_vertex_id(VertexHandle v);
	int s_edge_id(EdgeHandle e);
	int s_face_id(FaceHandle f);
	void compute_average_bface_area();
	std::vector<FaceHandle> f_3faces(FaceHandle f);
private:
	//incidence , transform between volume and surface id
	std::vector<Geometry::Vec3d>s_face_normals;
	std::vector<int> v2s_vertexid;
	std::vector<VertexHandle> s2v_vertex;
	std::vector<int> v2s_face_id;
	std::vector<FaceHandle> s2v_face_;
	std::map<int, EdgeHandle> maps2v_edge;
	std::map<int, int> mapv2s_edge;

	std::vector<std::vector<VertexHandle>> b_face2vertices;
	std::vector<std::vector<FaceHandle> > b_edge2faces_;
	std::vector<std::vector<EdgeHandle> > b_face2edges;
	std::vector < std::vector<int>> b_vertex2facesid;
	
	std::vector<std::vector<VertexHandle>>incidences_of_v;
	std::vector<OpenVolumeMesh::CellHandle>b_face2cell;
	
	//mesh info
private:
	VolumeMesh* mesh;
	SurfaceMesh boundaryMesh;
	//InteractiveWidget *VolumeMeshViewer;
	int n_boundary_faces;
	int n_boundary_edges;
	int n_boundary_vertices;
	double avg_bface_area;
	std::vector<OpenVolumeMesh::FaceHandle> SpecialFaces;
	std::vector<VertexHandle> specialPoints;

	//some needed input
	public:
		double qpsigma;//0.8
		int smooth_iteration;// 5
		int small_chart_threshold;//10;
		double area_threshold;//5.0

	//some needed output
	public:
		std::vector<double>chart_mean_value;
		std::vector<int>axis_align_tetLabel;
		std::vector<int>tetLabel;
		std::vector<int>tet_face_chart_id;
		std::vector<Eigen::Triplet<double> >ChartPairDist;
		std::vector<int>Init_EdgeBoudnaryTag;
		std::vector<int>EdgeBoundaryTag;
	//labeling consequence
	private:
		//bool checkStatus();
		bool filtered_extrema;//=true if all of extrema_area/avgbound_facearea <area_threshold
		bool swallowed_chart;//=true if all of charts has neighbors>=4
};

#endif