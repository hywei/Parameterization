#include "CrossParameter.h"
#include "Parameter.h"
#include "ParamPatch.h"
#include "ChartCreator.h"
#include "TransFunctor.h"
#include "Barycentric.h"
#include "../ModelMesh/MeshModel.h"
#include <fstream>

namespace PARAM
{

	CrossParameter::CrossParameter(const Parameter& parameter_1, const Parameter& parameter_2)
		: m_parameter_1(parameter_1), m_parameter_2(parameter_2){}
	CrossParameter::~CrossParameter(){}

	bool CrossParameter::LoadCorrespondingFile(const std::string& corresponding_file)
	{
		ifstream file(corresponding_file.c_str());
		if(file.fail()){
			std::cerr << "Error : Can't load corresponding file!" << std::endl;
			file.close();
			return false;
		}

		int type;
		file >> type;

		int conner_corresponding_num, edge_corresponding_num, patch_corresponding_num;
		file >> conner_corresponding_num;

		int conner_idx1, conner_idx2;
		for(int i=0; i<conner_corresponding_num; ++i)
		{
			file >> conner_idx1 >> conner_idx2;
			m_conner_correspondingAB[conner_idx1] = conner_idx2;
			m_conner_correspondingBA[conner_idx2] = conner_idx1;
		}

		file >> type;
		file >> edge_corresponding_num;
		int edge_idx1, edge_idx2;
		for(int i=0; i<edge_corresponding_num; ++i)
		{
			file >> edge_idx1 >> edge_idx2;
			m_edge_correspondingAB[edge_idx1] = edge_idx2;
			m_edge_correspondingBA[edge_idx2] = edge_idx1;
		}
		file >> type;

		file >> patch_corresponding_num;
		int patch_idx1, patch_idx2;
		for(int i=0; i<patch_corresponding_num; ++i)
		{
			file >> patch_idx1 >> patch_idx2;
			m_patch_correspondingAB[patch_idx1] = patch_idx2;
			m_patch_correspondingBA[patch_idx2] = patch_idx1;
		}

		file.close();

// 		ofstream fout("corresponding_BA.corr");
// 		fout << 1 << std::endl;
// 		fout << m_conner_correspondingBA.size() << std::endl;
// 		for(std::map<int, int>::const_iterator im= m_conner_correspondingBA.begin(); im!=m_conner_correspondingBA.end(); ++im)
// 		{
// 			fout << im->first << " " << im->second <<" ";
// 		}
// 		fout << std::endl;
// 
// 		fout << 1 << std::endl;
// 		fout << m_edge_correspondingBA.size() << std::endl;
// 		for(std::map<int, int>::const_iterator im= m_edge_correspondingBA.begin(); im!=m_edge_correspondingBA.end(); ++im)
// 		{
// 			fout << im->first << " " << im->second <<" ";
// 		}
// 		fout << std::endl;
// 
// 		fout << 1 << std::endl;
// 		fout << m_patch_correspondingBA.size() << std::endl;
// 		for(std::map<int, int>::const_iterator im=m_patch_correspondingBA.begin(); im!= m_patch_correspondingBA.end(); ++im)
// 		{
// 			fout << im->first <<" " << im->second <<" ";
// 		}
// 		fout << std::endl;
// 		fout.close();
		return true;
	}

	bool CrossParameter::GetSurfaceCoordOnA(const ChartParamCoord& chart_param_coord_onB, SurfaceCoord& surface_coord_onA) const
	{
		if(m_conner_correspondingAB.size() == 0 || m_edge_correspondingAB.size() == 0 || m_patch_correspondingBA.size() == 0)
			return false;

		ChartParamCoord chart_param_coord_onA = GetChartParamCoord4CorrespondingChartOnA(chart_param_coord_onB);
		if(!m_parameter_1.FindCorrespondingOnSurface(chart_param_coord_onA, surface_coord_onA))
		{
			printf("@@@Error: Can't find corresponding on surface A!\n");
			return false;
		}
		return true;
	}

	bool CrossParameter::GetSurfaceCoordOnB(const ChartParamCoord& chart_param_coord_onA, SurfaceCoord& surface_coord_onB) const
	{				
		if(m_conner_correspondingAB.size() == 0 || m_edge_correspondingAB.size() == 0 || m_patch_correspondingBA.size() == 0)
			return false;

		ChartParamCoord chart_param_coord_onB = GetChartParamCoord4CorrespondingChartOnB(chart_param_coord_onA);
		if(!m_parameter_2.FindCorrespondingOnSurface(chart_param_coord_onB, surface_coord_onB))
		{
			printf("@@@Error: Can't find corresponding on surface B!\n");
			return false;
		}
		return true;
	}


	void CrossParameter::FindCorrespondingAB()
	{
		const boost::shared_ptr<MeshModel> p_mesh_1 = m_parameter_1.GetMeshModel();
		int vert_num = p_mesh_1->m_Kernel.GetModelInfo().GetVertexNum();
	    m_corresponding_AB.clear();
		m_corresponding_AB.resize(vert_num);
		printf("Find Corresponding from Surface A to Surface B: ##############");

		int len = 0, progress = 0;
		m_uncorresponding_vert_array_A.clear();
		for(int vid=0; vid < vert_num; ++vid)
		{
			int chart_id = m_parameter_1.GetVertexChartID(vid);
			PARAM::ParamCoord param_coord = m_parameter_1.GetVertexParamCoord(vid);

			PARAM::SurfaceCoord surface_coord;
			if(!GetSurfaceCoordOnB(PARAM::ChartParamCoord(param_coord, chart_id), surface_coord))
			{
				m_uncorresponding_vert_array_A.push_back(vid);
				if(vid != 0){
					surface_coord = m_corresponding_AB[vid-1];
				}
			}
			m_corresponding_AB[vid] = surface_coord;

			if(len == (vert_num/100+1)){
				printf("%3d%%", ++progress);
				printf("\b\b\b\b");
				len=0;
			}
			len ++;
		}
		printf("\n");
	}

	void CrossParameter::FindCorrespondingBA()
	{
		const boost::shared_ptr<MeshModel> p_mesh_2 = m_parameter_2.GetMeshModel();
		int vert_num = p_mesh_2->m_Kernel.GetModelInfo().GetVertexNum();

		m_corresponding_BA.clear();
		m_corresponding_BA.resize(vert_num);
		printf("Find Corresponding from Surface B to Surface A: ##############");

		int len = 0, progress = 0;

		m_uncorresponding_vert_array_B.clear();
		for(int vid=0; vid < vert_num; ++vid)
		{
			int chart_id = m_parameter_2.GetVertexChartID(vid);
			PARAM::ParamCoord param_coord = m_parameter_2.GetVertexParamCoord(vid);

			PARAM::SurfaceCoord surface_coord;
			if(!GetSurfaceCoordOnA(PARAM::ChartParamCoord(param_coord, chart_id), surface_coord))
			{
				m_uncorresponding_vert_array_B.push_back(vid);
				if(vid != 0){
					surface_coord = m_corresponding_BA[vid-1];
				}
			}
			m_corresponding_BA[vid] = surface_coord;

			if(len == (vert_num/100+1)){
				printf("%3d%%", ++progress);
				printf("\b\b\b\b");
				len=0;
			}
			len ++;
		}
		printf("\n");
	}

	void CrossParameter::VertTextureTransferBA()
	{
		const boost::shared_ptr<MeshModel> p_mesh_1 = m_parameter_1.GetMeshModel();
		const boost::shared_ptr<MeshModel> p_mesh_2 = m_parameter_2.GetMeshModel();

		const PolyTexCoordArray& face_tex_array_1 = p_mesh_1->m_Kernel.GetFaceInfo().GetTexCoord();

		int vert_num_2 = p_mesh_2->m_Kernel.GetModelInfo().GetVertexNum();
	    
		m_transfer_vert_tex_array_B.clear();
		m_transfer_vert_tex_array_B.resize(vert_num_2);

		for(int vid=0; vid<vert_num_2; ++vid)
		{
			const SurfaceCoord& surface_A = m_corresponding_BA[vid];
			int fid_A = surface_A.face_index;
			const Barycentrc& bary_A = surface_A.barycentric;
			m_transfer_vert_tex_array_B[vid] = face_tex_array_1[fid_A][0] * bary_A[0] +
				face_tex_array_1[fid_A][1] * bary_A[1] + face_tex_array_1[fid_A][2] * bary_A[2];
		}
	}

	void CrossParameter::FaceTextureTransferBA()
	{
		const boost::shared_ptr<MeshModel> p_mesh_1 = m_parameter_1.GetMeshModel();
		const boost::shared_ptr<MeshModel> p_mesh_2 = m_parameter_2.GetMeshModel();

		const PolyTexCoordArray& face_tex_array_1 = p_mesh_1->m_Kernel.GetFaceInfo().GetTexCoord();

		int face_num_2 = p_mesh_2->m_Kernel.GetModelInfo().GetFaceNum();
		const PolyIndexArray& face_list_2 = p_mesh_2->m_Kernel.GetFaceInfo().GetIndex();

		m_transfer_face_tex_array_B.clear();
		m_transfer_face_tex_array_B.resize(face_num_2);

		for(int i=0; i<face_num_2; ++i)
		{
			const IndexArray& face = face_list_2[i];
			TexCoordArray& face_tex_vec = m_transfer_face_tex_array_B[i];
			face_tex_vec.clear(); face_tex_vec.resize(3);
			for(int k=0; k<3; ++k)
			{
				int vid = face[k];
				const SurfaceCoord& surface_A = m_corresponding_BA[vid];
				int fid_A = surface_A.face_index;
				const Barycentrc& bary_A = surface_A.barycentric;
				face_tex_vec[k] = face_tex_array_1[fid_A][0] * bary_A[0] +
					face_tex_array_1[fid_A][1] * bary_A[1] + face_tex_array_1[fid_A][2] *bary_A[2];
			}
		}
	}

	ChartParamCoord CrossParameter::GetChartParamCoord4CorrespondingChartOnA(const ChartParamCoord& chart_param_coord_onB) const
	{
		int chart_idx_B = chart_param_coord_onB.chart_id;
		ParamCoord param_coord_B = chart_param_coord_onB.param_coord;

		int chart_idx_A = m_patch_correspondingBA.find(chart_idx_B)->second;

		boost::shared_ptr<ChartCreator> p_chart_creator_A = m_parameter_1.GetChartCreator();
		boost::shared_ptr<ChartCreator> p_chart_creator_B = m_parameter_2.GetChartCreator();

		const ParamPatch& patch_A = (p_chart_creator_A->GetPatchArray())[chart_idx_A];
		const ParamPatch& patch_B = (p_chart_creator_B->GetPatchArray())[chart_idx_B];

		std::vector<int> _patch_conner_B(patch_B.m_conner_index_array.size());
		std::vector<int> _patch_edge_B(patch_B.m_edge_index_array.size());
		
		for(size_t k=0; k<patch_B.m_conner_index_array.size(); ++k)
		{
			int conner_idx = patch_B.m_conner_index_array[k];
			int corr_conner_idx = m_conner_correspondingBA.find(conner_idx)->second;
			_patch_conner_B[k] = corr_conner_idx;
		}

		for(size_t k=0; k<patch_B.m_edge_index_array.size(); ++k)
		{
			int edge_idx = patch_B.m_edge_index_array[k];
			int corr_edge_idx = m_edge_correspondingBA.find(edge_idx)->second;
			_patch_edge_B[k] = corr_edge_idx;
		}

		const std::vector<int>& patch_conner_A = patch_A.m_conner_index_array;
		const std::vector<int>& patch_edge_A = patch_A.m_edge_index_array;
		size_t o_pos = find(patch_conner_A.begin(), patch_conner_A.end(), _patch_conner_B[0]) - patch_conner_A.begin();
		size_t x_pos = find(patch_edge_A.begin(), patch_edge_A.end(), _patch_edge_B[0]) - patch_edge_A.begin();

		assert(o_pos != patch_conner_A.size());
		assert(x_pos != patch_edge_A.size());

		std::pair<int, int> old_x_axis(o_pos, x_pos);
		std::pair<int, int> new_x_axis(0, 0);

		TransFunctor tran_func(p_chart_creator_A);
		zjucad::matrix::matrix<double> tran_mat = tran_func.GetTransMatrixInOneChart(chart_idx_A, old_x_axis, new_x_axis);
		zjucad::matrix::matrix<double> pc_mat(3, 1);
		pc_mat(0, 0) = param_coord_B.s_coord; pc_mat(1, 0) = param_coord_B.t_coord; pc_mat(2, 0) = 1;
		zjucad::matrix::matrix<double> new_pc_mat = tran_mat * pc_mat;

		ChartParamCoord ret;
		ret.chart_id = chart_idx_A;
		ret.param_coord.s_coord = new_pc_mat(0, 0); 
		ret.param_coord.t_coord = new_pc_mat(1, 0);

		return ret;
 	}

	ChartParamCoord CrossParameter::GetChartParamCoord4CorrespondingChartOnB(const ChartParamCoord& chart_param_coord_onA) const
	{
		int chart_idx_A = chart_param_coord_onA.chart_id;
		ParamCoord param_coord_A = chart_param_coord_onA.param_coord;

// 		for(map<int, int>::const_iterator im = m_conner_correspondingAB.begin(); im != m_conner_correspondingAB.end(); ++im)
// 		{
// 			std::cout << "( " << im->first <<" " << im->second << "), ";
// 		}
// 		std::cout << std::endl;
// 
// 		for(map<int, int>::const_iterator im = m_edge_correspondingAB.begin(); im != m_edge_correspondingAB.end(); ++im)
// 		{
// 			std::cout << "( " << im->first <<" " << im->second << "), ";
// 		}
// 		std::cout << std::endl;
// 
// 		for(map<int, int>::const_iterator im = m_patch_correspondingAB.begin(); im != m_patch_correspondingAB.end(); ++im)
// 		{
// 			std::cout << "( " << im->first <<" " << im->second << "), ";
// 		}
// 		std::cout << std::endl;
// 
// 		std::cout << chart_idx_A << std::endl;

		int chart_idx_B = (m_patch_correspondingAB.find(chart_idx_A))->second;
		boost::shared_ptr<ChartCreator> p_chart_creator_A = m_parameter_1.GetChartCreator();
		boost::shared_ptr<ChartCreator> p_chart_creator_B = m_parameter_2.GetChartCreator();

		const ParamPatch& patch_A = (p_chart_creator_A->GetPatchArray())[chart_idx_A];
		const ParamPatch& patch_B = (p_chart_creator_B->GetPatchArray())[chart_idx_B];

		std::vector<int> _patch_conner_A(patch_A.m_conner_index_array.size());
		std::vector<int> _patch_edge_A(patch_A.m_edge_index_array.size());

		for(size_t k=0; k<patch_A.m_conner_index_array.size(); ++k)
		{			
			int conner_idx = patch_A.m_conner_index_array[k];
			int corr_conner_idx = m_conner_correspondingAB.find(conner_idx)->second;
			_patch_conner_A[k] = corr_conner_idx;
		}

		for(size_t k=0; k<patch_A.m_edge_index_array.size(); ++k)
		{
			int edge_idx = patch_A.m_edge_index_array[k];
			int corr_edge_idx = m_edge_correspondingAB.find(edge_idx)->second;
			_patch_edge_A[k] = corr_edge_idx;
		}


		const std::vector<int>& patch_conner_B = patch_B.m_conner_index_array;
		const std::vector<int>& patch_edge_B = patch_B.m_edge_index_array;
		size_t o_pos = find(patch_conner_B.begin(), patch_conner_B.end(), _patch_conner_A[0]) - patch_conner_B.begin();
		size_t x_pos = find(patch_edge_B.begin(), patch_edge_B.end(), _patch_edge_A[0]) - patch_edge_B.begin();

		assert(o_pos != patch_conner_B.size());
		assert(x_pos != patch_edge_B.size());

		std::pair<int, int> old_x_axis(o_pos, x_pos);
		std::pair<int, int> new_x_axis(0, 0);

		TransFunctor tran_func(p_chart_creator_B);
		zjucad::matrix::matrix<double> tran_mat = tran_func.GetTransMatrixInOneChart(chart_idx_B, old_x_axis, new_x_axis);
		zjucad::matrix::matrix<double> pc_mat(3, 1);
		pc_mat(0, 0) = param_coord_A.s_coord; pc_mat(1, 0) = param_coord_A.t_coord; pc_mat(2, 0) = 1;
		zjucad::matrix::matrix<double> new_pc_mat = tran_mat * pc_mat;

		ChartParamCoord ret;
		ret.chart_id = chart_idx_B;
		ret.param_coord.s_coord = new_pc_mat(0, 0); 
		ret.param_coord.t_coord = new_pc_mat(1, 0);

		return ret;
	}
}