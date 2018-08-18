#ifndef _ATOMICSTRUCTURE_H_
#define _ATOMICSTRUCTURE_H_

#include "tools.hpp"
#include "AtomicTypeTable.h"
#include "properties.hpp"





struct Atom
{
	int id;
	double X, Y, Z;
};

class AtomicStructure
{
private:
	pair<double, double> x_limits_;
	pair<double, double> y_limits_;
	pair<double, double> z_limits_;

	vector<int> sections_id;
public:
	vector<Atom> atoms_;
	AtomicTypeTable *atomic_type_table_;
	AtomicStructure(AtomicTypeTable &InputAtomTypeTable);
	void read(string filename);
	void create(string input, vector<int> size);
	void sort_by_sections(double section_length);
	void create_center(string input, vector<int> size);//������� �� ����� � �������� ����������
	int size() const;
	AtomicStructure& operator+= (AtomicStructure &B);
	AtomicStructure operator+ (AtomicStructure &B);
	void to_file(string output_file_name);
	void to_file(ostream &output_stream);
	void to_file(ostream &output_stream, string atom_type);
	void delete_structure();
	void transform(Eigen::Matrix3d &transform_matrix);
	void cut_by_plane(double coeff[4]);
	void replicate(double* Translation0, double* Translation1,
		double* Translation2, int* size);
	void add_shifted(const AtomicStructure &SecondStructure, double X, double Y, double Z);//��������� � ������������ ��������� ������, ������� ����������� ��������� � ������ ����� (�� ����������� ����������� ����, ������������� ����������)
	void add_shifted(const AtomicStructure &SecondStructure, double distX, double distY, double distZ, string AlongTheBorderXYZ); //��������� � ������ ��������� ������ � � ����������� �� ���������� AlongTheBorderXYZ (X, Y ��� Z) ������� �� � ������������ ������� ������������ ��������� �� ���������� ����� ��������; 
	void find_x_limits();
	void find_y_limits();
	void find_z_limits();
	void shift(double* vec);
	void shift(double X, double Y, double Z);
	void shift_min_limit(double X, double Y, double Z);
	void reproduce(int nX, int nY, int nZ, double distX, double distY, double distZ);
	void clear();
	void dist_analysis(vector<int> type1id, vector<int> type2id, double left, double right, double deltaDist, vector<int> &dhyst);
	void coordination_number_analysis(vector<int> type1id, vector<int> type2id, 
		const BinaryPropertyTable<double> &connection_table, vector<int> &CoordNumHyst);
	void coordination_number_analysis(vector<int> type1id, vector<int> type2id, double left, double right, 
		const BinaryPropertyTable<double> &connection_table, vector<int> &CoordNumHyst);
	void plus_wth_conflicts(AtomicStructure &NotmainStructure,
		BinaryPropertyTable<double> &property_hard, BinaryPropertyTable<double> &property_soft);
	void BVS(const BinaryPropertyTable<double>& r_c, const BinaryPropertyTable<double>& r_0,
		const BinaryPropertyTable<double>& b, vector<int> &BVS_hyst, vector<double> &BVS_result_for_all_atoms, double precision);
	void BVS_time_test(const BinaryPropertyTable<double>& r_c, const BinaryPropertyTable<double>& r_0,
		const BinaryPropertyTable<double>& b, vector<int> &BVS_hyst,
		vector<double> &BVS_result_for_all_atoms, double precision,
		int schedule_type, int streams_num, int schedule_sec_par);
	double GII(vector<double> &BVS_results, double precision);
	double GII(vector<double> &BVS_results, double precision, const UnaryPropertyTable<int>& valence);
	~AtomicStructure();
};

#endif