#include <ctime>
#include <iostream> 
#include <fstream> 
#include <cstdlib> 
#include <cstdio> 
#include <string> 
#include <clocale> 
#include <vector>
#include <eigen3/Eigen/Dense>
#include <omp.h>
#include "ConsoleInterface.hpp"


using namespace std;
vector<AtomicStructure> create_vec_of_struct(string fileoffiles, AtomicTypeTable AtomType) {
	ifstream in(fileoffiles);
	int NF;
	in >> NF;
	string name;

	vector<AtomicStructure> FinalAtoms;
	for (int n = 0; n < NF; n++)
	{
		in >> name;
		vector<int> size;
		for (int i = 0; i < 3; i++)
		{
			double s;
			in >> s;
			size.push_back(s);
		}
		AtomicStructure File(AtomType);
		File.create(name, size);
		FinalAtoms.push_back(File);
		size.clear();
	}

	in.close();
	return FinalAtoms;
}

void to_files(vector<AtomicStructure> AtomicStuctures) {
	string xyz;
	xyz = ".xyz";
	for (int i = 0; i < (int)AtomicStuctures.size(); i++)
	{
		string res;
		res = to_string(i + 1) + xyz;
		AtomicStuctures[i].to_file(res);
	}
}

void test_coord_num_analysis(AtomicStructure structure, BinaryPropertyTable<double> connection_table)
{
	ofstream time("timeCoordNumAnalysis.dat");
	double tCoordAnStart, tCoordAnEnd;

	time << "N_atoms:" << "\t" << structure.size() << endl;

	vector<string> type1symb;
	type1symb.push_back("Na1+");

	vector<string> type2symb;
	type2symb.push_back("Cl1-");

	vector<int> type1 = structure.atomic_type_table_->string_vector_to_id_vector_wth_exist_check(type1symb);
	vector<int> type2 = structure.atomic_type_table_->string_vector_to_id_vector_wth_exist_check(type2symb);

	vector<int> CoordNumHyst;


	tCoordAnStart = omp_get_wtime();

	structure.coordination_number_analysis(type1, type2, 0, 10, connection_table, CoordNumHyst);
	//structure.CoordinationNum(type1, type2, "ConnectionTable.txt", CoordNumHyst);

	tCoordAnEnd = omp_get_wtime();
	time << "Time of Coordination Number Analysis:" << "\t" << (tCoordAnEnd - tCoordAnStart) << endl;

	plot_hyst("CoordNumAnalysis_Hyst.dat", "Координационное Число Структуры", "Количество атомов", "КЧ", 
		1.0, CoordNumHyst, false);

	time.close();
}

void test_dist_analysis(AtomicStructure structure)
{
	ofstream out("DistAnalysis_Hyst.dat");
	ofstream time("timeDistAnalysis.dat");
	double tDistAnStart, tDistAnEnd;
	
	time << "N_atoms:" << "\t" << structure.size() << endl;

	vector<string> type1symb;
	type1symb.push_back("Na1+");

	vector<string> type2symb;
	type2symb.push_back("Cl-1");

	vector<int> type1 = structure.atomic_type_table_->string_vector_to_id_vector_wth_exist_check(type1symb);
	vector<int> type2 = structure.atomic_type_table_->string_vector_to_id_vector_wth_exist_check(type2symb);

	double deltaDist = 0.1;
	vector<int> dhyst;

	tDistAnStart = omp_get_wtime();
	structure.dist_analysis(type1, type2, 0, 10, deltaDist,dhyst);

	tDistAnEnd = omp_get_wtime();
	time << "Time of Distance Analysis:" << "\t" << (tDistAnEnd - tDistAnStart) << endl;

	plot_hyst("DistAnalysis_Hyst.dat","Анализ расстояний","Количество атомов", "Расстояние (  )",
		deltaDist, dhyst, true);
	out.close();
	time.close();
}

void time_test_coord_num_analysis(AtomicStructure structure, BinaryPropertyTable<double> connection_table, 
	vector<double> &time_coord)
{

	double tCoordAnStart, tCoordAnEnd;

	vector<string> type1symb;
	type1symb.push_back("Na1+");

	vector<string> type2symb;
	type2symb.push_back("Cl1-");

	vector<int> type1 = structure.atomic_type_table_->string_vector_to_id_vector_wth_exist_check(type1symb);
	vector<int> type2 = structure.atomic_type_table_->string_vector_to_id_vector_wth_exist_check(type2symb);

	vector<int> CoordNumHyst;


	tCoordAnStart = omp_get_wtime();
	structure.coordination_number_analysis(type1, type2, 0, 10, connection_table, CoordNumHyst);
	tCoordAnEnd = omp_get_wtime();

	time_coord[structure.size()] = (tCoordAnEnd - tCoordAnStart);
}

void time_test_dist_analysis(AtomicStructure structure, vector<double> &time_dist)
{
	double tDistAnStart, tDistAnEnd;

	vector<string> type1symb;
	type1symb.push_back("Na1+");

	vector<string> type2symb;
	type2symb.push_back("Cl-1");

	vector<int> type1 = structure.atomic_type_table_->string_vector_to_id_vector_wth_exist_check(type1symb);
	vector<int> type2 = structure.atomic_type_table_->string_vector_to_id_vector_wth_exist_check(type2symb);

	double deltaDist = 0.1;
	vector<int> dhyst;

	tDistAnStart = omp_get_wtime();
	structure.dist_analysis(type1, type2, 0, 10, deltaDist, dhyst);
	tDistAnEnd = omp_get_wtime();

	time_dist[structure.size()] = (tDistAnEnd - tDistAnStart);
}

void test_BVS(AtomicStructure structure, BinaryPropertyTable<double> r_c, 
	BinaryPropertyTable<double> r_0, BinaryPropertyTable<double> b, UnaryPropertyTable<int> valence)
{
	ofstream out("BVS_hyst.dat");
	ofstream time("time_BVS_GII.dat");

	double t_bvs_start, t_bvs_end,
		t_gii_start, t_gii_end;
	vector<int> bvs_hyst;
	vector<double> bvs_GII;
	double GII_; 
	double precision_ = 0.01;


	time << "N_atoms:" << "\t" << structure.size() << endl;

	t_bvs_start = omp_get_wtime();

	structure.BVS(r_c,r_0,b,bvs_hyst, bvs_GII, precision_);

	t_bvs_end = omp_get_wtime();
	time << "BVS analysis time:" << "\t" << (t_bvs_end - t_bvs_start) << endl;

	t_gii_start = omp_get_wtime();
	GII_ = structure.GII(bvs_GII, precision_, valence);
	t_gii_end = omp_get_wtime();

	time << "GII analysis time:" << "\t" << (t_gii_end - t_gii_start) << endl;

	cout << "structure GII:" << "\t" << GII_ << endl;

	plot_hyst("BVS_hyst.dat", "BVS Структуры","количество атомов", "значение BVS", 
		precision_, bvs_hyst, true);
	//plot_hyst("BVS_results.dat", "BVS Структуры", "атомы по очередности", "BVS",
	//	1.0, bvs_GII);
	out.close();
}

void time_test_BVS(AtomicStructure structure, BinaryPropertyTable<double> r_c,
	BinaryPropertyTable<double> r_0, BinaryPropertyTable<double> b, UnaryPropertyTable<int> valence, 
	vector<double> &time_BVS, vector<double> &time_GII)
{
	double t_bvs_start, t_bvs_end,
		t_gii_start, t_gii_end;
	vector<int> bvs_hyst;
	vector<double> bvs_GII;
	double GII_;
	double precision_ = 0.01;


	t_bvs_start = omp_get_wtime();
	structure.BVS(r_c, r_0, b, bvs_hyst, bvs_GII, precision_);
	t_bvs_end = omp_get_wtime();

	time_BVS[structure.size()] = (t_bvs_end - t_bvs_start);


	t_gii_start = omp_get_wtime();
	GII_ = structure.GII(bvs_GII, precision_, valence);
	t_gii_end = omp_get_wtime();

	cout << "GII\t" << GII_ << endl;
	time_GII[structure.size()] = (t_gii_end - t_gii_start);
}

void test_BVS(AtomicStructure structure, BinaryPropertyTable<double> r_c,
	BinaryPropertyTable<double> r_0, BinaryPropertyTable<double> b)
{
	ofstream out("BVS_hyst.dat");
	ofstream time("time_BVS_GII.dat");
	double t_bvs_start, t_bvs_end,
		t_gii_start, t_gii_end;
	vector<int> bvs_hyst;
	vector<double> bvs_GII;
	double GII_;
	double precision_ = 0.01;


	time << "N_atoms:" << "\t" << structure.size() << endl;

	t_bvs_start = omp_get_wtime();

	structure.BVS(r_c, r_0, b, bvs_hyst, bvs_GII, precision_);

	t_bvs_end = omp_get_wtime();
	time << "BVS analysis time:" << "\t" << (t_bvs_end - t_bvs_start) << endl;

	t_gii_start = omp_get_wtime();
	GII_ = structure.GII(bvs_GII, precision_);
	t_gii_end = omp_get_wtime();

	time << "GII analysis time:" << "\t" << (t_gii_end - t_gii_start) << endl;

	cout << "structure GII:" << "\t" << GII_ << endl;

	plot_hyst("BVS_hyst.dat", "BVS Структуры", "количество атомов", "значение BVS",
		precision_, bvs_hyst, true);

	//plot_hyst("BVS_results.dat", "BVS Структуры", "атомы по очередности", "BVS",
		//1.0, bvs_GII);

	out.close();
	time.close();
}



void big_time_test_BVS(string structure_filename, int N_atoms_struct, int N_lim, int h)
{
	AtomicTypeTable AtomType;
	vector<double> time_BVS, time_GII;
	time_BVS.resize(N_lim*N_lim*N_lim *h*h*h *N_atoms_struct);
	time_GII.resize(N_lim*N_lim*N_lim *h*h*h *N_atoms_struct);

	BinaryPropertyTable<double> r_c(AtomType);
	BinaryPropertyTable<double> r_0(AtomType);
	BinaryPropertyTable<double> b(AtomType);
	UnaryPropertyTable<int> valence(AtomType);


	ifstream r_c_table("r_c.dat");
	ifstream r_0_table("r_0.dat");
	ifstream b_table("b.dat");
	ifstream valence_table("valence.dat");

	r_c.read_property_table(r_c_table);
	r_0.read_property_table(r_0_table);
	b.read_property_table(b_table);
	valence.read_property_table(valence_table);

	for (int N = 1; N <= N_lim; ++N)
	{
		AtomicStructure structure(AtomType);
		structure.create(structure_filename, { h*N, h*N, h*N });

		r_c.update_property_table();
		r_0.update_property_table();
		b.update_property_table();
		valence.update_property_table();

		time_test_BVS(structure, r_c, r_0, b, valence, time_BVS, time_GII);
	}

	log_plot_("BVS_time_analysis.dat", "BVS TIME", "t", "N atoms",
		time_BVS);
	//plot_hyst("GII_time_analysis.dat", "GII_TIME", "t", "N atoms",
	//	1.0, time_GII);


	r_c_table.close();
	r_0_table.close();
	b_table.close();
}

void big_time_test_COORD(string structure_filename, int N_atoms_struct, int N_lim, int h)
{
	AtomicTypeTable AtomType;
	vector<double> time_COORD;
	time_COORD.resize(N_lim*N_lim*N_lim *h*h*h *N_atoms_struct);

	BinaryPropertyTable<double> connection_property(AtomType);
	ifstream connection_table("connection_table.txt");
	connection_property.read_property_table(connection_table);

	for (int N = 1; N < (N_lim); ++N)
	{
		AtomicStructure structure(AtomType);
		structure.create(structure_filename,{ h*N, h*N, h*N });

		connection_property.update_property_table();

		time_test_coord_num_analysis(structure, connection_property, time_COORD);
	}

	log_plot_("COORD_time_analysis.dat", "COORD TIME", "t", "N atoms",
		time_COORD);

	connection_table.close();
}

void big_time_test_DIST(string structure_filename, int N_atoms_struct, int N_lim, int h)
{
	AtomicTypeTable AtomType;
	vector<double> time_DIST;

	time_DIST.resize(N_lim*N_lim*N_lim *h*h*h *N_atoms_struct);

	for (int N = 1; N < (N_lim); ++N)
	{
		AtomicStructure structure(AtomType);
		structure.create(structure_filename,{ h*N, h*N, h*N });
			
		time_test_dist_analysis(structure, time_DIST);
	}

	log_plot_("DIST_time_analysis.dat", "DIST TIME", "t", "N atoms",
		time_DIST);
}

int main()
{
	omp_set_num_threads(16);
	#pragma omp parallel
	{
		cout << "check \n";
	} 
}

// int main(int argc, char* argv[])
// {
// 	AtomicTypeTable AtomType;
// 	console mainconsole(AtomType);
// 	if (argc > 1)
// 	{
// 		bool result;
// 		ifstream script_file(argv[1]);
// 		result = mainconsole.script_mode(script_file);

// 		if (result == false);
// 		cout << "dead_end" << '\n';
// 		getc(stdin);
// 	}
// 	else
// 	{
// 		mainconsole.screen_mode();
// 	}
// 	return 0;
// }



// int main(int argc, char* argv[])
// {
// 	AtomicTypeTable TT;
// 	AtomicStructure NaCl(TT);
// 	/*NaCl.read("NaCl.xyz");
// 	Rotation a;
// 	a.inite_Euler(30, 30, 30);
// 	Matrix3d b = a.build_matrix();
// 	NaCl.transform(b);
// 	NaCl.to_file("NaCl_rot.xyz");*/
// 	plot_structure("XYZ//NaCl.xyz");
// 	cin.get();
// 	return 0;
// }











/*int main() {//BIG_TIME
	//big_time_test_BVS("NaCl.xyz", 8, 7, 5);
	//big_time_test_COORD("NaCl.xyz", 8, 8, 5);
	big_time_test_DIST("NaCl.xyz", 8, 8, 5);
	return 0;
}*/

/*
int main() {

	ofstream time("timeOfCreation.dat");
	int tBuildStart, tBuildEnd;

	tBuildStart = clock();
	AtomicTypeTable AtomType;
	BinaryPropertyTable<double> connection_property(AtomType);
	AtomicStructure NaCl(AtomType);
	NaCl.create("NaCl.xyz", { 2,2,2 });
	tBuildEnd = clock();

	ifstream connection_table("connection_table.txt");
	connection_property.update_property_table();
	connection_property.read_property_table(connection_table);
	//test_dist_analysis(NaCl);
	test_coord_num_analysis(NaCl, connection_property);

	AtomType.show_table();

	ofstream result_structure_out_Na("NaCl_result_Na.xyz");
	ofstream result_structure_out_Cl("NaCl_result_Cl.xyz");
	NaCl.to_file(result_structure_out_Na, "Na1+");
	NaCl.to_file(result_structure_out_Cl, "Cl1-");

	time << "T(sec)" << "\t" << "N_atoms = " << "\t" << NaCl.size() << endl;
	time << "Time of creation:" << "\t" << (float(tBuildEnd - tBuildStart)) / CLOCKS_PER_SEC << endl;


	//plot_structure("NaCl_result.xyz");

	time.close();
	connection_table.close();
	return 0;
}
*/

/*int main() {
	AtomicTypeTable AtomType;
	
	AtomicStructure NaCl(AtomType);
	NaCl.create("NaCl.xyz", { 25, 25, 25});

	BinaryPropertyTable<double> r_c(AtomType);
	BinaryPropertyTable<double> r_0(AtomType);
	BinaryPropertyTable<double> b(AtomType);
	UnaryPropertyTable<int> valence(AtomType);

	ifstream r_c_table("r_c.dat");
	ifstream r_0_table("r_0.dat");
	ifstream b_table("b.dat");
	ifstream valence_table("valence.dat");

	r_c.read_property_table(r_c_table);
	r_0.read_property_table(r_0_table);
	b.read_property_table(b_table);
	valence.read_property_table(valence_table);


	test_BVS(NaCl, r_c, r_0, b, valence);

	AtomType.show_table();

	r_c_table.close();
	r_0_table.close();
	b_table.close();
	return 0;
}*/


/*vector<int> size;
size.push_back(1);
size.push_back(1);
size.push_back(1);
AtomicStructure B;
B.create("in.xyz", size);

B.toFile("B.xyz");

AtomicStructure C = B + FinalAtoms;

vector<double> n;
n.push_back(1);
n.push_back(1);
n.push_back(1);
C.rotateQuater(n, 45);

C.toFile("C.xyz");
*/