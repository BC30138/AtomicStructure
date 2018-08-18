#ifndef _TEST_HPP_
#define _TEST_HPP_
#include <omp.h>
#include "AtomicStructure.h"

struct stc{
		size_t size;
		double time;
		long long int compares;  
};

void test_dist_analysis(AtomicStructure structure, AtomicTypeTable *AT, long long int comp)
{
	ofstream out("DistAnalysis_Hyst.dat");
	ofstream time("timeDistAnalysis.dat");
	double tDistAnStart, tDistAnEnd;

	time << "N_atoms:" << "\t" << structure.size() << endl;

	vector<string> type1symb;
	type1symb.push_back("Na1+");

	vector<string> type2symb;
	type2symb.push_back("Cl-1");

	vector<int> type1 = AT->string_vector_to_id_vector_wth_exist_check(type1symb);
	vector<int> type2 = AT->string_vector_to_id_vector_wth_exist_check(type2symb);

	double deltaDist = 0.1;
	vector<int> dhyst;

	tDistAnStart = omp_get_wtime();
	structure.dist_analysis_old(type1, type2, 10, deltaDist, dhyst, comp);
	// structure.dist_analysis(type1, type2, 10, deltaDist, dhyst, comp);
	tDistAnEnd = omp_get_wtime();

	time << "Time of Distance Analysis:" << "\t" << (tDistAnEnd - tDistAnStart) << endl;
	cout << structure.size() << '\t' << comp << '\n';
	plot_hyst("DistAnalysis_Hyst.dat", "DistAnalyze", "Atoms number", "Distance",
		deltaDist, dhyst, true);
	out.close();
	time.close();
}

void time_test_dist_analysis(AtomicStructure structure, AtomicTypeTable *AT, double &time_dist,long long int &comp)
{
	double tDistAnStart, tDistAnEnd;

	vector<string> type1symb;
	type1symb.push_back("Na1+");

	vector<string> type2symb;
	type2symb.push_back("Cl-1");

	vector<int> type1 = AT->string_vector_to_id_vector_wth_exist_check(type1symb);
	vector<int> type2 = AT->string_vector_to_id_vector_wth_exist_check(type2symb);

	double deltaDist = 0.1;
	vector<int> dhyst;
	

	tDistAnStart = omp_get_wtime();
	// structure.dist_analysis_old(type1, type2, 0, 10, deltaDist, dhyst, comp);
	structure.dist_analysis(type1, type2, 10, deltaDist, dhyst, comp);
	tDistAnEnd = omp_get_wtime();
	
		// plot_hyst("DistAnalysis_Hyst.dat", "DistAnalyze", "Atoms number", "Distance",
		// deltaDist, dhyst, true);

	cout << (tDistAnEnd - tDistAnStart) << "\n\n";
	time_dist = (tDistAnEnd - tDistAnStart);
}

void big_time_test_DIST(string structure_filename, int N_lim, int h)
{
	AtomicTypeTable AtomType;
	vector<stc> size_time_comp; 
	for (int N = 1; N < (N_lim); ++N)
	{		 
		
		stc single_stc; 
		ifstream in(structure_filename);
		AtomicStructure structure(AtomType);
		int size[3] = {h*N,h*N,h*N};
		structure.create(in, size, 3);
		single_stc.size = structure.size();
		time_test_dist_analysis(structure, &AtomType, single_stc.time, single_stc.compares);
		in.close();
		size_time_comp.push_back(single_stc);
	}

	ofstream out("data/time_n_comp_DIST");
	for (size_t i = 0; i < size_time_comp.size(); i++)
	{
			out << size_time_comp[i].size << '\t' << size_time_comp[i].time << '\t' << size_time_comp[i].compares << endl;
	}
	out.close();
	// log_plot_("DIST_time_analysis.dat", "DIST TIME", "t", "N atoms",
	// 	time_DIST);
}

void conflict()
{
	ifstream in("XYZ/FeClconflict.xyz");
	ifstream in_right("XYZ/HOconflict.xyz");
	AtomicTypeTable AtomType;
	AtomicStructure structure(AtomType);
	AtomicStructure structure_right(AtomType);
	BinaryPropertyTable<double> hard_conf(AtomType);
	ifstream hard_in("properties/hard_conflict_table.txt");
	hard_conf.read_property_table(hard_in);
	hard_conf.show_property_table();
	BinaryPropertyTable<double> soft_conf(AtomType);
	ifstream soft_in("properties/soft_conflict_table.txt");
	hard_conf.read_property_table(soft_in);
	int size[3] = {1,1,1};
	// int size[3] = {3,3,3};
	structure.create(in, size, 2);
	structure_right.create(in_right, size, 2);
	structure.plus_conflicts_old(structure_right,hard_conf,soft_conf,3);

	structure.show();
	AtomType.show_table();
}


#endif