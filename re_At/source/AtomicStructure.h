#ifndef _ATOMICSTRUCTURE_H_
#define _ATOMICSTRUCTURE_H_
#define _USE_MATH_DEFINES
#include <omp.h>
#include "tools.hpp"
#include "AtomicTypeTable.h"
#include "properties.hpp"



struct atom
{
	Vector3d coordinates;
	int id;

	atom(){}

	atom(atom&& other)
	{
		coordinates = other.coordinates;
		id = other.id;
	}

	atom(const atom& other)
	{
		coordinates = other.coordinates;
		id = other.id;
	}

	atom& operator=(atom&& other)
	{
		if (this != &other)
		{
			coordinates = other.coordinates;
			id = other.id;
		}
		return *this;
	}

	double &X = coordinates(0);
	double &Y = coordinates(1);
	double &Z = coordinates(2);
};


class AtomicStructure
{
	bool ONE_METHOD_FLAG;

	pair<double, double> x_limits_;
	pair<double, double> y_limits_;
	pair<double, double> z_limits_;
	double Translations_[3][3];
	AtomicTypeTable *atomic_type_table_;
	double section_width_;
	vector<size_t> sections_;
	size_t columns_number_;
	size_t rows_number_;
	size_t matrix_number_;
	size_t matrix_size_;

	vector<atom> atoms_;
	void find_x_limits()
	{
		x_limits_ = { atoms_[0].X, atoms_[0].X };
		for (unsigned i = 1; i < atoms_.size(); i++)
		{
			if (x_limits_.first > atoms_[i].X) x_limits_.first = atoms_[i].X;
			else if (x_limits_.second < atoms_[i].X) x_limits_.second = atoms_[i].X;
		}
	}
	void find_y_limits() {
		y_limits_ = { atoms_[0].Y, atoms_[0].Y };

		for (unsigned i = 1; i < atoms_.size(); i++)
		{
			if (y_limits_.first > atoms_[i].Y) y_limits_.first = atoms_[i].Y;
			else if (y_limits_.second < atoms_[i].Y) y_limits_.second = atoms_[i].Y;
		}
	}
	void find_z_limits() {
		z_limits_ = { atoms_[0].Z, atoms_[0].Z };
		for (unsigned i = 1; i < atoms_.size(); i++)
		{
			if (z_limits_.first > atoms_[i].Z) z_limits_.first = atoms_[i].Z;
			else if (z_limits_.second < atoms_[i].Z) z_limits_.second = atoms_[i].Z;
		}
	}

public:
	size_t size()
	{
		return atoms_.size();
	}
	AtomicStructure(AtomicTypeTable &ATtable) {
		atomic_type_table_ = &ATtable;
		ONE_METHOD_FLAG = false;
		section_width_ = 1.5;
		columns_number_ = 0;
		matrix_number_ = 0;
		matrix_size_ = 0;
	}
	void read(istream &file_in);
	void read(istream &file_in, double section_width);
	void make_sections(double width);
	void make_sections_3d(double width);
	void dist_analysis_old(vector<int> type1id, vector<int> type2id, double right,
		double deltaDist, vector<int> &dhyst, long long int &comp_it);
	void dist_analysis(vector<int> type1id, vector<int> type2id, double neighbourhood_width,
		double deltaDist, vector<int> &dhyst, long long int &comp_it);
	void replicate(double* Translation0, double* Translation1, 
	double* Translation2, int* size, double section_width);
	void create(istream &file_in , int* size);
	void create(istream &file_in , int* size, double section_width);
	void info();
	void show();
	void to_file(string output_file_name);
	void to_file_types(string output_file_name);
	void shift(double* vec);
	void plus_conflicts_old(AtomicStructure &right_hand, BinaryPropertyTable<double> &property_hard,
	BinaryPropertyTable<double> &property_soft, double right_bound);
};

#endif