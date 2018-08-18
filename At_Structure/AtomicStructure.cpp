#include "AtomicStructure.h"
#include "properties.hpp"
#include <fcntl.h>
#include <omp.h>
using namespace std;
using namespace Eigen;

AtomicStructure::AtomicStructure(AtomicTypeTable &ATtable) {
	atomic_type_table_ = &ATtable;
}

void AtomicStructure::read(string file_name)
{
	ifstream in("XYZ\\" + file_name);
	bool N_atoms_error = false;
	string N_atoms_str;
	getline(in, N_atoms_str);
	
	if (!N_atoms_str.empty())
	{
		unsigned char_it = 0;
		while (!N_atoms_error && (char_it < (N_atoms_str.size())))
		{
			if (!isdigit(N_atoms_str[char_it]))
				N_atoms_error = true;
			else ++char_it;
		}
	}
	else
	{
		N_atoms_error = true;;
	}

	if (!N_atoms_error)
	{
		int N_atoms = atof(N_atoms_str.c_str());
		string Symbol;
		getline(in, Symbol);
		Atom t;
		for (int i = 0; i < N_atoms; i++)
		{
			in >> Symbol;
			t.id = atomic_type_table_->get_Id_with_adding(Symbol);
			in >> t.X;
			in >> t.Y;
			in >> t.Z;
			atoms_.push_back(t);
		}

		find_x_limits();
		find_y_limits();
		find_z_limits();
	}
	else
	{
		cout << "ERROR in file '" << file_name << "': first line for number of atoms." << '\n';
	}
}

void AtomicStructure::create(string filename , vector<int> size)
{
	cout << "Creating structure from file..." << endl;
	ifstream in("XYZ\\" + filename);
	
	int N_all, N_tr;
	in >> N_all >> N_tr;

	int N_atoms = N_all - 3;


	vector<vector<double>> Tr;
	Tr.resize(3);
	Tr[0].resize(3);
	Tr[1].resize(3);
	Tr[2].resize(3);

	string Symbol;
	Atom t;
	
	for (int i = 0; i < N_tr; i++)
	{
		in >> Symbol;
		in >> Tr[i][0];
		in >> Tr[i][1];
		in >> Tr[i][2];
	}

	for (int i = 0; i < N_atoms; i++)
	{
		in >> Symbol;
		t.id = atomic_type_table_->get_Id_with_adding(Symbol);
		in >> t.X;
		in >> t.Y;
		in >> t.Z;
		atoms_.push_back(t);
	}

	double translationX = Tr[0][0];
	double translationY = Tr[0][1];
	double translationZ = Tr[0][2];
	
	if (size[0] != 0)
	{
		if (is_equal(translationX, translationY) && is_equal(translationY, translationZ) &&
			is_equal(translationX, 0))
		{
			
		}
		else
		{
			int N_atoms = atoms_.size();
			for (int i = 1; i < size[0]; i++)
			{
				for (int j = 0; j < N_atoms; j++)
				{
					Atom t;
					t.id = atoms_[j].id;
					t.X = atoms_[j].X + i*translationX;
					t.Y = atoms_[j].Y + i*translationY;
					t.Z = atoms_[j].Z + i*translationZ;
					atoms_.push_back(t);
				}
			}
		}
	}
	
	translationX = Tr[1][0];
	translationY = Tr[1][1];
	translationZ = Tr[1][2];

	if (size[1] != 0)
	{
		if (is_equal(translationX, translationY) && is_equal(translationY, translationZ) &&
			is_equal(translationX, 0))
		{

		}
		else
		{
			int N_atoms = atoms_.size();
			for (int i = 1; i < size[1]; i++)
			{
					for (int j = 0; j < N_atoms; j++)
					{
						Atom t;
						t.id = atoms_[j].id;
						t.X = atoms_[j].X + i*translationX;
						t.Y = atoms_[j].Y + i*translationY;
						t.Z = atoms_[j].Z + i*translationZ;
						atoms_.push_back(t);
					}
			}
		}
	}

	translationX = Tr[2][0];
	translationY = Tr[2][1];
	translationZ = Tr[2][2];

	if (size[2] != 0)
	{
		if (is_equal(translationX, translationY) && is_equal(translationY, translationZ) &&
			is_equal(translationX, 0))
		{

		}
		else
		{
			int N_atoms = atoms_.size();
			for (int i = 1; i < size[2]; i++)
			{
					for (int j = 0; j < N_atoms; j++)
					{
						Atom t;
						t.id = atoms_[j].id;
						t.X = atoms_[j].X + i*translationX;
						t.Y = atoms_[j].Y + i*translationY;
						t.Z = atoms_[j].Z + i*translationZ;
						atoms_.push_back(t);
					}
			}
		}
	}

	find_x_limits();
	find_y_limits();
	find_z_limits();
	in.close();
}

void AtomicStructure::create_center(string filename, vector<int> size)
{
	cout << "Creating structure from file..." << endl;
	ifstream in("XYZ\\" + filename);
	int N_all, N_tr;
	in >> N_all >> N_tr;

	N_tr = 3;
	int N_atoms = N_all - N_tr;


	vector<vector<double>> Tr;
	Tr.resize(3);
	Tr[0].resize(3);
	Tr[1].resize(3);
	Tr[2].resize(3);

	string Symbol;
	for (int i = 0; i < N_tr; i++)
	{
		in >> Symbol;
		in >> Tr[i][0];
		in >> Tr[i][1];
		in >> Tr[i][2];
	}

	for (int i = 0; i < N_atoms; i++)
	{
		Atom t;
		in >> Symbol;
		t.id = atomic_type_table_->get_Id_with_adding(Symbol);
		in >> t.X;
		in >> t.Y;
		in >> t.Z;
		atoms_.push_back(t);
	}

	double translationX = Tr[0][0];
	double translationY = Tr[0][1];
	double translationZ = Tr[0][2];

	if (size[0] != 0)
	{
		if (is_equal(translationX, translationY) && is_equal(translationY, translationZ) &&
			is_equal(translationX, 0))
		{

		}
		else
		{
			int N_atoms = atoms_.size();
			for (int i = 1; i < size[0]; i++)
			{
				for (int j = 0; j < N_atoms; j++)
				{
					Atom t;
					t.id = atoms_[j].id;
					t.X = atoms_[j].X + i*translationX;
					t.Y = atoms_[j].Y + i*translationY;
					t.Z = atoms_[j].Z + i*translationZ;
					atoms_.push_back(t);
				}
			}
		}
	}

	translationX = Tr[1][0];
	translationY = Tr[1][1];
	translationZ = Tr[1][2];

	if (size[1] != 0)
	{
		if (is_equal(translationX, translationY) && is_equal(translationY, translationZ) &&
			is_equal(translationX, 0))
		{

		}
		else
		{
			int N_atoms = atoms_.size();
			for (int i = 1; i < size[1]; i++)
			{
				for (int j = 0; j < N_atoms; j++)
				{
					Atom t;
					t.id = atoms_[j].id;
					t.X = atoms_[j].X + i*translationX;
					t.Y = atoms_[j].Y + i*translationY;
					t.Z = atoms_[j].Z + i*translationZ;
					atoms_.push_back(t);
				}
			}
		}
	}

	translationX = Tr[2][0];
	translationY = Tr[2][1];
	translationZ = Tr[2][2];

	if (size[2] != 0)
	{
		if (is_equal(translationX, translationY) && is_equal(translationY, translationZ) &&
			is_equal(translationX, 0))
		{

		}
		else
		{
			int N_atoms = atoms_.size();
			for (int i = 1; i < size[2]; i++)
			{
				for (int j = 0; j < N_atoms; j++)
				{
					Atom t;
					t.id = atoms_[j].id;
					t.X = atoms_[j].X + i*translationX;
					t.Y = atoms_[j].Y + i*translationY;
					t.Z = atoms_[j].Z + i*translationZ;
					atoms_.push_back(t);
				}
			}
		}
	}

	find_x_limits();
	find_y_limits();
	find_z_limits();
	in.close();

	double centerX = x_limits_.second - fabs(x_limits_.second - x_limits_.first) / 2;
	double centerY = y_limits_.second - fabs(y_limits_.second - y_limits_.first) / 2;
	double centerZ = z_limits_.second - fabs(z_limits_.second - z_limits_.first) / 2;

	for (int i = 0; i < atoms_.size(); i++)
	{
		atoms_[i].X -= centerX;
		atoms_[i].Y -= centerY;
		atoms_[i].Z -= centerZ;
	}

}

void AtomicStructure::find_x_limits()
{
	x_limits_ = {atoms_[0].X, atoms_[0].X};
	for (int i = 1; i < atoms_.size(); i++)
	{		
		if (x_limits_.first > atoms_[i].X) x_limits_.first = atoms_[i].X;
		else if (x_limits_.second < atoms_[i].X) x_limits_.second = atoms_[i].X;
	}
}

void AtomicStructure::find_y_limits() {
	y_limits_ = { atoms_[0].Y, atoms_[0].Y };

	for (int i = 1; i < atoms_.size(); i++)
	{
		if (y_limits_.first > atoms_[i].Y) y_limits_.first = atoms_[i].Y;
		else if (y_limits_.second < atoms_[i].Y) y_limits_.second = atoms_[i].Y;
	}
}

void AtomicStructure::find_z_limits() {
	z_limits_ = {atoms_[0].Z, atoms_[0].Z };
	for (int i = 1; i < atoms_.size(); i++)
	{
		if (z_limits_.first > atoms_[i].Z) z_limits_.first = atoms_[i].Z;
		else if (z_limits_.second < atoms_[i].Z) z_limits_.second = atoms_[i].Z;
	}
}

int AtomicStructure::size() const
{
	return atoms_.size();
}

AtomicStructure& AtomicStructure::operator+= (AtomicStructure &B)
{
	cout << "Summing..." << endl; 

	for (int i = 0; i < B.size(); i++)
	{
		atoms_.push_back(B.atoms_[i]);
		if (i % 524288 == 0) cout << "*";
	}
	find_x_limits();
	find_y_limits();
	find_z_limits();
	cout << "\n";
	return *this;
}

AtomicStructure AtomicStructure::operator+ ( AtomicStructure &B)
{
	return *this += B;
}

void AtomicStructure::to_file(string output_file_name)
{
	cout << "Export to file..." << endl;
	ofstream out("XYZ\\" + output_file_name);
	out << /*"#" <<*/ atoms_.size() << "\n\n";

	for (int n = 0; n < atoms_.size(); n++)
	{
		out << longid_to_string(atomic_type_table_->get_element(atoms_[n].id).long_id) << "\t" << atoms_[n].X << "\t" << atoms_[n].Y << "\t" << atoms_[n].Z << "\n";
		if (n % 524288 == 0) cout << "*";
	}
	cout << "\n";
	out << "\n\n";
	out.close();
}

void AtomicStructure::to_file(ostream &out)
{
	cout << "Export to file..." << endl;
	out << /*"#" <<*/ atoms_.size() << "\n\n";

	for (int n = 0; n < atoms_.size(); n++)
	{
		out << longid_to_string(atomic_type_table_->get_element(atoms_[n].id).long_id) << "\t" << atoms_[n].X << "\t" << atoms_[n].Y << "\t" << atoms_[n].Z << "\n";
		if (n % 524288 == 0) cout << "*";
	}
	cout << "\n";
	out << "\n\n";
}

void AtomicStructure::to_file(ostream& out, string atom_type)
{
	cout << "Export to file..." << endl;
	int id_atom_type = atomic_type_table_->get_id(atom_type);

	vector<int> atom_type_vec;

	if (id_atom_type != -1)
	{
		for (int n = 0; n < atoms_.size(); n++)
		{
			if (atoms_[n].id == id_atom_type)
			{
				atom_type_vec.push_back(n);
			}
		}
	}
	else cout << "Export to file error: No such atom type in structure \n";

	out << /*"#" <<*/ atom_type_vec.size() << "\n\n";

	for (unsigned long i = 0; i < atom_type_vec.size(); ++i)
	{
		out << longid_to_string(atomic_type_table_->get_element(atoms_[atom_type_vec[i]].id).long_id)
			<< "\t" << atoms_[atom_type_vec[i]].X 
			<< "\t" << atoms_[atom_type_vec[i]].Y << "\t" << atoms_[atom_type_vec[i]].Z << "\n";
	}
	

	cout << "\n";
	out << "\n\n";
}

void AtomicStructure::delete_structure()
{
	atoms_.clear();
}

void AtomicStructure::replicate(double* Translation0, double* Translation1, 
	double* Translation2, int* size)
{
	double translationX = Translation0[0];
	double translationY = Translation0[1];
	double translationZ = Translation0[2];

	if (size[0] != 0)
	{
		if (is_equal(translationX, translationY) && is_equal(translationY, translationZ) &&
			is_equal(translationX, 0))
		{

		}
		else
		{
			int N_atoms = atoms_.size();
			for (int i = 1; i < size[0]; i++)
			{
				for (int j = 0; j < N_atoms; j++)
				{
					Atom t;
					t.id = atoms_[j].id;
					t.X = atoms_[j].X + i*translationX;
					t.Y = atoms_[j].Y + i*translationY;
					t.Z = atoms_[j].Z + i*translationZ;
					atoms_.push_back(t);
				}
			}
		}
	}

	translationX = Translation1[0];
	translationY = Translation1[1];
	translationZ = Translation1[2];

	if (size[1] != 0)
	{
		if (is_equal(translationX, translationY) && is_equal(translationY, translationZ) &&
			is_equal(translationX, 0))
		{

		}
		else
		{
			int N_atoms = atoms_.size();
			for (int i = 1; i < size[1]; i++)
			{
				for (int j = 0; j < N_atoms; j++)
				{
					Atom t;
					t.id = atoms_[j].id;
					t.X = atoms_[j].X + i*translationX;
					t.Y = atoms_[j].Y + i*translationY;
					t.Z = atoms_[j].Z + i*translationZ;
					atoms_.push_back(t);
				}
			}
		}
	}

	translationX = Translation2[0];
	translationY = Translation2[1];
	translationZ = Translation2[2];

	if (size[2] != 0)
	{
		if (!(is_equal(translationX, translationY) && is_equal(translationY, translationZ) &&
			is_equal(translationX, 0)))
		{
			int N_atoms = atoms_.size();
			for (int i = 1; i < size[2]; i++)
			{
				for (int j = 0; j < N_atoms; j++)
				{
					Atom t;
					t.id = atoms_[j].id;
					t.X = atoms_[j].X + i*translationX;
					t.Y = atoms_[j].Y + i*translationY;
					t.Z = atoms_[j].Z + i*translationZ;
					atoms_.push_back(t);
				}
			}
		}
	}

	find_x_limits();
	find_y_limits();
	find_z_limits();
}

void AtomicStructure::shift(double* vec) {
	cout << "Shifting..." << endl;
	double shiftX = vec[0];
	double shiftY = vec[1];
	double shiftZ = vec[2];
	for (int i = 0; i < atoms_.size(); i++)
	{
		atoms_[i].X += shiftX;
		atoms_[i].Y += shiftY;
		atoms_[i].Z += shiftZ;
		if (i % 524288 == 0) cout << "*";
	}
	cout << "\n";
	find_x_limits();
	find_y_limits();
	find_z_limits();
}

void AtomicStructure::shift(double X, double Y, double Z) {
	cout << "Shifting..." << endl;
	for (int i = 0; i < atoms_.size(); i++)
	{
		atoms_[i].X += X;
		atoms_[i].Y += Y;
		atoms_[i].Z += Z;
		if (i % 524288 == 0) cout << "*";
	}
	cout << "\n";
	find_x_limits();
	find_y_limits();
	find_z_limits();
}

void AtomicStructure::shift_min_limit(double X, double Y, double Z) {
	cout << "Shifting..." << endl;
	double shiftX = X - x_limits_.first;
	double shiftY = Y - y_limits_.first;
	double shiftZ = Z - z_limits_.first;

	for (int i = 0; i < atoms_.size(); i++)
	{
		atoms_[i].X += shiftX;
		atoms_[i].Y += shiftY;
		atoms_[i].Z += shiftZ;
		if (i % 524288 == 0) cout << "*";
	}
	cout << "\n";
	find_x_limits();
	find_y_limits();
	find_z_limits();
}

void AtomicStructure::clear() {
		atoms_.clear();
}

void AtomicStructure::reproduce(int nX, int nY, int nZ, double distX, double distY, double distZ)
{
	cout << "Reproducing..." << endl;
	if (nX != 0)
	{
		int N_atoms = atoms_.size();
		for (int i = 1; i < nX; i++)
		{
				for (int j = 0; j < N_atoms; j++)
				{
					Atom t;
					t.id = atoms_[j].id;
					t.X = atoms_[j].X + i*(distX);
					t.Y = atoms_[j].Y;
					t.Z = atoms_[j].Z;
					atoms_.push_back(t);
				}
		}
	}
	
	if (nY != 0)
	{
		int N_atoms = atoms_.size();
		for (int i = 1; i < nY; i++)
		{
				for (int j = 0; j < N_atoms; j++)
				{
					Atom t;
					t.id = atoms_[j].id;
					t.X = atoms_[j].X;
					t.Y = atoms_[j].Y + i*(distY);
					t.Z = atoms_[j].Z;
					atoms_.push_back(t);
				}
		}
	}

	if (nZ != 0)
	{
		int N_atoms = atoms_.size();
		for (int i = 1; i < nZ; i++)
		{
				for (int j = 0; j < N_atoms; j++)
				{
					Atom t;
					t.id = atoms_[j].id;
					t.X = atoms_[j].X;
					t.Y = atoms_[j].Y;
					t.Z = atoms_[j].Z + i*(distZ);
					atoms_.push_back(t);
				}
		}
	}
	find_x_limits();
	find_y_limits();
	find_z_limits();
}

void AtomicStructure::dist_analysis(vector<int> type1id, vector<int> type2id, double left, double right,
	double deltaDist, vector<int> &dhyst) {
	cout << "!!!Distance Analysis!!!" << endl;
	//ofstream plot("plotAnalysisDist.dat");
	//plot.precision(5)

	vector<int> type1atoms;
	vector<int> type2atoms;


	cout << "Search of type 1 atoms..." << endl;
	for (int i = 0; i < type1id.size(); i++)
	{
		for (int n = 0; n < atoms_.size(); n++)
		{
			if (atoms_[n].id == type1id[i]) type1atoms.push_back(n);
		}
	}

	cout << "Search of type 2 atoms..." << endl;
	for (int i = 0; i < type2id.size(); i++)
	{
		for (int n = 0; n < atoms_.size(); n++)
		{
			if (atoms_[n].id == type2id[i]) type2atoms.push_back(n);
		}
	}

	cout << "Calculating Process..." << endl;


	double squareleft, squareright;
	squareleft = Square(left);
	squareright = Square(right);

	dhyst.resize(int((right - left) / deltaDist) + 1);
	for (int i = 0; i < dhyst.size(); i++)
	{
		dhyst[i] = 0;
	}

	int block;
	//int i = 0;

	omp_set_num_threads(16);
#pragma omp parallel
	{
#pragma omp for schedule(static, 1000)
		for (int n1 = 0; n1 < type1atoms.size(); n1++)
		{
			for (int n2 = 0; n2 < type2atoms.size(); n2++)
			{
				double sqrdist = Square(atoms_[type1atoms[n1]].X - atoms_[type2atoms[n2]].X)
					+ Square(atoms_[type1atoms[n1]].Y - atoms_[type2atoms[n2]].Y)
					+ Square(atoms_[type1atoms[n1]].Z - atoms_[type2atoms[n2]].Z);
				if (sqrdist >= squareleft && sqrdist <= squareright)
				{
					block = int((sqrt(sqrdist) - left) / deltaDist);
					dhyst[block]++;
					//i++;
					//if (i % 524288 == 0) cout << "*";
				}
			}
		}
	}
	cout << "\n";


	
//	plot << "width = " << deltaDist << endl;
//	plot << "bin(x, s) = s*int(x / s)" << endl;
	//plot << "set xtics " << deltaDist << endl;
	/*plot << "set boxwidth 1" << endl;
	plot << "plot 'AnalysisDist.dat' f w boxes fs solid 0.5 title '�����������'" << endl;
	plot << "set yrange[0:GPVAL_Y_MAX + 1]" << endl;
	plot << "set xrange[GPVAL_X_MIN - width:GPVAL_X_MAX + width]" << endl;
	plot << "replot" << endl;
	plot << "pause - 1" << endl;
	*/
	//system("gnuplot plotAnalysisDist.dat");
}

void AtomicStructure::coordination_number_analysis(vector<int> type1id, vector<int> type2id, double left, double right, 
	const BinaryPropertyTable<double> &connection_table, vector<int> &CoordNumHyst) {
	cout << "!!!Coordination Number Analysis!!!" << endl;

	
	int size_ = atomic_type_table_->size();
	BinaryPropertyTable<double> square_connection_table([](double x) {
		return x*x;
	}, connection_table);

	vector<int> type1atoms;
	vector<int> type2atoms;

	cout << "Search of type 1 atoms..." << endl;
	for (int i = 0; i < type1id.size(); i++)
	{
		for (int n = 0; n < atoms_.size(); n++)
		{
			if (atoms_[n].id == type1id[i]) type1atoms.push_back(n);
		}
	}


	cout << "Search of type 2 atoms..." << endl;
	for (int i = 0; i < type2id.size(); i++)
	{
		for (int n = 0; n < atoms_.size(); n++)
		{
			if (atoms_[n].id == type2id[i]) type2atoms.push_back(n);
		}
	}


	double sqrleft, sqrright;
	sqrleft = Square(left);
	sqrright = Square(right);

	CoordNumHyst.clear();
	CoordNumHyst.resize(type2atoms.size());
	
	for (int i = 0; i < CoordNumHyst.size(); i++)
	{
		CoordNumHyst[i] = 0;
	}

	//int i = 0;
	double Dist;
	cout << "Calculating Process..." << endl;
	
	int CoordNumMax = 0;
	omp_set_num_threads(16);
#pragma omp parallel 
	{
#pragma omp for schedule(static, 1000)
		for (int n1 = 0; n1 < type1atoms.size(); n1++)
		{
			int CoordNum = 0;
			for (int n2 = 0; n2 < type2atoms.size(); n2++)
			{
				double sqrdist = Square(atoms_[type1atoms[n1]].X - atoms_[type2atoms[n2]].X)
					+ Square(atoms_[type1atoms[n1]].Y - atoms_[type2atoms[n2]].Y)
					+ Square(atoms_[type1atoms[n1]].Z - atoms_[type2atoms[n2]].Z);
				if (sqrdist >= sqrleft && sqrdist <= sqrright) {
					double TableVal = square_connection_table.get_property(atoms_[type2atoms[n2]].id, atoms_[type1atoms[n1]].id);
					if (sqrdist <= TableVal)
					{
						CoordNum++;
					}
				}

				//i++;
				//33554432
				//if (i % 1000000 == 0) { 
				//	cout << "*"; 
				//};
			}
			CoordNumHyst[CoordNum]++;
			if (CoordNumMax < CoordNum) CoordNumMax = CoordNum;
		}
	}
	cout << "\n";

	while (CoordNumHyst.size() > (CoordNumMax + 1)) CoordNumHyst.erase(CoordNumHyst.end() - 1);
}

void AtomicStructure::coordination_number_analysis(vector<int> type1id, vector<int> type2id,
	const BinaryPropertyTable<double> &connection_table, vector<int> &CoordNumHyst) {
	cout << "!!!Coordination Number Analysis!!!" << endl;
	vector<int> type1atoms;
	vector<int> type2atoms;

	int size_ = atomic_type_table_->size();
	BinaryPropertyTable<double> square_connection_table([](double x) {
		return x*x;
	}, connection_table);


	cout << "Search of type 1 atoms..." << endl;
	for (int i = 0; i < type1id.size(); i++)
	{
		for (int n = 0; n < atoms_.size(); n++)
		{
			if (atoms_[n].id == type1id[i]) type1atoms.push_back(n);
		}
	}

	cout << "Search of type 2 atoms..." << endl;
	for (int i = 0; i < type2id.size(); i++)
	{
		for (int n = 0; n < atoms_.size(); n++)
		{
			if (atoms_[n].id == type2id[i]) type2atoms.push_back(n);
		}
	}


	CoordNumHyst.clear();
	CoordNumHyst.resize(type2atoms.size());

	for (int i = 0; i < CoordNumHyst.size(); i++)
	{
		CoordNumHyst[i] = 0;
	}

	int i = 0;
	double Dist;
	cout << "Calculating Process..." << endl;
	for (int n1 = 0; n1 < type1atoms.size(); n1++)
	{
		int CoordNum = 0;
		CoordNumHyst.push_back(0);
		for (int n2 = 0; n2 < type2atoms.size(); n2++)
		{
			double sqrdist = Square(atoms_[type1atoms[n1]].X - atoms_[type2atoms[n2]].X)
				+ Square(atoms_[type1atoms[n1]].Y - atoms_[type2atoms[n2]].Y)
				+ Square(atoms_[type1atoms[n1]].Z - atoms_[type2atoms[n2]].Z);
				double TableVal = square_connection_table.get_property(atoms_[type2atoms[n2]].id, atoms_[type1atoms[n1]].id);
				if (sqrdist <= TableVal)
				{
					CoordNum++;
				}
			i++;
			if (i % 33554432 == 0) { cout << "*"; };
		}
		CoordNumHyst[CoordNum]++;
	}
	cout << "\n";
}

void AtomicStructure::add_shifted(const AtomicStructure &SecondStructure, double X, double Y, double Z) {
	AtomicStructure TemperaryStructure(*atomic_type_table_);
	TemperaryStructure = SecondStructure;
	TemperaryStructure.shift_min_limit(X, Y, Z);
	*this += TemperaryStructure;
}

void AtomicStructure::add_shifted(const AtomicStructure &SecondStructure, double distX, double distY, 
	double distZ, string AlongTheBorderXYZ) {
	string Xaxis = "X";
	string Yaxis = "Y";
	string Zaxis = "Z";

	AtomicStructure TemperaryStructure(*atomic_type_table_);
	TemperaryStructure = SecondStructure;

	if (AlongTheBorderXYZ == Xaxis)
	{
		double lengthX = fabs(x_limits_.second - x_limits_.first);
		TemperaryStructure.shift_min_limit(x_limits_.first + lengthX + distX, y_limits_.first + distY, z_limits_.first + distZ);
		*this += TemperaryStructure;
	}
	else
		if (AlongTheBorderXYZ == Yaxis)
		{
			double lengthY = fabs(y_limits_.second - y_limits_.first);
			TemperaryStructure.shift_min_limit(x_limits_.first + distX, y_limits_.first + lengthY + distY, z_limits_.first + distZ);
			*this += TemperaryStructure;
		}
		else
			if (AlongTheBorderXYZ == Zaxis)
			{
				double lengthZ = fabs(z_limits_.second - z_limits_.first);
				TemperaryStructure.shift_min_limit(x_limits_.first + distX, y_limits_.first + distY, z_limits_.first + lengthZ + distZ);
				*this += TemperaryStructure;
			}
			else
				cout << "ERROR in setting the boundary";

}

void AtomicStructure::plus_wth_conflicts(AtomicStructure &other,BinaryPropertyTable<double> &property_hard,
	BinaryPropertyTable<double> &property_soft)
{
	cout << "!!!Plus wth Conflicts!!!" << endl;

	int NotExcludedCounter = 0;
	double sqrDist;
	int sizeofMainStruct = atoms_.size();
	int sizeofAtTpTbl = atomic_type_table_->size();
	int sizeofother = other.size();

	Eigen::MatrixXi warningcounter(sizeofAtTpTbl, sizeofAtTpTbl);
	
	for (unsigned long i = 0; i < sizeofAtTpTbl; ++i) {
		for (unsigned long j = 0; j < sizeofAtTpTbl; ++j)
		{
			warningcounter(i, j) = 0;
		}
	}


	double hard_value;
	double soft_value;
	int time = 0;
	for (unsigned long j = 0; j < sizeofother; ++j)
		{
		for (unsigned long i = 0; i < sizeofMainStruct; ++i)
			{
			sqrDist = Square(atoms_[i].X - other.atoms_[j].X)
				+ Square(atoms_[i].Y - other.atoms_[j].Y)
				+ Square(atoms_[i].Z - other.atoms_[j].Z);
			hard_value = property_hard.get_property(other.atoms_[j].id, atoms_[i].id);
			soft_value = property_soft.get_property(other.atoms_[j].id, atoms_[i].id);
			time++;
			if (time % 54432 == 0) cout << '*';

			if (sqrDist > hard_value)
			{
				++NotExcludedCounter;
				if (sqrDist < soft_value)
				{
					++warningcounter(atoms_[i].id, other.atoms_[j].id);
				}
				if (NotExcludedCounter == sizeofMainStruct - 1)
				{
					atoms_.push_back(other.atoms_[j]);
				}
			}
			else
			{
				i = sizeofMainStruct;
				for (unsigned n = 0; n < sizeofAtTpTbl; ++n)
				{
					warningcounter(n, other.atoms_[j].id) = 0;
				}
			}
		}
	}

	for (unsigned long i = 0; i < sizeofAtTpTbl; ++i)
		for (unsigned long j = 0; j < sizeofAtTpTbl; ++j)
		{
			if (warningcounter(i, j) != 0)
			{
				cout << longid_to_string(atomic_type_table_->get_element(i).long_id) << "\t" << longid_to_string(atomic_type_table_->get_element(j).long_id) << "\t"
					<< warningcounter(i, j) << endl;
			}
		}
	system("pause");
	
	find_x_limits();
	find_y_limits();
	find_z_limits();
}

void AtomicStructure::BVS(const BinaryPropertyTable<double>& r_c,const BinaryPropertyTable<double>& r_0,
	const BinaryPropertyTable<double>& b, vector<int> &BVS_hyst, 
	vector<double> &BVS_result_for_all_atoms, double precision) {
	cout << "!!!BVS Analysis!!!" << endl;
	cout << "Calculating Process..." << endl;
	int size_ = atoms_.size();

	BinaryPropertyTable<double> square_r_c([](double x) {
	return x*x;
	}, r_c);
	
	BVS_hyst.resize(int(20/precision) + 1);
	BVS_result_for_all_atoms.resize(size_);
	for (int i = 0; i < size_; i++)
	{
		BVS_result_for_all_atoms[i] = 0;
	}
	for (int i = 0; i < BVS_hyst.size(); i++)
	{
		BVS_hyst[i] = 0;
	}
	
	

	int block;
	int n = 0;

omp_set_num_threads(16);
#pragma omp parallel 
{
#pragma omp for schedule(dynamic, 100)
	for (int i = 0; i < size_; i++)
	{
		for (int j = i; j < size_; j++)
		{
			double dist = Square(atoms_[i].X - atoms_[j].X) 
				+ Square(atoms_[i].Y - atoms_[j].Y)
				+ Square(atoms_[i].Z - atoms_[j].Z);
			if (dist < square_r_c.get_property(atoms_[i].id, atoms_[j].id))
			{
				double bv = exp(
					(r_0.get_property(atoms_[i].id, atoms_[j].id) - sqrt(dist)) /
					b.get_property(atoms_[i].id, atoms_[j].id)
				);
				BVS_result_for_all_atoms[i] += bv;
				BVS_result_for_all_atoms[j] += bv;
				n++;
				if (n % 10000 == 0) cout << "*";
			}
		}
		block = int(BVS_result_for_all_atoms[i] / precision);
		BVS_hyst[block]++;
	}
}
	cout << "\n";

}

double AtomicStructure::GII(vector<double> &BVS_result_for_all_atoms, double precision)
{
	double GII_ = 0;
	int size_ = atoms_.size();
	unsigned long i = 0;
	for (unsigned long i = 0; i < size_; ++i)
		{
				GII_ += Square(BVS_result_for_all_atoms[i] 
					- abs((*this).atomic_type_table_->get_element(atoms_[i].id).id_char[2]));
		}
	GII_ /= size_;
	return GII_;
}

double AtomicStructure::GII(vector<double> &BVS_result_for_all_atoms, double precision,const UnaryPropertyTable<int> &valence) //��������� �������� ����������� v
{
	double GII_ = 0;
	int size_ = atoms_.size();
	unsigned long i = 0;
	for (unsigned long i = 0; i < size_; ++i)
	{
		GII_ += Square(BVS_result_for_all_atoms[i]
			- abs(valence.get_property(atoms_[i].id)));
	}
	GII_ /= size_;
	return GII_;
}

AtomicStructure::~AtomicStructure() 
{
	atoms_.clear();
}

void AtomicStructure::BVS_time_test(const BinaryPropertyTable<double>& r_c, const BinaryPropertyTable<double>& r_0,
	const BinaryPropertyTable<double>& b, vector<int> &BVS_hyst,
	vector<double> &BVS_result_for_all_atoms, double precision, int schedule_type, int streams_num, int schedule_sec_par) {
	cout << "!!!BVS Analysis!!!" << endl;
	cout << "Calculating Process..." << endl;
	int size_ = atoms_.size();

	BinaryPropertyTable<double> square_r_c([](double x) {
		return x*x;
	}, r_c);

	BVS_hyst.resize(int(20 / precision) + 1);
	BVS_result_for_all_atoms.resize(size_);
	for (int i = 0; i < size_; i++)
	{
		BVS_result_for_all_atoms[i] = 0;
	}
	for (int i = 0; i < BVS_hyst.size(); i++)
	{
		BVS_hyst[i] = 0;
	}



	int block;
	int n = 0;

	omp_set_num_threads(16);
#pragma omp parallel 
	{
#pragma omp for schedule(dynamic, 100)
		for (int i = 0; i < size_; i++)
		{
			for (int j = i; j < size_; j++)
			{
				double dist = Square(atoms_[i].X - atoms_[j].X)
					+ Square(atoms_[i].Y - atoms_[j].Y)
					+ Square(atoms_[i].Z - atoms_[j].Z);
				if (dist < square_r_c.get_property(atoms_[i].id, atoms_[j].id))
				{
					double bv = exp(
						(r_0.get_property(atoms_[i].id, atoms_[j].id) - sqrt(dist)) /
						b.get_property(atoms_[i].id, atoms_[j].id)
					);
					BVS_result_for_all_atoms[i] += bv;
					BVS_result_for_all_atoms[j] += bv;
					n++;
					if (n % 10000 == 0) cout << "*";
				}
			}
			block = int(BVS_result_for_all_atoms[i] / precision);
			BVS_hyst[block]++;
		}
	}
	cout << "\n";

}

void AtomicStructure::transform(Matrix3d &transform_matrix)
{
	for (unsigned i = 0; i < atoms_.size(); ++i)
	{
		Vector3d temp = {atoms_[i].X, atoms_[i].Y, atoms_[i].Z};
		temp = transform_matrix * temp;
		atoms_[i] = {atoms_[i].id, temp(0), temp(1), temp(2)};
	}
}

void AtomicStructure::cut_by_plane(double coeff[4])
{
	for (unsigned i = 0; i < size(); ++i)
	{
		if ((atoms_[i].X * coeff[0] + atoms_[i].Y * coeff[1] + atoms_[i].Z * coeff[2] + coeff[3]) >= 0)
		{
			atoms_.erase(atoms_.begin() + i);
		}
	}
}

void AtomicStructure::sort_by_sections(double length)
{

}