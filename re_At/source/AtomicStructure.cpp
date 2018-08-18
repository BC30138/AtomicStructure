#include "AtomicStructure.h"



void AtomicStructure::read(istream &file_in, double section_width) {
	bool N_atoms_error = false;
	string N_atoms_str;

	getline(file_in, N_atoms_str);

	if (!N_atoms_str.empty())
	{
		unsigned int char_it = 0;
		while ((char_it < (N_atoms_str.size() - 1)) && !N_atoms_error)
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
		string comment;
		getline(file_in, comment);

		if (!(comment.empty()))
			{
				stringstream Translation_in(comment);
				string temp;
				Translation_in >> temp;
				if (temp == "translations")
				{
					for (int tr_count = 0; tr_count < 3; ++tr_count)
					{
						for (int arg_count = 0; arg_count < 3; ++arg_count)
						{
							Translation_in >> Translations_[tr_count][arg_count];
						}
					}
				}
			}

		int N_atoms = stoi(N_atoms_str);
		atom t;
		for (int i = 0; i < N_atoms; i++)
		{
			string Symbol;
			file_in >> Symbol;
			t.id = atomic_type_table_->get_Id_with_adding(Symbol);
			file_in >> t.X;
			file_in >> t.Y;
			file_in >> t.Z;
			getline(file_in, Symbol);
			if (file_in.bad())
			{
				cout << "ERROR in file: incorrect number of atoms." << endl;
				break;
			}
			atoms_.push_back(t);
		}

		if (!ONE_METHOD_FLAG) make_sections(section_width);
	}
	else
	{
		cout << "ERROR in file: first line for number of atoms." << '\n';
	}
};

void AtomicStructure::read(istream &file_in) {
	bool N_atoms_error = false;
	string N_atoms_str;
	
	//file_in >> N_atoms_str;

	std::getline(file_in, N_atoms_str);

	if (!N_atoms_str.empty())
	{
		unsigned int char_it = 0;
		while ((char_it < (N_atoms_str.size() - 1)) && !N_atoms_error)
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
		string comment;
		getline(file_in, comment);

		if (!(comment.empty()))
			{
				stringstream Translation_in(comment);
				string temp;
				Translation_in >> temp;
				if (temp == "translations")
				{
					for (int tr_count = 0; tr_count < 3; ++tr_count)
					{
						for (int arg_count = 0; arg_count < 3; ++arg_count)
						{
							Translation_in >> Translations_[tr_count][arg_count];
						}
					}
				}
			}

		int N_atoms = stoi(N_atoms_str);
		atom t;
		for (int i = 0; i < N_atoms; i++)
		{
			string Symbol;
			file_in >> Symbol;
			t.id = atomic_type_table_->get_Id_with_adding(Symbol);
			file_in >> t.X;
			file_in >> t.Y;
			file_in >> t.Z;
			getline(file_in, Symbol);
			atoms_.push_back(t);
		}

		if (!ONE_METHOD_FLAG) make_sections(section_width_);
	}
	else
	{
		cout << "ERROR in file: first line for number of atoms." << '\n';
	}
}

void AtomicStructure::make_sections(double width)
{
	section_width_ = width; 
	find_x_limits();
	find_y_limits();
	find_z_limits();

	columns_number_ = (int)((x_limits_.second - x_limits_.first) / width) + 1;
	rows_number_ = (int)((y_limits_.second - y_limits_.first) / width) + 1;
	matrix_size_ = columns_number_ * rows_number_;
	matrix_number_ = (int)((z_limits_.second - z_limits_.first) / width) + 1;
	

	vector<vector<size_t>> block_sections;
	block_sections.resize(matrix_size_  * matrix_number_);

	//quickSort_parallel(atoms_, size(),8);

	omp_set_num_threads(8);
	#pragma omp parallel
 	{
	 #pragma omp for schedule(static)
		for (size_t atoms_it = 0; atoms_it < size(); ++atoms_it)
		{
			block_sections[floor((atoms_[atoms_it].X - x_limits_.first) / width)
			+ columns_number_ * floor((atoms_[atoms_it].Y - y_limits_.first) / width)
			+ matrix_size_ * floor((atoms_[atoms_it].Z - z_limits_.first) / width)
			].push_back(atoms_it);
		}
	}
	vector<atom> atoms_temp;
	sections_.push_back(0); 
	for (size_t block_it = 0; block_it < block_sections.size(); ++block_it)
	{
		for (size_t atoms_it = 0; atoms_it < block_sections[block_it].size(); ++atoms_it)
		{
			atoms_temp.push_back(atoms_[block_sections[block_it][atoms_it]]);
		}
		sections_.push_back(sections_[block_it] + block_sections[block_it].size());
	}

	atoms_.clear();
	for (size_t it = 0; it < atoms_temp.size(); ++it)
	{
		atoms_.push_back(atoms_temp[it]);
	}
}

void AtomicStructure::dist_analysis_old(vector<int> type1id, vector<int> type2id, double right,
	double deltaDist, vector<int> &dhyst, long long int &comp_it) {
	cout << "!!!Distance Analysis!!!" << endl;
	//ofstream plot("plotAnalysisDist.dat");
	//plot.precision(5)

	vector<int> type1atoms;
	vector<int> type2atoms;


	cout << "Search of type 1 atoms..." << endl;
	for (int i = 0; i < (int)type1id.size(); i++)
	{
		for (int n = 0; n < (int)atoms_.size(); n++)
		{
			if (atoms_[n].id == type1id[i]) type1atoms.push_back(n);
		}
	}

	cout << "Search of type 2 atoms..." << endl;
	for (int i = 0; i < (int)type2id.size(); i++)
	{
		for (int n = 0; n < (int)atoms_.size(); n++)
		{
			if (atoms_[n].id == type2id[i]) type2atoms.push_back(n);
		}
	}

	cout << "Calculating Process..." << endl;


	double squareright;
	squareright = Square(right);

	dhyst.resize(int((right) / deltaDist) + 1);
	for (int i = 0; i < (int)dhyst.size(); i++)
	{
		dhyst[i] = 0;
	}

	int block;

	comp_it = 0;

	omp_set_num_threads(4);
	#pragma omp parallel
	{
	#pragma omp for schedule(static, 1000)
		for (int n1 = 0; n1 < (int)type1atoms.size(); n1++)
		{
			for (int n2 = 0; n2 < (int)type2atoms.size(); n2++)
			{
				double sqrdist = Square(atoms_[type1atoms[n1]].X - atoms_[type2atoms[n2]].X)
					+ Square(atoms_[type1atoms[n1]].Y - atoms_[type2atoms[n2]].Y)
					+ Square(atoms_[type1atoms[n1]].Z - atoms_[type2atoms[n2]].Z);
					++comp_it;
				if (sqrdist <= squareright)
				{
					block = (int)floor((sqrt(sqrdist)) / deltaDist);
					dhyst[block]++;
				}
			}
		}
	}
}

void AtomicStructure::dist_analysis(vector<int> type1id, vector<int> type2id, double neighbour,
	double deltaDist, vector<int> &dhyst, long long int &comp_it) {
	cout << "!!!Distance Analysis!!!" << endl;
	//ofstream plot("plotAnalysisDist.dat");
	//plot.precision(5)
	cout << "Calculating Process..." << endl;
	size_t sec_around = (size_t)ceil(neighbour / section_width_);
	dhyst.resize((size_t)ceil((sqrt(
		  Square(x_limits_.first - x_limits_.second) 
		+ Square(y_limits_.first - y_limits_.second)
		+ Square(z_limits_.first - z_limits_.second)
		)
		)/ deltaDist));
	for (size_t i = 0; i < dhyst.size(); i++)
	{
		dhyst[i] = 0;
	}
	comp_it = 0;
	//int block_per_thread = (int)(sections_.size()/8); 
	omp_set_num_threads(8);
	#pragma omp parallel
	{
	#pragma omp for schedule(static)
	//cout << columns_number_ << "\t" << rows_number_ << "\t" << matrix_number_ << "\t" << sec_around << endl;
		for (size_t sec_count = 0; sec_count < sections_.size() - 1; ++sec_count)
		{
			size_t to_left_z = (int)(sec_count / matrix_size_) - (int)(max(0, (int)(sec_count / matrix_size_) - (int)sec_around));
			size_t to_right_z = min((int)matrix_number_, (int)(sec_count / matrix_size_) + (int)sec_around) - 
								(int)(sec_count / matrix_size_);
			size_t to_left_y = (int)(sec_count / columns_number_) - max((int)((int)(sec_count / matrix_size_) * rows_number_), 
										 (int)(sec_count / columns_number_) - (int)sec_around);
			size_t to_right_y = min((int)((int)(sec_count / matrix_size_) * rows_number_ + rows_number_), 
										 (int)(sec_count / columns_number_) + (int)sec_around)
										  - (int)(sec_count / columns_number_);
			size_t to_left_x = sec_count - max((int)((int)(sec_count / columns_number_) * columns_number_), (int)(sec_count - sec_around));
			size_t to_right_x = min((int)((sec_count / columns_number_) * columns_number_ + columns_number_) , (int)(sec_count + sec_around)) 
						- sec_count;
			size_t start = sec_count - to_left_z * matrix_size_ - to_left_y * columns_number_ - to_left_x;
			//cout << to_left_x << "\t" << to_right_x << "\t" << to_left_y << "\t" << to_right_y << "\t" << to_left_z << "\t" << to_right_z << endl;  
			for (size_t n1 = sections_[sec_count]; n1 < sections_[sec_count + 1]; ++n1)
			{
				if (find(type1id.begin(), type1id.end(), atoms_[n1].id) != type1id.end())
				{
					for (size_t z_it = 0; z_it < to_left_z + to_right_z; ++z_it)
					{
						for (size_t y_it = 0; y_it < to_left_y + to_right_y; ++y_it)
						{
							for(size_t n2 = sections_[start + z_it * matrix_size_ + y_it * columns_number_];
									   n2 < sections_[start + z_it * matrix_size_ + y_it * columns_number_ + to_left_x + to_right_x]; ++n2) 
							{
								
								if (find(type2id.begin(), type2id.end(), atoms_[n2].id) != type2id.end())
								{
									++comp_it;
									double sqrdist = Square(atoms_[n1].X - atoms_[n2].X)
													+ Square(atoms_[n1].Y - atoms_[n2].Y)
													+ Square(atoms_[n1].Z - atoms_[n2].Z);
									if (sqrdist < Square(neighbour))
									{
										dhyst[(int)floor((sqrt(sqrdist)) / deltaDist)]++;
									}
								}
							}
						}	
					}
				}
			}
		}
	}
}


void AtomicStructure::replicate(double* Translation0, double* Translation1, 
	double* Translation2, int* size, double section_width)
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
					atom t;
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
					atom t;
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
					atom t;
					t.id = atoms_[j].id;
					t.X = atoms_[j].X + i*translationX;
					t.Y = atoms_[j].Y + i*translationY;
					t.Z = atoms_[j].Z + i*translationZ;
					atoms_.push_back(t);
				}
			}
		}
	}

	make_sections(section_width);
}

void AtomicStructure::create(istream &file_in , int* size)
{
	ONE_METHOD_FLAG = true;
	read(file_in);
	replicate(Translations_[0], Translations_[1], Translations_[2], size, section_width_);
	ONE_METHOD_FLAG = false;
}

void AtomicStructure::create(istream &file_in , int* size, double section_width)
{
	ONE_METHOD_FLAG = true;
	read(file_in,section_width);
	replicate(Translations_[0], Translations_[1], Translations_[2], size, section_width);
	ONE_METHOD_FLAG = false;
}

void AtomicStructure::info()
{
		cout << "Atoms number\t=\t" << size() << '\n';
		cout << "Sections number\t=\t" << sections_.size() << '\n';
		cout << "Section width\t=\t" << section_width_ << '\n'; 
		cout << "Types number\t=\t" << atomic_type_table_->size() << '\n';
		cout << "X limits\t=\t" << x_limits_.first << "\t:\t" << x_limits_.second << '\n';
		cout << "Y limits\t=\t" << y_limits_.first << "\t:\t" << y_limits_.second << '\n';
		cout << "Z limits\t=\t" << z_limits_.first << "\t:\t" << z_limits_.second << '\n';
}

void AtomicStructure::show()
{
	for (size_t i = 0; i < size(); ++i)
		{
		cout << i << '\t' << atoms_[i].X << '\t' << atoms_[i].Y << '\t' <<
		 atoms_[i].Z << endl;
		}
	cout << "\nsections number:\t" << sections_.size() << endl;
	for (size_t it = 0; it < sections_.size(); ++it)
	{
		cout << sections_[it] << endl;
	}
}

void AtomicStructure::to_file(string output_file_name)
{
	cout << "Export to file..." << endl;
	ofstream out("XYZ/" + output_file_name);
	out << /*"#" <<*/ atoms_.size() << "\n\n";

	for (int n = 0; n < (int)atoms_.size(); n++)
	{
		out << longid_to_string(atomic_type_table_->get_element(atoms_[n].id).long_id) << "\t" 
		<< atoms_[n].X << "\t" << atoms_[n].Y << "\t" << atoms_[n].Z << "\n";
	}
	cout << "\n";
	out << "\n\n";
	out.close();
}

void AtomicStructure::to_file_types(string output_file_name)
{
	cout << "Export to file..." << endl;

	vector<vector<size_t>> types;
	types.resize(atomic_type_table_->size());
	for (size_t n = 0; n < atoms_.size(); ++n)
	{
		types[atoms_[n].id].push_back(n);		
	}

	for (size_t it = 0; it < types.size(); ++it)
	{
		string symb = longid_to_string(atomic_type_table_->get_element((int)it).long_id);
		ofstream out("XYZ/" + output_file_name + symb[0] + symb[1] + ".xyz");
		out << types[it].size() << "\n\n";
		for (size_t jt = 0; jt < types[it].size(); ++jt)
		{
			out << symb << "\t" << atoms_[types[it][jt]].X << "\t" << atoms_[types[it][jt]].Y << "\t" << atoms_[types[it][jt]].Z << "\n";
		}
		out.close();
	}
	cout << "\n";
}

void AtomicStructure::shift(double* vec) {
	cout << "Shifting..." << endl;
	double shiftX = vec[0];
	double shiftY = vec[1];
	double shiftZ = vec[2];
	for (size_t i = 0; i < atoms_.size(); i++)
	{
		atoms_[i].X += shiftX;
		atoms_[i].Y += shiftY;
		atoms_[i].Z += shiftZ;
		//if (i % 524288 == 0) cout << "*";
	}
	cout << "\n";
	find_x_limits();
	find_y_limits();
	find_z_limits();
}

void AtomicStructure::plus_conflicts_old(AtomicStructure &right_hand, BinaryPropertyTable<double> &property_hard,
	BinaryPropertyTable<double> &property_soft, double right_bound)
{
	BinaryPropertyTable<double> square_hard([](double x) {
		return x*x;}, property_hard);
	BinaryPropertyTable<double> square_soft([](double x) {
		return x*x;}, property_soft);
	double square_bound = Square(right_bound);
	
	int warning_counter_[atomic_type_table_->size()][atomic_type_table_->size()];
	for (size_t it = 0; it < atomic_type_table_->size(); ++it)
		for (size_t jt = 0; jt < atomic_type_table_->size(); ++jt) warning_counter_[it][jt] = 0;
	
	size_t left_at_size = size(); 
	for(size_t right_at = 0; right_at < right_hand.size(); ++right_at)
	{
		bool hard_conflict_flag = false;
		int warning_counter[atomic_type_table_->size()];
		for (size_t it = 0; it < left_at_size; ++it) warning_counter[it] = 0;
		for (size_t left_at = 0; (left_at < size()) && !hard_conflict_flag; ++left_at)
		{
			double square_dist = Square(atoms_[left_at].X - right_hand.atoms_[right_at].X) + Square(atoms_[left_at].Y - right_hand.atoms_[right_at].Y)
				+ Square(atoms_[left_at].Z - right_hand.atoms_[right_at].Z); 
			cout << longid_to_string(atomic_type_table_->get_element(atoms_[left_at].id).long_id) << "\t" << 
				longid_to_string(atomic_type_table_->get_element(right_hand.atoms_[right_at].id).long_id) << "\t" << square_dist 
					//<< "\t"	<< square_hard.get_property(atoms_[left_at].id, right_hand.atoms_[right_at].id) 
					<< endl;
			if(square_dist > square_bound)
			{	
				if (square_dist < square_hard.get_property(atoms_[left_at].id, right_hand.atoms_[right_at].id))
					{
						hard_conflict_flag = true;
						left_at = size();
					}
				else
					{
						if(square_dist < square_soft.get_property(atoms_[left_at].id, right_hand.atoms_[right_at].id)) 
							++warning_counter[atoms_[right_at].id];
					} 
			}
		}
		if (!hard_conflict_flag)
		{
			atoms_.push_back(right_hand.atoms_[right_at]);
			for (size_t it = 0; it < atomic_type_table_->size(); ++it)
			{
				if (it < right_hand.atoms_[right_at].id)
				warning_counter_[right_hand.atoms_[right_at].id][it] += warning_counter[it];
				else
					warning_counter_[it][right_hand.atoms_[right_at].id] += warning_counter[it];
			}
		}
	}
	cout << "\n\n";
	for (size_t it = 0; it < atomic_type_table_->size(); ++it)
		for (size_t jt = it; jt < atomic_type_table_->size(); ++jt)
		{
			if (warning_counter_[it][jt] != 0)
			{
				cout << longid_to_string(atomic_type_table_->get_element(it).long_id) << "\t" << 
				longid_to_string(atomic_type_table_->get_element(jt).long_id) << "\t" << warning_counter_[it][jt] << endl;
			}
		}
}