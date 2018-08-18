#ifndef _CONSOLEINTERFACE_HPP_
#define _CONSOLEINTERFACE_HPP_

#include <sstream>
#include "AtomicStructure.h"

class console
{
enum cmds_id_ {CMD_HELP, CMD_EXIT, CMD_CREATE,CMD_READ,CMD_REPLICATE,CMD_SAVE,CMD_ROTATEBY,
				CMD_ROTATEBY_EULER, CMD_ROTATEBY_AIR,CMD_ROTATEBY_QUAT, CMD_MIN_MAX, CMD_SHIFT, CMD_SHIFTMIN,
				CMD_UPROP, CMD_BINPROP, CMD_DELETE, CMD_BVS, CMD_COORD, CMD_DIST, CMD_PLOT, CMD_SCRIPT, 
				CMD_EQUAL, CMD_PRINT, CMD_ERROR, CMD_SIZE, CMD_INT, CMD_DOUBLE, CMD_VECTOR, CMD_ATOMICSTRUCT,
				CMD_ATOM, CMD_ATOM_TR, CMD_2DARRAY};
enum type_id_ {TPI_ATOMICSTRUCT, TPI_VECTOR, TPI_UPROP_INT, TPI_UPROP_DOUBLE,TPI_BINPROP_INT, TPI_BINPROP_DOUBLE,
				TPI_ARRAY_OF_VECTORS};

	std::unordered_map<cmds_id_, string> infos_;
	std::unordered_map<string, cmds_id_> cmds_;
	std::unordered_map<string, type_id_> types_;
	std::unordered_map<string, vector<double*>> arrays2d_;
	std::unordered_map<string, double*> vectors_;
	std::unordered_map<string, AtomicStructure*> AtomicStructs_;
	std::unordered_map<string, UnaryPropertyTable<int>*> U_props_int_;
	std::unordered_map<string, UnaryPropertyTable<double>*> U_props_double_;
	std::unordered_map<string, BinaryPropertyTable<int>*> Bin_props_int_;
	std::unordered_map<string, BinaryPropertyTable<double>*> Bin_props_double_;

	AtomicTypeTable *atomic_type_table_;
public:
	console(AtomicTypeTable &atomic_type_table)
	{
		atomic_type_table_ = &atomic_type_table;

		cmds_["help"] = CMD_HELP;
		cmds_["exit"] = CMD_EXIT;
		cmds_["delete"] = CMD_DELETE;
		cmds_["script"] = CMD_SCRIPT;
		cmds_["create"] = CMD_CREATE;
		cmds_["print"] = CMD_PRINT;
		cmds_["="] = CMD_EQUAL;
		cmds_["replicate"] = CMD_REPLICATE;
		cmds_["save"] = CMD_SAVE;
		cmds_["bvs"] = CMD_BVS;
		cmds_["coord"] = CMD_COORD;
		cmds_["dist"] = CMD_DIST;
		cmds_["plot"] = CMD_PLOT;
		cmds_["rotateby"] = CMD_ROTATEBY;
		cmds_["read"] = CMD_READ;
		cmds_["shift"] = CMD_SHIFT;

		cmds_["tr_struct"] = CMD_ATOM_TR;
		cmds_["struct"] = CMD_ATOM;
		cmds_["vectors"] = CMD_2DARRAY;
		cmds_["vector"] = CMD_VECTOR;
		cmds_["atomicstruct"] = CMD_ATOMICSTRUCT;
		cmds_["uprop"] = CMD_UPROP;
		cmds_["binprop"] = CMD_BINPROP;
		cmds_["size"] = CMD_SIZE;
			cmds_["euler"] = CMD_ROTATEBY_EULER;
			cmds_["air"] = CMD_ROTATEBY_AIR;
			cmds_["quat"] = CMD_ROTATEBY_QUAT;
		cmds_["minmax"] = CMD_MIN_MAX;
		cmds_["shftmin"] = CMD_SHIFTMIN;
			cmds_["int"] = CMD_INT;
			cmds_["double"] = CMD_DOUBLE;
		
		ifstream info_in("data\\command_info.dat");
		string word;
		cmds_id_ _cmd;
		string _line;
		
		while (std::getline(info_in, _line))
		{
			stringstream _line_in(_line);
			_line_in >> word;
			auto word_it = cmds_.find(word);
			if (word_it == cmds_.end())
			{
				infos_[_cmd] += _line + '\n';
			}
			else
			{
				_cmd = word_it->second;
			}
			_line_in.clear();
		}

		info_in.close();
	}
	
private:


	cmds_id_ cmd_process_(istream &in, int program_counter)
	{
		string first_console_word;
		in >> first_console_word;
		auto cmd_ = cmds_.find(first_console_word);
		if (cmd_ != cmds_.end())
		{
			switch (cmd_->second)
			{
			case (CMD_HELP):
			{
				string parameter;
				in >> parameter;
				if ((parameter.empty()) || (parameter[0] == '#'))
				{
					cout << "type \"help <command_name>\" for more info." << '\n';
					cout << "command list:" << '\n';
					cout << "read" << '\t' << "create" << '\t' << "print" << '\t' << "delete" << '\n';
					cout << "exit" << '\t' << "replicate" << '\t' << "bvs" << '\n';
				}
				else
				{	
					string comment;
					in >> comment;
					if ((comment.empty()) || (comment[0] == '#'))
					{
						auto parameter_it = cmds_.find(parameter);
						auto info_it = infos_.find(parameter_it->second);
						if (info_it != infos_.end())
						{
							cout << infos_[parameter_it->second];
						}
						else
						{
							cout << "ERROR in line " << program_counter <<
								": incorrect command '" << parameter << "'.\n";
							return CMD_ERROR;
						}
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							": invalid value '" << comment << "'.\n";
						return CMD_ERROR;
					}
				}
				break;
			}

			case (CMD_READ):
			{
				string val_type;
				in >> val_type;
				auto val_type_it = cmds_.find(val_type);
				switch (val_type_it->second)
				{
				case (CMD_ATOM):
				{
					string var_name;
					string input_file_path;
					in >> var_name;

					char first_symb;
					char symb = 0;
					in.get(first_symb);
					in.get(first_symb);
					if ((first_symb == 34) || (first_symb == 39))
					{
						while (symb != first_symb) //eof ����� �������� ��� ����� ��� ���� ��� ������
						{
							if (in.eof())
							{
								cout << "ERROR in line " << program_counter <<
									": requires a closing quote after the filename." << '\n';
								return CMD_ERROR;
							}
							else
							{
								in.get(symb);
								if ((symb != first_symb))
								{
									input_file_path.push_back(symb);
								}
							}
						}
						string comment;
						in >> comment;
						if (!(comment[0] == '#') && !(comment.empty()))
						{
							cout << "ERROR in line " << program_counter <<
								": command 'read' must have three parameters." << '\n';
							return CMD_ERROR;
						}
						else
						{
							auto struct_it = AtomicStructs_.find(var_name);
							if (struct_it == AtomicStructs_.end())
							{
								if (isalpha(var_name[0]))
								{
									auto it = types_.find(var_name);
									if (it == types_.end())
									{
										AtomicStructure *structure = new AtomicStructure(*atomic_type_table_);
										AtomicStructs_[var_name] = structure;
										types_[var_name] = TPI_ATOMICSTRUCT;
										cout << "structure '" << var_name << "' successfully created." << '\n';
										struct_it = AtomicStructs_.find(var_name);
									}
									else
									{
										cout << "ERROR in line " << program_counter <<
											": variable's name already used for other type variable." << '\n';
										return CMD_ERROR;
									}
								}
								else
								{
									cout << "ERROR in line " << program_counter <<
										": variable's name must begin with letter." << '\n';
									return CMD_ERROR;
								}
							}
							struct_it->second->read(input_file_path);
							cout << "reading completed successfully." << '\n';
						}
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							":required filename in quotes." << '\n';
						return CMD_ERROR;
					}
					break;
				}

				case (CMD_2DARRAY):
				{
					string var_name;
					string input_file_path;
					in >> var_name;

					char first_symb;
					char symb = 0;
					in.get(first_symb);
					in.get(first_symb);
					if ((first_symb == 34) || (first_symb == 39))
					{
						while (symb != first_symb) //eof ����� �������� ��� ����� ��� ���� ��� ������
						{
							if (in.eof())
							{
								cout << "ERROR in line " << program_counter <<
									": requires a closing quote after the filename." << '\n';
								return CMD_ERROR;
							}
							else
							{
								in.get(symb);
								if ((symb != first_symb))
								{
									input_file_path.push_back(symb);
								}
							}
						}
						string comment;
						in >> comment;
						if (!(comment[0] == '#') && !(comment.empty()))
						{
							cout << "ERROR in line " << program_counter <<
								": command 'read' must have three parameters." << '\n';
							return CMD_ERROR;
						}
						else
						{
							ifstream ifs("XYZ\\" + input_file_path);
							bool N_vectors_error = false;
							string temp_str;
							getline(ifs, temp_str);
							if (!temp_str.empty())
							{
								unsigned char_it = 0;
								while (!N_vectors_error && (char_it < (temp_str.size())))
								{
									if (!isdigit(temp_str[char_it]))
										N_vectors_error = true;
									else ++char_it;
								}
							}
							else
							{
								N_vectors_error = true;;
							}

							if (!N_vectors_error)
							{
								int N_vectors = atof(temp_str.c_str());
								auto vectors_it = arrays2d_.find(var_name);
								if (vectors_it == arrays2d_.end())
								{
									if (isalpha(var_name[0]))
									{
										auto it = types_.find(var_name);
										if (it == types_.end())
										{

											vector<double*> _vecs;
											_vecs.resize(N_vectors);
											arrays2d_[var_name] = _vecs;
											types_[var_name] = TPI_ARRAY_OF_VECTORS;
											cout << "array '" << var_name << "' successfully created." << '\n';
											vectors_it = arrays2d_.find(var_name);
										}
										else
										{
											cout << "ERROR in line " << program_counter <<
												": variable's name already used for other type variable." << '\n';
											return CMD_ERROR;
										}
									}
									else
									{
										cout << "ERROR in line " << program_counter <<
											": variable's name must begin with letter." << '\n';
										return CMD_ERROR;
									}
								}
								else
								{
									vectors_it->second.resize(N_vectors);
								}

								getline(ifs, temp_str);

								for (unsigned i = 0; i < N_vectors; ++i)
								{
									getline(ifs, temp_str);
									stringstream parameters_stream(temp_str);
									if (isalpha(temp_str[0]))
									{
										string temp_word;
										parameters_stream >> temp_word;
									}
									vectors_it->second[i] = new double[3];
									parameters_stream >> vectors_it->second[i][0];
									parameters_stream >> vectors_it->second[i][1];
									parameters_stream >> vectors_it->second[i][2];
								}
								cout << "reading completed successfully." << '\n';
							}
							else
							{
								cout << "ERROR in file '" << input_file_path << "': first line for number of atoms\\vectors." << '\n';
							}
						}
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							":required filename in quotes." << '\n';
						return CMD_ERROR;
					}
					break;
				}

				case (CMD_ATOM_TR):
				{
					string var_name;
					string input_file_path;
					in >> var_name;

					char first_symb;
					char symb = 0;
					in.get(first_symb);
					in.get(first_symb);
					if ((first_symb == 34) || (first_symb == 39))
					{
						while (symb != first_symb) //eof ����� �������� ��� ����� ��� ���� ��� ������
						{
							if (in.eof())
							{
								cout << "ERROR in line " << program_counter <<
									": requires a closing quote after the filename." << '\n';
								return CMD_ERROR;
							}
							else
							{
								in.get(symb);
								if ((symb != first_symb))
								{
									input_file_path.push_back(symb);
								}
							}
						}
						string comment;
						in >> comment;
						if (!(comment[0] == '#') && !(comment.empty()))
						{
							cout << "ERROR in line " << program_counter <<
								": command 'read' must have three parameters." << '\n';
							return CMD_ERROR;
						}
						else
						{
							bool var_name_error = false;
							auto struct_it = AtomicStructs_.find(var_name);
							if (struct_it == AtomicStructs_.end())
							{
								if (isalpha(var_name[0]))
								{
									auto it = types_.find(var_name);
									if (it == types_.end())
									{
										AtomicStructure *structure = new AtomicStructure(*atomic_type_table_);
										AtomicStructs_[var_name] = structure;
										types_[var_name] = TPI_ATOMICSTRUCT;
										cout << "structure '" << var_name << "' successfully created." << '\n';
										struct_it = AtomicStructs_.find(var_name);
									}
									else
									{
										cout << "ERROR in line " << program_counter <<
											": variable's name already used for other type variable." << '\n';
										var_name_error = true;
										return CMD_ERROR;
									}
								}
								else
								{
									cout << "ERROR in line " << program_counter <<
										": variable's name must begin with letter." << '\n';
									var_name_error = true;
									return CMD_ERROR;
								}
							}

							if (!var_name_error)
							{
								ifstream ifs("XYZ\\" + input_file_path);
								bool N_vectors_error = false;
								string temp_str;
								getline(ifs, temp_str);
								if (!temp_str.empty())
								{
									unsigned char_it = 0;
									while (!N_vectors_error && (char_it < (temp_str.size())))
									{
										if (!isdigit(temp_str[char_it]))
											N_vectors_error = true;
										else ++char_it;
									}
								}
								else
								{
									N_vectors_error = true;;
								}

								if (!N_vectors_error)
								{
									int N_vectors = atof(temp_str.c_str());
									auto vectors_it = arrays2d_.find(var_name + "_tr");
									if (vectors_it == arrays2d_.end())
									{
										if (isalpha(var_name[0]))
										{
											auto it = types_.find(var_name + "_tr");
											if (it == types_.end())
											{

												vector<double*> _vecs;
												_vecs.resize(3);
												arrays2d_[var_name + "_tr"] = _vecs;
												types_[var_name + "_tr"] = TPI_ARRAY_OF_VECTORS;
												cout << "array '" << var_name + "_tr" << "' successfully created." << '\n';
												vectors_it = arrays2d_.find(var_name + "_tr");
											}
											else
											{
												cout << "ERROR in line " << program_counter <<
													": variable's name already used for other type variable." << '\n';
												return CMD_ERROR;
											}
										}
										else
										{
											cout << "ERROR in line " << program_counter <<
												": variable's name must begin with letter." << '\n';
											return CMD_ERROR;
										}
									}
									else
									{
										vectors_it->second.resize(3);
									}

									getline(ifs, temp_str);

									for (unsigned i = 0; i < 3; ++i)
									{
										getline(ifs, temp_str);
										stringstream parameters_stream(temp_str);
										if (isalpha(temp_str[0]))
										{
											string temp_word;
											parameters_stream >> temp_word;
										}
										vectors_it->second[i] = new double[3];
										parameters_stream >> vectors_it->second[i][0];
										parameters_stream >> vectors_it->second[i][1];
										parameters_stream >> vectors_it->second[i][2];
									}

									Atom t;
									string Symbol;
									for (unsigned i = 0; i < (N_vectors - 3); ++i)
									{
										ifs >> Symbol;
										t.id = atomic_type_table_->get_Id_with_adding(Symbol);
										ifs >> t.X;
										ifs >> t.Y;
										ifs >> t.Z;
										struct_it->second->atoms_.push_back(t);
									}
									struct_it->second->find_x_limits();
									struct_it->second->find_y_limits();
									struct_it->second->find_z_limits();

									cout << "reading completed successfully." << '\n';
								}
								else
								{
									cout << "ERROR in file '" << input_file_path << "': first line for number of atoms\\vectors." << '\n';
								}
							}
						}
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							":required filename in quotes." << '\n';
						return CMD_ERROR;
					}
					break;
				}

				default:
				{
					cout << "ERROR in line " << program_counter <<
						": incorrect value type '" << val_type << "'.\n";
					return CMD_ERROR;
				}
				}
				break;
			}

			case (CMD_REPLICATE):
			{
				string structure_name;
				in >> structure_name;

				auto it = AtomicStructs_.find(structure_name);
				if (it != AtomicStructs_.end())
				{
					bool parameter_input_error = false;
					string vector_param_name;
					in >> vector_param_name;
					auto param_it = types_.find(vector_param_name);
					vector<double*> parameter_vectors;
					parameter_vectors.resize(3);
					switch (param_it->second)
					{
					case (TPI_VECTOR):
					{
						vector<string> vector_parameters_names;
						vector_parameters_names.resize(3);
						vector_parameters_names[0] = vector_param_name;
						for (unsigned parameters_it = 1; parameters_it < 3; ++parameters_it)
						{
							
							in >> vector_parameters_names[parameters_it];
							if (vectors_.find(vector_parameters_names[parameters_it]) == vectors_.end())
							{
								parameter_input_error = true;
								parameters_it = 3;
							}
						}
						if (!parameter_input_error)
						{
							parameter_vectors[0] = vectors_[vector_parameters_names[0]];
							parameter_vectors[1] = vectors_[vector_parameters_names[1]];
							parameter_vectors[2] = vectors_[vector_parameters_names[2]];
						}
						break;
					}
					case (TPI_ARRAY_OF_VECTORS):
					{
						auto array2d_it = arrays2d_.find(vector_param_name);
						if (array2d_it->second.size() == 3)
						{
							parameter_vectors = array2d_it->second;
						}
						else
						{
							parameter_input_error = true;
						}
						break;
					}
					default:
					{
						parameter_input_error = true;
					}
					}
					
					
					if (!parameter_input_error)
					{
						vector<string> size_str;
						size_str.resize(3);
						for (unsigned parameters_it = 0; parameters_it < 3; ++parameters_it)
						{
							in >> size_str[parameters_it];
							if (!size_str[parameters_it].empty())
							{
								unsigned char_it = 0;
								while (!parameter_input_error && (char_it < (size_str[parameters_it].size())))
								{
									if (!isdigit(size_str[parameters_it][char_it]))
									{
										parameter_input_error = true;
										parameters_it = 3;
									}
									else ++char_it;
								}
							}
							else
							{
								while (parameters_it < 3)
								{
									size_str[parameters_it] = "0";
									++parameters_it;
								}
							}
						}

						if (!parameter_input_error)
						{
							string comment;
							in >> comment;
							if (!(comment[0] == '#') && !(comment.empty()))
							{
								cout << "ERROR in line " << program_counter <<
									": command 'replicate' must have five parameters." << '\n';
								return CMD_ERROR;
							}
							else
							{
								int* size = new int[3];
								size[0] = atof(size_str[0].c_str());
								size[1] = atof(size_str[1].c_str());
								size[2] = atof(size_str[2].c_str());
								AtomicStructs_[structure_name]->replicate(parameter_vectors[0],
									parameter_vectors[1], parameter_vectors[2], size);
								cout << "size of '" << structure_name << "':" << AtomicStructs_[structure_name]->size() << '\n';
							}
						}
						else
						{
							cout << "ERROR in line " << program_counter <<
								": you must type replication size in the form of 3 int numbers" << '\n';
							return CMD_ERROR;
						}
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							": you must type array of size 3 or 3 vectors." << '\n';
						return CMD_ERROR;
					}
				}
				else
				{
					cout << "ERROR in line " << program_counter <<
						":structure '" << structure_name << "' not exist." << '\n';
					return CMD_ERROR;
				}
				break;
			}

			case (CMD_EXIT):
			{
				string comment;
				in >> comment;
				if (!(comment[0] == '#') && !(comment.empty()))
				{
					cout << "ERROR in line " << program_counter <<
						": command 'exit' have no parameters." << '\n';
					return CMD_ERROR;
				}
				else return CMD_EXIT;
				break;
			}

			case (CMD_PRINT):
			{
				string var_name;
				in >> var_name;

				
				auto it = types_.find(var_name);
				switch (it->second)
				{
				case (TPI_VECTOR):
				{
					string comment;
					in >> comment;
					if (!(comment[0] == '#') && !(comment.empty()))
					{
						cout << "ERROR in line " << program_counter <<
							": incorrect parameter '" << comment << "'.\n";
						return CMD_ERROR;
					}
					else
					{
						cout << '{' << vectors_[var_name][0] << ','
							<< vectors_[var_name][1] << ',' << vectors_[var_name][2]
							<< '}' << '\n';
					}
					break;
				}

				case (TPI_ATOMICSTRUCT):
				{
					string parameter;
					in >> parameter;
					
					if (parameter.empty() || (parameter[0] == '#'))
					{
						for (unsigned i = 0; (i < 10) && (i < AtomicStructs_[var_name]->size()); ++i)
						{
							cout << longid_to_string(atomic_type_table_->
								get_element(AtomicStructs_[var_name]->atoms_[i].id).long_id) << '\t' <<
								AtomicStructs_[var_name]->atoms_[i].X << '\t' <<
								AtomicStructs_[var_name]->atoms_[i].Y << '\t' <<
								AtomicStructs_[var_name]->atoms_[i].Z << '\n';
						}
					}
					else
					{
						auto it = cmds_.find(parameter);
						if (it != cmds_.end())
						{
							if (it->second == CMD_SIZE)
							{
								string comment;
								in >> comment;
								if (!(comment[0] == '#') && !(comment.empty()))
								{
									cout << "ERROR in line " << program_counter <<
										": incorrect parameter '" << comment << "'.\n";
									return CMD_ERROR;
								}
								else
								{
									cout << "size of '" << var_name << "': " << AtomicStructs_[var_name]->size() << '\n';
								}
							}
							else
							{
								cout << "ERROR in line " << program_counter <<
									": incorrect command '" << parameter << "'.\n";
								return CMD_ERROR;
							}
						}
						else
						{
							bool parameter_error = false;
							unsigned char_it = 0;
							while (!parameter_error && (char_it < (parameter.size())))
							{
								if (!isdigit(parameter[char_it]))
									parameter_error = true;
								else ++char_it;
							}

							if (!parameter_error)
							{
								string second_parameter;
								in >> second_parameter;
								if (!second_parameter.empty() && (second_parameter[0] != '#'))
								{
									unsigned char_it = 0;
									while (!parameter_error && (char_it < (second_parameter.size())))
									{
										if (!isdigit(parameter[char_it]))
											parameter_error = true;
										else ++char_it;
									}
									if (!parameter_error)
									{
										string comment;
										in >> comment;
										if (!comment.empty() && (comment[0] != '#'))
										{
											cout << "ERROR in line " << program_counter <<
												": incorrect parameter '" << comment << "'.\n";
											return CMD_ERROR;
										}
										else
										{
											int left = atoi(parameter.c_str());
											int right = atoi(second_parameter.c_str());
											if ((right < AtomicStructs_[var_name]->size())
												&& (left < right))
											{
												for (int i = left; i < right; ++i)
												{
													cout << longid_to_string(atomic_type_table_->
														get_element(AtomicStructs_[var_name]->atoms_[i].id).long_id) << '\t' <<
														AtomicStructs_[var_name]->atoms_[i].X << '\t' <<
														AtomicStructs_[var_name]->atoms_[i].Y << '\t' <<
														AtomicStructs_[var_name]->atoms_[i].Z << '\n';
												}
											}
											else
											{
												cout << "ERROR in line " << program_counter <<
													": incorrect parameters.\n";
												return CMD_ERROR;
											}
										}
									}
								}
								else
								{
									int i = atoi(parameter.c_str());
									if (i < AtomicStructs_[var_name]->size())
									{
										cout << longid_to_string(atomic_type_table_->
											get_element(AtomicStructs_[var_name]->atoms_[i].id).long_id) << '\t' <<
											AtomicStructs_[var_name]->atoms_[i].X << '\t' <<
											AtomicStructs_[var_name]->atoms_[i].Y << '\t' <<
											AtomicStructs_[var_name]->atoms_[i].Z << '\n';
									}
									else
									{
										cout << "ERROR in line " << program_counter <<
											": incorrect parameter '" << parameter << "'.\n";
										return CMD_ERROR;
									}
								}
							}
							else
							{
								cout << "ERROR in line " << program_counter <<
									": incorrect parameter '" << parameter << "'.\n";
								return CMD_ERROR;
							}
						}
					}
					break;
				}

				case (TPI_BINPROP_DOUBLE):
				{
					string parameter;
					in >> parameter;

					if (parameter.empty() || (parameter[0] == '#'))
					{
						cout << "binary property\t;\t"
							<< "size:" << '\t' << Bin_props_double_[var_name]->size()
							<< "\t;\t type:" << '\t' << "double" << '\n';
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							": incorrect parameter '" << parameter << "'.\n";
						return CMD_ERROR;
					}
					break;
				}

				case (TPI_BINPROP_INT):
				{
					string parameter;
					in >> parameter;

					if (parameter.empty() || (parameter[0] == '#'))
					{
						cout << "binary property\t;\t"
							<< "size:" << '\t' << Bin_props_int_[var_name]->size()
							<< "\t;\t type:" << '\t' << "int" << '\n';
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							": incorrect parameter '" << parameter << "'.\n";
						return CMD_ERROR;
					}
					break;
				}

				case (TPI_UPROP_DOUBLE):
				{
					string parameter;
					in >> parameter;

					if (parameter.empty() || (parameter[0] == '#'))
					{
						cout << "unary property\t;\t"
							<< "size:" << '\t' << U_props_double_[var_name]->size()
							<< "\t;\t type:" << '\t' << "double" << '\n';
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							": incorrect parameter '" << parameter << "'.\n";
						return CMD_ERROR;
					}
					break;
				}

				case (TPI_UPROP_INT):
				{
					string parameter;
					in >> parameter;

					if (parameter.empty() || (parameter[0] == '#'))
					{
						cout << "unary property\t;\t"  
							<< "size:" << '\t' << U_props_int_[var_name]->size()
							<< "\t;\t type:" << '\t' << "int" << '\n';
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							": incorrect parameter '" << parameter << "'.\n";
						return CMD_ERROR;
					}
					break;
				}

				case (TPI_ARRAY_OF_VECTORS):
				{
					string parameter;
					in >> parameter;


					auto var_it = arrays2d_.find(var_name);
					if (parameter.empty() || (parameter[0] == '#'))
					{
						for (unsigned i = 0; (i < 10) && (i < var_it->second.size()); ++i)
						{
							cout << var_name << '[' << i << ']' << '\t' <<
									var_it->second[i][0] << '\t' <<
									var_it->second[i][1] << '\t' <<
									var_it->second[i][2] << '\n';
						}
					}
					else
					{
						auto it = cmds_.find(parameter);
						if (it != cmds_.end())
						{
							if (it->second == CMD_SIZE)
							{
								string comment;
								in >> comment;
								if (!(comment[0] == '#') && !(comment.empty()))
								{
									cout << "ERROR in line " << program_counter <<
										": incorrect parameter '" << comment << "'.\n";
									return CMD_ERROR;
								}
								else
								{
									cout << "size of '" << var_name << "': " << var_it->second.size() << '\n';
								}
							}
							else
							{
								cout << "ERROR in line " << program_counter <<
									": incorrect command '" << parameter << "'.\n";
								return CMD_ERROR;
							}
						}
						else
						{
							bool parameter_error = false;
							unsigned char_it = 0;
							while (!parameter_error && (char_it < (parameter.size())))
							{
								if (!isdigit(parameter[char_it]))
									parameter_error = true;
								else ++char_it;
							}

							if (!parameter_error)
							{
								string second_parameter;
								in >> second_parameter;
								if (!second_parameter.empty() && (second_parameter[0] != '#'))
								{
									unsigned char_it = 0;
									while (!parameter_error && (char_it < (second_parameter.size())))
									{
										if (!isdigit(parameter[char_it]))
											parameter_error = true;
										else ++char_it;
									}

									if (!parameter_error)
									{
										string comment;
										in >> comment;
										if (!comment.empty() && (comment[0] != '#'))
										{
											cout << "ERROR in line " << program_counter <<
												": incorrect parameter '" << comment << "'.\n";
											return CMD_ERROR;
										}
										else
										{
											int left = atoi(parameter.c_str());
											int right = atoi(second_parameter.c_str());
											if ((right < var_it->second.size())
												&& (left < right))
											{
												for (int i = left; i < right; ++i)
												{
													cout << var_name << '[' << i << ']' << '\t' <<
														var_it->second[i][0] << '\t' <<
														var_it->second[i][1] << '\t' <<
														var_it->second[i][2] << '\n';
												}
											}
											else
											{
												cout << "ERROR in line " << program_counter <<
													": incorrect parameters.\n";
												return CMD_ERROR;
											}
										}
									}
								}
								else
								{
									int i = atoi(parameter.c_str());
									if (i < var_it->second.size())
									{
										cout << var_name << '[' << i << ']' << '\t' <<
											var_it->second[i][0] << '\t' <<
											var_it->second[i][1] << '\t' <<
											var_it->second[i][2] << '\n';
									}
									else
									{
										cout << "ERROR in line " << program_counter <<
											": incorrect parameter '" << parameter << "'.\n";
										return CMD_ERROR;
									}
								}
							}
							else
							{
								cout << "ERROR in line " << program_counter <<
									": incorrect parameter '" << parameter << "'.\n";
								return CMD_ERROR;
							}
						}
					}
					break;
				}
							
				default:
					{
						cout << "ERROR in line " << program_counter << ":undeclared variable '" << var_name << "'.\n";
						return CMD_ERROR;
					}
				}
			break;
			}

			case (CMD_CREATE):
			{
				string _val_type;
				in >> _val_type;
				
				auto _val_type_id = cmds_.find(_val_type);
				switch (_val_type_id->second)
				{
				case (CMD_ATOMICSTRUCT):
					{
						string comment;
						string struct_name;
						in >> struct_name;
						in >> comment;

						if (!(comment[0] == '#') && !(comment.empty()))
						{
							cout << "ERROR in line " << program_counter <<
								": command 'createstruct' must have one parameter only." << '\n';
							return CMD_ERROR;
						}
						else
						{
							if (isalpha(struct_name[0]))
							{
								auto it = types_.find(struct_name);
								if ((it == types_.end()) || (it->second == TPI_ATOMICSTRUCT))
								{
									if (it != types_.end())
									{
										AtomicStructs_.erase(struct_name);
										AtomicStructure *structure = new AtomicStructure(*atomic_type_table_);
										AtomicStructs_[struct_name] = structure;
										cout << "structure '" << struct_name << "' successfully created." << '\n';
									}
									else
									{
										AtomicStructure *structure = new AtomicStructure(*atomic_type_table_);
										AtomicStructs_[struct_name] = structure;
										types_[struct_name] = TPI_ATOMICSTRUCT;
										cout << "structure '" << struct_name << "' successfully created." << '\n';
									}
								}
								else
								{
									cout << "ERROR in line " << program_counter <<
										": variable's name already used for other type variable." << '\n';
									return CMD_ERROR;
								}
							}
							else
							{
								cout << "ERROR in line " << program_counter <<
									": variable's name must begin with letter." << '\n';
								return CMD_ERROR;
							}
						}
						break;
					}

				case (CMD_VECTOR):
				{
					string vector_name;
					in >> vector_name;

					if (isalpha(vector_name[0]))
					{
						auto it = types_.find(vector_name);
						if ((it == types_.end()) || (it->second == TPI_VECTOR))
						{
							vector<double> parameters;
							parameters.resize(3);
							bool parameter_input_error = false;
							for (unsigned parameters_it = 0; parameters_it < 3; ++parameters_it)
							{
								in >> parameters[parameters_it];
								if (!in)
								{
									parameters_it = 3;
									parameter_input_error = true;
									in.clear();
								}
							}

							if (!parameter_input_error)
							{
								string comment;
								in >> comment;
								if (!(comment[0] == '#') && !(comment.empty()))
								{
									cout << "ERROR in line " << program_counter <<
										": vector must have three parameters." << '\n';
									return CMD_ERROR;
								}
								else
								{
									if (it != types_.end())
									{
										vectors_[vector_name][0] = parameters[0];
										vectors_[vector_name][1] = parameters[1];
										vectors_[vector_name][2] = parameters[2];
										cout << vector_name << " = {" << vectors_[vector_name][0] << ','
											<< vectors_[vector_name][1] << ',' << vectors_[vector_name][2]
											<< '}' << '\n';
									}
									else
									{
										double* tmp = new double[3];
										tmp[0] = parameters[0];
										tmp[1] = parameters[1];
										tmp[2] = parameters[2];
										vectors_[vector_name] = tmp;
										types_[vector_name] = TPI_VECTOR;
										cout << vector_name << " = {" << vectors_[vector_name][0] << ','
											<< vectors_[vector_name][1] << ',' << vectors_[vector_name][2]
											<< '}' << '\n';
									}
								}
							}
							else
							{
								cout << "ERROR in line " << program_counter <<
									": incorrect vector's parameters." << '\n';
								return CMD_ERROR;
							}

						}
						else
						{
							cout << "ERROR in line " << program_counter <<
								": variable's name already used for other type variable." << '\n';
							return CMD_ERROR;
						}
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							": variable's name must begin with letter." << '\n';
						return CMD_ERROR;
					}
					break;
				}

				case (CMD_UPROP):
				{
					string var_name;
					in >> var_name;

					if (isalpha(var_name[0]))
					{
						string type_of_prop;
						in >> type_of_prop;

						if (isalpha(var_name[0]))
						{
							auto type_it = cmds_.find(type_of_prop);
							switch (type_it->second)
							{
							case (CMD_INT):
							{
								string input_file_path;
								char first_symb;
								char symb = 0;
								in.get(first_symb);
								in.get(first_symb);
								if ((first_symb == 34) || (first_symb == 39))
								{
									while (symb != first_symb) //eof ����� �������� ��� ����� ��� ���� ��� ������
									{
										if (in.eof())
										{
											cout << "ERROR in line " << program_counter <<
												": requires a closing quote after the filename." << '\n';
											return CMD_ERROR;
										}
										else
										{
											in.get(symb);
											if ((symb != first_symb))
											{
												input_file_path.push_back(symb);
											}
										}
									}

									string comment;
									in >> comment;
									if (!(comment[0] == '#') && !(comment.empty()))
									{
										cout << "ERROR in line " << program_counter <<
											": command 'crtuprop' must have three parameters." << '\n';
										return CMD_ERROR;
									}
									else
									{
										auto var_it = types_.find(var_name);
										if ((var_it == types_.end()) || (var_it->second == TPI_UPROP_INT))
										{
											if (var_it != types_.end())
											{
												ifstream input_file(input_file_path);
												U_props_int_.erase(var_name);
												UnaryPropertyTable<int> *tmp_property = new UnaryPropertyTable<int>(*atomic_type_table_);
												tmp_property->read_property_table(input_file);
												U_props_int_[var_name] = tmp_property;
												cout << "property '" << var_name << "' successfully created, "
													<< "property size: " << U_props_int_[var_name]->size() << '\n';
											}
											else
											{
												ifstream input_file(input_file_path);
												UnaryPropertyTable<int> *tmp_property = new UnaryPropertyTable<int>(*atomic_type_table_);
												tmp_property->read_property_table(input_file);
												U_props_int_[var_name] = tmp_property;
												types_[var_name] = TPI_UPROP_INT;
												cout << "property '" << var_name << "' successfully created, "
													<< "property size: " << U_props_int_[var_name]->size() << '\n';
											}
										}
										else
										{
											cout << "ERROR in line " << program_counter <<
												": variable's name already used for other type variable." << '\n';
											return CMD_ERROR;
										}
									}
								}
								else
								{
									cout << "ERROR in line " << program_counter <<
										":required filename in quotes." << '\n';
									return CMD_ERROR;
								}
								break;
							}

							case (CMD_DOUBLE):
							{
								string input_file_path;
								char first_symb;
								char symb = 0;
								in.get(first_symb);
								in.get(first_symb);
								if ((first_symb == 34) || (first_symb == 39))
								{
									while (symb != first_symb) //eof ����� �������� ��� ����� ��� ���� ��� ������
									{
										if (in.eof())
										{
											cout << "ERROR in line " << program_counter <<
												": requires a closing quote after the filename." << '\n';
											return CMD_ERROR;
										}
										else
										{
											in.get(symb);
											if ((symb != first_symb))
											{
												input_file_path.push_back(symb);
											}
										}
									}

									string comment;
									in >> comment;
									if (!(comment[0] == '#') && !(comment.empty()))
									{
										cout << "ERROR in line " << program_counter <<
											": command 'crtuprop' must have three parameters." << '\n';
										return CMD_ERROR;
									}
									else
									{
										auto var_it = types_.find(var_name);
										if ((var_it == types_.end()) || (var_it->second == TPI_UPROP_DOUBLE))
										{
											if (var_it != types_.end())
											{
												ifstream input_file(input_file_path);
												U_props_double_.erase(var_name);
												UnaryPropertyTable<double> *tmp_property
													= new UnaryPropertyTable<double>(*atomic_type_table_);
												tmp_property->read_property_table(input_file);
												U_props_double_[var_name] = tmp_property;
												cout << "property '" << var_name << "' successfully created, "
													<< "property size: " << U_props_double_[var_name]->size() << '\n';
											}
											else
											{
												ifstream input_file(input_file_path);
												UnaryPropertyTable<double> *tmp_property
													= new UnaryPropertyTable<double>(*atomic_type_table_);
												tmp_property->read_property_table(input_file);
												U_props_double_[var_name] = tmp_property;
												types_[var_name] = TPI_UPROP_DOUBLE;
												cout << "property '" << var_name << "' successfully created, "
													<< "property size: " << U_props_double_[var_name]->size() << '\n';
											}
										}
										else
										{
											cout << "ERROR in line " << program_counter <<
												": variable's name already used for other type variable." << '\n';
											return CMD_ERROR;
										}
									}
								}
								else
								{
									cout << "ERROR in line " << program_counter <<
										":required filename in quotes." << '\n';
									return CMD_ERROR;
								}
								break;
							}

							default:
							{
								cout << "ERROR in line " << program_counter <<
									": incorrect type of property '" << type_of_prop << "'.\n";
								return CMD_ERROR;
							}
							}
						}
						else
						{
							cout << "ERROR in line " << program_counter <<
								": variable's name must begin with letter." << '\n';
							return CMD_ERROR;
						}

					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							": variable's name must begin with letter." << '\n';
						return CMD_ERROR;
					}
					break;
				}

				case (CMD_BINPROP):
				{
					string var_name;
					in >> var_name;

					if (isalpha(var_name[0]))
					{
						string type_of_prop;
						in >> type_of_prop;

						auto type_it = cmds_.find(type_of_prop);
						switch (type_it->second)
						{
						case (CMD_INT):
						{
							string input_file_path;
							char first_symb;
							char symb = 0;
							in.get(first_symb);
							in.get(first_symb);
							if ((first_symb == 34) || (first_symb == 39))
							{
								while (symb != first_symb) //eof ����� �������� ��� ����� ��� ���� ��� ������
								{
									if (in.eof())
									{
										cout << "ERROR in line " << program_counter <<
											": requires a closing quote after the filename." << '\n';
										return CMD_ERROR;
									}
									else
									{
										in.get(symb);
										if ((symb != first_symb))
										{
											input_file_path.push_back(symb);
										}
									}
								}

								string comment;
								in >> comment;
								if (!(comment[0] == '#') && !(comment.empty()))
								{
									cout << "ERROR in line " << program_counter <<
										": command 'crtbinprop' must have three parameters." << '\n';
									return CMD_ERROR;
								}
								else
								{
									auto var_it = types_.find(var_name);
									if ((var_it == types_.end()) || (var_it->second == TPI_BINPROP_INT))
									{
										if (var_it != types_.end())
										{
											ifstream input_file(input_file_path);
											Bin_props_int_.erase(var_name);
											BinaryPropertyTable<int> *tmp_property
												= new BinaryPropertyTable<int>(*atomic_type_table_);
											tmp_property->read_property_table(input_file);
											Bin_props_int_[var_name] = tmp_property;
											cout << "property '" << var_name << "' successfully created, "
												<< "property size: " << Bin_props_int_[var_name]->size() << '\n';
										}
										else
										{
											ifstream input_file(input_file_path);
											BinaryPropertyTable<int> *tmp_property
												= new BinaryPropertyTable<int>(*atomic_type_table_);
											tmp_property->read_property_table(input_file);
											Bin_props_int_[var_name] = tmp_property;
											types_[var_name] = TPI_BINPROP_INT;
											cout << "property '" << var_name << "' successfully created, "
												<< "property size: " << Bin_props_int_[var_name]->size() << '\n';
										}
									}
									else
									{
										cout << "ERROR in line " << program_counter <<
											": variable's name already used for other type variable." << '\n';
										return CMD_ERROR;
									}
								}
							}
							else
							{
								cout << "ERROR in line " << program_counter <<
									":required filename in quotes." << '\n';
								return CMD_ERROR;
							}
							break;
						}

						case (CMD_DOUBLE):
						{
							string input_file_path;
							char first_symb;
							char symb = 0;
							in.get(first_symb);
							in.get(first_symb);
							if ((first_symb == 34) || (first_symb == 39))
							{
								while (symb != first_symb) //eof ����� �������� ��� ����� ��� ���� ��� ������
								{
									if (in.eof())
									{
										cout << "ERROR in line " << program_counter <<
											": requires a closing quote after the filename." << '\n';
										return CMD_ERROR;
									}
									else
									{
										in.get(symb);
										if ((symb != first_symb))
										{
											input_file_path.push_back(symb);
										}
									}
								}

								string comment;
								in >> comment;
								if (!(comment[0] == '#') && !(comment.empty()))
								{
									cout << "ERROR in line " << program_counter <<
										": command 'crtbinprop' must have three parameters." << '\n';
									return CMD_ERROR;
								}
								else
								{
									auto var_it = types_.find(var_name);
									if ((var_it == types_.end()) || (var_it->second == TPI_BINPROP_DOUBLE))
									{
										if (var_it != types_.end())
										{
											ifstream input_file(input_file_path);
											Bin_props_double_.erase(var_name);
											BinaryPropertyTable<double> *tmp_property
												= new BinaryPropertyTable<double>(*atomic_type_table_);
											tmp_property->read_property_table(input_file);
											Bin_props_double_[var_name] = tmp_property;
											cout << "property '" << var_name << "' successfully created, "
												<< "property size: " << Bin_props_double_[var_name]->size() << '\n';
										}
										else
										{
											ifstream input_file(input_file_path);
											BinaryPropertyTable<double> *tmp_property
												= new BinaryPropertyTable<double>(*atomic_type_table_);
											tmp_property->read_property_table(input_file);
											Bin_props_double_[var_name] = tmp_property;
											types_[var_name] = TPI_BINPROP_DOUBLE;
											cout << "property '" << var_name << "' successfully created, "
												<< "property size: " << Bin_props_double_[var_name]->size() << '\n';
										}
									}
									else
									{
										cout << "ERROR in line " << program_counter <<
											": variable's name already used for other type variable." << '\n';
										return CMD_ERROR;
									}
								}
							}
							else
							{
								cout << "ERROR in line " << program_counter <<
									":required filename in quotes." << '\n';
								return CMD_ERROR;
							}
							break;
						}

						default:
						{
							cout << "ERROR in line " << program_counter <<
								": incorrect type of property '" << type_of_prop << "'.\n";
							return CMD_ERROR;
						}
						}

					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							": variable's name must begin with letter." << '\n';
						return CMD_ERROR;
					}
					break;
				}

				default:
				{
					cout << "ERROR in line " << program_counter <<
						": incorrect value type '" << _val_type << "'.\n";
					return CMD_ERROR;
				}
				
				}
				break;
			}

			case (CMD_DELETE):
			{
				string var_name;
				string comment;
				in >> var_name;
				in >> comment;
				if (!(comment[0] == '#') && !(comment.empty()))
				{
					cout << "ERROR in line " << program_counter <<
						": command 'delete' have one parameter only." << '\n';
					return CMD_ERROR;
				}
				else
				{
					auto it = types_.find(var_name);
					if (it != types_.end())
					{
						switch (it->second)
						{
						case (TPI_ATOMICSTRUCT):
						{
							AtomicStructs_.erase(var_name);
							cout << "structure '" << var_name << "' succesfully deleted." << '\n';
							break;
						}
						case (TPI_VECTOR):
						{
							vectors_.erase(var_name);
							cout << "vector '" << var_name << "' succesfully deleted." << '\n';
							break;
						}
						case (TPI_BINPROP_DOUBLE):
						{
							Bin_props_double_.erase(var_name);
							cout << "binary property(double) '" << var_name << "' succesfully deleted." << '\n';
							break;
						}
						case (TPI_BINPROP_INT):
						{
							Bin_props_int_.erase(var_name);
							cout << "binary property(int) '" << var_name << "' succesfully deleted." << '\n';
							break;
						}
						case (TPI_UPROP_DOUBLE):
						{
							U_props_double_.erase(var_name);
							cout << "unary property(double) '" << var_name << "' succesfully deleted." << '\n';
							break;
						}
						case (TPI_UPROP_INT):
						{
							U_props_int_.erase(var_name);
							cout << "unary property(int) '" << var_name << "' succesfully deleted." << '\n';
							break;
						}
						}
						types_.erase(it);
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							": variable '" << var_name << "' not exist." << '\n';
						return CMD_ERROR;
					}
				}
				break;
			}

			case (CMD_BVS):
			{
				string var_name;
				in >> var_name;

				auto var = AtomicStructs_.find(var_name);
				if (var != AtomicStructs_.end())
				{
					string r_c_str;
					in >> r_c_str;

					auto r_c = Bin_props_double_.find(r_c_str);
					if (r_c != Bin_props_double_.end())
					{
						string r_0_str;
						in >> r_0_str;

						auto r_0 = Bin_props_double_.find(r_0_str);
						if (r_0 != Bin_props_double_.end())
						{
							string b_str;
							in >> b_str;

							auto b = Bin_props_double_.find(b_str);
							if (b != Bin_props_double_.end())
							{
								string valence_str;
								in >> valence_str;

								auto valence = U_props_int_.find(valence_str);
								if (valence != U_props_int_.end())
								{
									double precision_;
									in >> precision_;

									if (!in)
									{
										cout << "ERROR in line " << program_counter <<
											": precision parameter incorrect." << '\n';
										return CMD_ERROR;
									}
									else
									{
										string comment;
										in >> comment;

										if (!(comment[0] == '#') && !(comment.empty()))
										{
											cout << "ERROR in line " << program_counter <<
												": command 'bvs' must have six parameters." << '\n';
											return CMD_ERROR;
										}
										else
										{
											ofstream out("BVS_hyst.dat");
											double t_bvs_start, t_bvs_end,
												t_gii_start, t_gii_end;
											vector<int> bvs_hyst;
											vector<double> bvs_GII;
											double GII_;

											cout << "N_atoms:" << "\t" << AtomicStructs_[var_name]->size() << endl;

											r_c->second->update_property_table();
											r_0->second->update_property_table();
											b->second->update_property_table();
											valence->second->update_property_table();

											t_bvs_start = omp_get_wtime();

											AtomicStructs_[var_name]->BVS(*r_c->second, *r_0->second, *b->second,
												bvs_hyst, bvs_GII, precision_);

											t_bvs_end = omp_get_wtime();
											cout << "BVS analysis time:" << "\t" << (t_bvs_end - t_bvs_start) << endl;

											t_gii_start = omp_get_wtime();
											GII_ = AtomicStructs_[var_name]->GII(bvs_GII, precision_, *valence->second);
											t_gii_end = omp_get_wtime();

											cout << "GII analysis time:" << "\t" << (t_gii_end - t_gii_start) << endl;

											cout << "structure GII:" << "\t" << GII_ << endl;

											plot_hyst("BVS_hyst.dat", "BVS ���������", "���������� ������", "�������� BVS",
												precision_, bvs_hyst, true);
											//plot_hyst("BVS_results.dat", "BVS ���������", "����� �� �����������", "BVS",
											//	1.0, bvs_GII);
											out.close();
										}
									}
								}
								else
								{
									cout << "ERROR in line " << program_counter <<
										": parameter '" << valence_str << "' incorrect." << '\n';
									return CMD_ERROR;
								}

							}
							else
							{
								cout << "ERROR in line " << program_counter <<
									": parameter '" << b_str << "' incorrect." << '\n';
								return CMD_ERROR;
							}
						}
						else
						{
							cout << "ERROR in line " << program_counter <<
								": parameter '" << r_0_str << "' incorrect." << '\n';
							return CMD_ERROR;
						}
					}
					else
					{
						cout << "ERROR in line " << program_counter <<
							": parameter '" << r_c_str << "' incorrect." << '\n';
						return CMD_ERROR;
					}
				}
				else
				{
					cout << "ERROR in line " << program_counter <<
						": variable '" << var_name << "' incorrect." << '\n';
					return CMD_ERROR;
				}
				break;
			}

			default:
			{
				cout << "ERROR in line " << program_counter << 
					": incorrect use of the command '" << first_console_word << "'.\n";
				return CMD_ERROR;
			}
			}
		}
		else
		{
			if (vectors_.find(first_console_word) != vectors_.end())
			{
				string second_var;
				string operator_;
				in >> operator_;
				
				if (!operator_.empty())
				{
					auto cmd_ = cmds_.find(operator_);
					switch (cmd_->second)
					{
					case (CMD_EQUAL):
					{
						in >> second_var;
						if (vectors_.find(second_var) != vectors_.end())
						{
							vectors_[first_console_word] = vectors_[second_var];
							cout << first_console_word << " = {" << vectors_[first_console_word][0] << ','
								<< vectors_[first_console_word][1] << ',' << vectors_[first_console_word][2]
								<< '}' << '\n';
						}
						else
						{
							cout << "ERROR in line " << program_counter << ":undeclared variable '" << second_var << "'.\n";
							return CMD_ERROR;
						}
						break;
					}
					default:
					{
						cout << "ERROR in line " << program_counter <<
							":no such operator'" << operator_ << "' for " << first_console_word << '\n';
						return CMD_ERROR;
					}
					}
				}
				else
				{
					cout << first_console_word << " = {" << vectors_[first_console_word][0] << ','
						<< vectors_[first_console_word][1] << ',' << vectors_[first_console_word][2]
						<< '}' << '\n';
				}
			}
			else
			{
				cout << "ERROR in line " << program_counter <<
					":'" << first_console_word << "' command not exist." << '\n';
				return CMD_ERROR;
			}
		}
	}


public: 
	
	void screen_mode()
	{
		int program_counter=0;
		bool stop = false;
		cmds_id_ result;
		string cmd_string;
		//std::cout << "(~>_<)~ \\( ^w^)/ AtomicStructure \\(^w^ )/ ~(>_<~)" << '\n';
		std::cout << "----AtomicStructure----" << '\n';

		while (!stop)
		{
			++program_counter;
			cout << program_counter << ": ->";
			std::getline(std::cin,cmd_string);
			if (!cmd_string.empty())
			{
				std::stringstream cmd_in(cmd_string);
				result = cmd_process_(cmd_in, program_counter);
				if (result == CMD_EXIT) stop = true;
				cmd_in.clear();
			}
		}
	}

	bool script_mode(istream &fin)
	{
		int program_counter = 0;
		bool stop = false;
		cmds_id_ result;
		string cmd_string;
		std::cout << "----AtomicStructure----" << '\n';

		while (!stop)
		{
			++program_counter;
			std::getline(fin, cmd_string);
			if (!cmd_string.empty())
			{
				std::stringstream cmd_in((string)cmd_string);
				result = cmd_process_(cmd_in, program_counter);
				if ((result == CMD_EXIT) || (result == CMD_ERROR)) stop = true;
				cmd_in.clear();
			}
			
		}

		if (result == CMD_ERROR) return false;
		else return true;
	}
};

#endif