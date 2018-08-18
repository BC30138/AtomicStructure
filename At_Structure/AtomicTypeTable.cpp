#include "AtomicTypeTable.h"

using namespace std;

void empty_till_end(char(&idch)[8], int firstempty)
{

	for (unsigned long i = firstempty; i < 8; ++i)
	{
		idch[i] = at_empty;
	}
}

void comment_read(string inputstring, char(&idch)[8], int firstletter)
{
	int size = inputstring.size();
	switch (firstletter)
	{
	case 2:
	{
		if (size < 8)
		{
			for (unsigned long i = 2; i < size; ++i)
			{
				idch[i + 1] = inputstring[i];
			}
			empty_till_end(idch, size + 1);
		}
		else
		{
			cout << "Error:comment max size = 5 symbols" << endl;
			system("pause");
		}
		break;
	}
	case 3:
	{

		if ((size < 9) && (size > 3))
		{
			for (unsigned long i = 3; i < size; ++i)
			{
				idch[i] = inputstring[i];
			}
			empty_till_end(idch, size);
		}
		else
		{
			if (size <= 3)
			{
				empty_till_end(idch, 3);
			}
			else
			{
				cout << "Error:comment max size = 5 symbols" << endl;
				system("pause");
			}
		}
		break;
	}
	case 4:
	{
		if (size < 10)
		{
			for (unsigned long i = 4; i < size; ++i)
			{
				idch[i - 1] = inputstring[i];
			}
			empty_till_end(idch, size - 1);
		}
		else
		{
			cout << "Error:comment max size = 5 symbols" << endl;
			system("pause");
		}
		break;
	}
	default:
	{
		cout << "wrong!" << endl;
	}
	}
}

long long int string_to_longid(string inputstring)
{
	element tmp;

	if (isalpha(inputstring[0]))
	{
		tmp.id_char[0] = toupper(inputstring[0]);
		if (isalpha(inputstring[1]))  // inputstring[1] - �����
		{
			tmp.id_char[1] = tolower(inputstring[1]);
			if (isdigit(inputstring[2]) && (inputstring[2] != '0')) //inputstring[2] - ����� �� ������ ����
			{
				switch (inputstring[3])
				{
				case '+':
				{
					tmp.id_char[2] = inputstring[2] - '0';
					comment_read(inputstring, tmp.id_char, 4);
					break;
				}

				case '-':
				{
					tmp.id_char[2] = (inputstring[2] - '0') * (-1);
					comment_read(inputstring, tmp.id_char, 4);
					break;
				}

				default:
				{
					cout << "Error:You must enter valence number's sign" << endl;
					system("pause");
				}
				}
			}

			else
			{
				if (inputstring.size() > 2)
				{

					switch (inputstring[2])
					{
					case '+':
					{
						if (isdigit(inputstring[3]))
						{
							tmp.id_char[2] = inputstring[3] - '0';
							comment_read(inputstring, tmp.id_char, 4);
						}
						else
						{
							cout << "Error:You must enter valence number";
							system("pause");
						}
						break;
					}

					case '-':
					{
						if (isdigit(inputstring[3]))
						{
							tmp.id_char[2] = (inputstring[3] - '0') * (-1);
							comment_read(inputstring, tmp.id_char, 4);
						}
						else
						{
							cout << "Error:You must enter valence number";
							system("pause");
						}
						break;
					}

					case '0':
					{
						tmp.id_char[2] = 0;
						comment_read(inputstring, tmp.id_char, 3);
						break;
					}

					case '*':
					{
						tmp.id_char[2] = at_undefined;
						comment_read(inputstring, tmp.id_char, 3);
						break;
					}

					default:
					{
						cout << "Error: 3rd symbol can be: Let,Num,'*','+' or '-'" << endl;
					}
					}
				}
				else {
					tmp.id_char[2] = at_undefined;
					comment_read(inputstring, tmp.id_char, 3);
				}
			}
		}

		else {
			tmp.id_char[1] = at_empty;

			if (isdigit(inputstring[1]))
			{
				switch (inputstring[2])
				{
				case '+':
				{
					tmp.id_char[2] = inputstring[1] - '0';
					comment_read(inputstring, tmp.id_char, 3);
					break;
				}

				case '-':
				{
					tmp.id_char[2] = (inputstring[1] - '0') * (-1);
					comment_read(inputstring, tmp.id_char, 3);
					break;
				}

				default:
				{
					cout << "Error:You must enter valence number's sign" << endl;
					system("pause");
				}
				}
			}
			else
			{
				if (inputstring.size() > 1)
				{
					switch (inputstring[1])
					{
					case '+':
					{
						if (isdigit(inputstring[2]))
						{
							tmp.id_char[2] = inputstring[2] - '0';
							comment_read(inputstring, tmp.id_char, 3);
						}
						else
						{
							cout << "Error:You must enter valence number" << endl;
							system("pause");
						}
						break;
					}

					case '-':
					{
						if (isdigit(inputstring[2]))
						{
							tmp.id_char[2] = (inputstring[2] - '0') * (-1);
							comment_read(inputstring, tmp.id_char, 3);
						}
						else
						{
							cout << "Error:You must enter valence number" << endl;
							system("pause");
						}
						break;
					}

					case '_':
					{
						if (isdigit(inputstring[2]))
						{
							switch (inputstring[3])
							{
							case '+':
							{
								tmp.id_char[2] = inputstring[2] - '0';
								comment_read(inputstring, tmp.id_char, 4);
								break;
							}

							case '-':
							{
								tmp.id_char[2] = (inputstring[2] - '0') * (-1);
								comment_read(inputstring, tmp.id_char, 4);
								break;
							}
							default:
							{
								cout << "Error:You must enter valence number's sign" << endl;
								system("pause");
							}
							}
						}
						else
						{
							switch (inputstring[2])
							{
							case '+':
							{
								if (isdigit(inputstring[3]))
								{
									tmp.id_char[2] = inputstring[3] - '0';
									comment_read(inputstring, tmp.id_char, 4);
								}
								else
								{
									cout << "Error:You must enter valence number" << endl;
									system("pause");
								}
								break;
							}

							case '-':
							{
								if (isdigit(inputstring[3]))
								{
									tmp.id_char[2] = (inputstring[3] - '0')*(-1);
									comment_read(inputstring, tmp.id_char, 4);
								}
								else
								{
									cout << "Error:You must enter valence number" << endl;
									system("pause");
								}
								break;
							}

							case '*':
							{
								tmp.id_char[2] = at_undefined;
								comment_read(inputstring, tmp.id_char, 3);
								break;
							}
							default:
							{
								cout << "Error: 3rd symbol can be: Let,Num,'*','+' or '-'" << endl;
								system("pause");
							}
							}
						}
						break;
					}

					case '*':
					{
						tmp.id_char[2] = at_undefined;
						comment_read(inputstring, tmp.id_char, 2);
						break;
					}
					default:
					{
						cout << "Error: 2nd symbol can be: Let,Num,'*','+' or '-'" << endl;
						system("pause");
					}
					}
				}
				else
				{
					tmp.id_char[2] = at_undefined;
					comment_read(inputstring, tmp.id_char, 2);
				}
			}
		}
	}
	else
	{
		cout << "Error: 1st symbol can be : Let, Num, '*', '+' or '-'" << endl;
		system("pause");
	}

	return tmp.long_id;
}

string longid_to_string(long long int inputlongid)
{
	string output;
	element tmp;
	tmp.long_id = inputlongid;

	output.push_back(tmp.id_char[0]);
	if (tmp.id_char[1] != at_empty)
	{
		output.push_back(tmp.id_char[1]);
	}
	if (tmp.id_char[2] != at_undefined)
	{
		int num;
		num = tmp.id_char[2];
		if (num > 0)
		{
			output.push_back('+');
			output.push_back(num + '0');
		}
		else
		{
			if (num == 0)
			{
				output.push_back('0');
			}
			else
			{
				output.push_back(abs(num) + '0');
				output.push_back('-');
			}
		}
	}
	int i = 3;
	while (tmp.id_char[i] != at_empty && i != 8)
	{
		output.push_back(tmp.id_char[i]);
		++i;
	}
	return output;
}





int AtomicTypeTable::size() {
	return elements_.size();
}

void AtomicTypeTable::plus(AtomicTypeTable &B) {
	for (int i = 0; i < (B).size(); i++)
	{
		(*this).get_Id_with_adding(B.elements_[i].long_id);
	}
}

void AtomicTypeTable::add_element(element el)
{
	if (longid_to_id_table.find(el.long_id) == longid_to_id_table.end())
	{
		elements_.push_back(el);
		longid_to_id_table[el.long_id] = elements_.size() - 1;
	}
	else
	{
		cout << "Error:element type object alredy exist in AtomicTypeTable" << endl;
		system("pause");
	}
}

void AtomicTypeTable::add_element(long long int elementlongid)
{
	if (longid_to_id_table.find(elementlongid) == longid_to_id_table.end())
	{
		element Elem;
		Elem.long_id = elementlongid;
		elements_.push_back(Elem);
		longid_to_id_table[elementlongid] = elements_.size() - 1;
	}
	else
	{
		cout << "Error:element type object alredy exist in AtomicTypeTable" << endl;
		system("pause");
	}
}

void AtomicTypeTable::read_type_table_from_file(string txtname)
{
	element Elem;
	string strid;
	int intid;

	ifstream in(txtname);

	while (!in.eof())
	{
		in >> intid;
		in >> strid;
		Elem.long_id = string_to_longid(strid);
		elements_.push_back(Elem);
		longid_to_id_table[Elem.long_id] = elements_.size() - 1;
	}
}

AtomicTypeTable::~AtomicTypeTable()
{
	elements_.resize(0);
}

element AtomicTypeTable::get_element(int Id)
{
	if (Id < 0 && Id >(elements_.size() + 1))
	{
		cout << "not found such Id in Table" << endl;
		system("pause");
	}
	return elements_[Id];
}

element AtomicTypeTable::get_element(unsigned long Id)
{
	if (Id < 0 && Id >(elements_.size() + 1))
	{
		cout << "not found such Id in Table" << endl;
		system("pause");
	}
	return elements_[Id];
}

element AtomicTypeTable::get_element(string Symb)
{
	return elements_[longid_to_id_table[string_to_longid(Symb)]];
}

element AtomicTypeTable::get_element(long long int longid)
{
	return elements_[longid_to_id_table[longid]];
}

int AtomicTypeTable::get_Id_with_adding(long long int elementlongid)
{
	auto it = longid_to_id_table.find(elementlongid);
	if (it == longid_to_id_table.end())
	{
		(*this).add_element(elementlongid);
		return longid_to_id_table[elementlongid];
	}
	else return it->second;
}

int AtomicTypeTable::get_Id_with_adding(string Symb)
{
	long long int elementlongid = string_to_longid(Symb);
		
	auto it = longid_to_id_table.find(elementlongid);
	if (it == longid_to_id_table.end())
	{
		(*this).add_element(elementlongid);
		return longid_to_id_table[elementlongid];
	}
	else return it->second;
}

int AtomicTypeTable::get_Id_with_adding(element temp)
{
	auto it = longid_to_id_table.find(temp.long_id);
	if (it == longid_to_id_table.end())
	{
		(*this).add_element(temp);
		return longid_to_id_table[temp.long_id];
	}
	else return it->second;
}

int AtomicTypeTable::get_id(string Symb)
{
	auto it = longid_to_id_table.find(string_to_longid(Symb));
	if (it == longid_to_id_table.end())
	{
		return -1;
	}
	else return it->second;
}

int AtomicTypeTable::get_id(element temp)
{
	auto it = longid_to_id_table.find(temp.long_id);
	if (it == longid_to_id_table.end())
	{
		return -1;
	}
	else return it->second;
}

int AtomicTypeTable::get_id(long long int longid)
{
	auto it = longid_to_id_table.find(longid);
	if (it == longid_to_id_table.end())
	{
		return -1;
	}
	else return it->second;
}

void AtomicTypeTable::show_table()
{
	cout << "AtomicTypeTable:" << endl;
	for (int i = 0; i < elements_.size(); i++)
	{
		cout << i << "\t" << longid_to_string(elements_[i].long_id) << endl;
	}
	system("pause");
}

vector<int> AtomicTypeTable::string_vector_to_id_vector_wth_exist_check(vector<string> type) {
	vector<int> TypeIdVector;
	for (int i = 0; i < type.size(); i++)
	{
		int notexistcheck;
		notexistcheck = (*this).get_id(type[i]);
		if (notexistcheck != -1)
		{
			TypeIdVector.push_back(notexistcheck);
		}
	}
	return TypeIdVector;
}