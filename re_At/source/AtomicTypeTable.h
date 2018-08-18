#ifndef _ATOMICTYPETABLE_H_
#define _ATOMICTYPETABLE_H_

#include <string> 
#include <fstream> 
#include <vector>
#include <unordered_map>
#include <map>
#include <iostream>
using namespace std;

#define at_empty -127
#define at_undefined -99

union element
{
	char id_char[8];
	long long int long_id;
};

class AtomicTypeTable //Dictionary 
{
private:
	void add_element(element);
	void add_element(long long int elementlongid);
public:
	vector<element> elements_;
	unordered_map<long long int, int> longid_to_id_table;

	AtomicTypeTable() {}
	size_t size();
	void plus(AtomicTypeTable &B);
	void read_type_table_from_file(string txtname);
	element get_element(size_t Id);
	element get_element(int Id);
	element get_element(string Symb);
	element get_element(long long int longid);
	int get_Id_with_adding(string);
	int get_Id_with_adding(element);
	int get_Id_with_adding(long long int);
	int get_id(string Symb);
	int get_id(element);
	int get_id(long long int);
	void show_table();
	vector<int> string_vector_to_id_vector_wth_exist_check(vector<string> type);
	~AtomicTypeTable();
};

void empty_till_end(char(&)[8], int);
void comment_read(string, char(&)[8], int);
long long int string_to_longid(string);
string longid_to_string(long long int);

#endif 