#ifndef _PROPERTIES_H_
#define _PROPERTIES_H_

#include <map>
#include <iostream>
#include <string> 
#include <utility>
#include<functional>
#include "AtomicTypeTable.h"

template<typename Typename_>
class UnaryPropertyTable
{
	AtomicTypeTable& atomic_type_table_;
	int work_size_;
public:
	std::map<long long int, Typename_> longid_property_table_;
	std::vector<Typename_> id_property_table_;
	inline UnaryPropertyTable(AtomicTypeTable& atomic_type_table) :
		atomic_type_table_(atomic_type_table) {
		work_size_ = atomic_type_table.size();
		id_property_table_.resize(work_size_);
	}
	UnaryPropertyTable(std::function<Typename_(Typename_)> func,
		const UnaryPropertyTable<Typename_>& origin_property_table) :
		atomic_type_table_(origin_property_table.atomic_type_table_)
	{
		for (auto it : origin_property_table.longid_property_table_)
		{
			longid_property_table_[it.first] = func(it.second);
		}


		work_size_ = origin_property_table.work_size_;
		int sizeof_origin_id_property_table = origin_property_table.id_property_table_.size();
		id_property_table_.resize(sizeof_origin_id_property_table);

		for (int i = 0; i < sizeof_origin_id_property_table; ++i)
		{
			id_property_table_[i] = func(origin_property_table.id_property_table_[i]);
		}
	}
public:
	int size()
	{
		return longid_property_table_.size();
	}

	inline void read_property_table(std::istream& input_stream)
	{
		int swap_buffer;
		std::string symbol;
		int id;
		long long int longid;
		Typename_ property_;

		while (input_stream.good())
		{
			input_stream >> symbol;
			longid = string_to_longid(symbol);
			id = atomic_type_table_.get_id(symbol);
			input_stream >> property_;
			input_stream >> std::ws;

			longid_property_table_[longid] = property_;

			if (id != -1)
			{
				id_property_table_[id] = property_;
			}
		}
	}

	inline double get_property(const int& id) const
	{
		if (0 <= id && id < atomic_type_table_.size())
		{
			return id_property_table_[id];
		}
		else
		{
			cout << "Error: Element not found by ids" << endl;
			system("pause");
		}
	}

	inline double get_property(const std::string& symbol) const
	{
		int id = atomic_type_table_.get_id(symbol);
		if (id != -1)
		{
			return id_property_table_[id];
		}
		else
		{
			std::cout << "Error:Elements not found by symbols" << std::endl;
			system("pause");
		}
	}

	inline int update_property_table()
	{
		work_size_ = atomic_type_table_.size();
		id_property_table_.resize(work_size_);
		for (unsigned long i = 0; i < atomic_type_table_.size(); ++i)
		{
			auto it = longid_property_table_.find(atomic_type_table_.elements_[i].long_id);
			if (it == longid_property_table_.end())
			{
				i = atomic_type_table_.size();
				return -1;
			}
			else
			{
				id_property_table_[i] = (*it).second;
				if (i == work_size_ - 1)
				{
					return 0;
				}
			}
		}
	}

	inline void show_property_table()
	{
		for (auto it = longid_property_table_.begin(); it != longid_property_table_.end(); ++it)
		{
			cout << longid_to_string(it->first) << "\t" << it->second << '\n';
		}
	}

};



// struct PairCompare {
// 	inline bool operator()(const std::pair<int, int>& p1, const std::pair<int, int>& p2) const
// 	{
// 		return (p1.first < p2.first) || (p1.first == p2.first && p1.second < p2.second);
// 	}
// };

template<typename Typename_>
class BinaryPropertyTable
{
	AtomicTypeTable& atomic_type_table_;
	int work_size_;
	int size_;
public:
	std::map< std::pair<long long int, long long int>, Typename_> longid_property_table_;
	std::vector<Typename_> id_property_table_;
	inline BinaryPropertyTable(AtomicTypeTable& atomic_type_table) :
		atomic_type_table_(atomic_type_table) {
		size_ = 0;
		work_size_ = atomic_type_table.size();
		id_property_table_.resize(work_size_ * work_size_);
	}
	BinaryPropertyTable(std::function<Typename_ (Typename_)> func,
		const BinaryPropertyTable<Typename_>& origin_property_table) :
		atomic_type_table_(origin_property_table.atomic_type_table_)
	{
		for (auto it : origin_property_table.longid_property_table_)
		{
			longid_property_table_[it.first] = func(it.second);
		}


		work_size_ = origin_property_table.work_size_;
		int sizeof_origin_id_property_table = origin_property_table.id_property_table_.size();
		id_property_table_.resize(sizeof_origin_id_property_table);


		for (int i = 0; i < sizeof_origin_id_property_table; ++i)
		{
			id_property_table_[i] = func(origin_property_table.id_property_table_[i]);
		}
		size_ = origin_property_table.size_;
	}
public:
	size_t size()
	{
		return size_;
	}

	inline void read_property_table(std::istream& input_stream)
	{
		int swap_buffer;
		std::string symbol;
		std::pair<int, int> id_pair;
		std::pair<long long int, long long int> longid_pair;
		Typename_ property;
		

		while (input_stream.good())
		{
			input_stream >> symbol;
			longid_pair.first = string_to_longid(symbol);
			id_pair.first = atomic_type_table_.get_id(symbol);
			input_stream >> symbol;
			longid_pair.second = string_to_longid(symbol);
			id_pair.second = atomic_type_table_.get_id(symbol);
			input_stream >> property;
			input_stream >> std::ws;

			longid_property_table_[longid_pair] = property;
			swap_buffer = longid_pair.second;
			longid_pair.second = longid_pair.first;
			longid_pair.first = swap_buffer;
			longid_property_table_[longid_pair] = property;

			++size_;

			if (id_pair.first != -1 && id_pair.second != -1)
			{
				id_property_table_[id_pair.first*work_size_ + id_pair.second] = property;
				if (id_pair.first != id_pair.second)
				{
					id_property_table_[id_pair.second*work_size_ + id_pair.first] = property;
				}
			}
		}
	}

	inline double get_property(const int& id_first, const int& id_second) const
	{
		if (0 <= id_first && id_first < atomic_type_table_.size() &&
			0 <= id_second && id_second < atomic_type_table_.size())
		{
			return id_property_table_[id_first * work_size_ + id_second];
		}
		else
		{
			cout << "Error: Element not found by ids" << endl;
			system("pause");
		}
	}

	inline double get_property(const std::string& symbol_first,
		const std::string& symbol_second) const
	{
		int id_first = atomic_type_table_.get_id(symbol_first);
		int id_second = atomic_type_table_.get_id(symbol_second);
		if (id_first != -1 && id_second != -1)
		{
			return id_property_table_[id_first * work_size_ + id_second];
		}
		else
		{
			std::cout << "Error:Elements not found by symbols" << std::endl;
			system("pause");
		}
	}

	inline int update_property_table()
	{
		work_size_ = atomic_type_table_.size();
		id_property_table_.resize(work_size_*work_size_);
		for (unsigned long i = 0; i < atomic_type_table_.size(); ++i)
			for (unsigned long j = 0; j < atomic_type_table_.size(); ++j)
			{
				auto it = longid_property_table_.find(make_pair(atomic_type_table_.elements_[i].long_id,
					atomic_type_table_.elements_[j].long_id));
				if (it == longid_property_table_.end())
				{
					i = atomic_type_table_.size();
					j = atomic_type_table_.size();
					return -1;
				}
				else
				{
					id_property_table_[i*work_size_ + j] = (*it).second;
					if (i != j)
					{
						id_property_table_[j*work_size_ + i] = (*it).second;
					}
					if (i == work_size_ - 1 &&
						j == work_size_ - 1)
					{
						return 0;
					}
				}
			}
	}

	inline void show_property_table()
	{
		for (//std::map<std::pair<long long int,long long int>, Typename_, PairCompare>::iterator it 
			auto it
		= longid_property_table_.begin(); it != longid_property_table_.end(); ++it)
		{
			cout << longid_to_string(it->first.first) << "\t" << longid_to_string(it->first.second) << "\t" 
				<< it->second << '\n';
		}
	}
};

#endif


