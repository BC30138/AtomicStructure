#include <iostream>

#include "AtomicStructure.h"
#include "tests.hpp"

// int main()
// {
// 	AtomicTypeTable AT;
// 	ifstream in("XYZ//NaCl_cool.xyz"); 
// 	AtomicStructure a(AT);
// 	int size[3] = {3,3,3};
// 	a.create(in, size);
// 	a.show();
// 	a.info();
// 	//test_dist_analysis(a, &AT);
// 	return 0;
// }

int main()
{
	AtomicTypeTable AtomType;
	ifstream in("XYZ/NaCl_cool.xyz");
	AtomicStructure structure(AtomType);
	int size[3] = {6*5,6*5,6*5};
	// int size[3] = {1,1,1};
	// int size[3] = {3,3,3};
	structure.create(in, size, 2);
	plot_structure("XYZ/NaClpresblock.xyz");
	// structure.show();
	structure.info();
	long long int comp;	
	test_dist_analysis(structure, &AtomType, comp); 
	// big_time_test_DIST("XYZ/NaCl_cool.xyz", 25, 5);
	// conflict();
	return 0;
}