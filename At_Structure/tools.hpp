#ifndef _TOOLS_H_
#define _TOOLS_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream> 
#include <vector>
#include <cstdio> 
#include <string> 
#include <Eigen/Dense>
#include <iostream>
#include <map>
using namespace std;
using namespace Eigen;


inline bool is_equal(double x, double y) {
	return fabs(x - y) < numeric_limits<double>::epsilon();
}

inline double Square(double x) {
	return x * x;
}

inline void plot_hyst(string filename, string HystName, string YaxisName, string XaxisName, double width,
	const vector<int> &hyst, bool is_frequency)
{
	ofstream out("data\\" + filename);
	ofstream plot("data\\plotHyst.dat");
	if (is_frequency)
	{
		for (int i = 0; i < hyst.size(); i++)
		{
			if (hyst[i] != 0)
				out << width * i << "\t" << hyst[i] << endl;
		}
	}
	else
	{
		for (int i = 0; i < hyst.size(); i++)
		{
			if (hyst[i] != 0)
				out << i << "\t" << hyst[i] << endl;
		}
	}

	//plot << "set encoding iso_8859_1" << endl;
	//plot << L"set terminal pdfcairo enhanced font 'Helvetica,5'" << endl;
	//plot << "set encoding utf8" << endl;

	//plot << "set terminal png size 800,600" << endl;
	//plot << "set output '"<< << '"
	plot << "width = " << width << endl;
	plot << "bin(x, s) = s*int(x / s)" << endl;
	plot << "set ylabel '" << YaxisName << "'" << endl;//font 'Times - Roman, 20'" << endl;
	plot << "set xlabel '" << XaxisName << "'" << endl;//font 'Times - Roman, 20'" << endl;
	plot << "set boxwidth width" << endl;
	plot << "plot '" << filename << "' u (bin($1, width) + width/2.0):($2) w boxes fs solid 0.5 title '" << HystName << "'"
		<< endl;
	plot << "set yrange[0:GPVAL_Y_MAX + 1]" << endl;
	plot << "set xrange[GPVAL_X_MIN - width:GPVAL_X_MAX + width]" << endl;
	plot << "replot" << endl;
	plot << "pause - 1" << endl;


	system("gnuplot data\\plotHyst.dat");
	plot.close();
	out.close();
}

inline void plot_hyst(string filename, string HystName, string YaxisName, string XaxisName, double width,
	const vector<double> &hyst)
{
	ofstream out("data\\" + filename);
	ofstream plot("data\\plotHyst.dat");
	for (int i = 0; i < hyst.size(); i++)
	{
		if (hyst[i] != 0)
			out << i << "\t" << hyst[i] << endl;
	}

	//plot << "set encoding iso_8859_1" << endl;
	//plot << L"set terminal pdfcairo enhanced font 'Helvetica,5'" << endl;
	//plot << "set encoding utf8" << endl;
	plot << "set size 0.6, 0.6" << endl;
	plot << "width = " << width << endl;
	plot << "bin(x, s) = s*int(x / s)" << endl;
	plot << "set ylabel '" << YaxisName << "'" << endl;//font 'Times - Roman, 20'" << endl;
	plot << "set xlabel '" << XaxisName << "'" << endl;//font 'Times - Roman, 20'" << endl;
	plot << "set boxwidth width" << endl;
	plot << "plot '" << filename << "' u (bin($1, width) + width/2.0):($2) w boxes fs solid 0.5 title '" << HystName << "'"
		<< endl;
	plot << "set yrange[0:GPVAL_Y_MAX + 1]" << endl;
	plot << "set xrange[GPVAL_X_MIN - width:GPVAL_X_MAX + width]" << endl;
	plot << "replot" << endl;
	plot << "pause - 1" << endl;


	system("gnuplot data\\plotHyst.dat");
	plot.close();
	out.close();
}

inline void log_plot_(string filename, string PlotName, string YaxisName, string XaxisName,
	const vector<double> &hyst)
{
	ofstream out("data\\" + filename);
	ofstream plot("data\\plotHyst.dat");
	for (int i = 0; i < hyst.size(); i++)
	{
		if (hyst[i] != 0)
			out << i << "\t" << hyst[i] << endl;
	}

	plot << "set logscale y" << endl;
	plot << "set ylabel '" << YaxisName << "'" << endl;
	plot << "set xlabel '" << XaxisName << "'" << endl;
	plot << "plot '" << filename << "' w lines lw 3 title '" << PlotName << "'"
		<< endl;
	plot << "pause - 1" << endl;


	system("gnuplot data\\plotHyst.dat");
	plot.close();
	out.close();
}

inline void plot_structure(string outputfilename) {
	ofstream plot("data\\plot.dat");
	plot << "set xlabel 'X axis'" << endl;
	plot << "set ylabel 'Y axis'" << endl;
	plot << "set zlabel 'Z axis'" << endl;
	plot << "splot '" << outputfilename << "' using 2:3 : 4 w points lt rgb 'red' pt 7 title '" << outputfilename << "'" << endl;
	plot << "pause - 1" << endl;

	system("gnuplot data\\plot.dat");
}

inline vector<string> input_file_names_list(string listname)
{
	ifstream in(listname);
	vector<string> filelist;
	string filename;
	while (in.peek() != EOF)
	{
		in >> filename;
		filelist.push_back(filename);
	}
	filelist.erase(filelist.end() - 1);
	return filelist;
}

inline double degree_to_rad(double degree)
{
	return degree * M_PI / 180;
}

class Rotation
{
	Matrix3d this_rot;
public:
	Rotation() 
	{
		clear();
	}

	void inite_around_X(double degree)
	{
		this_rot = AngleAxisd(degree_to_rad(degree), Vector3d::UnitX()) * this_rot;
	}

	void inite_around_Y(double degree)
	{
		this_rot = AngleAxisd(degree_to_rad(degree), Vector3d::UnitY()) * this_rot;
	}

	void inite_around_Z(double degree)
	{
		this_rot = AngleAxisd(degree_to_rad(degree), Vector3d::UnitZ()) * this_rot;
	}

	void inite_Euler(double alpha, double beta, double gamma)
	{
		inite_around_Z(alpha);
		inite_around_X(beta);
		inite_around_Z(gamma);
	}

	void inite_Air(double alpha, double beta, double gamma)
	{
		inite_around_Z(alpha);
		inite_around_X(beta);
		inite_around_Y(gamma);
	}

	void clear()
	{
		this_rot = Matrix3d::Identity();
	}

	Matrix3d build_matrix()
	{
		return this_rot;
	}
};

#endif