/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#pragma once
#include <io/io.hpp>
#include <core/matrix.h>

DGE_MATRICE_BEGIN
namespace example {
/// <summary>
/// \brief Example to show the IO use of CSV files with Matrice IO system.
/// </summary>
/// <param name="apath">Path to the folder that holds CSV files.</param>
void read_and_write_csv_eng(fs::path&& apath) try
{
	// \c Read data with a float type...
	auto raw_shape = dgelom::IO::read<float>(dgelom::concat(apath, "sphere_air.txt"));  // from file pcloud_0mm.txt
	auto true_shape = dgelom::IO::read<float>(dgelom::concat(apath, "sphere_our.txt")); // from file pcloud_1mm.txt

	auto diff = (raw_shape - true_shape).eval();

	{
		dgelom::IO::CSV csv(dgelom::concat(apath, "sphere_dx.csv"));
		csv.open(dgelom::IO::out | dgelom::IO::app);
		csv.append("GL_POINTS");
		for (auto i = 0; i < diff.rows(); ++i) {
			dgelom::auto_vector_t<float_t, 3> tmp;
			tmp.x = raw_shape[i][0];
			tmp.y = raw_shape[i][1];
			tmp.z = diff[i][0];
			csv.append(tmp.begin(), tmp.end());
		}
		csv.close();

		csv.reset(dgelom::concat(apath, "sphere_dy.csv"));
		csv.open(dgelom::IO::out | dgelom::IO::app);
		csv.append("GL_POINTS");
		for (auto i = 0; i < diff.rows(); ++i) {
			dgelom::auto_vector_t<float_t, 3> tmp;
			tmp.x = raw_shape[i][0];
			tmp.y = raw_shape[i][1];
			tmp.z = diff[i][1];
			csv.append(tmp.begin(), tmp.end());
		}
		csv.close();

		csv.reset(dgelom::concat(apath, "sphere_dz.csv"));
		csv.open(dgelom::IO::out | dgelom::IO::app);
		csv.append("GL_POINTS");
		for (auto i = 0; i < diff.rows(); ++i) {
			dgelom::auto_vector_t<float_t, 3> tmp;
			tmp.x = raw_shape[i][0];
			tmp.y = raw_shape[i][1];
			tmp.z = (diff[i][2]+0.1)/5;
			csv.append(tmp.begin(), tmp.end());
		}
		csv.close();
	}

	// \c Do some processing...
	std::vector<dgelom::Vec4_<float>> raw_dmap, true_dmap;
	for (auto row = raw_shape.rwbegin(); row != raw_shape.rwend(); ++row) {
		const auto x = row[0], y = row[1], z = row[2];
		const auto dz = std::sqrt(400 - dgelom::sqsum(x, y)) - std::abs(z)+1.5f;
		raw_dmap.push_back({ x, y, std::abs(z), dz });
	}
	for (auto row = true_shape.rwbegin(); row != true_shape.rwend(); ++row) {
		const auto x = row[0], y = row[1], z = row[2];
		const auto dz = std::sqrt(400 - dgelom::sqsum(x, y)) - std::abs(z)+1.5f;
		true_dmap.push_back({ x, y, std::abs(z), dz });
	}

	// \c Write data to files...
	dgelom::IO::CSV csv_raw(apath.append("pcloud_0mm.csv").string());  // to file pcloud_0mm.csv
	csv_raw.open(dgelom::IO::out | dgelom::IO::app);
	for (auto row = raw_shape.rwbegin(); row != raw_shape.rwend(); ++row) {
		csv_raw.append(row.begin(), row.end());
	}
	csv_raw.close();

	dgelom::IO::CSV csv_true(apath.append("pcloud_1mm.csv").string()); // to file pcloud_1mm.csv
	csv_true.open(dgelom::IO::out | dgelom::IO::app);
	for (auto row = true_shape.rwbegin(); row != true_shape.rwend(); ++row) {
		csv_true.append(row.begin(), row.end());
	}
	csv_true.close();
}
catch (dgelom::exception::error& e) {
	std::cout << "Err: " << e.what()
		<< "in function: " << e.location()._func
		<< ". (See line " << e.location()._line
		<< " in file '" << e.location()._file << "')\n";
}
catch (std::exception& e) {
	std::cout << "Unknown exception with info: [" << e.what() << "]" << std::endl;
}
catch (...) {
	std::cout << "Unknown exception." << std::endl;
}
}
DGE_MATRICE_END