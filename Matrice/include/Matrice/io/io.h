#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <functional>
#include "../core/matrix.h"
#include "../core/vector.h"

namespace dgelom {
class IO 
{
public:
	IO() {};

	///<summary>
	//@brief: Template function to cvt. number to std::string
	//@author: Zhilong Su - Jan.10.2017 @SEU
	///</summary>
	template<typename _Ty> static inline  std::string strf(_Ty _val)
	{
		std::ostringstream strs; strs << _val; return strs.str();
	}

	template<typename _Ty, int _M, int _N = 0> static inline 
	types::Matrix_<_Ty, (_N == 0 ? 0 : _M), _N> read(const std::string& _path)
	{
		types::Vec_<_Ty, _M> cell;
		std::vector<types::Vec_<_Ty, _M>> data;
		std::ifstream fin(_path); assert(fin);
		if constexpr (_M == 3)
			while (fin >> cell(0) >> cell(1) >> cell(2))
				data.push_back(std::move(cell));
		types::Matrix_<_Ty, _N == 0 ? 0 : _M, _N> ret(_M, data.size());
		std::size_t i = 0;
		for (const auto& p : data) {
			for (int j = 0; j < _M; ++j)
				ret[j][i] = p(j);
			i++;
		}
		fin.close();
		return (ret);
	}
	template<typename _Op> static
	MATRICE_HOST_FINL auto read(const std::string& _path, _Op _op)
	{
		std::ifstream _Fin(_path);
		if (!_Fin.is_open()) std::runtime_error("Cannot open file.");
		return _op(std::forward<std::ifstream>(_Fin));
	}
	template<typename _Op> static 
	MATRICE_HOST_FINL void write(const std::string& _path, _Op _op)
	{
		std::ofstream _Fout(_path);
		if (!_Fout.is_open()) std::runtime_error("Cannot open file.");
		_op(std::forward<std::ofstream>(_Fout));
		_Fout.close();
	}
};
}

