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
		return (ret);
	}

	static inline bool write(std::function<void(void)> func) { func; }
};
}

