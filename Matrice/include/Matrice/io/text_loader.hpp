/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/
#pragma once

#if _MSC_VER >= 1926
#include <filesystem>
#else
#include <experimental/filesystem>
#endif
#include <array>
#include <vector>
#include <algorithm>

#include "_dir.hpp"

DGE_MATRICE_BEGIN
namespace io {
template<typename _Ty>
class TextLoader {
public:
	using dir_t = Dir<leaf_tag>;

	MATRICE_HOST_INL auto num_files()const noexcept {
		return _Mynum;
	}

protected:
	explicit TextLoader(dir_t dir) noexcept
		:_Mydir{ dir } {
		_Mynum = _Count_files();
	}
	explicit TextLoader(dir_t::path_t path) noexcept
		:TextLoader(dir_t(path)) {
	}

private:
	MATRICE_HOST_INL size_t _Count_files() {
		return std::count_if(_Mydir.begin(), _Mydir.end(),
			[](const auto& _Ent) {
				return std::filesystem::is_regular_file(_Ent);
			}
			);
	}

	dir_t _Mydir;
	size_t _Mynum;
};

template<typename _Ty>
class csv_loader : public TextLoader<_Ty> {
	using _Mybase = TextLoader<_Ty>;
public:
	using typename _Mybase::dir_t;
	csv_loader(const dir_t& dir) : _Mybase(dir) {}
	csv_loader(const char* dir) : _Mybase(dir) {}

	MATRICE_HOST_INL auto operator()() const {

	}
};

}
DGE_MATRICE_END