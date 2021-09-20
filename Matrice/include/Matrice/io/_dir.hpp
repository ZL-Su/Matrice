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
#include "util/_macros.h"
#if _MSC_VER >= 1926
#include <filesystem>
#else
#include <experimental/filesystem>
#endif

DGE_MATRICE_BEGIN
namespace io {

struct leaf_tag{};

/// <summary>
/// \brief CLASS TEMPLATE "Dir" defines a file directory instance.
/// </summary>
/// <typeparam name="_Tag"></typeparam>
template<typename _Tag = leaf_tag>
class Dir {
public:
	using path_t = std::filesystem::path;
	using entry_t = std::filesystem::directory_entry;

	Dir(path_t&& _Pval) noexcept 
		: _Myval{ _Seek_entry(std::move(_Pval)) } {
	}
	Dir(const path_t& _Pval) noexcept 
		: Dir{ path_t{_Pval} } {
	}

	/// <summary>
	/// \brief Entry of the Dir.
	/// </summary>
	/// <returns>reference of the dir::entry</returns>
	MATRICE_HOST_INL decltype(auto) entry() const noexcept {
		return (_Myval);
	}
	MATRICE_HOST_INL decltype(auto) entry() noexcept {
		return (_Myval);
	}

	MATRICE_HOST_INL auto begin() const noexcept {
		return std::filesystem::directory_iterator{ _Myval };
	}
	MATRICE_HOST_INL auto begin() noexcept {
		return std::filesystem::directory_iterator{ _Myval };
	}
	MATRICE_HOST_INL auto end() const noexcept {
		return std::filesystem::directory_iterator{};
	}
	MATRICE_HOST_INL auto end() noexcept {
		return std::filesystem::directory_iterator{};
	}

private:
	auto _Seek_entry(path_t&& _Pval) const noexcept {
		const auto _Txt = _Pval.string();
		const auto _Ret = std::find(_Txt.begin(), _Txt.end(), ':');
		auto _Path = path_t{};
		if (_Ret != _Txt.end()) {
			_Path = _Pval;
		}
		else {
			_Path = std::filesystem::current_path();
			_Path += (_Pval);
		}
		return std::forward<entry_t>(entry_t{ _Path });
	}

	entry_t _Myval;
};

MATRICE_HOST_FINL Dir<> make_dir(Dir<>::path_t&& _Pval) noexcept {
	return { std::forward<decltype(_Pval)>(_Pval) };
}
MATRICE_HOST_FINL Dir<> make_dir(const Dir<>::path_t& _Pval) noexcept {
	return { _Pval };
}
}
DGE_MATRICE_END