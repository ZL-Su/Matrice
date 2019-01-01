/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
**************************************************************************/
#pragma once
#include<experimental/filesystem>
#include <unordered_map>
#include <vector>
#include "../util/_macros.h"

DGE_MATRICE_BEGIN namespace io {
/**
 * \for filesystem namespace in STD library
 */
namespace fs = std::experimental::filesystem;

_DETAIL_BEGIN

struct folder_tag {};
struct loader_tag {};

template<typename _Tag = folder_tag> class _Dir_impl {};

template<> class _Dir_impl<folder_tag> {
public:
	using category = folder_tag;
	using value_type = std::string;
	using path_type = fs::path;
	using container = std::vector<value_type>;
	/**
	 * \_Root must has form "./folder_name" or "/folder_name" 
	 */
	_Dir_impl(value_type _Root, path_type&& _Workdir = fs::current_path())
		:_Mypath(_Workdir) {
		if (_Root[0] == '.') _Mypath.concat(_Root.begin() + 1, _Root.end());
		else _Mypath.concat(_Root.begin(), _Root.end());

		_Ext_subfolder();
	}

	/**
	 * \return full path of _Idx-th folder
	 */
	MATRICE_HOST_INL auto operator[] (std::size_t _Idx) const {
		return (_Mypath.string() + folder(_Idx));
	}

	/**
	 * \return _Idx-th folder name under work-path _Mypath
	 */
	MATRICE_HOST_INL value_type folder(std::size_t _Idx) const {
#ifdef _DEBUG
		if (_Idx >= _Mysubfolders.size()) throw
			std::exception("_Idx over range of _Mysubfolders");
#endif
		return (_Mysubfolders[_Idx]);
	}

	/**
	 * \return size of current directory
	 */
	MATRICE_HOST_INL const auto size() const {
		return (_Mysubfolders.size());
	}

private:
	MATRICE_HOST_INL void _Ext_subfolder() {
		_Mysubfolders.clear();
		auto _Path = _Mypath.string();
		fs::directory_iterator _End;
		for (decltype(_End) _Begin(_Mypath); _Begin != _End; ++_Begin) {
			if (fs::is_directory(_Begin->status())) {
				auto _Sub = _Begin->path().string().substr(_Path.size());
				_Mysubfolders.emplace_back(_Sub);
			}
		}
		_Mysubfolders.shrink_to_fit();
	}

	path_type _Mypath;
	container _Mysubfolders;
};

template<typename _Tag = loader_tag> class _Data_loader_impl{};

template<> class _Data_loader_impl<loader_tag> {
	using _Mydir_type = _Dir_impl<folder_tag>;
	using _Myt = _Data_loader_impl;
	struct _Loader_iterator {
	public:
		_Loader_iterator(const std::add_pointer_t<_Myt> _This)
			:_Myptr(_This){
		}


	private:
		std::size_t _Mypos = 0;
		std::add_pointer_t<_Myt> _Myptr = nullptr;
	};
public:
	using category = loader_tag;
	using dir_type = _Mydir_type;
	using iterator = _Loader_iterator;

	_Data_loader_impl(const _Mydir_type& _Dir) 
		: _Mydir(_Dir) {
		_Mynames.resize(_Mydir.size());
	}

	MATRICE_HOST_INL auto begin() const {

	}
	MATRICE_HOST_INL auto end() const {

	}
private:
	_Mydir_type _Mydir;
	std::vector<_Mydir_type::container> _Mynames;
};

_DETAIL_END

using directory = detail::_Dir_impl<detail::folder_tag>;
using data_loader = detail::_Data_loader_impl<detail::loader_tag>;

} DGE_MATRICE_END