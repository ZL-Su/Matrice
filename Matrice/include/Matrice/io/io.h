/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include <iosfwd>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <functional>
#include <type_traits>
#include <experimental/filesystem>
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../util/utils.h"
#include "../util/genalgs.h"
#include "../util/_macro_conditions.h"
#include "../private/_range.h"

DGE_MATRICE_BEGIN
/**
 * \for filesystem namespace in STD library
 */
namespace fs = std::experimental::filesystem;
enum is_skip_folder { N = 0, Y = 1 };

class IO {
public:
	using FileMap = std::map<std::string, std::string>;
	using FileInfo = std::vector<FileMap>;

	/**
	 * \extracts folder(s) for a given absolute path.
	 */
	template<size_t _Nsubfolder = 0> struct Dir_ 
	{
		Dir_() : m_path(fs::current_path().string()) { _Init(); }
		Dir_(const std::string& _path) : m_path(_path) { _Init(); }

		MATRICE_HOST_INL const auto operator()(std::size_t i, std::size_t j) const {
			return std::string(m_path + m_subpaths[i] + m_names[i][j]);
		}
		MATRICE_HOST_INL const auto operator()(std::size_t i) const {
			return std::string(m_path + m_names[0][i]);
		}
		MATRICE_HOST_FINL auto& names(std::size_t i = 0) { return (m_names[i]); }
		MATRICE_HOST_FINL const auto& names(std::size_t i = 0) const { return (m_names[i]); }
		MATRICE_HOST_FINL const auto count(std::size_t _Idx = 0) const { return m_names[_Idx].size(); }

		MATRICE_HOST_FINL auto& path() { return m_path; }
		MATRICE_HOST_FINL const auto& path() const { return m_path; }
	private:
		MATRICE_HOST_INL void _Init() {
			m_subpaths.clear();
			fs::directory_iterator _End; size_t _Cnt = 0;
			for (decltype(_End) _Begin(m_path); _Begin != _End; ++_Begin) {
				if (fs::is_directory(_Begin->status())) {
					auto _Path = _Begin->path().string().substr(m_path.size());
					m_subpaths.emplace_back(_Path);
					++_Cnt;
				}
			}
			if constexpr (_Nsubfolder == 0) m_subpaths.emplace_back("");
			if (_Nsubfolder != 0 && _Nsubfolder > m_subpaths.size())
				std::runtime_error("No enough subfolders");
			
			_Cnt = 0; 
			m_names.resize(_Nsubfolder == 0 ? 1 : min(_Nsubfolder, m_subpaths.size()));
			for_each(m_names, [&](auto& _List) {
				_List = filenames<Y>(m_path + m_subpaths[_Cnt++]);
			});
		}
		std::string m_path;
		std::vector<std::string> m_subpaths;
		std::vector<std::vector<std::string>> m_names;
	};

	/**
	 * \csv writer
	 */
	class CSV {
	public:
		CSV(const std::string& _Fname, const std::string& _Delm = ",")
			:_My_filename(_Fname), _My_delimeter(_Delm){}

		MATRICE_HOST_FINL void open(int _Mode) {
			_My_file.open(_My_filename, _Mode);
		}
		MATRICE_HOST_FINL void close() { 
			_My_file.close(); 
		}

		template<typename _It> 
		MATRICE_HOST_FINL auto append(_It _First, _It _Last) {
			for (; _First != _Last;) {
				_My_file << *_First;
				if (++_First != _Last) _My_file << _My_delimeter;
			}
			_My_file << "\n";
			return (_My_linecnt++);
		}

	private:
		std::string _My_filename, _My_delimeter;
		std::size_t _My_linecnt = 0;
		std::fstream _My_file;
	};

	IO() {};

	static MATRICE_HOST_FINL auto workspace() { return fs::current_path(); }

	template<is_skip_folder _Fopt = is_skip_folder::Y>
	static MATRICE_HOST_INL auto filenames(const std::string& dir)
	{
		std::vector<std::string> _Ret;

		fs::path _Path(dir);
		if (!fs::exists(_Path)) return std::move(_Ret);

		fs::directory_iterator _End;
		for (fs::directory_iterator _Begin(_Path); _Begin != _End; ++_Begin) {
			if (fs::is_regular_file(_Begin->status()))
				_Ret.emplace_back(_Begin->path().string().substr(dir.size()));

			if constexpr (_Fopt == is_skip_folder::N)
				if (fs::is_directory(_Begin->status())) {
					auto _Inner = filenames<_Fopt>(_Begin->path().string());
					_Ret.insert(_Ret.end(), _Inner.begin(), _Inner.end());
				}
		}

		return std::move(_Ret);
	}

	///<summary>
	//@brief: Template function to cvt. any value to std::string 
	//@author: Zhilong Su - Jan.10.2017 @SEU
	///</summary>
	template<typename _Ty> static MATRICE_HOST_FINL auto strf(_Ty _val) {
		std::ostringstream strs; strs << _val; return strs.str();
	}
	// \TEMPLATE function to cvt. numeric value to std::string
	template<typename _Ty> static MATRICE_HOST_FINL auto strfn(_Ty _val) {
#ifdef __CXX11__
		return (std::to_string(_val));
#elif
		return strf(_val);
#endif
	}

	template<typename _Ty, int _M, int _N = 0> 
	static MATRICE_HOST_INL auto read(std::string _path) {
		types::Vec_<_Ty, _M> _Cell;
		std::vector<decltype(_Cell)> data;
		std::ifstream fin(_path); assert(fin);
		if constexpr (_M == 3)
			while (fin >> _Cell(0) >> _Cell(1) >> _Cell(2))
				data.push_back(std::move(_Cell));
		types::Matrix_<_Ty, _N == 0 ? 0 : _M, _N> _Ret(_M, data.size());
		std::size_t i = 0;
		for (const auto& p : data) {
			for (int j = 0; j < _M; ++j)
				_Ret[j][i] = p(j);
			i++;
		}
		fin.close();
		return (_Ret);
	}
	template<typename _Op> 
	static MATRICE_HOST_INL auto read(std::string _path, _Op _op) {
		if (!fs::exists(fs::path(_path))) 
			throw std::runtime_error("'" + _path + "' does not exist!");

		std::ifstream _Fin(_path);
		if (!_Fin.is_open()) 
			throw std::runtime_error("Cannot open file.");
		return _op(_Fin);
	}

	template<typename _Ty, std::size_t _N0, std::size_t _N1 = _N0>
	static MATRICE_HOST_INL auto read(const std::string& _Path, std::size_t _Skips = 0) {
		using value_type = _Ty;

		if (!fs::exists(fs::path(_Path)))
			throw std::runtime_error("'" + _Path + "' does not exist!");
		std::ifstream _Fin(_Path);
		if (!_Fin.is_open())
			throw std::runtime_error("Cannot open file in " + _Path);

		std::string _Line;
		for (const auto _Idx : range(0, _Skips)) std::getline(_Fin, _Line);

		std::array<value_type, _N1 - _N0 + 1> _Myline;
		std::vector<decltype(_Myline)> _Data;
		while (std::getline(_Fin, _Line)) {
			auto _Res = split<value_type>(_Line, ',');
			for (auto _Idx = _N0-1; _Idx < _N1; ++_Idx) {
				_Myline[_Idx - _N0 + 1] = _Res[_Idx];
			}
			_Data.emplace_back(_Myline);
		}

		return std::forward<decltype(_Data)>(_Data);
	}

	//
	template<typename _Op> 
	static MATRICE_HOST_FINL void write(std::string _path, _Op _op) {
		std::ofstream _Fout(_path);
		if (!_Fout.is_open()) 
			throw std::runtime_error("Cannot open file.");

		_op(std::forward<std::ofstream>(_Fout));
		_Fout.close();
	}
	// \multi-path write
	template<typename _Op> static
	MATRICE_HOST_FINL void write(const std::initializer_list<std::string> _paths, _Op _op)
	{
		for_each(_paths, [&_op](const auto& path) {write(path, _op); });
	}

	// \split a string with the given token
	template<typename _Ty = std::string> static 
	MATRICE_HOST_FINL auto split(const std::string& _string, char _token) {
		std::vector<_Ty> _Res;

		const auto _String = std::string(1, _token) + _string;
		auto _Pos = _String.begin(), _End = _String.end();
		if (_String.back() == _token) --_End;

		while (_Pos != _End) {
			auto _Pos_last = _Pos+1;
			_Pos = std::find(_Pos_last, _End, _token);
			_Res.push_back(stonv<_Ty>(_String.substr(std::distance(_String.begin(), _Pos_last), std::distance(_Pos_last, _Pos))));
		}

		return (_Res);
	}

};

// \read interface
template<typename... _Args> HOST_INL_CXPR_T
read(_Args... _args) { try { return IO::read(_args...); } catch (std::exception& _E) { throw _E; } };
// \write interface
template<typename... _Args> HOST_INL_CXPR_T
write(_Args... _args) { try { return IO::write(_args...); } catch (std::exception& _E) { throw _E; } };
// \definite a path under current folder
template<typename T> HOST_INL_CXPR_T
defpath(const T local) { return std::forward<std::string>(IO::workspace().string() + "\\" + IO::strf(local)); };

// \Class: std::string helper  
// \Coded by: dgelom su
class string_helper final {

	using basic_value_type = char;
	using value_type = std::string;
	using const_value_type = std::add_const_t<value_type>;
	using value_reference = std::add_lvalue_reference_t<value_type>;
	using const_value_reference = std::add_const_t<value_reference>;

public:

	// \Given a string and a token, return data items with type of "_Ty"
	template<typename _Ty = value_type> static 
	_INLINE_VAR auto split(const_value_reference _Str, basic_value_type _Token) {
		std::vector<_Ty> _Res;
		const_value_type _String = value_type(1, _Token) + _Str;

		auto _Pos = _String.begin();
		while (_Pos != _String.end()) {
			auto _Pos_last = _Pos + 1;

			_Pos = std::find(_Pos_last, _String.end(), _Token);

			_Res.push_back(dgelom::stonv<_Ty>(_String.substr(
				std::distance(_String.begin(), _Pos_last), 
				std::distance(_Pos_last, _Pos)
			)));
		}
		return std::forward<decltype(_Res)>(_Res);
	}

};

namespace io { _DETAIL_BEGIN

struct folder_tag {};
struct loader_tag {};
struct file_tag {};

/**
 * \file or folder name collector
 */
template<typename _Tag> struct _Collector {};
/**
 * \build file directory according to a given path
 */
template<typename _Tag = folder_tag> class _Dir_impl {};
/**
 * \load data file from given directory
 */
template<typename _Ty, typename _Tag> class _Data_loader_impl{};

_DETAIL_END

} DGE_MATRICE_END
#include "inline\_directory.inl"
#include "inline\_data_loader.inl"

DGE_MATRICE_BEGIN namespace io {
using directory = detail::_Dir_impl<detail::folder_tag>;
using folder_collector = detail::_Collector<detail::folder_tag>;
using file_collector = detail::_Collector<detail::file_tag>;
template<typename _Ty = std::float_t>
using data_loader = detail::_Data_loader_impl<_Ty, detail::loader_tag>;
using default_data_loader = data_loader<std::float_t>;

template<std::size_t _N, typename _Cont>
MATRICE_HOST_FINL auto serial(const _Cont& _L) {
	return tuple_n<_N - 1>::_(_L.data());
}

} DGE_MATRICE_END