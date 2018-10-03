#pragma once
#include <iosfwd>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <functional>
#include <type_traits>
#include <experimental\filesystem>
#include "../core/matrix.h"
#include "../core/vector.h"
namespace dgelom {
namespace fs = std::experimental::filesystem;
enum is_skip_folder { N = 0, Y = 1 };
class IO 
{
public:
	using FileMap = std::map<std::string, std::string>;
	using FileInfo = std::vector<FileMap>;
	template<size_t _Nsubfolder = 0> struct Dir_ 
	{
		Dir_() : m_path(fs::current_path().string()) { _Init(); }
		Dir_(const std::string& _path) : m_path(_path) { _Init(); }

		MATRICE_HOST_INL const auto operator()(size_t i, size_t j) const {
			return std::string(m_path + m_subpaths[i] + m_names[i][j]);
		}
		MATRICE_HOST_INL const auto operator()(size_t i) const {
			return std::string(m_path + m_names[0][i]);
		}
		MATRICE_HOST_FINL auto& names(size_t i = 0) { return (m_names[i]); }
		MATRICE_HOST_FINL const auto& names(size_t i = 0) const { return (m_names[i]); }
		MATRICE_HOST_FINL const auto count() const { return m_names.front().size(); }
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
	template<typename _Ty>
	static MATRICE_HOST_FINL std::string strf(_Ty _val)
	{
		std::ostringstream strs; strs << _val; return strs.str();
	}
	// \TEMPLATE function to cvt. numeric value to std::string
	template<typename _Ty>
	static MATRICE_HOST_FINL std::string strfn(_Ty _val)
	{
#ifdef __CXX11__
		return (std::to_string(_val));
#elif
		return strf(_val);
#endif
	}

	template<typename _Ty, int _M, int _N = 0> static 
	MATRICE_HOST_INL auto read(std::string _path)
	{
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
	template<typename _Op> static
	MATRICE_HOST_INL auto read(std::string _path, _Op _op)
	{
		if (!fs::exists(fs::path(_path))) 
			throw std::runtime_error("'" + _path + "' does not exist!");

		std::ifstream _Fin(_path);
		if (!_Fin.is_open()) 
			throw std::runtime_error("Cannot open file.");
		return _op(_Fin);
	}
	//
	template<typename _Op> static 
	MATRICE_HOST_FINL void write(std::string _path, _Op _op)
	{
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
		auto _Pos = _String.begin();
		while (_Pos != _String.end()) {
			auto _Pos_last = _Pos+1;
			_Pos = std::find(_Pos_last, _String.end(), _token);
			_Res.push_back(stonv<_Ty>(_String.substr(std::distance(_String.begin(), _Pos_last), std::distance(_Pos_last, _Pos))));
		}

		return (_Res);
	}

};

// \read interface
template<typename... _Args> MATRICE_HOST_ICEA
read(_Args... _args) { try { return IO::read(_args...); } catch (std::exception& _E) { throw _E; } };
// \write interface
template<typename... _Args> MATRICE_HOST_ICEA
write(_Args... _args) { try { return IO::write(_args...); } catch (std::exception& _E) { throw _E; } };
// \definite a path under current folder
template<typename T> MATRICE_HOST_ICEA 
defpath(const T local) { return std::forward<std::string>(IO::workspace().string() + "\\" + IO::strf(local)); };
}

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