/*  *************************************************************************
	This file is an Open Source C++ High Resolution Clock.
	Copyright(C) 2017, Zhilong Su, all rights reserved.

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
*	*************************************************************************/
#pragma once

#include <iostream>
#include <chrono>
#include <iomanip>
#include <string>

/* 
 * @Name: _HRC -- High Resolution Clock
 * @Calling pipeline: _HRC.start -> "target block" -> _HRC_stop -> elapsed_time()
 * @Copyright(c): Zhilong Su (su-zl@seu.edu.cn) 2017
 */
template<typename _Clock = std::chrono::high_resolution_clock> class HRC_ final
{
	using time_point_type =  std::chrono::time_point<_Clock>;
public:
	HRC_() { }
	~HRC_() { }

public:
	// \start-timer of high resolution clock
	__declspec(property(get = _prop_time_getter)) time_point_type start;
	// \stop-timer of high resolution clock
	__declspec(property(get = _prop_time_getter)) time_point_type stop;
	constexpr time_point_type _prop_time_getter() { return (m_start = _Clock::now()); }

	// \return elapsed time
	inline auto elapsed_time(const time_point_type& _start)
	{
		using namespace std;
		const auto now = _Clock::now();
		auto interval = chrono::duration_cast<chrono::nanoseconds>(
			now - _start).count();
		m_time = static_cast<double>(interval) / m_ns2msconst;
		return m_time;
	}

	// \output elasped time of runing target represented by _txt
	inline auto elapsed_time(const time_point_type& _start, const std::string& _txt)
	{
		using namespace std;
		const auto now = _Clock::now();
		auto interval = chrono::duration_cast<chrono::nanoseconds>(
			now - _start).count();
		m_time = static_cast<double>(interval) / m_ns2msconst;
		cout <<" >> [Timer Message] Elapsed time of " 
			<< _txt << setprecision(9)
			<< " " << m_time << "ms" << endl;
		return m_time;
	}

	// \output elasped time of runing target represented by _txt
	inline auto elapsed_time(const std::string& _txt)
	{
		using namespace std;
		auto interval = chrono::duration_cast<chrono::nanoseconds>(_Clock::now() - m_start).count();
		cout << " >> [Timer Message] Elapsed time of "
			<< _txt << setprecision(9)
			<< " " << (m_time = static_cast<double>(interval) / m_ns2msconst) << "ms" << endl;
		return m_time;
	}

private:
	double m_time;
	time_point_type m_start;
	const double m_ns2msconst = 1.0e6;
};

namespace dgelom {
	using delicate_clock_t = HRC_< std::chrono::high_resolution_clock>;
}