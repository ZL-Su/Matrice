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

#include <iostream>
#include <iomanip>
#include <map>
#include <list>
#include <utility>
#include <algorithm>
#include <numeric>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/type_traits.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/function.hpp>

#include "../util/genalgs.h"
#include "../arch/ixpacket.h"
#include "../core/vector.h"

namespace dgelom {
	template<typename T, size_t N>
	using plane_array = std::array<T, N>;
	namespace detail {

		template <size_t N, size_t M>
		struct power : std::integral_constant<size_t, N * power<N, M - 1>::value> {};
		template <size_t N>
		struct power<N, 0> : std::integral_constant<size_t, 1> {};

		/// N-dimensional grid iterator (nested loop with variable depth).
		template <unsigned _Nod> class grid_iterator 
		{
		public:
			typedef plane_array<size_t, _Nod> index;

			explicit grid_iterator(const plane_array<size_t, _Nod> &dims)
				: N(dims), idx(0)
			{
				dgelom::fill(i, 0);
				done = (i == N);
			}

			explicit grid_iterator(size_t dim) : idx(0) {
				dgelom::fill(N, dim);
				dgelom::fill(i, 0);
				done = (0 == dim);
			}

			MATRICE_HOST_FINL size_t operator[](size_t d) const { return i[d]; }
			MATRICE_HOST_FINL const index& operator*() const { return i; }
			MATRICE_HOST_FINL size_t position() const { return idx; }

			grid_iterator& operator++() {
				done = true;
				for (size_t d = _Nod; d--; ) {
					if (++i[d] < N[d]) {
						done = false;
						break;
					}
					i[d] = 0;
				}
				++idx;
				return *this;
			}

			MATRICE_HOST_FINL operator bool() const { return !done; }

		private:
			index N, i;
			bool  done;
			size_t idx;
		};

		template <typename T, size_t N> inline constexpr
		plane_array<T, N> operator+(plane_array<T, N> a, plane_array<T, N> b) {
			static_assert(N <= 4, "Dimensions cannot be greater than 4.");
			(dgelom::simd::Packet_<T, 4>(b.data()) + dgelom::simd::Packet_<T, 4>(a.data())).unpack(a.data());
			return a;
		}

		template <typename T, size_t N, typename C> inline constexpr
		plane_array<T, N> operator-(plane_array<T, N> a, C b) {
			static_assert(N <= 4, "Dimensions cannot be greater than 4.");
			(dgelom::simd::Packet_<T, 4>(a.data()) - dgelom::simd::Packet_<T, 4>(b)).unpack(a.data());
			return a;
		}

		template <typename T, size_t N, typename C> inline constexpr
		plane_array<T, N> operator*(plane_array<T, N> a, C b) {
			static_assert(N <= 4, "Dimensions cannot be greater than 4.");
			(dgelom::simd::Packet_<T, 4>(a.data()) * dgelom::simd::Packet_<T, 4>(b)).unpack(a.data());
			return a;
		}

		// Value of k-th B-Spline basic function at t.
		template<typename _T, typename = std::enable_if_t<std::is_arithmetic_v<_T>>>
		MATRICE_HOST_FINL constexpr _T _Bspline_kernel(size_t k, _T t) {
			assert(0 <= t && t < 1);
			assert(k < 4);

			switch (k) {
			case 0: return (t * (t * (-t + 3) - 3) + 1) / 6;
			case 1: return (t * t * (3 * t - 6) + 4) / 6;
			case 2: return (t * (t * (-3 * t + 3) + 3) + 1) / 6;
			case 3: return t * t * t / 6;
			default: return 0;
			}
		}

		// Checks if p is between lo and hi
		template <typename T, size_t N> MATRICE_HOST_FINL
		bool boxed(const plane_array<T, N> &lo, const plane_array<T, N> &p, const plane_array<T, N> &hi) {
			for (unsigned i = 0; i < N; ++i) {
				if (p[i] < lo[i] || p[i] > hi[i]) return false;
			}
			return true;
		}
		template<typename _T, typename = std::enable_if_t<std::is_arithmetic_v<_T>>>
		MATRICE_HOST_FINL _T safe_divide(_T a, _T b) {
			return b == 0.0 ? 0.0 : a / b;
		}

		template<typename _T, unsigned _Nod> class control_lattice 
		{
		public:
			using value_t = _T;
			typedef plane_array<size_t, _Nod> index;
			typedef plane_array<value_t, _Nod> point;

			virtual ~control_lattice() {}

			virtual value_t operator()(const point &p) const = 0;

			virtual void report(std::ostream& os) const = 0;

			template <class CooIter, class ValIter>
			value_t residual(CooIter coo_begin, CooIter coo_end, ValIter val_begin) const {
				value_t res = 0.0;

				CooIter p = coo_begin;
				ValIter v = val_begin;

				for (; p != coo_end; ++p, ++v) {
					(*v) -= (*this)(*p);
					res = std::max(res, std::abs(*v));
				}

				return res;
			}
		};

		template <typename _T, unsigned _Nod>
		class initial_approximation : public control_lattice<_T, _Nod> {
		public:
			using value_t = _T;
			typedef typename control_lattice<value_t, _Nod>::point point;

			initial_approximation(std::function<value_t(const point&)> f)
				: f(f) {}

			MATRICE_HOST_FINL value_t operator()(const point &p) const { return f(p); }

			void report(std::ostream &os) const {os << "initial approximation";}

		private:
			std::function<value_t(const point&)> f;
		};

		template <typename _T, unsigned _Nod>
		class control_lattice_dense : public control_lattice<_T, _Nod> {
		public:
			using value_t = _T;
			typedef typename control_lattice<value_t, _Nod>::index index;
			typedef typename control_lattice<value_t, _Nod>::point point;

			template <class CooIter, class ValIter>
			control_lattice_dense(
				const point &coo_min, const point &coo_max, index grid_size,
				CooIter coo_begin, CooIter coo_end, ValIter val_begin
			) : cmin(coo_min), cmax(coo_max), grid(grid_size)
			{
				for (unsigned i = 0; i < _Nod; ++i) {
					hinv[i] = (grid[i] - 1) / (cmax[i] - cmin[i]);
					cmin[i] -= 1 / hinv[i];
					grid[i] += 2;
				}

				boost::multi_array<value_t, _Nod> delta(grid);
				boost::multi_array<value_t, _Nod> omega(grid);

				std::fill(delta.data(), delta.data() + delta.num_elements(), 0.0);
				std::fill(omega.data(), omega.data() + omega.num_elements(), 0.0);

				CooIter p = coo_begin;
				ValIter v = val_begin;

				for (; p != coo_end; ++p, ++v) {
					if (!boxed(coo_min, *p, coo_max)) continue;

					index i;
					point s;

					for (unsigned d = 0; d < _Nod; ++d) {
						value_t u = ((*p)[d] - cmin[d]) * hinv[d];
						i[d] = floor(u) - 1;
						s[d] = u - floor(u);
					}

					plane_array< value_t, power<4, _Nod>::value > w;
					value_t sum_w2 = 0.0;

					for (grid_iterator<_Nod> d(4); d; ++d) {
						value_t prod = 1.0;
						for (unsigned k = 0; k < _Nod; ++k)
							prod *= _Bspline_kernel(d[k], s[k]);

						w[d.position()] = prod;
						sum_w2 += prod * prod;
					}

					for (grid_iterator<_Nod> d(4); d; ++d) {
						value_t w1 = w[d.position()];
						value_t w2 = w1 * w1;
						value_t phi = (*v) * w1 / sum_w2;

						index j = i + (*d);

						delta(j) += w2 * phi;
						omega(j) += w2;
					}
				}

				phi.resize(grid);

				std::transform(
					delta.data(), delta.data() + delta.num_elements(),
					omega.data(), phi.data(), safe_divide<value_t>);
			}

			value_t operator()(const point &p) const {
				index i;
				point s;

				for (unsigned d = 0; d < _Nod; ++d) {
					value_t u = (p[d] - cmin[d]) * hinv[d];
					i[d] = floor(u) - 1;
					s[d] = u - floor(u);
				}

				value_t f = 0;

				for (grid_iterator<_Nod> d(4); d; ++d) {
					value_t w = 1.0;
					for (unsigned k = 0; k < _Nod; ++k) w *= _Bspline_kernel(d[k], s[k]);

					f += w * phi(i + (*d));
				}
				return f;
			}

			void report(std::ostream &os) const {
				boost::io::ios_all_saver stream_state(os);

				os << "dense  [" << grid[0];
				for (unsigned i = 1; i < _Nod; ++i)
					os << ", " << grid[i];
				os << "] (" << phi.num_elements() * sizeof(value_t) << " bytes)";
			}

			void append_refined(const control_lattice_dense &r) {
				static const plane_array<value_t, 5> s = {
					0.125, 0.500, 0.750, 0.500, 0.125
				};

				for (grid_iterator<_Nod> i(r.grid); i; ++i) {
					value_t f = r.phi(*i);

					if (f == 0.0) continue;

					for (grid_iterator<_Nod> d(5); d; ++d) {
						index j;
						bool skip = false;
						for (unsigned k = 0; k < _Nod; ++k) {
							j[k] = 2 * i[k] + d[k] - 3;
							if (j[k] >= grid[k]) {
								skip = true;
								break;
							}
						}

						if (skip) continue;

						value_t c = 1.0;
						for (unsigned k = 0; k < _Nod; ++k) c *= s[d[k]];

						phi(j) += f * c;
					}
				}
			}

			value_t fill_ratio() const {
				size_t total = phi.num_elements();
				size_t nonzeros = total - std::count(phi.data(), phi.data() + total, 0.0);

				return static_cast<value_t>(nonzeros) / total;
			}

		private:
			point cmin, cmax, hinv;
			index grid;

			boost::multi_array<value_t, _Nod> phi;

		};

		template <typename _T, unsigned _Nod>
		class control_lattice_sparse : public control_lattice<_T, _Nod> 
		{
		public:
			using value_t = _T;
			typedef typename control_lattice<value_t, _Nod>::index index;
			typedef typename control_lattice<value_t, _Nod>::point point;

			template <class CooIter, class ValIter>
			control_lattice_sparse(
				const point &coo_min, const point &coo_max, index grid_size,
				CooIter coo_begin, CooIter coo_end, ValIter val_begin
			) : cmin(coo_min), cmax(coo_max), grid(grid_size)
			{
				for (unsigned i = 0; i < _Nod; ++i) {
					hinv[i] = (grid[i] - 1) / (cmax[i] - cmin[i]);
					cmin[i] -= 1 / hinv[i];
					grid[i] += 2;
				}

				std::map<index, two_doubles> dw;

				CooIter p = coo_begin;
				ValIter v = val_begin;

				for (; p != coo_end; ++p, ++v) {
					if (!boxed(coo_min, *p, coo_max)) continue;

					index i;
					point s;

					for (unsigned d = 0; d < _Nod; ++d) {
						value_t u = ((*p)[d] - cmin[d]) * hinv[d];
						i[d] = floor(u) - 1;
						s[d] = u - floor(u);
					}

					plane_array< value_t, power<4, _Nod>::value > w;
					value_t sum_w2 = 0.0;

					for (grid_iterator<_Nod> d(4); d; ++d) {
						value_t prod = 1.0;
						for (unsigned k = 0; k < _Nod; ++k) 
							prod *= _Bspline_kernel(d[k], s[k]);

						w[d.position()] = prod;
						sum_w2 += prod * prod;
					}

					for (grid_iterator<_Nod> d(4); d; ++d) {
						value_t w1 = w[d.position()];
						value_t w2 = w1 * w1;
						value_t phi = (*v) * w1 / sum_w2;

						two_doubles delta_omega = { w2 * phi, w2 };

						_Ewise_append(dw[i + (*d)], delta_omega);
					}
				}

				phi.insert(//boost::container::ordered_unique_range,
					boost::make_transform_iterator(dw.begin(), _Delta_over_omega),
					boost::make_transform_iterator(dw.end(), _Delta_over_omega)
				);
			}

			value_t operator()(const point &p) const {
				index i;
				point s;

				for (unsigned d = 0; d < _Nod; ++d) {
					value_t u = (p[d] - cmin[d]) * hinv[d];
					i[d] = floor(u) - 1;
					s[d] = u - floor(u);
				}

				value_t f = 0;

				for (grid_iterator<_Nod> d(4); d; ++d) {
					value_t w = 1.0;
					for (unsigned k = 0; k < _Nod; ++k) 
						w *= _Bspline_kernel(d[k], s[k]);

					f += w * _Get_phi(i + (*d));
				}

				return f;
			}

			void report(std::ostream &os) const {
				boost::io::ios_all_saver stream_state(os);

				size_t grid_size = grid[0];

				os << "sparse [" << grid[0];
				for (unsigned i = 1; i < _Nod; ++i) {
					os << ", " << grid[i];
					grid_size *= grid[i];
				}

				size_t bytes = phi.size() * sizeof(std::pair<index, value_t>);
				size_t dense_bytes = grid_size * sizeof(value_t);

				value_t compression = static_cast<value_t>(bytes) / dense_bytes;
				os << "] (" << bytes << " bytes, compression: "
					<< std::fixed << std::setprecision(2) << compression << ")";
			}
		private:
			point cmin, cmax, hinv;
			index grid;

			typedef boost::container::flat_map<index, value_t> sparse_grid;
			sparse_grid phi;

			typedef plane_array<value_t, 2> two_doubles;

			static std::pair<index, value_t> _Delta_over_omega(const std::pair<index, two_doubles> &dw) {
				return std::make_pair(dw.first, safe_divide(dw.second[0], dw.second[1]));
			}

			static void _Ewise_append(two_doubles &a, const two_doubles &b) {
				boost::transform(a, b, boost::begin(a), std::plus<value_t>());
			}

			value_t _Get_phi(const index &i) const {
				typename sparse_grid::const_iterator c = phi.find(i);
				return c == phi.end() ? 0.0 : c->second;
			}
		};

	} // namespace detail

	template <typename _T, size_t _Nod> class linear_approximation 
	{
	public:
		using value_t = _T;
		typedef typename detail::control_lattice<value_t, _Nod>::point point;

		template <class CooIter, class ValIter> inline
		linear_approximation(CooIter coo_begin, CooIter coo_end, ValIter val_begin)
		{
			namespace ublas = boost::numeric::ublas;

			size_t n = std::distance(coo_begin, coo_end);

			if (n <= _Nod) {
				// Not enough points to get a unique plane
				boost::fill(C, 0.0);
				C[_Nod] = std::accumulate(val_begin, val_begin + n, 0.0) / n;
				return;
			}

			ublas::matrix<value_t> A(_Nod + 1, _Nod + 1); A.clear();
			ublas::vector<value_t> f(_Nod + 1);           f.clear();

			CooIter p = coo_begin;
			ValIter v = val_begin;

			value_t sum_val = 0.0;

			// Solve least-squares problem to get approximation with a plane.
			for (; p != coo_end; ++p, ++v, ++n) {
				plane_array<value_t, _Nod + 1> x;
				boost::copy(*p, boost::begin(x));
				x[_Nod] = 1.0;

				for (unsigned i = 0; i <= _Nod; ++i) {
					for (unsigned j = 0; j <= _Nod; ++j) {
						A(i, j) += x[i] * x[j];
					}
					f(i) += x[i] * (*v);
				}

				sum_val += (*v);
			}

			ublas::permutation_matrix<size_t> pm(_Nod + 1);
			ublas::lu_factorize(A, pm);

			bool singular = false;
			for (unsigned i = 0; i <= _Nod; ++i) {
				if (A(i, i) == 0.0) {
					singular = true;
					break;
				}
			}

			if (singular) {
				boost::fill(C, 0.0);
				C[_Nod] = sum_val / n;
			}
			else {
				ublas::lu_substitute(A, pm, f);
				for (unsigned i = 0; i <= _Nod; ++i) C[i] = f(i);
			}
		}

		inline value_t operator()(const point &p) const {
			value_t f = C[_Nod];

			for (unsigned i = 0; i < _Nod; ++i)
				f += C[i] * p[i];

			return f;
		}
	private:
		plane_array<value_t, _Nod + 1> C;
	};

	template <typename _T, size_t _Nod, size_t _Max_levels = 8> class MBA {
	public:
		using value_t = _T;
		typedef plane_array<size_t, _Nod> index;
		typedef plane_array<value_t, _Nod> point;
		using point_t = types::Vec_<value_t, _Nod>;
		struct configration {
			value_t tolerance = std::numeric_limits<value_t>::epsilon();
			value_t min_fill = 0.5;
			types::Vec4_<value_t> range; //{min_x, min_y, max_x, max_y}
			plane_array<size_t, _Nod> grid = { 3, 3 }; //
			
		};
		template <class CooIter, class ValIter>
		MBA(
			const point &coo_min, const point &coo_max, index grid,
			CooIter coo_begin, CooIter coo_end, ValIter val_begin,
			unsigned max_levels = 8, value_t tol = 1e-8, value_t min_fill = 0.5,
			std::function<value_t(point)> initial = std::function<value_t(point)>()
		)
		{
			_Init(
				coo_min, coo_max, grid,
				coo_begin, coo_end, val_begin,
				max_levels, tol, min_fill, initial
			);
		}

		template <class CooRange, class ValRange>
		MBA(
			const point &coo_min, const point &coo_max, index grid,
			CooRange coo, ValRange val,
			unsigned max_levels = 8, value_t tol = 1e-8, value_t min_fill = 0.5,
			std::function<value_t(point)> initial = std::function<value_t(point)>()
		)
		{
			_Init(
				coo_min, coo_max, grid,
				boost::begin(coo), boost::end(coo), boost::begin(val),
				max_levels, tol, min_fill, initial
			);
		}

		template <class CooRange, class ValRange>
		MBA(const configration& config, CooRange coo, ValRange val, std::function<value_t(point)> initial = std::function<value_t(point)>())
		{
			_Init(point{config.range.x, config.range.y}, 
				   point{config.range.z, config.range.w}, config.grid,
				   boost::begin(coo), boost::end(coo), boost::begin(val),
				   _Max_levels, config.tolerance, config.min_fill, initial);
		}
		// \param: p - the interpolation point
		inline value_t operator()(const types::Vec_<value_t, _Nod> &p) const 
		{
			value_t f = 0.0; point _P = p;
			dgelom::for_each(cl, [&](const auto& psi) {f += (*psi)(_P); });
			return f;
		}

		friend std::ostream& operator<<(std::ostream &os, const MBA &h) {
			size_t level = 0;
			dgelom::for_each(h.cl, [&](const auto& psi) {
				os << "level " << ++level << ": ";
				psi->report(os);
				os << std::endl;
			});
			return os;
		}

	private:
		typedef detail::control_lattice<value_t, _Nod>        lattice;
		typedef detail::initial_approximation<value_t, _Nod>  initial_approximation;
		typedef detail::control_lattice_dense<value_t, _Nod>  dense_lattice;
		typedef detail::control_lattice_sparse<value_t, _Nod> sparse_lattice;


		std::list<std::shared_ptr<lattice>> cl;

		template <class CooIter, class ValIter>
		void _Init(
			const point &cmin, const point &cmax, index grid,
			CooIter coo_begin, CooIter coo_end, ValIter val_begin,
			unsigned max_levels, value_t tol, value_t min_fill,
			std::function<value_t(point)> initial
		)
		{
			using namespace dgelom::detail;

			const ptrdiff_t n = std::distance(coo_begin, coo_end);
			std::vector<value_t> val(val_begin, val_begin + n);

			value_t res = value_t(0);
			auto minmax = std::minmax_element(val.begin(), val.end());
			auto eps = std::max(std::abs(*minmax.first), std::abs(*minmax.second));
			eps *= tol;

			if (initial) {
				// Start with the given approximation.
				cl.push_back(std::make_shared<initial_approximation>(initial));
				res = cl.back()->residual(coo_begin, coo_end, val.begin());
				if (res <= eps) return;
			}

			size_t lev = 1;
			// Create dense head of the hierarchy.
			{
				std::shared_ptr<dense_lattice> psi = std::make_shared<dense_lattice>(
					cmin, cmax, grid, coo_begin, coo_end, val.begin());

				res = psi->residual(coo_begin, coo_end, val.begin());
				value_t fill = psi->fill_ratio();

				for (; (lev < max_levels) && (res > eps) && (fill > min_fill); ++lev) {
					grid = grid * 2ul - 1ul;

					std::shared_ptr<dense_lattice> f = std::make_shared<dense_lattice>(
						cmin, cmax, grid, coo_begin, coo_end, val.begin());

					res = f->residual(coo_begin, coo_end, val.begin());
					fill = f->fill_ratio();

					f->append_refined(*psi);
					psi.swap(f);
				}

				cl.push_back(psi);
			}

			// Create sparse tail of the hierrchy.
			for (; (lev < max_levels) && (res > eps); ++lev) {
				grid = grid * 2ul - 1ul;

				cl.push_back(std::make_shared<sparse_lattice>(
					cmin, cmax, grid, coo_begin, coo_end, val.begin()));

				res = cl.back()->residual(coo_begin, coo_end, val.begin());
			}
		}
	};

} // namespace mba