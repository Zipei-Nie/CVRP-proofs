#pragma once

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

typedef kv::interval<double> itv;
struct Area;
struct IntegralDistance;
template<typename F> itv RightTriangle(const itv& x, const itv& y);
template<typename F> itv DiskSector(const itv& theta);
template <typename F> itv DiskSegment(const itv& h);
template <typename F> itv DiskDoubleSegment(const itv& h1, const itv& h2);
itv AverageDistance(const itv& x, const itv& y);
template <typename F> itv OverlapRegion(const itv& x, const itv& y, const itv& R);
