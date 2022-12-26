//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#ifndef SMARTMPW_CHECKER_VISUALIZER_HPP
#define SMARTMPW_CHECKER_VISUALIZER_HPP

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/geometries.hpp>

using coord_t = double;

/// �������ݽṹ����
namespace utils_visualize {

	namespace bg = boost::geometry;

	using T = coord_t;

	// n ά������
	template<size_t dimension = 2>
	using bg_point_base = bg::model::point<T, dimension, bg::cs::cartesian>;

	/*****************************
	* ���¶���ȫ����Զ�ά������ *
	******************************/

	// ��ά�����
	using bg_point_t = bg::model::d2::point_xy<T>;
	const bg_point_t origin_point(0, 0);  // ����ԭ��

	// ����
	using bg_linestring_t = bg::model::linestring<bg_point_t>;

	// ����Σ�����0�����߶���ڻ� inner rings����ʱ�룬���=�յ㣩
	using bg_polygon_t = bg::model::polygon<bg_point_t, false, true>;

	// ����������ߣ�˳ʱ�룬�����յ㣩
	using bg_ring_t = bg::model::ring<bg_point_t, true, false>;
	//using bg_ring_t = bg_polygon_t::ring_type;

	// �㼯��
	using bg_multi_point_t = bg::model::multi_point<bg_point_t>;

	// ���߼���
	using bg_multi_linestring_t = bg::model::multi_linestring<bg_linestring_t>;

	// ����μ���
	using bg_multi_polygon_t = bg::model::multi_polygon<bg_polygon_t>;

	// ����
	using bg_box_t = bg::model::box<bg_point_t>;

	// �߶Σ������ԣ�
	using bg_segment_t = bg::model::segment<bg_point_t>;
}

#endif // SMARTMPW_CHECKER_VISUALIZER_HPP
