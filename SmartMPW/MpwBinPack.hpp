//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
//todo:见307
#ifndef SMARTMPW_MPWBINPACK_HPP
#define SMARTMPW_MPWBINPACK_HPP

#include <list>
#include <string>
#include <unordered_set>
#include <numeric>
#include <algorithm>

#include "Data.hpp"

namespace mbp {

	using namespace std;

	class MpwBinPack {

		/// 排序规则定义
		struct SortRule {
			vector<size_t> sequence;
			coord_t target_area;

			string tabu_key_str() const {
				string key = to_string(sequence.front());
				for (size_t i = 1; i < sequence.size(); ++i) { key += "," + to_string(sequence[i]); }
				return key;
			}

			size_t tabu_key_hash() const { return hash<string>{}(tabu_key_str()); }
		};

		struct intermediate_solu
		{
			skyline_t skyline;    //中间解当前形成的skyline
			vector<polygon_ptr> _dst;      //当前做出的_dst
			double area_percentage;        //已放置的图形的总面积占包络矩形面积的百分比
			list<size_t> polygons;        //剩余要放置的形状
			size_t current_height;           //当前达到的最大高度
			size_t future_height;         //按多边形输入顺序把剩余的形状放置完毕后能达到的最大高度
			bool operator < (const intermediate_solu& x) const  //current_height越小越好，以此为依据时把这个小的排在前面
			{
				return current_height < x.current_height;
			}

		};

		using TabuTable = unordered_set<string>;
		using TabuFunc = string(SortRule::*)()const;
		TabuFunc tabu_key = &SortRule::tabu_key_str;

		//using TabuTable = unordered_set<size_t>;
		//using TabuFunc = size_t(SortRule::*)()const;
		//TabuFunc tabu_key = &SortRule::tabu_key_hash;

	public:

		MpwBinPack() = delete;

		MpwBinPack(const vector<polygon_ptr> &src, coord_t width, coord_t height, default_random_engine &gen) :
			_src(src), _bin_width(width), _bin_height(height), _obj_area(numeric_limits<coord_t>::max()),
			_gen(gen), _uniform_dist(0, _src.size() - 1) {
			reset();
			init_sort_rules();
		}

		const vector<polygon_ptr> &get_dst() const { return _dst; }

		coord_t get_obj_area() const { return _obj_area; }

		void set_obj_area(coord_t area) { _obj_area = area; }

		void set_bin_height(coord_t height) { _bin_height = height; } // 上界

		coord_t get_skyline_height() const { // 排版后的上边界
			return max_element(_skyline.begin(), _skyline.end(),
				[](const skylinenode_t &lhs, const skylinenode_t &rhs) { return lhs.y < rhs.y; })->y;
		}

		/// 基于bin_width进行RLS
		void random_local_search(int iter) {
			// the first time to call RLS on W_k
			if (iter == 1) {
				for (auto &rule : _sort_rules) {
					_polygons.assign(rule.sequence.begin(), rule.sequence.end());
					vector<polygon_ptr> target_dst;
					bool first_insert = insert_bottom_left_score(target_dst);
					assert(first_insert); // 第一次必能全部放下
					rule.target_area = _bin_width * get_skyline_height();
					if (rule.target_area < _obj_area) {
						_obj_area = rule.target_area;
						_dst = target_dst;
						
					}
				}
				// 降序排列，越后面的目标函数值越小选中概率越大
				sort(_sort_rules.begin(), _sort_rules.end(), [](const SortRule &lhs, const SortRule &rhs) {
					return lhs.target_area > rhs.target_area; });
			}
			// 迭代优化
			SortRule &picked_rule = _sort_rules[_discrete_dist(_gen)];
			for (int i = 1; i <= iter; ++i) {
				SortRule new_rule = picked_rule;
				if (iter % 4) { swap_sort_rule(new_rule); }
				else { rotate_sort_rule(new_rule); }

				//if (_tabu_table.count((new_rule.*tabu_key)())) { continue; } // 已禁忌
				//_tabu_table.insert((new_rule.*tabu_key)());

				_polygons.assign(new_rule.sequence.begin(), new_rule.sequence.end());
				vector<polygon_ptr> target_dst;
				if (!insert_bottom_left_score(target_dst)) { continue; } // 放不下
				coord_t target_height = get_skyline_height();
				new_rule.target_area = _bin_width * target_height;
				if (new_rule.target_area < picked_rule.target_area) {
					picked_rule = new_rule;
					if (picked_rule.target_area < _obj_area) {
						_obj_area = picked_rule.target_area;
						cout << _obj_area << endl;
						_dst = target_dst;
						set_bin_height(target_height);
					}
				}
			}
			// 更新排序规则列表
			sort(_sort_rules.begin(), _sort_rules.end(), [](const SortRule &lhs, const SortRule &rhs) {
				return lhs.target_area > rhs.target_area; });
		}
		static bool compare(shared_ptr<intermediate_solu> a, shared_ptr<intermediate_solu> b)  //area_percentage越大越好，以此为依据排序时保证这个大的solution排在前面
		{
			return a->area_percentage > b->area_percentage;
		}
		static bool compare_byfuture(shared_ptr<intermediate_solu> a, shared_ptr<intermediate_solu> b) 
		//future_height越小越好，以此为依据排序时保证这个小的solution排在前面		
		{
			return a->future_height < b->future_height;
		}
		void insert_bottom_left_score_beam_search(shared_ptr<intermediate_solu> child_solu)
		{
			list<size_t> tmp_polygons = child_solu->polygons;
			_skyline = child_solu->skyline;


			while (!tmp_polygons.empty()) {
				//找到最低的skyline，iter是指向这个skyline的指针，index是在vector中的序号
				auto bottom_skyline_iter = min_element(_skyline.begin(), _skyline.end(), [](skylinenode_t& lhs, skylinenode_t& rhs) { return lhs.y < rhs.y; });
				auto best_skyline_index = distance(_skyline.begin(), bottom_skyline_iter);

				polygon_ptr best_dst_node;
				size_t best_polygon_index;
				coord_t best_skyline_height;
				if (find_polygon_for_skyline_bottom_left_all(best_skyline_index, tmp_polygons, best_dst_node, best_polygon_index, best_skyline_height)) {
					tmp_polygons.remove(best_polygon_index);

				}
				else { // 填坑 将skyline补到与左右两个skyline中最小的同样的高度
					if (best_skyline_index == 0) { _skyline[best_skyline_index].y = _skyline[best_skyline_index + 1].y; }
					else if (best_skyline_index == _skyline.size() - 1) { _skyline[best_skyline_index].y = _skyline[best_skyline_index - 1].y; }
					else { _skyline[best_skyline_index].y = min(_skyline[best_skyline_index - 1].y, _skyline[best_skyline_index + 1].y); }
					merge_skylines(_skyline);
				}
			}
			child_solu->future_height =  get_skyline_height();
		}

		void check_child_solution(shared_ptr<intermediate_solu>& parent_solu, vector<shared_ptr<intermediate_solu>>& child_solu)
		{
			child_solu.clear();
			child_solu.reserve(cfg.filter_width);
			
			if (parent_solu->polygons.size() == 0)  //这个已经放完了
			{
				child_solu.push_back(parent_solu);
				return;
			}


			intermediate_solu test_solu;
			
			std::minstd_rand0 generator(cfg.random_seed);
			list<size_t>::iterator it = parent_solu->polygons.begin();

			int cnt = 0;
		
			//逐个图形尝试往里放
			while (it != parent_solu->polygons.end())
			{
				
				test_solu.polygons = parent_solu->polygons;  //初始化startnode的polygon list
				test_solu.skyline = parent_solu->skyline;
				test_solu._dst = parent_solu->_dst;
				test_solu.current_height = parent_solu->current_height;
				test_solu.area_percentage = parent_solu->area_percentage;

				size_t no = *it;
				
				

				_skyline = test_solu.skyline;

				//当前要放的skyline
				auto bottom_skyline_iter = min_element(_skyline.begin(), _skyline.end(), [](skylinenode_t& lhs, skylinenode_t& rhs) { return lhs.y < rhs.y; });
				auto best_skyline_index = distance(_skyline.begin(), bottom_skyline_iter);

				switch (_src.at(no)->shape())
				{
					
					
					case Shape::R: {
				
						int best_rect_score = -1; // R使用打分策略
						auto rect = dynamic_pointer_cast<rect_t>(_src.at(no));
						coord_t x; int rect_score;
						for (int rotate = 0; rotate <= 1; ++rotate) {            //矩形可以旋转1次
							coord_t w = rect->width, h = rect->height;
							if (rotate) { swap(w, h); }
							if (score_rect_for_skyline_bottom_left(best_skyline_index, w, h, x, rect_score)) {          //打分
								if (best_rect_score < rect_score) {       //分越大越好
									best_rect_score = rect_score;
									rect->lb_point.x = x;
									rect->lb_point.y = _skyline[best_skyline_index].y;
									rect->rotation = rotate ? Rotation::_90_ : Rotation::_0_;
								}
							}
						}
						if (best_rect_score != -1)  //矩形能放下
						{
							auto rect = dynamic_pointer_cast<rect_t>(_src.at(no));
							coord_t w = rect->width, h = rect->height;
							if (rect->rotation == Rotation::_90_) { swap(w, h); }
							skylinenode_t new_skyline_node{ rect->lb_point.x, rect->lb_point.y + h, w };
							if (rect->lb_point.x == _skyline[best_skyline_index].x) { // 靠左
								test_solu.skyline.insert(test_solu.skyline.begin() + best_skyline_index, new_skyline_node);
								test_solu.skyline[best_skyline_index + 1].x += w;
								test_solu.skyline[best_skyline_index + 1].width -= w;
								merge_skylines(test_solu.skyline);
							}
							else { // 靠右
								test_solu.skyline.insert(test_solu.skyline.begin() + best_skyline_index + 1, new_skyline_node);
								test_solu.skyline[best_skyline_index].width -= w;
								merge_skylines(test_solu.skyline);
							}
							
							test_solu.current_height = new_skyline_node.y < test_solu.current_height ? test_solu.current_height : new_skyline_node.y;
							test_solu._dst.push_back(make_shared<rect_t>(*dynamic_pointer_cast<rect_t>(_src.at(no))));
							test_solu.polygons.remove(no);
							//修改area
							double total_area = 0;
							for (int i = 0; i < test_solu._dst.size(); ++i)
							{
								total_area += test_solu._dst[i]->area;
							}
							test_solu.area_percentage = total_area / (test_solu.current_height * _bin_width);
						}
						else                        //这个skyline放不进去了
						{
							//填充skyline , 填坑
							if (best_skyline_index == 0) { test_solu.skyline[best_skyline_index].y = test_solu.skyline[best_skyline_index + 1].y; }
							else if (best_skyline_index == test_solu.skyline.size() - 1) { test_solu.skyline[best_skyline_index].y = test_solu.skyline[best_skyline_index - 1].y; }
							else { test_solu.skyline[best_skyline_index].y = min(test_solu.skyline[best_skyline_index - 1].y, test_solu.skyline[best_skyline_index + 1].y); }
							merge_skylines(test_solu.skyline);
						}
						break;
					}

					
					case Shape::L: {
				
						auto lshape = dynamic_pointer_cast<lshape_t>(_src.at(no));
						skyline_t score_skyline; coord_t score_height; coord_t score_waste;
						if (score_lshape_for_skyline_bottom_left(0, lshape, score_skyline, score_height, score_waste)) {  //L型能放下
							test_solu.skyline = score_skyline;
							//best_skyline_height = best_ltc_height;
							test_solu.current_height = score_height > test_solu.current_height ? score_height : test_solu.current_height;
							test_solu._dst.push_back(make_shared<lshape_t>(*dynamic_pointer_cast<lshape_t>(_src.at(no))));
							test_solu.polygons.remove(no);
							//修改area
							double total_area = 0;
							for (int i = 0; i < test_solu._dst.size(); ++i)
							{
								total_area += test_solu._dst[i]->area;
							}
							test_solu.area_percentage = total_area / (test_solu.current_height * _bin_width);
						}
						else                  //填充skyline
						{
							if (best_skyline_index == 0) { test_solu.skyline[best_skyline_index].y = test_solu.skyline[best_skyline_index + 1].y; }
							else if (best_skyline_index == test_solu.skyline.size() - 1) { test_solu.skyline[best_skyline_index].y = test_solu.skyline[best_skyline_index - 1].y; }
							else { test_solu.skyline[best_skyline_index].y = min(test_solu.skyline[best_skyline_index - 1].y, test_solu.skyline[best_skyline_index + 1].y); }
							merge_skylines(test_solu.skyline);
						}
						break;
					}
					case Shape::T: {
					
						auto tshape = dynamic_pointer_cast<tshape_t>(_src.at(no));
						skyline_t score_skyline; coord_t score_height;
						if (score_tshape_for_skyline_bottom_left(0, tshape, score_skyline, score_height)) {            //T型能放下
							test_solu.skyline = score_skyline;
							//best_skyline_height = best_ltc_height;
							test_solu.current_height = score_height > test_solu.current_height ? score_height : test_solu.current_height;
							test_solu._dst.push_back(make_shared<tshape_t>(*dynamic_pointer_cast<tshape_t>(_src.at(no))));
							test_solu.polygons.remove(no);
							//修改area
							double total_area = 0;
							for (int i = 0; i < test_solu._dst.size(); ++i)
							{
								total_area += test_solu._dst[i]->area;
							}
							test_solu.area_percentage = total_area / (test_solu.current_height * _bin_width);
						}
						else                               //填充skyline
						{
							if (best_skyline_index == 0) { test_solu.skyline[best_skyline_index].y = test_solu.skyline[best_skyline_index + 1].y; }
							else if (best_skyline_index == test_solu.skyline.size() - 1) { test_solu.skyline[best_skyline_index].y = test_solu.skyline[best_skyline_index - 1].y; }
							else { test_solu.skyline[best_skyline_index].y = min(test_solu.skyline[best_skyline_index - 1].y, test_solu.skyline[best_skyline_index + 1].y); }
							merge_skylines(test_solu.skyline);
						}
						break;
					}
					case Shape::C: {
					
						auto concave = dynamic_pointer_cast<concave_t>(_src.at(no));
						skyline_t score_skyline; coord_t score_height;
						if (score_concave_for_skyline_bottom_left(0, concave, score_skyline, score_height)) {           //U型能放下
							test_solu.skyline = score_skyline;
							//best_skyline_height = best_ltc_height;
							test_solu.current_height = score_height > test_solu.current_height ? score_height : test_solu.current_height;
							test_solu._dst.push_back(make_shared<concave_t>(*dynamic_pointer_cast<concave_t>(_src.at(no))));
							test_solu.polygons.remove(no);
							//修改area
							double total_area = 0;
							for (int i = 0; i < test_solu._dst.size(); ++i)
							{
								total_area += test_solu._dst[i]->area;
							}
							test_solu.area_percentage = total_area / (test_solu.current_height * _bin_width);
						}
						else                                  //填充skyline
						{
							if (best_skyline_index == 0) { test_solu.skyline[best_skyline_index].y = test_solu.skyline[best_skyline_index + 1].y; }
							else if (best_skyline_index == test_solu.skyline.size() - 1) { test_solu.skyline[best_skyline_index].y = test_solu.skyline[best_skyline_index - 1].y; }
							else { test_solu.skyline[best_skyline_index].y = min(test_solu.skyline[best_skyline_index - 1].y, test_solu.skyline[best_skyline_index + 1].y); }
							merge_skylines(test_solu.skyline);
						}
						break;
					}
					default: { assert(false); break; }
				}

				//cout << "percentage" << endl;
				//cout << test_solu.area_percentage << endl;
				//以area_percentage为依据在child_solu里面替换
				if (cnt < cfg.filter_width)
				{
					child_solu.push_back(make_shared<intermediate_solu> (test_solu));

					if (cnt == cfg.filter_width - 1)
					{
						sort(child_solu.begin(), child_solu.end(), compare);
					
					}
				}
				else          //检查是否能够替换parent_solu
				{
					int j = cfg.filter_width - 1;  //从后往前找
					while (j >= 0 && (child_solu[j]->area_percentage < test_solu.area_percentage || child_solu[j]->area_percentage == test_solu.area_percentage && generator()%2 == 0))
					{
						j--;
					}
					if (j < cfg.filter_width - 1)  //可以替换
					{
						for (int k = cfg.filter_width - 1; k > j + 1; k--)
						{
							child_solu[k] = child_solu[k - 1];
						}

						child_solu[j + 1] = make_shared<intermediate_solu>(test_solu);

					}

				}

				++cnt;
				++it;
			}
			/*cout << "res" << endl;
			for (int i = 0; i < cfg.filter_width; ++i)
			{
				cout << child_solu[i]->area_percentage << endl;
				cout << "dst::" << endl;
				for (int k = 0; k < child_solu[i]->_dst.size(); k++)
					cout << child_solu[i]->_dst[k]->id << "  ";
				cout << endl;
				list<size_t>::iterator tp = child_solu[i]->polygons.begin();
				cout << "polygon::" << endl;
				while (tp != child_solu[i]->polygons.end())
				{
					cout << *tp << "  ";
					++tp;
				}
				cout << endl;
			}*/
		}

		int mbp_based_beamsearch()      //mbp_solver有当前对应的特定的宽度
		{
			vector<shared_ptr<intermediate_solu>> parent_solu;    //父辈解，每代parent的个数由beam_width决定
			parent_solu.reserve(cfg.beam_width);
			vector<vector<shared_ptr<intermediate_solu>>> child_solu(cfg.beam_width);     //父辈解生成的经过filter后的子代解，每代子代解最多有beam_width*filter_width个
			child_solu.reserve(cfg.beam_width);
			for (auto child : child_solu)
			{
				child.reserve(cfg.filter_width);
			}
			shared_ptr<intermediate_solu> start_node;

			vector<size_t> seq(_src.size());
			iota(seq.begin(), seq.end(), 0);        //seq为输入顺序 0...N
			std::minstd_rand0 generator(cfg.random_seed);
			shuffle(seq.begin(), seq.end(), default_random_engine(cfg.random_seed));
			int total_area = 0;
			
			/*list<size_t>::iterator it = start_node.polygons.begin();
			while (it != start_node.polygons.end())
			{
				
				cout << *it << endl;
				++it;
			}*/
			

			//第一轮
			for (int i = 0; i < seq.size(); ++i)    //最开始一步，先把每个形状都放一遍试试，找出filter_width个空间利用率最大的放置方式
			{
				start_node = make_shared<intermediate_solu>();
				start_node->polygons.assign(seq.begin(), seq.end());  //初始化startnode的polygon list
				start_node->skyline.clear();
				start_node->skyline.push_back({ 0,0,_bin_width });
				start_node->_dst.clear();
				start_node->_dst.reserve(seq.size());

				
				

				_skyline = start_node->skyline;
				
				switch (_src.at(i)->shape()) {

			
				case Shape::R: {
					int best_rect_score = -1; // R使用打分策略
					auto rect = dynamic_pointer_cast<rect_t>(_src.at(i));
					coord_t x; int rect_score;
					for (int rotate = 0; rotate <= 1; ++rotate) {            //矩形可以旋转1次
						coord_t w = rect->width, h = rect->height;
						if (rotate) { swap(w, h); }
						if (score_rect_for_skyline_bottom_left(0, w, h, x, rect_score)) {          //打分
							if (best_rect_score < rect_score) {       //分越大越好
								best_rect_score = rect_score;
								rect->lb_point.x = x;
								rect->lb_point.y = _skyline[0].y;
								rect->rotation = rotate ? Rotation::_90_ : Rotation::_0_;
							}
						}
					}
					if (best_rect_score != -1)  //矩形能放下
					{
						auto rect = dynamic_pointer_cast<rect_t>(_src.at(i));
						coord_t w = rect->width, h = rect->height;
						if (rect->rotation == Rotation::_90_) { swap(w, h); }
						skylinenode_t new_skyline_node{ rect->lb_point.x, rect->lb_point.y + h, w };
						if (rect->lb_point.x == start_node->skyline[0].x) { // 靠左
							start_node->skyline.insert(start_node->skyline.begin(), new_skyline_node);
							start_node->skyline[1].x += w;
							start_node->skyline[1].width -= w;
							
							merge_skylines(start_node->skyline);
						}
						else { // 靠右
							start_node->skyline.insert(start_node->skyline.begin() + 1, new_skyline_node);
							start_node->skyline[0].width -= w;
							merge_skylines(start_node->skyline);
						}
						//best_rect_height = new_skyline_node.y;
						//best_skyline_height = best_rect_height;
						start_node->current_height = new_skyline_node.y;         //更新现有的高度
						start_node->_dst.push_back(make_shared<rect_t>(*dynamic_pointer_cast<rect_t>(_src.at(i))));
					}
					else                        //填充skyline
					{
						assert(false);
					}
					break;
				}


				case Shape::L: {
					auto lshape = dynamic_pointer_cast<lshape_t>(_src.at(i));
					skyline_t score_skyline; coord_t score_height; coord_t score_waste;
					if (score_lshape_for_skyline_bottom_left(0, lshape, score_skyline, score_height, score_waste)) {  //L型能放下
						start_node->skyline = score_skyline;
						//best_skyline_height = best_ltc_height;
						start_node->current_height = score_height;
						start_node->_dst.push_back(make_shared<lshape_t>(*dynamic_pointer_cast<lshape_t>(_src.at(i))));
					}
					else                  //填充skyline
					{
						assert(false);
					}
					break;
				}
				case Shape::T: {
					auto tshape = dynamic_pointer_cast<tshape_t>(_src.at(i));
					skyline_t score_skyline; coord_t score_height;
					if (score_tshape_for_skyline_bottom_left(0, tshape, score_skyline, score_height)) {            //T型能放下
						start_node->skyline = score_skyline;
						//best_skyline_height = best_ltc_height;
						start_node->current_height = score_height;
						start_node->_dst.push_back(make_shared<tshape_t>(*dynamic_pointer_cast<tshape_t>(_src.at(i))));
					}
					else                               //填充skyline
					{
						assert(false);
					}
					break;
				}
				case Shape::C: {
					auto concave = dynamic_pointer_cast<concave_t>(_src.at(i));
					skyline_t score_skyline; coord_t score_height;
					if (score_concave_for_skyline_bottom_left(0, concave, score_skyline, score_height)) {           //U型能放下
						start_node->skyline = score_skyline;
						//best_skyline_height = best_ltc_height;
						start_node->current_height = score_height;
						start_node->_dst.push_back(make_shared<concave_t>(*dynamic_pointer_cast<concave_t>(_src.at(i))));
					}
					else                                  //填充skyline
					{
						assert(false);
					}
					break;
				}
				default: { assert(false); break; }
				}
				start_node->polygons.remove(i);
				start_node->area_percentage = (double)_src.at(i)->area / (_bin_width*start_node->current_height);
				if (i < cfg.beam_width)
				{
					parent_solu.push_back(start_node);
					if (i == cfg.beam_width - 1)  //parent_solu排满了，进行一次排序，按有效面积占比从大到小的顺序从前往后排
					{
						sort(parent_solu.begin(), parent_solu.end(), compare);
					}
				}

				else          //检查是否能够替换parent_solu
				{
					int j = cfg.beam_width - 1;  //从后往前找
					while (j>=0 && (parent_solu[j]->area_percentage < start_node->area_percentage || parent_solu[j]->area_percentage == start_node->area_percentage && generator() % 2 == 0))
					{
						j--;
					}
					if (j < cfg.beam_width - 1)  //可以替换
					{
						for (int k = cfg.beam_width - 1; k > j + 1; k--)
						{
							parent_solu[k] = parent_solu[k - 1];
						}
						
						parent_solu[j + 1] = start_node;
						
					}

				}
				


				
			}
			

			bool endflag = true;
			int cn = 0;
			int total_best = INT_MAX;
			//后续，beamsearch流程
			while (true)
			{
				++cn;
				int bestheight = INT_MAX;
				for (int i = 0; i < cfg.beam_width; ++i)
				{
					check_child_solution(parent_solu[i],child_solu[i]);
					//system("pause");
				}

				int parentnum = 0;
				
				for (int i = 0; i < cfg.beam_width; ++i)
				{
					//cout << child_solu[i].size() << endl;
					for (int j = 0; j < child_solu[i].size(); ++j)
					{
						
						insert_bottom_left_score_beam_search(child_solu[i][j]);//计算按剩余顺序放满后的高度
						if (child_solu[i][j]->future_height < bestheight)
							bestheight = child_solu[i][j]->future_height;
						if (child_solu[i][j]->future_height < total_best)
						{
							total_best = child_solu[i][j]->future_height;
						}
						//cout << child_solu[i][j]->future_height << endl;
						if (parentnum < cfg.beam_width) //还没放满
						{
							//按顺序往parent_solu里面放
							parent_solu[parentnum] = child_solu[i][j];
							++parentnum;
							if (parentnum == cfg.beam_width)
							{
								sort(parent_solu.begin(), parent_solu.end(), compare_byfuture);
							}
						}
						else
						{
							//比较future_height，以此为依据排除
							int k = cfg.beam_width - 1;  //从后往前找
							while (k >= 0 && (parent_solu[k]->future_height > child_solu[i][j]->future_height || parent_solu[k]->future_height == child_solu[i][j]->future_height && generator() % 2 == 0))
							{
								k--;
							}
							if (k < cfg.beam_width - 1)  //可以替换
							{
								for (int m = cfg.beam_width - 1; m > k + 1; m--)
								{
									parent_solu[m] = parent_solu[m - 1];
								}

								parent_solu[k + 1] = child_solu[i][j];

							}
						}
					}
					
				}
				
				//cout << cn << "  loop best: " << bestheight << endl;
				//system("pause");
				endflag = true;   //默认结束搜索
				//选出来的三个parent_solu都放完了，即为结束
				for (int i = 0; i < cfg.beam_width; ++i)
				{
					if (parent_solu[i]->polygons.size() != 0) //还有一个没放完的
					{
						endflag = false;
						break;
					}
				}


				if (endflag)
					break;
			}
			
			//cout << "best: " << total_best << endl;
			return total_best * _bin_width;
			//system("pause");
		}

		/// 基于最下最左和打分策略，贪心构造一个完整解
		bool insert_bottom_left_score(vector<polygon_ptr> &dst) {
			reset();                  //skyline从0开始
			dst.clear(); dst.reserve(_polygons.size());

			while (!_polygons.empty()) {
				//找到最低的skyline，iter是指向这个skyline的指针，index是在vector中的序号
				auto bottom_skyline_iter = min_element(_skyline.begin(), _skyline.end(), [](skylinenode_t &lhs, skylinenode_t &rhs) { return lhs.y < rhs.y; });
				auto best_skyline_index = distance(_skyline.begin(), bottom_skyline_iter);

				polygon_ptr best_dst_node;
				size_t best_polygon_index;
				coord_t best_skyline_height;
				if (find_polygon_for_skyline_bottom_left_all(best_skyline_index, _polygons, best_dst_node, best_polygon_index, best_skyline_height)) {
					_polygons.remove(best_polygon_index);
					dst.push_back(best_dst_node);
					if (best_skyline_height > _bin_height) { return false; } // 超出_bin_height
				}
				else { // 填坑 将skyline补到与左右两个skyline中最小的同样的高度
					if (best_skyline_index == 0) { _skyline[best_skyline_index].y = _skyline[best_skyline_index + 1].y; }
					else if (best_skyline_index == _skyline.size() - 1) { _skyline[best_skyline_index].y = _skyline[best_skyline_index - 1].y; }
					else { _skyline[best_skyline_index].y = min(_skyline[best_skyline_index - 1].y, _skyline[best_skyline_index + 1].y); }
					merge_skylines(_skyline);
				}
			}

			return true;
		}

	private:
		void reset() {
			_skyline.clear();
			_skyline.push_back({ 0,0,_bin_width });
		}

		void init_sort_rules() {
			// 0_输入顺序
			vector<size_t> seq(_src.size());
			iota(seq.begin(), seq.end(), 0);        //seq为输入顺序 0...N
			_sort_rules.reserve(4);
			for (size_t i = 0; i < 4; ++i) { _sort_rules.push_back({ seq, numeric_limits<coord_t>::max() }); }   //初始target_area置为无穷大
			//_tabu_table.insert((_sort_rules[0].*tabu_key)());
			// 1_面积递减
			sort(_sort_rules[1].sequence.begin(), _sort_rules[1].sequence.end(), [this](size_t lhs, size_t rhs) {
				return _src.at(lhs)->area > _src.at(rhs)->area; });
			//_tabu_table.insert((_sort_rules[1].*tabu_key)());
			// 2_最长边递减
			sort(_sort_rules[2].sequence.begin(), _sort_rules[2].sequence.end(), [this](size_t lhs, size_t rhs) {
				return _src.at(lhs)->max_length > _src.at(rhs)->max_length; });
			//_tabu_table.insert((_sort_rules[2].*tabu_key)());
			// 3_随机排序
			shuffle(_sort_rules[3].sequence.begin(), _sort_rules[3].sequence.end(), _gen);
			//_tabu_table.insert((_sort_rules[3].*tabu_key)());

			// 默认输入顺序
			_polygons.assign(_sort_rules[0].sequence.begin(), _sort_rules[0].sequence.end());

			// 离散概率分布初始化
			vector<int> probs; probs.reserve(_sort_rules.size());
			for (int i = 1; i <= _sort_rules.size(); ++i) { probs.push_back(2 * i); }
			_discrete_dist = discrete_distribution<>(probs.begin(), probs.end());
		}

		//LINE 181-190:更改reference sequence

		/// 邻域动作1：交换两个块的顺序
		void swap_sort_rule(SortRule &rule) {
			size_t a = _uniform_dist(_gen);
			size_t b = _uniform_dist(_gen);
			while (a == b) { b = _uniform_dist(_gen); }
			swap(rule.sequence[a], rule.sequence[b]);
		}

		/// 邻域动作2：连续多个块移动
		void rotate_sort_rule(SortRule &rule) {
			size_t a = _uniform_dist(_gen);
			rotate(rule.sequence.begin(), rule.sequence.begin() + a, rule.sequence.end());
		}

		/// 给最下最左的角选择最好的块
		bool find_polygon_for_skyline_bottom_left_partial(size_t skyline_index, const list<size_t> &polygons,
			polygon_ptr &best_dst_node, size_t &best_polygon_index, coord_t &best_skyline_height) {

			int best_score = -1;
			for (size_t p : polygons) {
				switch (_src.at(p)->shape()) {
				case Shape::R: {
					auto rect = dynamic_pointer_cast<rect_t>(_src.at(p));
					coord_t x; int score;
					for (int rotate = 0; rotate <= 1; ++rotate) {
						coord_t w = rect->width, h = rect->height;
						if (rotate) { swap(w, h); }
						if (score_rect_for_skyline_bottom_left(skyline_index, w, h, x, score)) {
							if (best_score < score) {
								best_score = score;
								rect->lb_point.x = x;
								rect->lb_point.y = _skyline[skyline_index].y;
								rect->rotation = rotate ? Rotation::_90_ : Rotation::_0_;
								best_polygon_index = p;
							}
						}
					}
					break;
				}
				case Shape::L: {
					auto lshape = dynamic_pointer_cast<lshape_t>(_src.at(p));
					coord_t waste; // no use
					if (score_lshape_for_skyline_bottom_left(skyline_index, lshape, _skyline, best_skyline_height, waste)) {
						best_polygon_index = p;
						best_dst_node = make_shared<lshape_t>(*lshape);
						return true; // _skyline已被更新
					}
					break;
				}
				case Shape::T: {
					auto tshape = dynamic_pointer_cast<tshape_t>(_src.at(p));
					if (score_tshape_for_skyline_bottom_left(skyline_index, tshape, _skyline, best_skyline_height)) {
						best_polygon_index = p;
						best_dst_node = make_shared<tshape_t>(*tshape);
						return true; // _skyline已被更新
					}
					break;
				}
				case Shape::C: {
					auto concave = dynamic_pointer_cast<concave_t>(_src.at(p));
					if (score_concave_for_skyline_bottom_left(skyline_index, concave, _skyline, best_skyline_height)) {
						best_polygon_index = p;
						best_dst_node = make_shared<concave_t>(*concave);
						return true; // _skyline已被更新
					}
					break;
				}
				default: { assert(false); break; }
				}
			}

			if (best_score == -1) { return false; }

			// 运行到此处一定是矩形，更新`_skyline`和`best_skyline_height`
			assert(_src.at(best_polygon_index)->shape() == Shape::R);
			auto rect = dynamic_pointer_cast<rect_t>(_src.at(best_polygon_index));
			coord_t w = rect->width, h = rect->height;
			if (rect->rotation == Rotation::_90_) { swap(w, h); }
			skylinenode_t new_skyline_node{ rect->lb_point.x, rect->lb_point.y + h, w };
			if (rect->lb_point.x == _skyline[skyline_index].x) { // 靠左
				_skyline.insert(_skyline.begin() + skyline_index, new_skyline_node);
				_skyline[skyline_index + 1].x += w;
				_skyline[skyline_index + 1].width -= w;
				merge_skylines(_skyline);
			}
			else { // 靠右
				_skyline.insert(_skyline.begin() + skyline_index + 1, new_skyline_node);
				_skyline[skyline_index].width -= w;
				merge_skylines(_skyline);
			}
			best_skyline_height = new_skyline_node.y;
			best_dst_node = make_shared<rect_t>(*rect);
			return true;
		}

		/// 给最下最左的角选择最好的块
		bool find_polygon_for_skyline_bottom_left_all(size_t skyline_index, const list<size_t> &polygons,
			polygon_ptr &best_dst_node, size_t &best_polygon_index, coord_t &best_skyline_height) {

			int best_rect_score = -1; // R使用打分策略
			int best_ltc_delta = numeric_limits<int>::max(); // LTC使用skyline.size()变化量delta  numeric_limits<int>::max() int类型的最大值
			coord_t best_l_waste = numeric_limits<coord_t>::max(); // L同时使用最小浪费

			size_t best_rect_index, best_ltc_index;
			skyline_t best_rect_skyline, best_ltc_skyline;
			coord_t best_rect_height, best_ltc_height;

			for (size_t p : polygons) {
				switch (_src.at(p)->shape()) {
				case Shape::R: {
					auto rect = dynamic_pointer_cast<rect_t>(_src.at(p));
					coord_t x; int rect_score;
					for (int rotate = 0; rotate <= 1; ++rotate) {            //矩形可以旋转1次
						coord_t w = rect->width, h = rect->height;
						if (rotate) { swap(w, h); }
						if (score_rect_for_skyline_bottom_left(skyline_index, w, h, x, rect_score)) {          //打分
							if (best_rect_score < rect_score) {       //分越大越好
								best_rect_score = rect_score;
								rect->lb_point.x = x;
								rect->lb_point.y = _skyline[skyline_index].y;
								rect->rotation = rotate ? Rotation::_90_ : Rotation::_0_;
								best_rect_index = p;
							}
						}
					}
					break;
				}
				case Shape::L: {
					auto lshape = dynamic_pointer_cast<lshape_t>(_src.at(p));
					skyline_t score_skyline; coord_t score_height; coord_t score_waste;
					if (score_lshape_for_skyline_bottom_left(skyline_index, lshape, score_skyline, score_height, score_waste)) {  //旋转在打分策略里面
						if (best_l_waste > score_waste ||
							best_l_waste == score_waste && best_ltc_delta > score_skyline.size() - _skyline.size()) {
							best_l_waste = score_waste;
							best_ltc_delta = score_skyline.size() - _skyline.size();
							best_ltc_index = p;
							best_ltc_skyline = score_skyline;
							best_ltc_height = score_height;
						}
					}
					break;
				}
				case Shape::T: {
					auto tshape = dynamic_pointer_cast<tshape_t>(_src.at(p));
					skyline_t score_skyline; coord_t score_height;
					if (score_tshape_for_skyline_bottom_left(skyline_index, tshape, score_skyline, score_height)) {
						if (best_ltc_delta > score_skyline.size() - _skyline.size()) {
							best_ltc_delta = score_skyline.size() - _skyline.size();
							best_ltc_index = p;
							best_ltc_skyline = score_skyline;
							best_ltc_height = score_height;
						}
					}
					break;
				}
				case Shape::C: {
					auto concave = dynamic_pointer_cast<concave_t>(_src.at(p));
					skyline_t score_skyline; coord_t score_height;
					if (score_concave_for_skyline_bottom_left(skyline_index, concave, score_skyline, score_height)) {
						if (best_ltc_delta > score_skyline.size() - _skyline.size()) {
							best_ltc_delta = score_skyline.size() - _skyline.size();
							best_ltc_index = p;
							best_ltc_skyline = score_skyline;
							best_ltc_height = score_height;
						}
					}
					break;
				}
				default: { assert(false); break; }
				}
			}

			if (best_rect_score == -1) { // R放不下
				if (best_ltc_delta == numeric_limits<int>::max())  // LTC放不下
					return false;
				else  // LTC能放下
					best_polygon_index = best_ltc_index;
			}
			else { // R能放下，计算`best_rect_skyline`和`best_rect_height`
				auto rect = dynamic_pointer_cast<rect_t>(_src.at(best_rect_index));
				coord_t w = rect->width, h = rect->height;
				if (rect->rotation == Rotation::_90_) { swap(w, h); }
				best_rect_skyline = _skyline;
				skylinenode_t new_skyline_node{ rect->lb_point.x, rect->lb_point.y + h, w };
				if (rect->lb_point.x == best_rect_skyline[skyline_index].x) { // 靠左
					best_rect_skyline.insert(best_rect_skyline.begin() + skyline_index, new_skyline_node);
					best_rect_skyline[skyline_index + 1].x += w;  
					best_rect_skyline[skyline_index + 1].width -= w;
					merge_skylines(best_rect_skyline);
				}
				else { // 靠右
					best_rect_skyline.insert(best_rect_skyline.begin() + skyline_index + 1, new_skyline_node);
					best_rect_skyline[skyline_index].width -= w;
					merge_skylines(best_rect_skyline);
				}
				best_rect_height = new_skyline_node.y;

				if (best_ltc_delta == numeric_limits<int>::max())  // LTC放不下
					best_polygon_index = best_rect_index;
				else // LTC能放下
					best_polygon_index = _src.at(best_rect_index)->area > _src.at(best_ltc_index)->area ? best_rect_index : best_ltc_index;
			}

			switch (_src.at(best_polygon_index)->shape()) {
			case Shape::R:
				_skyline = best_rect_skyline;
				best_skyline_height = best_rect_height;
				best_dst_node = make_shared<rect_t>(*dynamic_pointer_cast<rect_t>(_src.at(best_polygon_index)));
				break;
			case Shape::L:
				_skyline = best_ltc_skyline;
				best_skyline_height = best_ltc_height;
				best_dst_node = make_shared<lshape_t>(*dynamic_pointer_cast<lshape_t>(_src.at(best_polygon_index)));
				break;
			case Shape::T:
				_skyline = best_ltc_skyline;
				best_skyline_height = best_ltc_height;
				best_dst_node = make_shared<tshape_t>(*dynamic_pointer_cast<tshape_t>(_src.at(best_polygon_index)));
				break;
			case Shape::C:
				_skyline = best_ltc_skyline;
				best_skyline_height = best_ltc_height;
				best_dst_node = make_shared<concave_t>(*dynamic_pointer_cast<concave_t>(_src.at(best_polygon_index)));
				break;
			default:
				assert(false);
				break;
			}
			return true;
		}

		/// Space定义
		struct SkylineSpace {
			coord_t x;
			coord_t y;
			coord_t width;
			coord_t hl;
			coord_t hr;
		};

		SkylineSpace skyline_nodo_to_space(size_t skyline_index) {          //返回一个space，这个space的x，y是对应skyline的x和y，宽度是space的宽度
			coord_t hl, hr;                                                 //hl和hr是左右墙的相对高度
			if (_skyline.size() == 1) {
				hl = hr = INF - _skyline[skyline_index].y;
			}
			else if (skyline_index == 0) {
				hl = INF - _skyline[skyline_index].y;
				hr = _skyline[skyline_index + 1].y - _skyline[skyline_index].y;
			}
			else if (skyline_index == _skyline.size() - 1) {
				hl = _skyline[skyline_index - 1].y - _skyline[skyline_index].y;
				hr = INF - _skyline[skyline_index].y;
			}
			else {
				hl = _skyline[skyline_index - 1].y - _skyline[skyline_index].y;
				hr = _skyline[skyline_index + 1].y - _skyline[skyline_index].y;
			}
			return { _skyline[skyline_index].x, _skyline[skyline_index].y, _skyline[skyline_index].width, hl, hr };
		}

		/// R打分策略
		bool score_rect_for_skyline_bottom_left(size_t skyline_index, coord_t width, coord_t height, coord_t &x, int &score) {
			if (width > _skyline[skyline_index].width) { return false; }

			SkylineSpace space = skyline_nodo_to_space(skyline_index);                 //确定这个skyline对应的space
			if (space.hl >= space.hr) {                                                //左墙更高
				if (width == space.width && height == space.hl) { score = 7; }
				else if (width == space.width && height == space.hr) { score = 6; }
				else if (width == space.width && height > space.hl) { score = 5; }
				else if (width < space.width && height == space.hl) { score = 4; }
				else if (width == space.width && height < space.hl && height > space.hr) { score = 3; }
				else if (width < space.width && height == space.hr) { score = 2; } // 靠右
				else if (width == space.width && height < space.hr) { score = 1; }
				else if (width < space.width && height != space.hl) { score = 0; }
				else { return false; }

				if (score == 2) { x = _skyline[skyline_index].x + _skyline[skyline_index].width - width; }  //靠右放置
				else { x = _skyline[skyline_index].x; }
			}
			else { // hl < hr，右墙更高
				if (width == space.width && height == space.hr) { score = 7; }
				else if (width == space.width && height == space.hl) { score = 6; }
				else if (width == space.width && height > space.hr) { score = 5; }
				else if (width < space.width && height == space.hr) { score = 4; } // 靠右
				else if (width == space.width && height < space.hr && height > space.hl) { score = 3; }
				else if (width < space.width && height == space.hl) { score = 2; }
				else if (width == space.width && height < space.hl) { score = 1; }
				else if (width < space.width && height != space.hr) { score = 0; } // 靠右
				else { return false; }

				if (score == 4 || score == 0) { x = _skyline[skyline_index].x + _skyline[skyline_index].width - width; } //靠右放置
				else { x = _skyline[skyline_index].x; }
			}
			if (x + width > _bin_width) { return false; }

			return true;
		}

		/// L打分策略
		bool score_lshape_for_skyline_bottom_left(size_t skyline_index, lshape_ptr &lshape, skyline_t &skyline, coord_t &skyline_height, coord_t &min_waste) {
			SkylineSpace space = skyline_nodo_to_space(skyline_index);
			skyline_t skyline_0l = _skyline, skyline_0r = _skyline,
				skyline_90 = _skyline, skyline_180 = _skyline,
				skyline_270l = _skyline, skyline_270r = _skyline;    //几个skyline的副本，计算每种放置方式的delta用
			int min_delta = numeric_limits<int>::max(); 
			min_waste = numeric_limits<coord_t>::max(); // 考虑浪费最小

			if (lshape->hd <= space.width) { // 0&靠左
				// lb_point
				point_t lb_point_0l = { skyline_0l[skyline_index].x, skyline_0l[skyline_index].y };    //lshape放在skyline的开始位置（最左端）
				// update
				skyline_0l[skyline_index].y += lshape->vl;           //lshape靠左放置之后，当前的skyline就分成了三部分，一部分是hu长度，高度增加了lshape的高度
				skyline_0l[skyline_index].width = lshape->hu;        //第二部分是L的hm部分，第三部分是lshape右侧剩余的原来的skyline部分
				// add
				skyline_0l.insert(skyline_0l.begin() + skyline_index + 1, {              //hm部分
					skyline_0l[skyline_index].x + skyline_0l[skyline_index].width,
					skyline_0l[skyline_index].y - lshape->vm,
					lshape->hm });
				skyline_0l.insert(skyline_0l.begin() + skyline_index + 2, {              //skyline剩余部分
					skyline_0l[skyline_index + 1].x + skyline_0l[skyline_index + 1].width,
					skyline_0l[skyline_index + 1].y - lshape->vr,
					space.width - lshape->hd });
				// merge
				coord_t skyline_height_bk = skyline_0l[skyline_index].y;  // 防止merge后skyline_index失效
				merge_skylines(skyline_0l);
				// delta
				if (min_delta > skyline_0l.size() - _skyline.size()) {       //delta是增加的skyline，所以越小越好
					min_delta = skyline_0l.size() - _skyline.size();
					min_waste = 0;
					lshape->rotation = Rotation::_0_;
					lshape->lb_point = lb_point_0l;
					skyline = skyline_0l;
					skyline_height = skyline_height_bk;       //这个skyline_height的意思是放下这个lshape后，这个lshape在全局的最高高度，不一定是所有区域中的最高高度
				}
			}

			if (lshape->hd < space.width) { // 0&靠右
				// lb_point
				point_t lb_point_0r = { skyline_0r[skyline_index].x + space.width - lshape->hd, skyline_0r[skyline_index].y };  //靠右放置的参考坐标
				// update
				skyline_0r[skyline_index].width -= lshape->hd;         //更新原skyline的这一部分
				// add
				skyline_0r.insert(skyline_0r.begin() + skyline_index + 1, {
					skyline_0r[skyline_index].x + skyline_0r[skyline_index].width,   //lshape的顶部，hu
					skyline_0r[skyline_index].y + lshape->vl,
					lshape->hu });
				skyline_0r.insert(skyline_0r.begin() + skyline_index + 2, {
					skyline_0r[skyline_index + 1].x + skyline_0r[skyline_index + 1].width,  //lshape的中间，hm
					skyline_0r[skyline_index + 1].y - lshape->vm,
					lshape->hm });
				// merge
				coord_t skyline_height_bk = skyline_0r[skyline_index + 1].y;  //放置后lshape的最高点
				merge_skylines(skyline_0r);                                    //整合skyline
				// delta
				if (min_delta > skyline_0r.size() - _skyline.size()) {
					min_delta = skyline_0r.size() - _skyline.size();
					min_waste = 0;
					lshape->rotation = Rotation::_0_;
					lshape->lb_point = lb_point_0r;
					skyline = skyline_0r;
					skyline_height = skyline_height_bk;
				}
			}

			if (lshape->vl <= space.width) { // 270&靠左
				// lb_point
				point_t lb_point_270l = { skyline_270l[skyline_index].x + lshape->vl, skyline_270l[skyline_index].y }; //参考点的位置
				// update
				skyline_270l[skyline_index].y += lshape->hu;
				skyline_270l[skyline_index].width = lshape->vm;
				// add
				skyline_270l.insert(skyline_270l.begin() + skyline_index + 1, {
					skyline_270l[skyline_index].x + skyline_270l[skyline_index].width,
					skyline_270l[skyline_index].y + lshape->hm,
					lshape->vr });
				skyline_270l.insert(skyline_270l.begin() + skyline_index + 2, {
					skyline_270l[skyline_index + 1].x + skyline_270l[skyline_index + 1].width,
					skyline_270l[skyline_index + 1].y - lshape->hd,
					space.width - lshape->vl });
				// merge
				coord_t skyline_height_bk = skyline_270l[skyline_index + 1].y;   //放置后lshape的最高点
				merge_skylines(skyline_270l);
				// delta
				if (min_delta > skyline_270l.size() - _skyline.size()) {
					min_delta = skyline_270l.size() - _skyline.size();
					min_waste = 0;                                           //不会产生waste
					lshape->rotation = Rotation::_270_;
					lshape->lb_point = lb_point_270l;
					skyline = skyline_270l;
					skyline_height = skyline_height_bk;
				}
			}

			if (lshape->vl < space.width) { // 270&靠右
				// lb_point
				point_t lb_point_270r = { skyline_270r[skyline_index].x + space.width, skyline_270r[skyline_index].y }; //参考点在这个skyline的最右端
				// update
				skyline_270r[skyline_index].width -= lshape->vl;
				// add
				skyline_270r.insert(skyline_270r.begin() + skyline_index + 1, {
					skyline_270r[skyline_index].x + skyline_270r[skyline_index].width,
					skyline_270r[skyline_index].y + lshape->hu,
					lshape->vm });
				skyline_270r.insert(skyline_270r.begin() + skyline_index + 2, {
					skyline_270r[skyline_index + 1].x + skyline_270r[skyline_index + 1].width,
					skyline_270r[skyline_index + 1].y + lshape->hm,
					lshape->vr });
				// merge
				coord_t skyline_height_bk = skyline_270r[skyline_index + 2].y;
				merge_skylines(skyline_270r);
				// delta
				if (min_delta > skyline_270r.size() - _skyline.size()) {
					min_delta = skyline_270r.size() - _skyline.size();
					min_waste = 0;
					lshape->rotation = Rotation::_270_;
					lshape->lb_point = lb_point_270r;
					skyline = skyline_270r;
					skyline_height = skyline_height_bk;
				}
			}
			//旋转90度只能紧靠右侧放置
			if (skyline_index + 1 < skyline_90.size()  //右侧还有skyline
				&& lshape->vm <= skyline_90[skyline_index + 1].width   //右侧的skyline的长度能够让lshape搭上去
				//&& lshape->hm == space.hr // 注释：允许浪费
				&& lshape->vr <= space.width) {                        //中间的空间能让lshape嵌进去
				coord_t new_skyline_height = max(lshape->hd, lshape->hu + space.hr);   //右侧不一定完全紧密相接，中间可能存在空隙 这个height是相对space对应的skyline的高度
				// lb_point  紧靠右侧放置
				point_t lb_point_90 = { skyline_90[skyline_index].x + space.width - lshape->vr, skyline_90[skyline_index].y + new_skyline_height };
				// update
				skyline_90[skyline_index].width -= lshape->vr;
				// add
				skyline_90.insert(skyline_90.begin() + skyline_index + 1, {
					skyline_90[skyline_index].x + skyline_90[skyline_index].width,
					skyline_90[skyline_index].y + new_skyline_height,
					lshape->vl });
				// delete right  修改原先space右侧的skyline
				skyline_90[skyline_index + 2].x += lshape->vm;
				skyline_90[skyline_index + 2].width -= lshape->vm;
				// merge
				coord_t skyline_height_bk = skyline_90[skyline_index + 1].y;
				merge_skylines(skyline_90);
				// waste
				//coord_t old_space = min(space.hl, space.hr) * space.width;
				//coord_t new_space = min(space.hl, new_skyline_height) * skyline_90[skyline_index].width;
				coord_t waste_90 = lshape->hd > lshape->hu + space.hr ?
					(lshape->hd - lshape->hu - space.hr) * lshape->vm : // 上方浪费  原先skyline右侧的浪费
					(lshape->hu + space.hr - lshape->hd) * lshape->vr;  // 下方浪费  lshape的vr边够不到skyline的边
				// delta
				if (min_waste > waste_90 ||         //越小越好
					min_waste == waste_90 && min_delta > skyline_90.size() - _skyline.size()) {
					min_waste = waste_90;
					min_delta = skyline_90.size() - _skyline.size();
					lshape->rotation = Rotation::_90_;
					lshape->lb_point = lb_point_90;
					skyline = skyline_90;
					skyline_height = skyline_height_bk;
				}
			}

			if (skyline_index >= 1        //左边要有
				&& lshape->hm <= skyline_180[skyline_index - 1].width   //旋转后lshape的左边的部分不超过左边的skyline
				//&& lshape->vm == space.hl // 注释：允许浪费
				&& lshape->hu <= space.width) {
				coord_t new_skyline_height = max(lshape->vl, lshape->vr + space.hl);   //相对于底下skyline的高度
				// lb_point
				point_t lb_point_180 = { skyline_180[skyline_index].x + lshape->hu, skyline_180[skyline_index].y + new_skyline_height };
				// update
				skyline_180[skyline_index].x -= lshape->hm;
				skyline_180[skyline_index].y += new_skyline_height;  //lshape旋转180度后最顶端的边
				skyline_180[skyline_index].width = lshape->hd;
				// add
				skyline_180.insert(skyline_180.begin() + skyline_index + 1, {
					skyline_180[skyline_index].x + skyline_180[skyline_index].width,         //原先skyline剩余的部分
					skyline_180[skyline_index].y - new_skyline_height,
					space.width - lshape->hu });
				// delete left
				skyline_180[skyline_index - 1].width -= lshape->hm;         //左侧skyline修改
				// merge
				coord_t skyline_height_bk = skyline_180[skyline_index].y;
				merge_skylines(skyline_180);
				// waste
				//coord_t old_space = min(space.hl, space.hr) * space.width;
				//coord_t new_space = min(space.hr, new_skyline_height) * skyline_180[skyline_index + 1].width;
				coord_t waste_180 = lshape->vl > lshape->vr + space.hl ?
					(lshape->vl - lshape->vr - space.hl) * lshape->hm :   //底部或左侧浪费
					(lshape->vr + space.hl - lshape->vl) * lshape->hu;
				// delta
				if (min_waste > waste_180 ||
					min_waste == waste_180 && min_delta > skyline_180.size() - _skyline.size()) {
					min_waste = waste_180;
					min_delta = skyline_180.size() - _skyline.size();
					lshape->rotation = Rotation::_180_;
					lshape->lb_point = lb_point_180;
					skyline = skyline_180;
					skyline_height = skyline_height_bk;
				}
			}

			if (min_delta == numeric_limits<int>::max()) { return false; }
			return true;
		}

		/// T打分策略   不旋转的时候是倒T型，参考点在左下角
		bool score_tshape_for_skyline_bottom_left(size_t skyline_index, tshape_ptr &tshape, skyline_t &skyline, coord_t &skyline_height) {
			SkylineSpace space = skyline_nodo_to_space(skyline_index);
			skyline_t skyline_0l = _skyline, skyline_0r = _skyline,
				skyline_90 = _skyline, skyline_180 = _skyline, skyline_270 = _skyline;
			int min_delta = numeric_limits<int>::max();

			if (tshape->hd <= space.width) { // 0&靠左
				// lb_point
				point_t lb_point_0l = { skyline_0l[skyline_index].x, skyline_0l[skyline_index].y };  //参考点和skyline的起始点重合了
				// update
				skyline_0l[skyline_index].y += tshape->vld;
				skyline_0l[skyline_index].width = tshape->hl;
				// add
				skyline_0l.insert(skyline_0l.begin() + skyline_index + 1, {
					skyline_0l[skyline_index].x + skyline_0l[skyline_index].width,   //倒T的最顶端
					skyline_0l[skyline_index].y + tshape->vlu,
					tshape->hu });
				skyline_0l.insert(skyline_0l.begin() + skyline_index + 2, {
					skyline_0l[skyline_index + 1].x + skyline_0l[skyline_index + 1].width,  //属于倒T的右侧平面
					skyline_0l[skyline_index + 1].y - tshape->vru,
					tshape->hr });
				skyline_0l.insert(skyline_0l.begin() + skyline_index + 3, {
					skyline_0l[skyline_index + 2].x + skyline_0l[skyline_index + 2].width,  //原skyline的剩余
					skyline_0l[skyline_index + 2].y - tshape->vrd,
					space.width - tshape->hd });
				// merge
				coord_t skyline_height_bk = skyline_0l[skyline_index + 1].y;
				merge_skylines(skyline_0l);
				// delta
				if (min_delta > skyline_0l.size() - _skyline.size()) {
					min_delta = skyline_0l.size() - _skyline.size();
					tshape->rotation = Rotation::_0_;
					tshape->lb_point = lb_point_0l;
					skyline = skyline_0l;
					skyline_height = skyline_height_bk;
				}
			}

			if (tshape->hd < space.width) { // 0&靠右
				// lb_point
				point_t lb_point_0r = { skyline_0r[skyline_index].x + space.width - tshape->hd, skyline_0r[skyline_index].y };
				// update
				skyline_0r[skyline_index].width -= tshape->hd;
				// add  新加三个属于tshape的skyline
				skyline_0r.insert(skyline_0r.begin() + skyline_index + 1, {
					skyline_0r[skyline_index].x + skyline_0r[skyline_index].width,
					skyline_0r[skyline_index].y + tshape->vld,
					tshape->hl });
				skyline_0r.insert(skyline_0r.begin() + skyline_index + 2, {
					skyline_0r[skyline_index + 1].x + skyline_0r[skyline_index + 1].width,
					skyline_0r[skyline_index + 1].y + tshape->vlu,
					tshape->hu });
				skyline_0r.insert(skyline_0r.begin() + skyline_index + 3, {
					skyline_0r[skyline_index + 2].x + skyline_0r[skyline_index + 2].width,
					skyline_0r[skyline_index + 2].y - tshape->vru,
					tshape->hr });
				// merge
				coord_t skyline_height_bk = skyline_0r[skyline_index + 2].y;
				merge_skylines(skyline_0r);
				// delta
				if (min_delta > skyline_0r.size() - _skyline.size()) {
					min_delta = skyline_0r.size() - _skyline.size();
					tshape->rotation = Rotation::_0_;
					tshape->lb_point = lb_point_0r;
					skyline = skyline_0r;
					skyline_height = skyline_height_bk;
				}
			}

			if (skyline_index + 1 < skyline_90.size()  //顺时针转90，右侧要有
				&& tshape->vru <= skyline_90[skyline_index + 1].width   //vru不超出右侧
				&& tshape->hr == space.hr    //右侧正好能搭上
				&& tshape->vrd <= space.width) {          //底部空间足够
				// lb_point
				point_t lb_point_90 = { skyline_90[skyline_index].x + space.width - tshape->vrd, skyline_90[skyline_index].y + tshape->hd };
				// update
				skyline_90[skyline_index].width -= tshape->vrd;
				// add
				skyline_90.insert(skyline_90.begin() + skyline_index + 1, {
					skyline_90[skyline_index].x + skyline_90[skyline_index].width,
					skyline_90[skyline_index].y + tshape->hd,
					tshape->vld });
				skyline_90.insert(skyline_90.begin() + skyline_index + 2, {
					skyline_90[skyline_index + 1].x + skyline_90[skyline_index + 1].width,
					skyline_90[skyline_index + 1].y - tshape->hl,
					tshape->vlu });
				// delete right
				skyline_90[skyline_index + 3].x += tshape->vru;
				skyline_90[skyline_index + 3].width -= tshape->vru;
				// merge
				coord_t skyline_height_bk = skyline_90[skyline_index + 1].y;
				merge_skylines(skyline_90);
				// delta
				if (min_delta > skyline_90.size() - _skyline.size()) {
					min_delta = skyline_90.size() - _skyline.size();
					tshape->rotation = Rotation::_90_;
					tshape->lb_point = lb_point_90;
					skyline = skyline_90;
					skyline_height = skyline_height_bk;
				}
			}

			if (skyline_index + 1 < skyline_180.size() && skyline_index >= 1  //旋转180度是正T型，左右都要有skyline
				&& tshape->hl <= skyline_180[skyline_index + 1].width
				&& tshape->hr <= skyline_180[skyline_index - 1].width   //左右都不超出
				&& tshape->vlu == space.hr
				&& tshape->vru == space.hl                                //左右正好贴合
				&& tshape->hu == space.width) {                         //底部正好嵌入
				// lb_point
				point_t lb_point_180 = { skyline_180[skyline_index].x + tshape->hu + tshape->hl, skyline_180[skyline_index].y + tshape->vlu + tshape->vld };
				// update
				skyline_180[skyline_index].x -= tshape->hr;
				skyline_180[skyline_index].y += (tshape->vru + tshape->vrd);
				skyline_180[skyline_index].width = tshape->hd;
				// delete left
				skyline_180[skyline_index - 1].width -= tshape->hr;
				// delete right
				skyline_180[skyline_index + 1].x += tshape->hl;
				skyline_180[skyline_index + 1].width -= tshape->hl;
				// merge
				coord_t skyline_height_bk = skyline_180[skyline_index].y;
				merge_skylines(skyline_180);
				// delta
				if (min_delta > skyline_180.size() - _skyline.size()) {
					min_delta = skyline_180.size() - _skyline.size();
					tshape->rotation = Rotation::_180_;
					tshape->lb_point = lb_point_180;
					skyline = skyline_180;
					skyline_height = skyline_height_bk;
				}
			}

			if (skyline_index >= 1                               //顺时针270，类似90
				&& tshape->vlu <= skyline_270[skyline_index - 1].width
				&& tshape->hl == space.hl
				&& tshape->vld <= space.width) {    //底部足够
				// lb_point
				point_t lb_point_270 = { skyline_270[skyline_index].x + tshape->vld, skyline_270[skyline_index].y };
				// update
				skyline_270[skyline_index].x -= tshape->vlu;
				skyline_270[skyline_index].y += (tshape->hl + tshape->hu);
				skyline_270[skyline_index].width = tshape->vru;
				// add
				skyline_270.insert(skyline_270.begin() + skyline_index + 1, {
					skyline_270[skyline_index].x + skyline_270[skyline_index].width,
					skyline_270[skyline_index].y + tshape->hr,
					tshape->vrd });
				skyline_270.insert(skyline_270.begin() + skyline_index + 2, {
					skyline_270[skyline_index + 1].x + skyline_270[skyline_index + 1].width,
					skyline_270[skyline_index + 1].y - tshape->hd,
					space.width - tshape->vld });
				// delete left
				skyline_270[skyline_index - 1].width -= tshape->vlu;
				// merge
				coord_t skyline_height_bk = skyline_270[skyline_index + 1].y;
				merge_skylines(skyline_270);
				// delta
				if (min_delta > skyline_270.size() - _skyline.size()) {
					min_delta = skyline_270.size() - _skyline.size();
					tshape->rotation = Rotation::_270_;
					tshape->lb_point = lb_point_270;
					skyline = skyline_270;
					skyline_height = skyline_height_bk;
				}
			}

			if (min_delta == numeric_limits<int>::max()) { return false; }
			return true;
		}

		/// C打分策略
		bool score_concave_for_skyline_bottom_left(size_t skyline_index, concave_ptr &concave, skyline_t &skyline, coord_t &skyline_height) {
			SkylineSpace space = skyline_nodo_to_space(skyline_index);
			skyline_t skyline_0l = _skyline, skyline_0r = _skyline;
			int min_delta = numeric_limits<int>::max();

			if (concave->hd <= space.width) { // 0&靠左
				// lb_point
				point_t lb_point_0l = { skyline_0l[skyline_index].x, skyline_0l[skyline_index].y };  //参考点在U型的左下角
				// update
				skyline_0l[skyline_index].y += concave->vld;
				skyline_0l[skyline_index].width = concave->hl;
				// add
				skyline_0l.insert(skyline_0l.begin() + skyline_index + 1, {
					skyline_0l[skyline_index].x + skyline_0l[skyline_index].width,
					skyline_0l[skyline_index].y - concave->vlu,
					concave->hu });
				skyline_0l.insert(skyline_0l.begin() + skyline_index + 2, {
					skyline_0l[skyline_index + 1].x + skyline_0l[skyline_index + 1].width,
					skyline_0l[skyline_index + 1].y + concave->vru,
					concave->hr });
				skyline_0l.insert(skyline_0l.begin() + skyline_index + 3, {
					skyline_0l[skyline_index + 2].x + skyline_0l[skyline_index + 2].width,
					skyline_0l[skyline_index + 2].y - concave->vrd,
					space.width - concave->hd });
				// merge
				coord_t skyline_height_bk = max(skyline_0l[skyline_index].y, skyline_0l[skyline_index + 2].y);
				merge_skylines(skyline_0l);
				// delta
				if (min_delta > skyline_0l.size() - _skyline.size()) {
					min_delta = skyline_0l.size() - _skyline.size();
					concave->rotation = Rotation::_0_;
					concave->lb_point = lb_point_0l;
					skyline = skyline_0l;
					skyline_height = skyline_height_bk;
				}
			}

			if (concave->hd < space.width) { // 0&靠右
				// lb_point
				point_t lb_point_0r = { skyline_0r[skyline_index].x + space.width - concave->hd, skyline_0r[skyline_index].y };
				// update
				skyline_0r[skyline_index].width -= concave->hd;
				// add
				skyline_0r.insert(skyline_0r.begin() + skyline_index + 1, {
					skyline_0r[skyline_index].x + skyline_0r[skyline_index].width,
					skyline_0r[skyline_index].y + concave->vld,
					concave->hl });
				skyline_0r.insert(skyline_0r.begin() + skyline_index + 2, {
					skyline_0r[skyline_index + 1].x + skyline_0r[skyline_index + 1].width,
					skyline_0r[skyline_index + 1].y - concave->vlu,
					concave->hu });
				skyline_0r.insert(skyline_0r.begin() + skyline_index + 3, {
					skyline_0r[skyline_index + 2].x + skyline_0r[skyline_index + 2].width,
					skyline_0r[skyline_index + 2].y + concave->vru,
					concave->hr });
				// merge
				coord_t skyline_height_bk = max(skyline_0r[skyline_index + 1].y, skyline_0r[skyline_index + 3].y);
				merge_skylines(skyline_0r);
				// delta
				if (min_delta > skyline_0r.size() - _skyline.size()) {
					min_delta = skyline_0r.size() - _skyline.size();
					concave->rotation = Rotation::_0_;
					concave->lb_point = lb_point_0r;
					skyline = skyline_0r;
					skyline_height = skyline_height_bk;
				}
			}

			if (min_delta == numeric_limits<int>::max()) { return false; }
			return true;
		}

		/// 合并同一level的skyline节点.
		static void merge_skylines(skyline_t &skyline) {
			skyline.erase(
				remove_if(skyline.begin(), skyline.end(), [](skylinenode_t &lhs) { return lhs.width <= 0; }),
				skyline.end()
			);
			for (size_t i = 0; i < skyline.size() - 1; ++i) {
				if (skyline[i].y == skyline[i + 1].y) {
					skyline[i].width += skyline[i + 1].width;
					skyline.erase(skyline.begin() + i + 1);
					--i;
				}
			}
		}

	private:
		// 输入
		const vector<polygon_ptr> &_src;
		coord_t _bin_width;
		coord_t _bin_height;

		// 输出
		vector<polygon_ptr> _dst;
		coord_t _obj_area;

		skyline_t _skyline;
		vector<SortRule> _sort_rules; // 排序规则列表，用于RLS
		list<size_t> _polygons;		  // SortRule的sequence，相当于指针，使用list快速删除，放置完毕为空
		//TabuTable _tabu_table;        // 禁忌表
		discrete_distribution<> _discrete_dist;   // 离散概率分布，用于挑选规则(即挑选sequence赋给_polygons)
		uniform_int_distribution<> _uniform_dist; // 均匀分布，用于交换sequence顺序
		default_random_engine &_gen;
	};
	
}

#endif // SMARTMPW_MPWBINPACK_HPP
