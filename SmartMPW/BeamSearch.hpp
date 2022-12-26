#ifndef SMARTMPW_BEAMSEARCH_HPP
#define SMARTMPW_BEAMSEARCH_HPP


#include "Instance.hpp"
#include "MpwBinPack.hpp"


using namespace mbp;

class BeamSearch 
{
	/// 候选宽度定义
	struct CandidateWidth 
	{
		coord_t value;
		int iter;
		unique_ptr<MpwBinPack> mbp_solver; // 存指针，减少排序造成的开销
	};


public:
	BeamSearch() = delete;


	BeamSearch(const Environment& env, const Config& cfg) :
		_env(env), _cfg(cfg), _ins(env), _gen(_cfg.random_seed),
		_obj_area(numeric_limits<coord_t>::max()) {}

	void run()
	{
		_start = clock();

		int total = _ins.get_total_area();
		cout << "total: " << total << endl;
		//vector<coord_t> candidate_widths = cal_candidate_widths_on_interval();
		vector<coord_t> candidate_widths = cal_candidate_widths_on_sqrt();
		vector<CandidateWidth> cw_objs; cw_objs.reserve(candidate_widths.size());
		int best = -1;
		int best_width = 0;
		for (coord_t bin_width : candidate_widths) {

			cw_objs.push_back({ bin_width, 1, unique_ptr<MpwBinPack>(
				new MpwBinPack(_ins.get_polygon_ptrs(), bin_width, INF, _gen)) });
			//cout << " width " << bin_width << endl;
			int res = cw_objs.back().mbp_solver->mbp_based_beamsearch();
			//cout << height << endl; //对每个宽度运行一次完整的beamsearch构建一个完整解，这个完整解就是这个宽度下的最优
			if (best == -1 || best > res)
			{
				best = res;
				best_width = bin_width;
			}
		}
		clock_t _end = clock();
		double sec = (double)(_end - _start) / CLOCKS_PER_SEC;
		cout << sec << endl;
		cout << "best: " << best << endl;
		cout << "fill ratio: " << (double)total / best << endl;
		int height = best / best_width;
		cout << "AR: " << (double)height / best_width << endl;
		//system("pause");
	}


	void record_sol(const string& sol_path) const {
		ofstream sol_file(sol_path);
		for (auto& dst_node : _dst) {
			sol_file << "In Polygon:" << endl;
			for (auto& point : *dst_node->in_points) { sol_file << "(" << point.x << "," << point.y << ")"; }
			dst_node->to_out_points();
			sol_file << endl << "Out Polygon:" << endl;
			for_each(dst_node->out_points.begin(), dst_node->out_points.end(),
				[&](point_t& point) { sol_file << "(" << point.x << "," << point.y << ")"; });
			sol_file << endl;
		}
	}

	void draw_sol(const string& html_path) const {
		utils_visualize_drawer::Drawer html_drawer(html_path, _cfg.ub_width, _cfg.ub_height);
		for (auto& dst_node : _dst) {
			string polygon_str;
			for_each(dst_node->out_points.begin(), dst_node->out_points.end(),
				[&](point_t& point) { polygon_str += to_string(point.x * 0.01) + "," + to_string(point.y * 0.01) + " "; });
			html_drawer.polygon(polygon_str);
		}
	}

#ifndef SUBMIT
	void draw_ins() const {
		ifstream ifs(_env.ins_html_path());
		if (ifs.good()) { return; }
		utils_visualize_drawer::Drawer html_drawer(_env.ins_html_path(), _cfg.ub_width, _cfg.ub_height);
		for (auto& src_node : _ins.get_polygon_ptrs()) {
			string polygon_str;
			for_each(src_node->in_points->begin(), src_node->in_points->end(),
				[&](const point_t& point) { polygon_str += to_string(point.x) + "," + to_string(point.y) + " "; });
			html_drawer.polygon(polygon_str);
		}
	}

	void record_log() const {
		ofstream log_file(_env.log_path(), ios::app);
		log_file.seekp(0, ios::end);
		if (log_file.tellp() <= 0) {
			log_file << "Instance,"
				"InsArea,ObjArea,FillRatio,"
				"Width,Height,WHRatio,"
				"Iteration,Duration,TotalDuration,RandomSeed" << endl;
		}
		log_file << _env.instance_name() << ","
			<< _ins.get_total_area() << "," << _obj_area << "," << _fill_ratio << ","
			<< _width << "," << _height << "," << _wh_ratio << ","
			<< _iteration << "," << _duration << ","
			<< static_cast<double>(clock() - _start) / CLOCKS_PER_SEC << "," << _cfg.random_seed << endl;
	}
#endif // !SUBMIT



private:
	/// 在区间[lb_width, ub_width]内，等距地生成候选宽度
	vector<coord_t> cal_candidate_widths_on_interval(coord_t interval = 1) {
		vector<coord_t> candidate_widths;
		coord_t min_width = 0, max_width = 0;
		for (auto& ptr : _ins.get_polygon_ptrs()) {
			min_width = max(min_width, ptr->max_length);
			max_width += ptr->max_length;
		}
		//min_width = ceil(max(min_width, _cfg.lb_width));
		//max_width = ceil(min(max_width, _cfg.ub_width));

		candidate_widths.reserve(max_width - min_width + 1);
		for (coord_t cw = min_width; cw <= max_width; cw += interval) {
			//if (cw * _cfg.ub_height < _ins.get_total_area()) { continue; }
			candidate_widths.push_back(cw);
		}
		return candidate_widths;
	}

	/// 开平方限制长宽比，削减分支数目
	vector<coord_t> cal_candidate_widths_on_sqrt(coord_t interval = 1) {
		vector<coord_t> candidate_widths;
		coord_t min_width = floor(_cfg.lb_scale * sqrt(_ins.get_total_area()));
		coord_t max_width = ceil(_cfg.ub_scale * sqrt(_ins.get_total_area()));
		for_each(_ins.get_polygon_ptrs().begin(), _ins.get_polygon_ptrs().end(),
			[&](const polygon_ptr& ptr) { min_width = max(min_width, ptr->max_length); });
		max_width = max(max_width, min_width);
		candidate_widths.reserve(max_width - min_width + 1);
		for (coord_t cw = min_width; cw <= max_width; cw += interval) {
			//if (cw * _cfg.ub_height < _ins.get_total_area()) { continue; }
			candidate_widths.push_back(cw);
		}
		return candidate_widths;
	}

	/// 检查cw_obj的RLS结果
	void check_cwobj(const CandidateWidth& cw_obj, int curr_iter = 0) {
		coord_t cw_height = cw_obj.mbp_solver->get_obj_area() / cw_obj.value;
		//if (cw_height > _cfg.ub_height) { // 解高度超出上界，不合法
		//	cw_obj.mbp_solver->set_obj_area(numeric_limits<coord_t>::max());
		//}
		//if (cw_height < _cfg.lb_height) { // 解高度不足下界，按下界计算
		//	cw_obj.mbp_solver->set_obj_area(cw_obj.value * _cfg.lb_height);
		//}
		if (cw_obj.mbp_solver->get_obj_area() < _obj_area) {
			_obj_area = cw_obj.mbp_solver->get_obj_area();
			_fill_ratio = 1.0 * _ins.get_total_area() / _obj_area;
			_width = cw_obj.value;
			_height = cw_height;
			_wh_ratio = 1.0 * max(_width, _height) / min(_width, _height);
			_dst = cw_obj.mbp_solver->get_dst();
			/*for (auto& dst_node : _dst) {
				cout << "In Polygon:" << endl;
				for (auto& point : *dst_node->in_points) { cout << "(" << point.x << "," << point.y << ")"; }
				dst_node->to_out_points();
				cout << endl << "Out Polygon:" << endl;
				for_each(dst_node->out_points.begin(), dst_node->out_points.end(),
					[&](point_t& point) { cout << "(" << point.x << "," << point.y << ")"; });
				cout << endl;
			}*/
			_duration = static_cast<double>(clock() - _start) / CLOCKS_PER_SEC;
			_iteration = curr_iter;
		}
	}



private:
	const Environment& _env;
	const Config& _cfg;

	const Instance _ins;
	default_random_engine _gen;
	clock_t _start;
	double _duration; // 最优解出现时间
	int _iteration;   // 最优解出现迭代次数

	coord_t _obj_area;
	double _fill_ratio;
	coord_t _width;
	coord_t _height;
	double _wh_ratio;
	vector<polygon_ptr> _dst;
};















#endif // SMARTMPW_BEAMSEARCH_HPP

