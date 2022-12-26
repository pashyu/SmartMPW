//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#ifndef SMARTMPW_CONFIG_HPP
#define SMARTMPW_CONFIG_HPP

#ifdef NDEBUG
#define debug_run(x)
#else
#define debug_run(x) x
#endif // !NDEBUG

#include <random>

using coord_t = int;

static constexpr int INF = 0x3f3f3f3f;

static const char* ins_list[]{
	"polygon_area_etc_input_1",
	"polygon_area_etc_input_2",
	"polygon_area_etc_input_3",
	"polygon_area_etc_input_4",
	"polygon_area_etc_input_5",
	"polygon_area_etc_input_6",
	"polygon_area_etc_input_7",
    "polygon_area_etc_input_8",
	"polygon_area_etc_input_9",
	"polygon_area_etc_input_10",
	//"ami33_lt_Ma",
	//"ami49_lt_Xu_1",
	//"ami49_lt_Xu_2"
};

struct Config {
	//unsigned int random_seed = std::random_device{}();
	//unsigned int random_seed = 1196562269;
	unsigned int random_seed = 10101;
	int ub_rls_iter = 9999;  // RLS����������
	int ub_asa_iter = 9999;  // ASA����������
	int ub_asa_time = 180;  // ASA��ʱʱ��
	int filter_width = 2;   //BeamSearch��ÿ��parent node�����µ���node����
	int beam_width = 3;     //BeamSearch��ÿһ��������Ϊparent node������
	coord_t lb_width = 50, ub_width = 400;
	coord_t lb_height = 50, ub_height = 300;
	double lb_scale = 0.9, ub_scale = 1.1;    //���ɺ�ѡ��ȵ�CWGA�����еĦ��ͦ�
} cfg;

#endif // SMARTMPW_CONFIG_HPP
