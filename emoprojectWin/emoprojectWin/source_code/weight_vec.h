#pragma once
#include<iostream>
#include<vector>
#include<algorithm>

using namespace std;

struct Euclid {
	int no;
	double distance;

	bool operator<(const Euclid& another) const {
		if (distance != another.distance) {
			return distance < another.distance;
		}
		return no < another.no;
	}
};



class weight_vec {
public:
	int select_size;//選択近傍
	int update_size;//更新近傍
	int H1, H2;//分割数
	int pop_m;
	int num_obj;
	vector<vector<double>> vec;//重みベクトル
	vector < vector<Euclid> > neighbor;

	void resize(int tmp_pop, int tmp_num_obj, int tmp_H1, int tmp_H2, int tmp_selection_size, int tmp_update_size);//resize


	void initialize();//weight vectorの設定


	void reselect_nighbor();//同じ距離のベクトルから近傍を決定


	void shuffle();


	void cout_weight();
};


