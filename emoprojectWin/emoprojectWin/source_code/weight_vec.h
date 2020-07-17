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
	int select_size;//�I���ߖT
	int update_size;//�X�V�ߖT
	int H1, H2;//������
	int pop_m;
	int num_obj;
	vector<vector<double>> vec;//�d�݃x�N�g��
	vector < vector<Euclid> > neighbor;

	void resize(int tmp_pop, int tmp_num_obj, int tmp_H1, int tmp_H2, int tmp_selection_size, int tmp_update_size);//resize


	void initialize();//weight vector�̐ݒ�


	void reselect_nighbor();//���������̃x�N�g������ߖT������


	void shuffle();


	void cout_weight();
};


