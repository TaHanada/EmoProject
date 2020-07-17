#include"weight_vec.h"
#include"math.h"
#include <iomanip>


void weight_vec::resize(int tmp_pop, int tmp_num_obj, int tmp_H1, int tmp_H2, int tmp_selection_size, int tmp_update_size) {
	H1 = tmp_H1;
	H2 = tmp_H2;
	select_size = tmp_selection_size;
	update_size = tmp_update_size;
	pop_m = tmp_pop;
	num_obj = tmp_num_obj;

	//重みベクトルのリサイズ
	vec.resize(pop_m);
	for (int i = 0; i < pop_m; ++i) {
		vec[i].resize(num_obj);
	}


	//距離と近傍のリサイズ
	neighbor.resize(pop_m);
	for (int i = 0; i < pop_m; ++i) {
		neighbor[i].resize(pop_m);
	}


}//コンストラクタ




//weight vectorの設定
void weight_vec::initialize() {
	vector<double> numerator1;
	vector<double> numerator2;
	numerator1.resize(num_obj);
	numerator2.resize(num_obj);
	numerator1[0] = -1;
	numerator2[0] = -1;


	//H1//////////////////////////////////////////////////
	int count = 0;
	int no = 0;
	for (int i = 0; i < pow(H1+1, num_obj); ++i) {
		numerator1[0]++;
		//cout << i << endl;
		//cout << no << endl;
		for (int j = 0; j < num_obj; ++j) {
			if (numerator1[j] > H1) {
				numerator1[j] = 0;
				numerator1[j + 1]++;
			}
		}
		
		count = 0;
		for (int j = 0; j < num_obj; ++j) {
			count += numerator1[j];
		}
		
		if (count == H1) {
			for (int j = 0; j < num_obj; ++j) {

				vec[no][j] = numerator1[j] / (double)H1;
				//cout << vec[no][j] << endl;

				//specification
				if (vec[no][j] == 0.0) {
//					vec[no][j] = pow(10, -6);
				}
			}

			no++;
		}
	}

	//H2///////////////////////////////////////////////////
	if (H2 != 0) {
		count = 0;
		for (int i = 0; i < pow(H2 + 1, num_obj); ++i) {
			numerator2[0]++;

			for (int j = 0; j < num_obj; ++j) {
				if (numerator2[j] > H2) {
					numerator2[j] = 0;
					numerator2[j + 1]++;
				}
			}

			count = 0;
			for (int j = 0; j < num_obj; ++j) {
				count += numerator2[j];
			}
			if (count == H2) {
				for (int j = 0; j < num_obj; ++j) {
					vec[no][j] = numerator2[j] / (double)H2;

					//modification
					vec[no][j] = (vec[no][j] + 1.0 / (double)num_obj) / 2.0;

					//specification
					if (vec[no][j] == 0.0) {
//						vec[no][j] = pow(10, -6);
					}
				}
				no++;
			}

		}
	}


	//近傍の計算
	for (int i = 0; i < pop_m; ++i) {
		for (int j = 0; j < pop_m; ++j) {//ユークリッド距離の計算
			neighbor[i][j].distance = 0;
			neighbor[i][j].no = j;
			for (int k = 0; k < num_obj; ++k) {
				neighbor[i][j].distance += pow(vec[i][k] - vec[j][k], 2);
			}
			neighbor[i][j].distance = sqrt(neighbor[i][j].distance);
		}
		sort(neighbor[i].begin(), neighbor[i].end());//小さい順にソート

	}
}

//距離が同じベクトルが複数ある場合, 各試行の最初にランダムに選びなおす.
void weight_vec::reselect_nighbor() {
	int count;

	int size = update_size;

	if (select_size > update_size) {
		size = select_size;
	}

	for (int i = 0; i < pop_m; ++i) {
		for (int j = 0; j < size - 1; ++j) {

			if (neighbor[i][j].distance - neighbor[i][j + 1].distance == 0) {//隣の距離が同じ

				vector<Euclid> tmp_same;
				for (int k = j; k < pop_m; ++k) {//どこまで同じかチェック

					if (neighbor[i][j].distance - neighbor[i][k].distance == 0) {
						tmp_same.push_back(neighbor[i][k]);//保存しておく
					}
					else {
						random_shuffle(tmp_same.begin(), tmp_same.end());//同じだったやつをシャッフル


						for (int r = j; r < k; ++r) {
							neighbor[i][r] = tmp_same[r - j];//代入
						}
						j = k - 1;
						break;
					}
				}

			}

		}
	}
}


void weight_vec::shuffle() {
	vector<vector<double>> tmp_vec;//重みベクトル
	vector < vector<Euclid> > tmp_neighbor;
	tmp_vec = vec;
	tmp_neighbor = neighbor;
	
	vector<int> rnd;
	vector<int> rnd2;


	for (int i = 0; i < pop_m; i++) {
		rnd.push_back(i);
	}
	random_shuffle(rnd.begin(), rnd.end());//shuffle
	rnd2 = rnd;

	for (int i = 0; i < pop_m; ++i) {
		vec[i] = tmp_vec[rnd[i]];
		neighbor[i] = tmp_neighbor[rnd[i]];
	}

	for (int i = 0; i < pop_m; ++i) {
		rnd2[rnd[i]] = i;
	}

	for (int i = 0; i < pop_m; ++i) {
		for (int j = 0; j < pop_m; ++j) {
			neighbor[i][j].no = rnd2[neighbor[i][j].no];
		}
	}

}


void weight_vec::cout_weight() {
	
	for (int i = 0; i < pop_m; ++i) {
		cout << i << " (";
		for (int j = 0; j < num_obj; ++j) {
			cout << vec[i][j] << " ";
		}
		cout << ") (";
		for (int j = 0; j < update_size; ++j) {
			cout << neighbor[i][j].no << " ";
		}
		cout << ")" << endl;
	}





}