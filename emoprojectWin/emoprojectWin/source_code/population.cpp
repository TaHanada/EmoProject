#include"population.h"

using namespace std;

//解情報の表示
void population::cout_property() {
	for (int i = 0; i < ind.size(); ++i) {
		cout << i << ". (" << ind[i].cost[0] << ", " << ind[i].cost[1] << ") (";
		for (int j = 0; j < ind[i].constraint.size(); ++j) {
			cout << ind[i].constraint[j] << " ";
		}
		cout <<") ("<< ind[i].weight_no << " " << ind[i].fitness << endl;
	}
}

//個体群の表示
void population::cout_population() {

	for (int i = 0; i < ind.size(); ++i) {
		cout << i << "  (";
		for (int j = 0; j < ind[i].var.size(); ++j) {
			cout << ind[i].var[j] << ", ";
		}
		cout << ")" << endl;
	}
}

//erase
void population::erase() {
	ind.erase(ind.begin(), ind.end());
}

//等号演算子
void population::operator=(const population &a) {
	ind.resize(a.ind.size());
	ind = a.ind;
};

//+=のオーバーロード
void population::operator+=(const population &a) {
	for (int i = 0; i < a.ind.size();i++) {
		ind.push_back(a.ind[i]);
	}
};