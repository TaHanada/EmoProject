#include"population.h"

using namespace std;

//�����̕\��
void population::cout_property() {
	for (int i = 0; i < ind.size(); ++i) {
		cout << i << ". (" << ind[i].cost[0] << ", " << ind[i].cost[1] << ") (";
		for (int j = 0; j < ind[i].constraint.size(); ++j) {
			cout << ind[i].constraint[j] << " ";
		}
		cout <<") ("<< ind[i].weight_no << " " << ind[i].fitness << endl;
	}
}

//�̌Q�̕\��
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

//�������Z�q
void population::operator=(const population &a) {
	ind.resize(a.ind.size());
	ind = a.ind;
};

//+=�̃I�[�o�[���[�h
void population::operator+=(const population &a) {
	for (int i = 0; i < a.ind.size();i++) {
		ind.push_back(a.ind[i]);
	}
};