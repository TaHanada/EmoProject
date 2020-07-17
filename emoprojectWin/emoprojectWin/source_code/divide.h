#ifndef DIVIDE_H
#define DIVIDE_H

class Divide {
public:
	int file_divide;
	static Divide* getInstance();
	void change_file_divide(int tmp) {
		file_divide = tmp;
	}
	int get_file_divide() {
		return file_divide;
	}
private:
	Divide();
};

#endif