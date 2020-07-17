#include "divide.h"

Divide* Divide::getInstance() {
	static Divide divide;
	return &divide;
}

Divide::Divide() {

}

