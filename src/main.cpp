#include <algorithm> // std::swap
#include <cstddef>   // size_t
#include <cstring>   // strlen, strcpy

using namespace std;

struct String {
	String(const char *str = "");
	String(size_t n, char c);
	String(const String &other);
	~String();
	void swap(String &a){
		std::swap(size, a.size);
		std::swap(str, a.str);
	}

    /* Реализуйте оператор присваивания */
	String &operator=(const String &other){
		if(this != &other){
			String(other).swap(*this);
		}
		return *this;
	}

	void append(const String &other);

	size_t size;
	char *str;
};

int main()
{

	return 0;
}
