#ifndef HEMELBSETUPTOOL_INDEX_H
#define HEMELBSETUPTOOL_INDEX_H

#include <cmath>
#include <iterator>
#include <ostream>

#define BOUNDS_CHECK

#ifdef BOUNDS_CHECK
#include <exception>

class IndexError: public std::exception {
	virtual const char* what() const throw () {
		return "IndexError";
	}
};
#endif

template<typename T> class Iterator;

template<typename T>
class Vec3 {
public:
	T x, y, z;

	// Ctors
	Vec3(T InX, T InY, T InZ) :
		x(InX), y(InY), z(InZ) {
	}

	Vec3() :
		x(0), y(0), z(0) {
	}

	Vec3(const Vec3& other) :
		x(other.x), y(other.y), z(other.z) {
	}

	template<typename U>
	Vec3(const Vec3<U>& other) :
		x(T(other.x)), y(T(other.y)), z(T(other.z)) {

	}
	// Operator Overloads
	bool operator==(const Vec3& other) const {
		return (x == other.x && y == other.y && z == other.z);
	}

	T& operator[](const int i) {
#ifdef BOUNDS_CHECK
		if (i < 0)
			throw IndexError();
#endif
		if (i == 0)
			return x;
		else if (i == 1)
			return y;
		else if (i == 2)
			return z;

#ifdef BOUNDS_CHECK
		throw IndexError();
#endif
	}

	const T& operator[](const int i) const {
#ifdef BOUNDS_CHECK
		if (i < 0)
			throw IndexError();
#endif
		if (i == 0)
			return x;
		else if (i == 1)
			return y;
		else if (i == 2)
			return z;
#ifdef BOUNDS_CHECK
		throw IndexError();
#endif
	}

	//#else
	//	T& operator[](const int i) {
	//		if (i == 0)
	//		return x;
	//		else if (i == 1)
	//		return y;
	//		else
	//		return z;
	//	}
	//	const T& operator[](const int i) const {
	//		if (i == 0)
	//		return x;
	//		else if (i == 1)
	//		return y;
	//		else
	//		return z;
	//	}
	//#endif

	Vec3 operator+(const Vec3& V2) const {
		return Vec3(x + V2.x, y + V2.y, z + V2.z);
	}

	Vec3 operator-(const Vec3& V2) const {
		return Vec3(x - V2.x, y - V2.y, z - V2.z);
	}

	Vec3 operator+(const T& val) const {
		return Vec3(x + val, y + val, z + val);
	}

	Vec3 operator-(const T& val) const {
		return Vec3(x - val, y - val, z - val);
	}

	Vec3 operator-() const {
		return Vec3(-x, -y, -z);
	}

	Vec3 operator/(T S) const {
		return Vec3(x / S, y / S, z / S);
	}

	Vec3 operator%(T S) const {
		return Vec3(x % S, y % S, z % S);
	}

	Vec3 operator*(T S) const {
		return Vec3(x * S, y * S, z * S);
	}

	Vec3& operator+=(const Vec3& V2) {
		x += V2.x;
		y += V2.y;
		z += V2.z;
		return *this;
	}

	Vec3& operator-=(const Vec3& V2) {
		x -= V2.x;
		y -= V2.y;
		z -= V2.z;
		return *this;
	}

	Vec3& operator/=(T S) {
		x /= S;
		y /= S;
		z /= S;
		return *this;
	}

	Vec3 operator*=(T S) {
		x *= S;
		y *= S;
		z *= S;
		return *this;
	}

	// Functions
	static T Dot(const Vec3 &V1, const Vec3 &V2) {
		return V1.x * V2.x + V1.y * V2.y + V1.z * V2.z;
	}

	T Dot(const Vec3 &V1) const {
		return V1.x * x + V1.y * y + V1.z * z;
	}

	static Vec3 Cross(const Vec3& V1, const Vec3 &V2) {
		return Vec3(V1.y * V2.z - V1.z * V2.y, V1.z * V2.x - V1.x * V2.z,
				V1.x * V2.y - V1.y * V2.x);
	}

	Vec3 Cross(const Vec3 &V2) const {
		return Vec3(y * V2.z - z * V2.y, z * V2.x - x * V2.z,
				x * V2.y - y * V2.x);
	}

	template<typename U>
	U Magnitude() const {
		return U(std::sqrt(x * x + y * y + z * z));
	}

	static T Distance(const Vec3& V1, const Vec3& V2) {
		return (V1 - V2).Magnitude<T> ();
	}

	T Distance(const Vec3 &V1) const {
		return (*this - V1).Magnitude<T> ();
	}

	void Normalize() {
		T mag = this->Magnitude<T> ();
		if (mag == 0) {
			return;
		}
		(*this) * (1. / mag);
		return;
	}
	typedef Iterator<T> iterator;
	iterator begin() {
		return iterator(*this, 0U);
	}

	iterator end() {
		return iterator(*this, 3U);
	}

	//friend std::ostream& operator<< (std::ostream& o, Vec3<T> const& v3);
};

template<typename T>
std::ostream& operator<<(std::ostream& o, Vec3<T> const& v3) {
	return o << "x: " << v3.x << "; y: " << v3.y << "; z: " << v3.z;
}

template<typename T>
class Iterator: public std::iterator<std::forward_iterator_tag, T> {
protected:
	Vec3<T>* vec;
	unsigned int i;

public:
	Iterator() :
		vec(NULL), i(0) {
	}

	Iterator(Vec3<T>& vec, unsigned int i = 0) :
		vec(&vec), i(i) {
	}

	Iterator(const Iterator& other) :
		vec(other.vec), i(other.i) {
	}

	Iterator& operator=(const Iterator& other) {
		if (this == &other) {
			return (*this);
		}
		this->vec = other.vec;
		this->i = other.i;

		return (*this);
	}

	Iterator& operator++() {
		this->i++;
		return *this;
	}

	bool operator==(const Iterator& other) const {
		return (this->vec == other.vec) && (this->i == other.i);
	}

	bool operator!=(const Iterator& other) const {
		return !(*this == other);
	}

	T& operator*() {
		return (*this->vec)[this->i];
	}

	T* operator->() {
		return &(*(*this));
	}
};

typedef Vec3<int> Index;
typedef Vec3<double> Vector;

//class BlockIndex: public Index {
//};
//class GlobalSiteIndex: public Index {
//public:
//	BlockIndex ToBlock(int n) const {
//		BlockIndex ans (*reinterpret_cast<const BlockIndex*> (this));
//		ans /= n;
//		return ans;
//	}
//};
//class LocalSiteIndex: public Index {
//};

//template<typename T, typename U>
//Vec3<T> convert(const Vec3<U> source) {
//	return Vec3<T> (T(source.x), T(source.y), T(source.z));
//}
#endif // HEMELBSETUPTOOL_INDEX_H
