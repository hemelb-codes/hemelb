#ifndef HEMELBSETUPTOOL_INDEX_H
#define HEMELBSETUPTOOL_INDEX_H

#include <cmath>
#include <iterator>
#include <ostream>

class IndexError: public std::exception {
	virtual const char* what() const throw () {
		return "IndexError";
	}
};

// Fully defined below.
template<typename T> class Iterator;

template<typename T>
class Vec3 {
protected:
	// Change this to true for debugging of Vec3 problems
	enum {
		CHECKBOUNDS = false
	};
	T impl[3];
public:
	// Ctors
	Vec3(T InX, T InY, T InZ) {
		this->impl[0] = InX;
		this->impl[1] = InY;
		this->impl[2] = InZ;
	}

	Vec3() {
		this->impl[0] = 0;
		this->impl[1] = 0;
		this->impl[2] = 0;
	}

	Vec3(const Vec3& other) {
		this->impl[0] = other.impl[0];
		this->impl[1] = other.impl[1];
		this->impl[2] = other.impl[2];
	}

	template<typename U>
	Vec3(const Vec3<U>& other) {
		this->impl[0] = T(other.impl[0]);
		this->impl[1] = T(other.impl[1]);
		this->impl[2] = T(other.impl[2]);
	}

	// Operator Overloads
	bool operator==(const Vec3& other) const {
		return (this->impl[0] == other.impl[0] && this->impl[1]
				== other.impl[1] && this->impl[2] == other.impl[2]);
	}

	T& operator[](const int i) {
		if (CHECKBOUNDS) {
			if (i < 0 || i > 2)
				throw IndexError();
		}
		return this->impl[i];
	}

	const T& operator[](const int i) const {
		if (CHECKBOUNDS) {
			if (i < 0 || i > 2)
				throw IndexError();
		}
		return this->impl[i];
	}

	Vec3 operator+(const Vec3& V2) const {
		return Vec3(this->impl[0] + V2.impl[0], this->impl[1] + V2.impl[1],
				this->impl[2] + V2.impl[2]);
	}

	Vec3 operator-(const Vec3& V2) const {
		return Vec3(this->impl[0] - V2.impl[0], this->impl[1] - V2.impl[1],
				this->impl[2] - V2.impl[2]);
	}

	Vec3 operator+(const T& val) const {
		return Vec3(this->impl[0] + val, this->impl[1] + val,
				this->impl[2] + val);
	}

	Vec3 operator-(const T& val) const {
		return Vec3(this->impl[0] - val, this->impl[1] - val,
				this->impl[2] - val);
	}

	Vec3 operator-() const {
		return Vec3(-this->impl[0], -this->impl[1], -this->impl[2]);
	}

	Vec3 operator/(T S) const {
		return Vec3(this->impl[0] / S, this->impl[1] / S, this->impl[2] / S);
	}

	Vec3 operator%(T S) const {
		return Vec3(this->impl[0] % S, this->impl[1] % S, this->impl[2] % S);
	}

	Vec3 operator*(T S) const {
		return Vec3(this->impl[0] * S, this->impl[1] * S, this->impl[2] * S);
	}

	Vec3& operator+=(const Vec3& V2) {
		this->impl[0] += V2.impl[0];
		this->impl[1] += V2.impl[1];
		this->impl[2] += V2.impl[2];
		return *this;
	}

	Vec3& operator-=(const Vec3& V2) {
		this->impl[0] -= V2.impl[0];
		this->impl[1] -= V2.impl[1];
		this->impl[2] -= V2.impl[2];
		return *this;
	}

	Vec3& operator/=(T S) {
		this->impl[0] /= S;
		this->impl[1] /= S;
		this->impl[2] /= S;
		return *this;
	}

	Vec3 operator*=(T S) {
		this->impl[0] *= S;
		this->impl[1] *= S;
		this->impl[2] *= S;
		return *this;
	}

	// Functions
	static T Dot(const Vec3 &V1, const Vec3 &V2) {
		return V1.impl[0] * V2.impl[0] + V1.impl[1] * V2.impl[1] + V1.impl[2]
				* V2.impl[2];
	}

	T Dot(const Vec3 &V1) const {
		return V1.impl[0] * this->impl[0] + V1.impl[1] * this->impl[1]
				+ V1.impl[2] * this->impl[2];
	}

	static Vec3 Cross(const Vec3& V1, const Vec3 &V2) {
		return Vec3(V1.impl[1] * V2.impl[2] - V1.impl[2] * V2.impl[1],
				V1.impl[2] * V2.impl[0] - V1.impl[0] * V2.impl[2],
				V1.impl[0] * V2.impl[1] - V1.impl[1] * V2.impl[0]);
	}

	Vec3 Cross(const Vec3 &V2) const {
		return Vec3(this->impl[1] * V2.impl[2] - this->impl[2] * V2.impl[1],
				this->impl[2] * V2.impl[0] - this->impl[0] * V2.impl[2],
				this->impl[0] * V2.impl[1] - this->impl[1] * V2.impl[0]);
	}

	template<typename U>
	U Magnitude() const {
		return U(
				std::sqrt(
						this->impl[0] * this->impl[0] + this->impl[1]
								* this->impl[1] + this->impl[2] * this->impl[2]));
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

	typedef T value_type;
	typedef Iterator<T> iterator;
	iterator begin() {
		return iterator(*this, 0U);
	}

	iterator end() {
		return iterator(*this, 3U);
	}
	// Be friends with all other Vec3 instantiations.
	template<class > friend class Vec3;
};

template<typename T>
std::ostream& operator<<(std::ostream& o, Vec3<T> const& v3) {
	return o << "x: " << v3[0] << "; y: " << v3[1] << "; z: " << v3[2];
}

template<typename T>
class Iterator: public std::iterator<std::forward_iterator_tag, T> {
public:
	typedef Vec3<T> vector;
protected:
	vector* vec;
	unsigned int i;

public:
	Iterator() :
		vec(NULL), i(0) {
	}

	Iterator(vector& vec, unsigned int i = 0) :
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

#endif // HEMELBSETUPTOOL_INDEX_H
