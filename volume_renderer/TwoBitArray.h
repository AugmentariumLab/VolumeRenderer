#ifndef TWOBITARRAY_H
#define TWOBITARRAY_H

#include <iostream>
#include <vector>

/**
	2-bit array class.
	Source: https://stackoverflow.com/questions/25381981/c-2-bit-bitfield-arrays-possible
*/
class TwoBitArray {
public:
	typedef unsigned char byte;

	TwoBitArray(){}

	TwoBitArray(int64_t size) {
		bits.resize((size + 3) / 4, 0);
	}

	~TwoBitArray() {
		std::vector<byte>().swap(bits);
	}

	class tbproxy {
	public:
		tbproxy(byte& b, int pos) : b(b), pos(pos) {}

		// getter
		operator int() const {
			return (b >> (pos * 2)) & 3;
		}

		// setter
		tbproxy operator=(int value) {
			const byte mask = ~(3 << (pos * 2));
			b = (b & mask) | (value << (pos * 2));
			return *this;
		}

	private:
		byte& b;
		int pos;
	};

	// create proxy to manipulate object at index
	tbproxy operator[](int64_t index) {
		return tbproxy(bits[index / 4], index & 3);
	}

	void resize(int64_t size) {
		bits.resize((size + 3) / 4, 0);
	}
	void resize(int64_t size, byte byteFill) {
		bits.resize((size + 3) / 4, byteFill);
	}

	// Return number of 2-bit elemnts
	int64_t size() {
		return bits.size() * 4;
	}

	// Return number of bytes
	int64_t bytes() {
		return bits.size();
	}

	TwoBitArray& operator=(std::vector<byte> &v) {
		bits = v;
		return *this;
	}

	void shrink_to_fit() {
		bits.shrink_to_fit();
	}

	void swap(TwoBitArray swapArr) {
		bits.swap(swapArr.bits);
	}

//private:
	std::vector<byte> bits;

};

#endif TWOBITARRAY_H