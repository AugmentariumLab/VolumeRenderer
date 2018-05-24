#ifndef FOURBITARRAY_H
#define FOURBITARRAY_H

#include <iostream>
#include <vector>

/**
4-bit array class.
Source: https://stackoverflow.com/questions/25381981/c-2-bit-bitfield-arrays-possible
*/
class FourBitArray {
public:
	typedef unsigned char byte;

	FourBitArray() {}

	FourBitArray(int64_t size) {
		bits.resize((size + 1) / 2, 0);
	}

	~FourBitArray() {
		std::vector<byte>().swap(bits);
	}

	class tbproxy {
	public:
		tbproxy(byte& b, int pos) : b(b), pos(pos) {}

		// getter
		operator int() const {
			return (b >> (pos * 4)) & 1;
		}

		// setter
		tbproxy operator=(int value) {
			const byte mask = ~(1 << (pos * 4));
			b = (b & mask) | (value << (pos * 4));
			return *this;
		}

	private:
		byte& b;
		int pos;
	};

	// create proxy to manipulate object at index
	tbproxy operator[](int64_t index) {
		return tbproxy(bits[index / 2], index & 1);
	}

	void resize(int64_t size) {
		bits.resize((size + 1) / 2, 0);
	}
	void resize(int64_t size, byte byteFill) {
		bits.resize((size + 1) / 2, byteFill);
	}

	// Return number of 2-bit elemnts
	int64_t size() {
		return bits.size() * 2;
	}

	// Return number of bytes
	int64_t bytes() {
		return bits.size();
	}

	FourBitArray& operator=(std::vector<byte> &v) {
		bits = v;
		return *this;
	}

	void shrink_to_fit() {
		bits.shrink_to_fit();
	}

	void swap(FourBitArray swapArr) {
		bits.swap(swapArr.bits);
	}

	void clear() {
		bits.clear();
	}

	//private:
	std::vector<byte> bits;

};

#endif FOURBITARRAY_H