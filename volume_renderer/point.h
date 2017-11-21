#include <Eigen/Core>

//template <typename Scalar, int Dimension>  struct TVector;
template <typename Scalar, int Dimension>  struct TPoint;


/**
* \brief Generic N-dimensional point data structure based on Eigen::Matrix
*/
template <typename _Scalar, int _Dimension> struct TPoint : public Eigen::Matrix<_Scalar, _Dimension, 1> {
public:
	enum {
		Dimension = _Dimension
	};

	typedef _Scalar                             Scalar;
	typedef Eigen::Matrix<Scalar, Dimension, 1> Base;
	//typedef TVector<Scalar, Dimension>          VectorType;
	typedef TPoint<Scalar, Dimension>           PointType;

	/// Create a new point with constant component vlaues
	TPoint(Scalar value = (Scalar)0) { Base::setConstant(value); }

	/// Create a new 2D point (type error if \c Dimension != 2)
	TPoint(Scalar x, Scalar y) : Base(x, y) { }

	/// Create a new 3D point (type error if \c Dimension != 3)
	TPoint(Scalar x, Scalar y, Scalar z) : Base(x, y, z) { }

	/// Create a new 4D point (type error if \c Dimension != 4)
	TPoint(Scalar x, Scalar y, Scalar z, Scalar w) : Base(x, y, z, w) { }

	/// Construct a point from MatrixBase (needed to play nice with Eigen)
	template <typename Derived> TPoint(const Eigen::MatrixBase<Derived>& p)
		: Base(p) { }

	/// Assign a point from MatrixBase (needed to play nice with Eigen)
	template <typename Derived> TPoint &operator=(const Eigen::MatrixBase<Derived>& p) {
		this->Base::operator=(p);
		return *this;
	}

	/// Return a human-readable string summary
	std::string toString() const {
		std::string result;
		for (size_t i = 0; i<Dimension; ++i) {
			result += std::to_string(this->coeff(i));
			if (i + 1 < Dimension)
				result += ", ";
		}
		return "[" + result + "]";
	}
};