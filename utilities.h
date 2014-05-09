#ifndef UTILITIES_H
#define UTILITIES_H

#include <sstream>

#include <boost/array.hpp>
#include <cmath>

template <typename Iter>
inline void PrintToStream(std::ostream & stream,Iter it, Iter end) {
    stream << "{";
	bool first = true;
	for (; it!=end; ++it) 
	{
		stream << (first?"":",") << *it;
		first = false;
	}
	stream << "}";
}

template <typename Iter>
inline void PrintToStream2D(std::ostream & stream,Iter it, Iter end) {
    stream << "{";
	bool first = true;
	for (; it!=end; ++it) 
	{
		stream << (first?"":",");
		PrintToStream( stream, it->begin(), it->end() );
		first = false;
	}
	stream << "}";
}


class ParameterStream {
public:
	ParameterStream(int argc, char** argv) : argc_(argc), argv_(argv), current_(1) {}
	template<class T> T Read( std::string name )
	{
		T t;
		if( current_ < argc_ )
		{
			std::istringstream is(argv_[current_]);
			is >> t;
			std::cout << name << " = " << t << "\n";
			cmdline_ << argv_[current_] << " ";
		} else
		{
			std::cout << name << " = ";
			std::cin >> t;
			cmdline_ << t << " ";
		}
		current_++;
		return t;
	}
	bool UserInput() {
		return current_ > argc_;
	}
	std::string getCommandLine() {
		return cmdline_.str();
	}
private: 
	int current_;
	int argc_;
	char** argv_;
	std::ostringstream cmdline_;
};

typedef boost::array<double,2> Vector2D;

inline Vector2D AddVectors2D( const Vector2D & v1, const Vector2D & v2 )
{
	Vector2D sum(v1);
	sum[0] += v2[0];
	sum[1] += v2[1];
	return sum;
}

inline Vector2D SubtractVectors2D( const Vector2D & v1, const Vector2D & v2 )
{
	Vector2D diff(v1);
	diff[0] -= v2[0];
	diff[1] -= v2[1];
	return diff;
}
inline Vector2D NegateVector2D( const Vector2D & v )
{
	Vector2D neg;
	neg[0] = -v[0];
	neg[1] = -v[1];
	return neg;
}
inline double NormSquared2D( const Vector2D & v )
{
	return v[0]*v[0] + v[1]*v[1];
}
inline Vector2D InterpolateVectors2D( const Vector2D & v1, const Vector2D & v2, double t)
{
	Vector2D intp;
	intp[0] = (1.0-t)*v1[0] + t*v2[0];
	intp[1] = (1.0-t)*v1[1] + t*v2[1];
	return intp;
}

inline double VectorAngle( const Vector2D & from, const Vector2D & to )
{
	// Returns the angle (between -PI and PI) that one has to rotate "from" 
	// in clockwise direction to align it with "to".
	return atan2( from[0] * to[1] - from[1] * to[0], from[0] * to[0] + from[1] * to[1] );
}

inline int ProperMod(int x, int a)
{
	return (x%a +a)%a;
}

inline int logtwo(int x)
{
	int log=0;
	while(x >>= 1)
	{
		log++;
	}
	return log;
}

#endif