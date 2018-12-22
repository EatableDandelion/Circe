#include "MathObj.h"
#include <math.h>

namespace Circe{
	
	float min(const float& a, const float& b)
	{
		return std::min(a,b);
	}
	
	float max(const float& a, const float& b)
	{
		return std::max(a,b);
	}
	
	float operator ""_deg(long double angle)
	{
		return angle*3.141596f/180.0f;
	}
	
	float operator ""_rad(long double angle)
	{
		return angle;
	}
	
	float cross(const Vec<2>& a, const Vec<2>& b)
	{
		return a(0)*b(1)-a(1)*b(0);
	}
	
	Vec<3> cross(const Vec<3>& a, const Vec<3>& b)
	{
		return Vec<3>(a(1)*b(2)-a(2)*b(1),-a(0)*b(2)+a(2)*b(0),a(0)*b(1)-a(1)*b(0));
	}
	
}