#include "Spinor.h"
#include <math.h>

namespace Circe
{
	Complex::Complex(const float& angle):c(cos(angle)), s(sin(angle))
	{}
	
	Complex::Complex(const float& real, const float& imaginary):c(real), s(imaginary)
	{}
	
	float Complex::getAngle() const
	{
		return atan2(s,c);
	}
	
	Complex Complex::operator*(const Complex& q) const
	{
		float a=c;float b=s;
		float c=q.c;float d=q.s;
		return Complex(a*c-b*d, b*c+a*d);
	}
	
	Complex Complex::operator=(const Complex& q)
	{
		c=q.c;
		s=q.s;
		return *this;
	}
	
	float Complex::length() const
	{
		return sqrt(c*c+s*s);
	}
	
	void Complex::normalize()
	{
		float length = (*this).length();
		c/=length;
		s/=length;
	}
	
	void Complex::conjugate()
	{
		s*=-1.0f;
	}
	
	float Complex::getReal() const
	{
		return c;
	}
	
	float Complex::getImaginary() const
	{
		return s;
	}
	
	void Complex::addAngle(const float& dtheta)
	{
		Complex q(dtheta);
		float ctemp = c*q.c - s*q.s;
		s = c*q.s + s*q.c;
		c=ctemp;
	}
	
	Quaternion::Quaternion():w(0.0f), x(0.0f), y(0.0f), z(1.0f)
	{}
	
	/*Quaternion::Quaternion(float angle, Vec<3> axis):Spinor<3>(cos(angle/2.0f), axis(0)*sin(angle/2.0f), axis(1)*sin(angle/2.0f), axis(2)*sin(angle/2.0f))
	{}*/
	
	Quaternion::Quaternion(const float& w, const float& x, const float& y, const float& z):w(w), x(x), y(y), z(z)
	{}
	
	Quaternion Quaternion::operator=(const Quaternion& q)
	{
		w=q.w;
		x=q.x;
		y=q.y;
		z=q.z;
		return *this;
	}
	
	Quaternion Quaternion::operator*(const Quaternion& q) const
	{
		return Quaternion(
		w*q.w - x*q.x - y*q.y - z*q.z, //w
		w*q.x + x*q.w + y*q.z - z*q.y, //x
		w*q.y - x*q.z + y*q.w + z*q.x, //y
		w*q.z + x*q.y - y*q.x + z*q.w //z
		);
	}
	
	void Quaternion::addAngle(const float& roll, const float& pitch, const float& yaw)
	{
		Quaternion q = (*this)*Quaternion(0.0f, roll*0.5f, pitch*0.5f, yaw*0.5f);
		w+=q.w;
		x+=q.x;
		y+=q.y;
		z+=q.z;
		normalize();
	}
	
	float Quaternion::length() const
	{
		return sqrt(w*w+x*x+y*y+z*z);
	}
	
	void Quaternion::normalize()
	{
		float length = (*this).length();
		w/=length;
		x/=length;
		y/=length;
		z/=length;
	}
	
	void Quaternion::conjugate()
	{
		x*=-1.0f;
		y*=-1.0f;
		z*=-1.0f;
	}
	
	Quaternion Quaternion::getConjugate() const
	{
		return Quaternion(w, -x, -y, -z);
	}
	
	float Quaternion::getW() const
	{
		return w;
	}
	
	float Quaternion::getX() const
	{
		return x;
	}
	
	float Quaternion::getY() const
	{
		return y;
	}
	
	float Quaternion::getZ() const
	{
		return z;
	}
	
	std::ostream& operator<<(std::ostream &strm, const Complex &q)
	{
		strm << "[" << q.getReal() << " + i" << q.getImaginary();
		return strm << "]";
	}
	
	
	std::ostream& operator<<(std::ostream &strm, const Quaternion &q)
	{
		strm << "[" << q.getW() << ", " << q.getX() << ", " << q.getY() << ", " << q.getZ();
		return strm << "]";
	}
}