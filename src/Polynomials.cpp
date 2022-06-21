#include "Polynomials.h"

namespace Circe
{	

	Polynomial::Polynomial()
	{}

	Polynomial::Polynomial(const std::initializer_list<float>& values)
	{
		int i = 0;
		for(float f : values)
		{
			data.insert(std::pair<int,float>(i,f));
			degree = i;
			i++;
		}
	}
	
	float Polynomial::operator()(const float& x) const
	{
		
		float xTot = 1.0f;
		float res = 0.0f;
		for(int i = 0; i<=degree; i++)
		{
			if(hasDegree(i))
			{
				res += data.at(i)*xTot;
			}
			xTot = x*xTot;
		}
		return res;
	}
	
	Polynomial Polynomial::operator+(const Polynomial& other) const
	{
		Polynomial res;
		for(int i = 0; i<=std::max(degree, other.degree); i++)
		{		
			res.setCoefficient(i, getCoefficient(i)+other.getCoefficient(i));
		}
		return res;
	}
	
	void Polynomial::operator+=(const Polynomial& other)
	{
		for(int i = 0; i<=std::max(degree, other.degree); i++)
		{
			setCoefficient(i, getCoefficient(i)+other.getCoefficient(i));
		}
	}
	
	Polynomial Polynomial::operator-(const Polynomial& other) const
	{
		Polynomial res;
		for(int i = 0; i<=std::max(degree, other.degree); i++)
		{
			res.setCoefficient(i, getCoefficient(i)-other.getCoefficient(i));
		}
		return res;
	}
	
	Polynomial Polynomial::operator*(const Polynomial& other) const
	{
		Polynomial res;
		for(int i = 0; i<=degree; i++)
		{
			for(int j = 0; j<=other.degree; j++)
			{	
				res.setCoefficient(i+j, res.getCoefficient(i+j)+getCoefficient(i)*other.getCoefficient(j));
			}
		}
		return res;
	}
	
	Polynomial Polynomial::operator*(const float& f) const
	{
		Polynomial res;
		for(int i = 0; i<=degree; i++)
		{
			res.setCoefficient(i, getCoefficient(i)*f);
		}
		return res;
	}
	
	void Polynomial::operator*=(const float& f)
	{
		for(int i = 0; i<=degree; i++)
		{
			setCoefficient(i, getCoefficient(i)*f);
		}
	}
	
	void Polynomial::setCoefficient(const int& deg, const float& value)
	{	
		if(value != 0.0f)
		{
			data[deg] = value;
			degree = (int)std::max(deg, degree);
		}
		else
		{
			eraseCoefficient(deg);
		}
	}
	
	void Polynomial::eraseCoefficient(const int& deg)
	{
		if(hasDegree(deg))
		{
			data.erase(deg);
			
			if(deg == degree)
			{
				/** Look for the new degree **/
				int i = deg-1;
				while(!hasDegree(deg) && i>=0)
				{	
					i--;
				}
				degree = i;
			}
		}
	}
	
	float Polynomial::getCoefficient(const int& deg) const
	{
		if(hasDegree(deg))
		{
			return data.at(deg);
		}
		return 0.0f;
	}
	
	int Polynomial::getDegree() const
	{
		return degree;
	}
	
	bool Polynomial::hasDegree(const int& deg) const
	{
		return data.find(deg) != data.end();
	}
	
	void Polynomial::raisePower(const int& power)
	{
		for(int i = degree; i>=0; i--)
		{
			if(hasDegree(i))
			{
				setCoefficient(i+power, getCoefficient(i));
				
				eraseCoefficient(i);
			}
		}
	}
	
	std::ostream& operator<<(std::ostream &strm, const Polynomial& P)
	{
		strm << "[";
		for(int i = 0; i<P.getDegree()+1; i++)
		{
			if(P.hasDegree(i))
			{
				if(i==0)
				{
					strm << P.getCoefficient(i);
				}
				else
				{
					if(P.getCoefficient(i) > 0)
					{
						strm << " + " << P.getCoefficient(i) << "x^" << i;
					}
					else
					{
						strm << " " << P.getCoefficient(i) << " x^" << i;
					}
				}
			}
		}
		return strm << "]";
	}
	
	Polynomial LegendrePolynomial(const int& n)
	{
		if(n == 0)
		{
			return Polynomial({1});
		}
		else if(n == 1)
		{
			return Polynomial({0.0f,1.0f});
		}
		else{
			Polynomial Pnmin1 = LegendrePolynomial(n-1);
			Pnmin1.raisePower(1);;
			return Pnmin1*((2.0f*n-1.0f)/(float)(n))-LegendrePolynomial(n-2)*((float)(n-1.0f)/(float)(n));
		}	
	}
	
}