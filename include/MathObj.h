#pragma once
#include <assert.h>
#include <iostream>
#include <ostream>
#include <math.h>
#include <memory>
#include <vector>
#include "Spinor.h"

namespace Circe
{	
	//Fwd declaration
	template<std::size_t N, typename Real> struct Vec;
	template<typename Real> struct Vec<2,Real>;
	template<typename Real> struct Vec<3,Real>;
	template<typename Real> struct Vec<4,Real>;
	template<std::size_t N, std::size_t M, typename Real> struct Mat;
	struct Mat44;
	//struct ITransform;
	//template<std::size_t N> struct Transform;
	struct Transform3;
	struct Position2;
	struct Position3;
	struct Direction2;
	struct Direction3;
	/*template<std::size_t N> struct Direction;
	template<> struct Direction<3>;
	template<std::size_t N> struct Position;
	template<> struct Position<3>;*/
	
	enum REF_FRAME
	{
		LOCAL,
		GLOBAL
	};
	
	//A few useful functions
	template<typename Real>
	extern Real min(const Real& a, const Real& b)
	{
		return std::min(a,b);
	}
	
	template<typename Real>
	extern Real max(const Real& a, const Real& b)
	{
		return std::max(a,b);
	}
	
	template<std::size_t M, typename Real>
	Real max(const Vec<M,Real>& a)
	{
		Real maximum = a(0);
		for(int i = 1; i<M; i++)
		{
			maximum = std::max(maximum, a(i));
		}
		return maximum;
	}
	
	template<std::size_t M, typename Real>
	Real min(const Vec<M,Real>& a)
	{
		Real minimum = a(0);
		for(int i = 1; i<M; i++)
		{
			minimum = std::min(minimum, a(i));
		}
		return minimum;
	}
	
	/*template<typename Real>
	extern Real operator ""_deg(const Real& angle)
	{
		return angle*3.141596/180.0;
	}
	
	template<typename Real>
	extern Real operator ""_rad(long double angle)
	{
		return angle;
	}*/
	
	template<typename Real>
	Real cross(const Vec<2, Real>& a, const Vec<2, Real>& b)
	{
		return a(0)*b(1)-a(1)*b(0);
	}
	
	template<typename Real>
	Vec<3, Real> cross(const Vec<3, Real>& a, const Vec<3, Real>& b)
	{
		return Vec<3, Real>(a(1)*b(2)-a(2)*b(1),-a(0)*b(2)+a(2)*b(0),a(0)*b(1)-a(1)*b(0));
	}
	
	//Scalar product
	template<std::size_t N, typename Real>
	Real dot(const Vec<N, Real>& a, const Vec<N, Real>& b)
	{
		return a.dot(b);
	}
	
	//Inequality operator
	template<std::size_t N, typename Real>
	bool operator!=(const Vec<N, Real>& a, const Vec<N, Real>& b)
	{
		return !(a==b);
	}
	
	//Add two vectors
	template<std::size_t N, typename Real>
	Vec<N, Real> operator+(const Vec<N, Real>& a, const Vec<N, Real>& b)
	{
		Vec<N, Real> res;
		for(int i=0; i<N;++i)
		{
			res(i)=a(i)+b(i);
		}
		return res;
	}
				
	//Subtract two vectors
	template<std::size_t N, typename Real>
	Vec<N, Real> operator-(const Vec<N, Real>& a, const Vec<N, Real>& b)
	{
		Vec<N, Real> res;
		for(int i=0; i<N;++i)
		{
			res(i)=a(i)-b(i);
		}
		return res;
	}
	
	//Negative vector
	template<std::size_t N, typename Real>
	Vec<N, Real> operator-(const Vec<N, Real>& a)
	{
		Vec<N, Real> res;
		for(int i=0; i<N;++i)
		{
			res(i)=-a(i);
		}
		return res;
	}
	
	//Multiply vector term by term
	template<std::size_t N, typename Real>
	Vec<N, Real> operator*(const Vec<N, Real>& a, const Vec<N, Real>& b)
	{
		Vec<N, Real> res;
		for(int i=0; i<N;++i)
		{
			res(i)=a(i)*b(i);
		}
		return res;
	}
	
	//Divide vector term by term
	template<std::size_t N, typename Real>
	Vec<N, Real> operator/(const Vec<N, Real>& a, const Vec<N, Real>& b)
	{
		Vec<N, Real> res;
		for(int i=0; i<N;++i)
		{
			res(i)=a(i)/b(i);
		}
		return res;
	}
	
	//Add float to a vector
	template<std::size_t N, typename Real>
	Vec<N, Real> operator+(const Vec<N, Real>& a, const Real& b)
	{
		Vec<N, Real> res;
		for(int i=0; i<N;++i)
		{
			res(i)=a(i)+b;
		}
		return res;
	}
	
	template<std::size_t N, typename Real>
	Vec<N, Real> operator+(const Real& b, const Vec<N, Real>& a)
	{
		return a+b;
	}
	
	//Multiply float to this vector
	template<std::size_t N, typename Real>
	Vec<N, Real> operator*(const Vec<N, Real>& a, const Real& b)
	{
		Vec<N, Real> res;
		for(int i=0; i<N;++i)
		{
			res(i)=a(i)*b;
		}
		return res;
	}
	
	template<std::size_t N, typename Real>
	Vec<N, Real> operator*(const Real& b, const Vec<N, Real>& a)
	{
		return a*b;
	}
	
	//Divide float to this vector
	template<std::size_t N, typename Real>
	Vec<N, Real> operator/(const Vec<N, Real>& a, const Real& b)
	{
		Vec<N, Real> res;
		for(int i=0; i<N;++i)
		{
			res(i)=a(i)/b;
		}
		return res;
	}
	
	//Get the norm of the vector
	template<std::size_t N, typename Real>
	Real length(const Vec<N, Real>& a)
	{
		return sqrt(dot(a, a));
	}
	
	//Normalize the vector
	template<std::size_t N, typename Real>
	void normalize(Vec<N, Real>& v)
	{
		v=v/Circe::length(v);
	}
	
	//Get distance square
	template<std::size_t N, typename Real>
	Real distanceSquare(const Vec<N, Real>& a, const Vec<N, Real>& b)
	{
		return (a-b).dot(a-b);
	}
	
	//Get distance
	template<std::size_t N, typename Real>
	Real distance(const Vec<N, Real>& a, const Vec<N, Real>& b)
	{
		return sqrt(distanceSquare(a, b));
	}
	
	//Matrix multiplication, returns v*A
	template<std::size_t N, std::size_t M, typename Real>
	Vec<M, Real> operator*(const Vec<N, Real>& a, const Mat<N, M, Real>& m)
	{
		Vec<M, Real> result;
		for(int j = 0; j<M; ++j)
		{
			Real value = (Real)0.0;
			for(int i=0; i<N;++i)
			{
				value+=a(i)*m(i,j);
			}
			result(j)=value;
		}
		return result;
	}

	
	template<std::size_t N, typename Real = float>
	struct Vec
	{
		public:
			Vec()
			{
				for(int i=0;i<N;++i)
				{
					data[i]=(Real)0.0;
				}
			}
			
			Vec(const std::initializer_list<Real>& values)
			{
				size_t i = 0;
				for(Real f : values)
				{
					data[i]=f;
					i++;
				}
			}
			
			//Copy constructor
			Vec(const Vec<N, Real>& v)
			{
				for(int i=0;i<N;++i)
				{
					data[i]=v(i);
				}
			}
			
			Vec(const Real& f)
			{
				for(int i=0;i<N;++i)
				{
					data[i]=f;
				}
			}
			
			//Copy assignment
			Vec<N>& operator=(const Vec<N, Real>& b)
			{
				for(int i=0;i<N;++i)
				{
					data[i]=b(i);
				}
				return *this;
			}
			
			//Add a vector to this one
			void operator+=(const Vec<N, Real>& b)
			{
				for(int i=0; i<N;++i)
				{
					data[i]+=b(i);
				}
			}
			
			//Subtract a vector to this one
			void operator-=(const Vec<N, Real>& b)
			{
				for(int i=0; i<N;++i)
				{
					data[i]-=b(i);
				}
			}

			//Multiply vector term by term
			void operator*=(const Vec<N, Real>& b)
			{
				for(int i=0; i<N;++i)
				{
					data[i]*=b(i);
				}
			}
			
			//Divide vector term by term
			void operator/=(const Vec<N, Real>& b)
			{
				for(int i=0; i<N;++i)
				{
					data[i]/=b(i);
				}
			}

			//Add float to this vector
			void operator+=(const Real& b)
			{
				for(int i=0; i<N;++i)
				{
					data[i]+=b;
				}
			}
			
			//Subtract float to this vector
			void operator-=(const Real& b)
			{
				for(int i=0; i<N;++i)
				{
					data[i]-=b;
				}
			}

			//Multiply float to this vector
			void operator*=(const Real& b)
			{
				for(int i=0; i<N;++i)
				{
					data[i]*=b;
				}
			}
			
			//Divide float to this vector
			void operator/=(const Real& b)
			{
				for(int i=0; i<N;++i)
				{
					data[i]/=b;
				}
			}
			
			//Equality operator
			bool operator==(const Vec<N, Real>& b) const
			{
				for(int i=0;i<N;++i)
				{
					if(data[i]!=b(i))
					{
						return false;
					}
				}
				return true;
			}
	
			//Dot product
			Real dot(const Vec<N, Real>& b) const
			{
				Real result(0.0);
				for(int i=0; i<N;++i)
				{
					result+=data[i]*b(i);
				}
				return result;
			}
		
			//Get one element of this vector
			Real operator()(const unsigned int& index) const
			{
				assert(index<N);
				return data[index];
			}
	
			//Set one element of this vector
			Real &operator()(const unsigned int& index)
			{
				assert(index<N);
				return data[index];
			}
			
			void reset()
			{
				for(int i =0; i<N; ++i)
				{
					data[i]=0.0;
				}
			}
			
		private:
			Real data[N];
	};
	
	template<typename Real = float>
	struct Vec<2, Real>
	{
		public:
			Vec():x((Real)0.0), y((Real)0.0)
			{}
		
			Vec(const Real& x, const Real& y):x(x), y(y)
			{}
			
			//Copy constructor
			Vec(const Vec<2, Real>& v):x(v.x), y(v.y)
			{}
			
			Vec(const Real& f)
			{
				x=f;
				y=f;
			}
			
			//Copy assignment
			Vec<2, Real>& operator=(const Vec<2, Real>& b)
			{
				x=b.x;
				y=b.y;
				return *this;
			}
			
			//Add a vector to this one
			void operator+=(const Vec<2, Real>& b)
			{
				x+=b.x;
				y+=b.y;
			}
			
			//Subtract a vector to this one
			void operator-=(const Vec<2, Real>& b)
			{
				x-=b.x;
				y-=b.y;
			}

			//Multiply vector term by term
			void operator*=(const Vec<2, Real>& b)
			{
				x*=b.x;
				y*=b.y;
			}
			
			//Divide vector term by term
			void operator/=(const Vec<2, Real>& b)
			{
				x/=b.x;
				y/=b.y;
			}

			//Add float to this vector
			void operator+=(const Real& b)
			{
				x+=b;
				y+=b;
			}
			
			//Subtract float to this vector
			void operator-=(const Real& b)
			{
				x-=b;
				y-=b;
			}

			//Multiply float to this vector
			void operator*=(const Real& b)
			{
				x*=b;
				y*=b;
			}
			
			//Divide float to this vector
			void operator/=(const Real& b)
			{
				x/=b;
				y/=b;
			}
			
			//Equality operator
			bool operator==(const Vec<2, Real>& b) const
			{
				return (x==b.x && y==b.y);
			}
	
			//Dot product
			Real dot(const Vec<2, Real>& b) const
			{
				return x*b.x + y*b.y;
			}
	
			//Get one element of this vector
			Real operator()(const unsigned int& index) const
			{
				assert(index<2);
				if(index==0)return x;
				if(index==1)return y;
				return 0.0f;
			}
	
			//Set one element of this vector
			Real &operator()(const unsigned int& index)
			{
				assert(index<2);
				if(index==0)return x;
				if(index==1)return y;
				return x;
			}
			
			void reset()
			{
				x=(Real)0.0;
				y=(Real)0.0;
			}
			
			Vec<2, Real> rotate(const Complex& q)
			{
				float c=q.getReal();
				float s=q.getImaginary();
				float xtemp=c*x-s*y;
				y=s*x+c*y;
				x=xtemp;
				return *this;
			}
			
			Vec<2, Real> rotateInv(const Complex& q)
			{
				float c=q.getReal();
				float s=q.getImaginary();
				float xtemp=c*x+s*y;
				y=-s*x+c*y;
				x=xtemp;
				return *this;
			}
			
		private:
			Real x, y;
	};
	
	template<typename Real = float>
	struct Vec<3, Real>
	{
		public:
			Vec():x((Real)0.0), y((Real)0.0f), z((Real)0.0f)
			{}
		
			Vec(const Real& x, const Real& y, const Real& z):x(x), y(y), z(z)
			{}
			
			//Copy constructor
			Vec(const Vec<3, Real>& v):x(v.x), y(v.y), z(v.z)
			{}
			
			Vec(const Real& f)
			{
				x=f;
				y=f;
				z=f;
			}
			
			//Copy assignment
			Vec<3, Real>& operator=(const Vec<3, Real>& b)
			{
				x=b.x;
				y=b.y;
				z=b.z;
				return *this;
			}
			
			//Add a vector to this one
			void operator+=(const Vec<3, Real>& b)
			{
				x+=b.x;
				y+=b.y;
				z+=b.z;
			}
			
			//Subtract a vector to this one
			void operator-=(const Vec<3, Real>& b)
			{
				x-=b.x;
				y-=b.y;
				z-=b.z;
			}

			//Multiply vector term by term
			void operator*=(const Vec<3, Real>& b)
			{
				x*=b.x;
				y*=b.y;
				z*=b.z;
			}
			
			//Divide vector term by term
			void operator/=(const Vec<3, Real>& b)
			{
				x/=b.x;
				y/=b.y;
				z/=b.z;
			}

			//Add float to this vector
			void operator+=(const Real& b)
			{
				x+=b;
				y+=b;
				z+=b;
			}
			
			//Subtract float to this vector
			void operator-=(const Real& b)
			{
				x-=b;
				y-=b;
				z-=b;
			}

			//Multiply float to this vector
			void operator*=(const Real& b)
			{
				x*=b;
				y*=b;
				z*=b;
			}
			
			//Divide float to this vector
			void operator/=(const Real& b)
			{
				x/=b;
				y/=b;
				z/=b;
			}
			
			//Equality operator
			bool operator==(const Vec<3, Real>& b) const
			{
				return (x==b.x && y==b.y && z==b.z);
			}
	
			//Dot product
			Real dot(const Vec<3, Real>& b) const
			{
				return x*b.x + y*b.y + z*b.z;
			}
			
			//Get one element of this vector
			Real operator()(const unsigned int& index) const
			{
				assert(index<3);
				if(index==0)return x;
				if(index==1)return y;
				if(index==2)return z;
				return (Real)0.0;
			}
	
			//Set one element of this vector
			Real& operator()(const unsigned int& index)
			{
				assert(index<3);
				if(index==0)return x;
				if(index==1)return y;
				if(index==2)return z;
				return x;
			}
			
			void reset()
			{
				x=(Real)0.0;
				y=(Real)0.0;
				z=(Real)0.0;
			}
			
			Vec<3, Real> rotate(const Quaternion& q)
			{
				Quaternion v((Real)0.0, x, y, z);
				Quaternion p=q;
				p.normalize();
				Quaternion pConj=p.getConjugate();
				Quaternion result = p*v*pConj;
				x=result.getX();
				y=result.getY();
				z=result.getZ();
				return *this;
			}
			
			Vec<3, Real> rotateInv(const Quaternion& p)
			{
				return rotate(p.getConjugate());
			}
			
		private:
			Real x, y, z, padding;
	};
	
	template<typename Real = float>
	struct Vec<4, Real>
	{
		public:
			Vec():x((Real)0.0), y((Real)0.0), z((Real)0.0), w((Real)0.0)
			{}
		
			Vec(const Real& x, const Real& y, const Real& z, const Real& w):x(x), y(y), z(z), w(w)
			{}
			
			//Copy constructor
			Vec(const Vec<4, Real>& v):x(v.x), y(v.y), z(v.z), w(v.w)
			{}
			
			Vec(const Real& f)
			{
				x=f;
				y=f;
				z=f;
				w=f;
			}
			
			//Copy assignment
			Vec<4, Real>& operator=(const Vec<4, Real>& b)
			{
				x=b.x;
				y=b.y;
				z=b.z;
				w=b.w;
				return *this;
			}
			
			//Add a vector to this one
			void operator+=(const Vec<4, Real>& b)
			{
				x+=b.x;
				y+=b.y;
				z+=b.z;
				w+=b.w;
			}
			
			//Subtract a vector to this one
			void operator-=(const Vec<4, Real>& b)
			{
				x-=b.x;
				y-=b.y;
				z-=b.z;
				w-=b.w;
			}

			//Multiply vector term by term
			void operator*=(const Vec<4, Real>& b)
			{
				x*=b.x;
				y*=b.y;
				z*=b.z;
				w*=b.w;
			}
			
			//Divide vector term by term
			void operator/=(const Vec<4, Real>& b)
			{
				x/=b.x;
				y/=b.y;
				z/=b.z;
				z/=b.w;
			}

			//Add float to this vector
			void operator+=(const Real& b)
			{
				x+=b;
				y+=b;
				z+=b;
				w+=b;
			}
			
			//Subtract float to this vector
			void operator-=(const Real& b)
			{
				x-=b;
				y-=b;
				z-=b;
				w-=b;
			}

			//Multiply float to this vector
			void operator*=(const Real& b)
			{
				x*=b;
				y*=b;
				z*=b;
				w*=b;
			}
			
			//Divide float to this vector
			void operator/=(const Real& b)
			{
				x/=b;
				y/=b;
				z/=b;
				z/=b;
			}
			
			//Equality operator
			bool operator==(const Vec<4, Real>& b) const
			{
				return (x==b.x && y==b.y && z==b.z && w==b.w);
			}
	
			//Dot product
			Real dot(const Vec<4, Real>& b) const
			{
				return x*b.x + y*b.y + z*b.z + w*b.w;
			}
	
			//Get one element of this vector
			Real operator()(const unsigned int& index) const
			{
				assert(index<4);
				if(index==0)return x;
				if(index==1)return y;
				if(index==2)return z;
				if(index==3)return w;
				return (Real)0.0f;
			}
	
			//Set one element of this vector
			Real &operator()(const unsigned int& index)
			{
				assert(index<4);
				if(index==0)return x;
				if(index==1)return y;
				if(index==2)return z;
				if(index==3)return w;
				return x;
			}
			
			void reset()
			{
				x=(Real)0.0f;
				y=(Real)0.0f;
				z=(Real)0.0f;
				w=(Real)0.0f;
			}
			
		private:
			Real x, y, z, w;
	};
	
	using Vec2=Vec<2, float>;
	using Vec3=Vec<3, float>;
	using Vec4=Vec<4, float>;
	
	/*
	template<std::size_t N>
	struct Position : public Vec<N>
	{
		public:
			template<typename... Args>
			Position(const REF_FRAME& frame, Args... args):frame(frame), vector(std::make_shared<Vec<N>>(std::forward<Args>(args)...))
			{}
			
			Position(const REF_FRAME& frame, const Vec<N> vector):frame(frame), vector(std::make_shared<Vec<N>>(vector))
			{}
			
			Position(const Position<N>& other):frame(other.frame), vector(std::make_shared<Vec<N>>(*(other.vector)))
			{}
			
			Direction<N> operator-(const Position<N>& p2)
			{
				return Direction<N>(frame, vector-p2.vector);
			}
			
			Position<N> operator+(const Direction<N>& d)
			{
				return Position<N>(frame, vector+d.getValue());
			}
			
			Position<N> operator*(const float& f)
			{
				return Position<N>(frame, (*vector)*f);
			}
			
			void operator+=(const Direction<N>& d)
			{
				(*vector)+=d.getValue();
			}
			
			Position<N> operator-(const Direction<N>& d)
			{
				return Position<N>(frame, vector-d.getValue());
			}
			
			void operator-=(const Direction<N>& d)
			{
				(*vector)-=d.getValue();
			}
			
			float operator()(const unsigned int& index) const
			{
				float res = (*vector)(index);
				return res;
			}
			
			REF_FRAME getFrame() const
			{
				return frame;
			}
			
			Vec<N> getValue() const
			{
				return *vector;
			}
			
		private:
			std::shared_ptr<Vec<N>> vector;
			REF_FRAME frame;
	};
	
	template<std::size_t N>
	struct Direction : public Vec<N>
	{
		public:
			template<typename... Args>
			Direction(const REF_FRAME& frame, Args... args):frame(frame), vector(std::make_shared<Vec<N>>(std::forward<Args>(args)...))
			{}
			
			Direction(const REF_FRAME& frame, const Vec<N> vector):frame(frame), vector(std::make_shared<Vec<N>>(vector))
			{}
			
			Direction(const Direction<N>& other):frame(other.frame), vector(std::make_shared<Vec<N>>(*(other.vector)))
			{}
			
			Direction<N> operator+(const Direction<N>& d2)
			{
				return Direction<N>(frame, vector+d2.vector);
			}
			
			Direction<N> operator-(const Direction<N>& d2)
			{
				return Direction<N>(frame, vector-d2.vector);
			}
			
			void operator*=(const float& f)
			{
				(*vector) *= f;
			}
			
			Position<N> operator+(const Position<N>& p2)
			{
				return Position<N>(frame, vector+p2.getValue());
			}
			
			Position<N> operator-(const Position<N>& p2)
			{
				return Position<N>(frame, vector-p2.getValue());
			}
			
			float operator()(const unsigned int& index) const
			{
				float res = (*vector)(index);
				return res;
			}
			
			REF_FRAME getFrame() const
			{
				return frame;
			}
			
			Vec<N> getValue() const
			{
				return *vector;
			}
			
		private:
			std::shared_ptr<Vec<N>> vector;
			REF_FRAME frame;
	};
	*/

	struct Position2
	{
		public:
			Position2(const REF_FRAME& frame, const float& x, const float& y);
			
			Position2(const REF_FRAME& frame, const Vec2 vector);
			
			Position2(const Position2& other);
			
			Position2& operator=(const Position2& other);
			
			Direction2 operator-(const Position2& p2);
			
			Position2 operator+(const Direction2& d);
			
			Position2 operator*(const float& f);
			
			void operator+=(const Direction2& d);
			
			Position2 operator-(const Direction2& d);
			
			void operator-=(const Direction2& d);
			
			float& operator()(const unsigned int& index);
			
			float operator()(const unsigned int& index) const;
			
			float distance2(const Position2& other);
			
			float distance(const Position2& other);
			
			REF_FRAME getFrame() const;
			
			Vec2 getValue() const;
			
		private:
			float x, y;
			REF_FRAME frame;
	};

	struct Position3
	{
		public:
			Position3(const REF_FRAME& frame, const float& x, const float& y, const float& z);
			
			Position3(const REF_FRAME& frame, const Vec3 vector);
			
			Position3(const Position3& other);
			
			Position3& operator=(const Position3& other);
			
			Direction3 operator-(const Position3& p2);
			
			Position3 operator+(const Direction3& d);
			
			Position3 operator*(const float& f);
			
			void operator+=(const Direction3& d);
			
			Position3 operator-(const Direction3& d);
			
			void operator-=(const Direction3& d);
			
			float& operator()(const unsigned int& index);
			
			float operator()(const unsigned int& index) const;
			
			float distance2(const Position3& other);
			
			float distance(const Position3& other);
			
			REF_FRAME getFrame() const;
			
			Vec3 getValue() const;
			
		private:
			float x, y, z;
			REF_FRAME frame;
	};
	
	struct Direction2
	{
		public:
			Direction2(const REF_FRAME& frame, const float& x, const float& y);
			
			Direction2(const REF_FRAME& frame, const Vec2 vector);
			
			Direction2(const Direction2& other);
			
			Direction2& operator=(const Direction2& other);
			
			Direction2 operator+(const Direction2& d2);
			
			Direction2 operator-(const Direction2& d2);
			
			void operator*=(const float& f);
			
			Position2 operator+(const Position2& d2);
			
			Position2 operator-(const Position2& d2);
			
			float& operator()(const unsigned int& index);
			
			float operator()(const unsigned int& index) const;
			
			REF_FRAME getFrame() const;
			
			Vec2 getValue() const;
			
		private:
			float x, y;
			REF_FRAME frame;
	};
	
	struct Direction3
	{
		public:
			Direction3(const REF_FRAME& frame, const float& x, const float& y, const float& z);
			
			Direction3(const REF_FRAME& frame, const Vec3 vector);
			
			Direction3(const Direction3& other);
			
			Direction3& operator=(const Direction3& other);
			
			Direction3 operator+(const Direction3& d2);
			
			Direction3 operator-(const Direction3& d2);
			
			void operator*=(const float& f);
			
			Position3 operator+(const Position3& d2);
			
			Position3 operator-(const Position3& d2);
			
			float& operator()(const unsigned int& index);
			
			float operator()(const unsigned int& index) const;
			
			REF_FRAME getFrame() const;
			
			Vec3 getValue() const;
			
		private:
			float x, y, z;
			REF_FRAME frame;
	};
	
	
	
	template<std::size_t N, std::size_t M=N, typename Real = float>
	struct Mat
	{	
		public:
			Mat()
			{	
				for(int i = 0; i<N; ++i)
				{
					for(int j = 0; j<M; ++j)
					{
						data.push_back((Real)0.0);
					}
				}	
			}
			
			Mat(const std::initializer_list<Real>& values)
			{
				
				for(int i = 0; i<N; ++i)
				{
					for(int j = 0; j<M; ++j)
					{
						data.push_back((Real)0.0);
					}
				}	
				
				size_t i = 0;
				size_t j = 0;
				for(Real f : values)
				{
					if(j==M)
					{
						i++;
						j=0;
					}
					data[i*M+j] =f;
					j++;
				}
			}
			
			
			Mat<N,M,Real> &operator=(const Mat<N,M,Real>& m)
			{
				for(int i=0;i<N;++i)
				{
					for(int j = 0; j<M; ++j)
					{
						data[i*M+j] = m(i,j);
					}
				}
				return *this;
			}
			
			//Turn the matrix to I
			void setIdentity()
			{
				for(int i = 0; i<min(N, M); ++i)
				{
					data[i*M+i]=(Real)1.0;
				}	
			}
			
			//Set all the coeffs to 0
			void reset()
			{
				for(int i = 0; i<N; ++i)
				{
					for(int j = 0; j<M; ++j)
					{
						data[i*M+j]=(Real)0.0;
					}
				}	
			}
			
			Mat<M,N,Real> getTranspose() const
			{
				Mat<M,N> result;
				for(int i = 0; i<N; ++i)
				{
					for(int j = 0; j<M; ++j)
					{
						result(j,i)=get(i,j);
					}
				}	
				return result;
			}
			
			//Matrix multiplication to matrix
			template<std::size_t O>
			Mat<N,O,Real> operator*(const Mat<M,O,Real>& m) const
			{
				Mat<N,O> result;
				for(int i = 0; i<N; ++i)
				{
					for(int j = 0; j<O; ++j)
					{
						Real value = (Real)0.0;
						for(int k=0;k<M; ++k)
						{
							value+=get(i,k)*m(k,j);
						}
						result(i,j)=value;
					}
				}	
				return result;
			}
			
			//Matrix multiplication to vector, returns A*v
			Vec<N,Real> operator*(const Vec<M,Real>& v) const
			{
				Vec<N,Real> result;
				for(int i=0; i<N;++i)
				{
					Real value = (Real)0.0;
					for(int j = 0; j<M; ++j)
					{
						value+=v(j)*get(i,j);
					}
					result(i)=value;
				}
				return result;
			}
			
			//Get the matrix multiplied by a scalar
			Mat<N,M,Real> operator*(const Real& f) const
			{
				Mat<N,M,Real> result;
				for(int i = 0; i<N; ++i)
				{
					for(int j = 0; j<M; ++j)
					{
						result(i,j)=get(i,j)*f;
					}
				}	
				return result;	
			}
			
			Mat<N,M,Real> operator/(const Real& f) const
			{
				Mat<N,M,Real> result;
				for(int i = 0; i<N; ++i)
				{
					for(int j = 0; j<M; ++j)
					{
						result(i,j)=get(i,j)/f;
					}
				}	
				return result;	
			}
			
			void operator+=(const Mat<N,M,Real>& m)
			{
				for(int i = 0; i<N; ++i)
				{
					for(int j = 0; j<M; ++j)
					{
						data[i*M+j] += m(i,j);
					}
				}	
			}
			
			Mat<N,M,Real> operator+(const Mat<N,M,Real>& m) const
			{
				Mat<N,M,Real> res;
				for(int i = 0; i<N; ++i)
				{
					for(int j = 0; j<M; ++j)
					{
						res(i,j)=get(i,j)+m(i,j);
					}
				}
				return res;
			}
			
			void operator-=(const Mat<N,M,Real>& m)
			{
				for(int i = 0; i<N; ++i)
				{
					for(int j = 0; j<M; ++j)
					{
						data[i*M+j]-=m(i,j);
					}
				}	
			}
			
			//Multiply by a scalar
			void operator*=(const Real& f)
			{
				for(int i = 0; i<N; ++i)
				{
					for(int j = 0; j<M; ++j)
					{
						data[i*M+j]*=f;
					}
				}
			}
			
			//Get the coefficient at row i and column j
			Real operator()(const int& i, const int& j) const
			{
				return data[i*M+j];
			}
			
			//Set the coefficient at row i and column j
			Real &operator()(const int& i, const int& j)
			{
				return data[i*M+j];;
			}
			
			Real* getData()
			{
				static Real res[N*M];
				int k=0;
				for(int i=0; i<N; i++)
				{
					for(int j=0; j<M; j++)
					{
						res[k] = get(i,j);
						k++;
					}
				}
				return res;
			}
			
			Real get(const int& i, const int& j) const
			{
				return data[i*M+j];
			}
			
			Vec<M, Real> getRow(const int& i) const
			{
				Vec<N, Real> row;
				for(int j=0; j<M; j++)
				{
					row(j) = get(i,j);
				}
				return row;
			}
			
			Vec<N, Real> getColumn(const int& j) const
			{
				Vec<N, Real> column;
				for(int i=0; i<N; i++)
				{
					column(i) = get(i,j);
				}
				return column;
			}
			
			void setRow(const int& i, const Vec<M, Real>& row)
			{
				for(int j=0; j<M; j++)
				{
					data[i*M+j] = row(j);
				}
			}
			
			void setColumn(const int& j, const Vec<N, Real>& column)
			{
				for(int i=0; i<N; i++)
				{
					data[i*M+j] = column(i);
				}
			}
				
		private:
			std::vector<Real> data;
			//Real data[N*M];
	};
	
	
	struct Mat44 : public Mat<4>
	{
		public:
			static Mat44 positionMatrix(const Position2& v)
			{
				return Mat44::positionMatrix(v(0), v(1), 0.0f);
			}
			
			static Mat44 positionMatrix(const Position3& v)
			{
				return Mat44::positionMatrix(v(0), v(1), v(2));
			}
		
			static Mat44 positionMatrix(const float& x, const float& y, const float& z)
			{
				Mat44 m;
				m.setIdentity();
				m(0,3)=x;
				m(1,3)=y;
				m(2,3)=z;
				return m;
			}
			
			static Mat44 scaleMatrix(const Vec<2>& v)
			{
				return Mat44::scaleMatrix(v(0), v(1), 1.0f);
			}
			
			static Mat44 scaleMatrix(const Vec<3>& v)
			{
				return Mat44::scaleMatrix(v(0), v(1), v(2));
			}
			
			static Mat44 scaleMatrix(const float& x, const float& y, const float& z)
			{
				Mat44 m;
				m.setIdentity();
				m(0,0)=x;
				m(1,1)=y;
				m(2,2)=z;
				return m;
			}
			
			//Get the rotation matrix
			static Mat44 rotationMatrix(const Complex& q)
			{
				float c = q.getReal();
				float s = q.getImaginary();
				
				Mat44 m;
				
				m(0,0)=c;		m(0,1)=s;		m(0,2)=0.0f;	m(0,3)=0.0f;
				m(1,0)=-s;		m(1,1)=c;		m(1,2)=0.0f;	m(1,3)=0.0f;
				m(2,0)=0.0f;	m(2,1)=0.0f;	m(2,2)=1.0f;	m(2,3)=0.0f;
				m(3,0)=0.0f;	m(3,1)=0.0f;	m(3,2)=0.0f;	m(3,3)=1.0f;
				
				return m;
			}
			
			static Mat44 rotationMatrix(const Quaternion& q)
			{
				float x = q.getX();
				float y = q.getY();
				float z = q.getZ();
				float w = q.getW();
				
				Mat44 m;
				
				m(0,0)=1.0f-(2.0f*y*y+2.0f*z*z);	m(0,1)=2.0f*x*y+2.0f*z*w;			m(0,2)=2.0f*x*z-2.0f*y*w;			m(0,3)=0.0f;
				m(1,0)=2.0f*x*y-2.0f*z*w;			m(1,1)=1.0f-(2.0f*x*x+2.0f*z*z);	m(1,2)=2.0f*y*z+2.0f*x*w;			m(1,3)=0.0f;
				m(2,0)=2.0f*x*z+2.0f*y*w;			m(2,1)=2.0f*y*z-2.0f*x*w;			m(2,2)=1.0f-(2.0f*x*x+2.0f*y*y);	m(2,3)=0.0f;
				m(3,0)=0.0f;						m(3,1)=0.0f;						m(3,2)=0.0f;						m(3,3)=1.0f;
				
				return m;
			}
			
			static Mat44 rotationMatrix(const Vec<2>& fwdAxis)
			{
				Vec<2> f = fwdAxis;
				normalize(f);
				
				Mat44 rotationMatrix;

				rotationMatrix(0,0) = f(0);
				rotationMatrix(0,1) = f(1);
				rotationMatrix(1,0) = -f(0);
				rotationMatrix(0,1) = f(1);
				
				rotationMatrix(2,2) = 1.0f;
				rotationMatrix(3,3) = 1.0f;
				return rotationMatrix;
			}
			
			static Mat44 rotationMatrix(const Vec<3>& leftAxis, const Vec<3>& fwdAxis){
				Vec<3> l = leftAxis;
				normalize(l);
				Vec<3> f = fwdAxis;
				if(l.dot(fwdAxis) != 0)
				{
					f = f - (f.dot(l))*l;
				}
				normalize(f);
				Vec<3> u = cross(l, f);
				
				Mat44 rotationMatrix;
				for(int j = 0; j<3; j++)
				{
					rotationMatrix(0,j) = l(j);
					rotationMatrix(1,j) = f(j);
					rotationMatrix(2,j) = u(j);
				}
				rotationMatrix(3,3) = 1.0f;
				return rotationMatrix;
			}
			
			Complex toComplex() const
			{
				return Complex(get(0,0), get(0,1));//Much simpler than the Quaternion!
			}
			
			Quaternion toQuaternion() const
			{
				float trace = get(0, 0) + get(1, 1) + get(2, 2);
				
				float x,y,z,w;
				
				if(trace > 0)
				{
					float s = 0.5f / (float)std::sqrt(trace+ 1.0f);
					w = 0.25f / s;
					x = (get(1, 2) - get(2, 1)) * s;
					y = (get(2, 0) - get(0, 2)) * s;
					z = (get(0, 1) - get(1, 0)) * s;
				}
				else
				{
					if(get(0, 0) > get(1, 1) && get(0, 0) > get(2, 2))
					{
						float s = 2.0f * (float)std::sqrt(1.0f + get(0, 0) - get(1, 1) - get(2, 2));
						w = (get(1, 2) - get(2, 1)) / s;
						x = 0.25f * s;
						y = (get(1, 0) + get(0, 1)) / s;
						z = (get(2, 0) + get(0, 2)) / s;
					}
					else if(get(1, 1) > get(2, 2))
					{
						float s = 2.0f * (float)std::sqrt(1.0f + get(1, 1) - get(0, 0) - get(2, 2));
						w = (get(2, 0) - get(0, 2)) / s;
						x = (get(1, 0) + get(0, 1)) / s;
						y = 0.25f * s;
						z = (get(2, 1) + get(1, 2)) / s;
					}
					else
					{
						float s = 2.0f * (float)std::sqrt(1.0f + get(2, 2) - get(0, 0) - get(1, 1));
						w = (get(0, 1) - get(1, 0) ) / s;
						x = (get(2, 0) + get(0, 2) ) / s;
						y = (get(1, 2) + get(2, 1) ) / s;
						z = 0.25f * s;
					}
				}

				float length = (float)std::sqrt(x * x + y * y + z * z + w * w);
				x /= length;
				y /= length;
				z /= length;
				w /= length;
				
				return Quaternion(w,x,y,z);
			}
			
			static Mat44 identity()
			{
				Mat44 m;
				m.setIdentity();
				return m;
			}
			
			static Mat44 orthoProjection(const float& width, const float& height, const float& nearField, const float& farField)
			{
				float left=-width/2.0f;
				float right=width/2.0f;
				float bottom=-height/2.0f;
				float top=height/2.0f;
				Mat44 m;
				m(0,0)=2.0f/width;
				m(1,1)=2.0f/height;
				m(2,2)=2.0f/(farField-nearField);
				m(3,3)=1.0f;
				
				m(0,3)=-(right+left)/(right-left);
				m(1,3)=-(top+bottom)/(top-bottom);
				m(2,3)=-(farField+nearField)/(farField-nearField);
				
				return m;
			}
			
			static Mat44 perspectiveProjection(const float& fieldOfView, const float& aspectRatio, const float& nearField, const float& farField)
			{
				Mat44 m;
				float tanHalfFov=std::tan(fieldOfView*3.141593/(2.0*180));
				float range = farField-nearField;
				m(0,0) = 1/(aspectRatio*tanHalfFov);
				m(1,1) = 1/(tanHalfFov);
				m(2,2) = (nearField+farField)/range;
				m(2,3) = -2*farField*nearField/range;
				m(3,2) = 1;
				return m;
			}
	};
	
	struct Transform3
	{
		public:
			Transform3(const Position3& position0 = Position3(REF_FRAME::LOCAL, Vec3(0.0f, 0.0f, 0.0f)), const Quaternion& rotation = Quaternion(), const Direction3& scale = Direction3(REF_FRAME::LOCAL, Vec3(1.0f, 1.0f, 1.0f)));
			
			/*Transform3():Transform3(Position3(REF_FRAME::LOCAL, Vec3(0.0f, 0.0f, 0.0f)), Quaternion(), Direction3(REF_FRAME::LOCAL, Vec3(1.0f, 1.0f, 1.0f)))
			{}*/		
			
			Transform3(const Transform3& transform);
			
			Transform3& operator=(const Transform3& transform);
			
			Position3 getFramePosition(const REF_FRAME& reference = REF_FRAME::LOCAL) const;
			
			Quaternion getFrameRotation() const;
			
			Direction3 getFrameScale() const;
			
			void setFramePosition(const Position3& newPosition);
			
			void setFramePosition(const float& x, const float& y, const float& z);
			
			void setFrameRotation(const Vec2& fwdAxis);
			
			void setFrameRotation(const Vec3& leftAxis, const Vec3& fwdAxis);
			
			void setFrameScale(const Direction3& newScale);
			
			void setFrameScale(const float& x, const float& y, const float& z);
			
			void rotate(const Direction3& eulerAngles);
			
			void translate(const Direction3& v);
			
			void resize(const float& scaleRatio);
			
			Mat<4> getTransformMatrix() const;
			
			void attachTo(const std::shared_ptr<Transform3>& parent);
			
			void detachFrom(const std::shared_ptr<Transform3>& parent);
			
			Direction3 toLocalFrame(const Direction3& direction);
			
			Direction3 toGlobalFrame(const Direction3 direction);
			
			Position3 toLocalFrame(const Position3& pos);
			
			Position3 toGlobalFrame(const Position3 pos);
			
			std::weak_ptr<Transform3> getParent();

		private:
			/** Attitude stored in LOCAL reference frame */
			Position3 position;
			Quaternion rotation;
			Direction3 scale;
			std::weak_ptr<Transform3> m_parent;
			
			Direction3 inThisFrame(const Direction3& direction);
			
			Direction3 inParentFrame(const Direction3& direction);
			
			Position3 inThisFrame(const Position3& pos);
			
			Position3 inParentFrame(const Position3& pos);
	};
	
	
	/*struct ITransform
	{
		virtual Mat<4> getTransformMatrix() const = 0;
	};
	
	template<std::size_t N>
	struct Transform : public ITransform
	{
		public:
			Transform(const Position<N>& position0, const Quaternion& rotation, const Direction<N>& scale):position(position0), rotation(rotation), scale(scale)
			{
				setFramePosition(position0);
			}
			
			Transform():Transform(Position<N>(REF_FRAME::LOCAL, Vec<N>(0.0f)), Quaternion(), Direction<N>(REF_FRAME::LOCAL, Vec<N>(1.0f)))
			{}		
			
			Transform(const Transform<N>& transform):position(transform.position), rotation(transform.rotation), scale(transform.scale)
			{}
			
			Transform<N> &operator=(const Transform<N>& transform)
			{
				position=transform.position;
				rotation=transform.rotation;
				scale=transform.scale;
				return *this;
			}
			
			Position<N> getFramePosition(const REF_FRAME& reference = REF_FRAME::LOCAL) const
			{
				if(reference == REF_FRAME::LOCAL)
				{
					if(std::shared_ptr<Transform<N>>  parent = m_parent.lock())
					{
						return parent->toGlobalFrame(position);
					}
				}
				return position;
			}
			
			Quaternion getFrameRotation() const
			{
				return rotation;
			}
			
			Direction<N> getFrameScale() const
			{
				return scale;
			}
			
			void setFramePosition(const Position<N>& newPosition)
			{				
				if(newPosition.getFrame() == REF_FRAME::LOCAL)
				{
					if(std::shared_ptr<Transform<N>>  parent = m_parent.lock())
					{
						position = parent->toLocalFrame(newPosition);
						return;
					}
				}
				position = newPosition;			
			}
			
			void setFrameRotation(const Vec<2>& fwdAxis)
			{
				rotation = (Mat44::rotationMatrix(fwdAxis)).toQuaternion();
			}
			
			void setFrameRotation(const Vec<3>& leftAxis, const Vec<3>& fwdAxis)
			{
				rotation = (Mat44::rotationMatrix(leftAxis, fwdAxis)).toQuaternion();
			}
			
			void setFrameScale(const Direction<N>& newScale)
			{
				scale = newScale;
			}			
			
			void rotate(const Direction<N>& eulerAngles)
			{
				if(eulerAngles.getFrame() == REF_FRAME::LOCAL)
				{
					rotation.addAngle(eulerAngles(0), eulerAngles(1), eulerAngles(2), true);
				}
				else
				{
					rotation.addAngle(eulerAngles(0), eulerAngles(1), eulerAngles(2), false);
				}
			}
			
			void translate(const Direction<N>& v)
			{
				Direction<N> v2 = v;
				if(v.getFrame() == REF_FRAME::LOCAL)
				{
					v2 = toLocalFrame(v);
				}
				position += v2;
			}
			
			void resize(const float& scaleRatio)
			{
				scale*=scaleRatio;
			}
			
			Mat<4> getTransformMatrix() const
			{
				Mat<4> res = (Mat44::positionMatrix(position)) * (Mat44::rotationMatrix(rotation)) * (Mat44::scaleMatrix(scale.getValue()));
				if(std::shared_ptr<Transform<N>> parent = m_parent.lock())
				{
					return parent->getTransformMatrix()*res;
				}
				else
				{
					return res;
				}
			}
			
			void attachTo(const std::shared_ptr<Transform<N>>& parent)
			{
				position = parent->toLocalFrame(position);
				m_parent = parent;
			}
			
			void detachFrom(const std::shared_ptr<Transform<N>>& parent)
			{
				position = parent->toGlobalFrame(position);
				m_parent.reset();
			}
			
			Direction<N> toLocalFrame(const Direction<N>& direction)
			{
				if(std::shared_ptr<Transform<N>>  parent = m_parent.lock())
				{
					return inThisFrame(parent->toLocalFrame(direction));
				}
				else
				{
					return inThisFrame(direction);
				}
			}
			
			Direction<N> toGlobalFrame(const Direction<N> direction)
			{
				if(std::shared_ptr<Transform<N>>  parent = m_parent.lock())
				{
					return parent->toGlobalFrame(inParentFrame(direction));
				}
				else
				{
					return inParentFrame(direction);
				}
			}
			
			Position<N> toLocalFrame(const Position<N>& pos)
			{
				if(std::shared_ptr<Transform<N>>  parent = m_parent.lock())
				{
					return inThisFrame(parent->toLocalFrame(pos));
				}
				else
				{
					return inThisFrame(pos);
				}
			}
			
			Position<N> toGlobalFrame(const Position<N> pos)
			{
				if(std::shared_ptr<Transform<N>>  parent = m_parent.lock())
				{
					return parent->toGlobalFrame(inParentFrame(pos));
				}
				else
				{
					return inParentFrame(pos);
				}
			}
			
			std::weak_ptr<Transform<N>> getParent()
			{
				return m_parent;
			}

		private:
			/** Attitude stored in LOCAL reference frame */
		/*	Position3 position;
			Quaternion rotation;
			Direction3 scale;
			std::weak_ptr<Transform3> m_parent;
			
			Direction3 inThisFrame(const Direction3& direction)
			{
				Direction3 newDirection(REF_FRAME::LOCAL, direction.getValue().rotateInv(rotation));
				return newDirection;
			}
			
			Direction3 inParentFrame(const Direction3& direction)
			{
				Direction3 newDirection = Direction3(REF_FRAME::GLOBAL, direction.getValue().rotate(rotation));
				return newDirection;
			}
			
			Position3 inThisFrame(const Position3& pos)
			{
				Position3 newPosition(REF_FRAME::LOCAL, (pos.getValue()-position.getValue()).rotate(rotation));
				return newPosition;
			}
			
			Position3 inParentFrame(const Position3& pos)
			{
				
				Position3 newPosition = Position3(REF_FRAME::GLOBAL, pos.getValue().rotateInv(rotation) + position.getValue());
				return newPosition;
			}
	};*/
	
	
	
	
	template<std::size_t N, typename Real>
	std::ostream& operator<<(std::ostream &strm, const Vec<N, Real> &v)
	{
		strm << "[" << v(0);
		for(int i=1; i<N;++i){
			strm << ", " << v(i);
		}
		return strm << "]";
	}
	
	std::ostream& operator<<(std::ostream &strm, const Position2 &v);
	
	std::ostream& operator<<(std::ostream &strm, const Direction2 &v);

	std::ostream& operator<<(std::ostream &strm, const Position3 &v);
	
	std::ostream& operator<<(std::ostream &strm, const Direction3 &v);
	
	template<std::size_t N, std::size_t M, typename Real>
	std::ostream& operator<<(std::ostream &strm, const Mat<N, M, Real> &m)
	{
		std::cout.precision(3);
		for(int i=0; i<N;++i){
			strm << "|" << m(i,0);
			for(int j=1; j<M; ++j)
			{
				strm << ",\t" << m(i,j);
			}
			strm << "|" << std::endl;
		}
		return strm;
	}
	
}