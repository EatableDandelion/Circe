#include "Math.h"

namespace Circe
{
	
	Mat::Mat()
	{}

	Mat::Mat(const std::size_t newM, const std::size_t newN)
	{	
		init(newM, newN);	
	}

	void Mat::init(const std::size_t newM, const std::size_t newN)
	{
		M = newM;
		N = newN;
		data.clear();
		for(int k = 0; k<M*N; k++)
			data.push_back(0.0);

		hasQR = false;
	}

	Mat::Mat(const std::size_t M, const std::size_t N,
		const std::initializer_list<Real>& values)
	{
		init(M,N);
		std::size_t k = 0;
		for(Real f : values)
		{
			data[k] = f;
			k++;
		}
		hasQR = false;
	}

	Mat::Mat(const std::initializer_list<Mat>& values)
	{
		for(Mat m : values)
		{
			M  = max(M,m.M);
			N += m.N;
		}

		init(M,N);

		int j0 = 0;
		for(Mat m : values)
		{
			for(int j = 0; j<m.N; j++)
			{
				for(int i = 0; i<m.M; i++)
				{
					(*this)(i,j0) = m(i,j);
				}
				j0++;
			}
		}
		hasQR = false;
	}
	
	
	Mat& Mat::operator=(const Mat& m)
	{
		init(m.M, m.N);
		data = m.data;	
		return *this;
	}
	

	//Matrix multiplication to matrix
	Mat Mat::operator*(const Mat& m) const
	{
		assert(N == m.M);

		Mat result(M,m.N);
		for(int i = 0; i<M; ++i)
		{
			for(int j = 0; j<m.N; ++j)
			{
				Real value = (Real)0.0;
				for(int k=0;k<N; ++k)
				{
					value+=(*this)(i,k)*m(k,j);
				}
				result(i,j)=value;
			}
		}	
		return result;
	}
		
	Mat Mat::transposeMult(const Mat& m) const
	{
		assert(M == m.M);

		Mat result(N,m.N);
		for(int i = 0; i<N; i++)
		{
			for(int j = 0; j<m.N; j++)
			{
				Real value = (Real)0.0;
				for(int k=0; k<M; k++)
				{
					value+=(*this)(k,i)*m(k,j);
				}
				result(i,j)=value;
			}
		}	
		return result;
	}
	
	Mat Mat::multTranspose(const Mat& m) const
	{
		assert(N == m.N);

		Mat result(M,m.M);
		for(int i = 0; i<M; ++i)
		{
			for(int j = 0; j<m.M; ++j)
			{
				Real value = (Real)0.0;
				for(int k=0;k<N; ++k)
				{
					value+=(*this)(i,k)*m(j,k);
				}
				result(i,j)=value;
			}
		}	
		return result;
	}	

	//Get the matrix multiplied by a scalar
	Mat Mat::operator*(const Real& f) const
	{
		Mat result(M,N);
		for(int k = 0; k<M*N; k++)
		{
			result.data[k] = data[k]*f;
		}	
		return result;	
	}

	Mat Mat::operator/(const Real& f) const
	{
		Mat result(M,N);
		for(int k = 0; k<M*N; k++)
		{
			result.data[k] = data[k]/f;
		}	
		return result;	
	}
	
	void Mat::operator+=(const Mat& m)
	{
		assert(M == m.M && N == m.N);
		hasQR = false;
		for(int k = 0; k<M*N; k++)
		{
			data[k]+=m.data[k];
		}	
	}
	
	Mat Mat::operator+(const Mat& m) const
	{
		assert(M == m.M && N == m.N);
		Mat result(M,N);
		for(int k = 0; k<M*N; k++)
		{
			result.data[k] = data[k] + m.data[k];
		}	
		return result;	
	}
		
	Mat Mat::operator-(const Mat& m) const
	{
		assert(M == m.M && N == m.N);
		Mat result(M,N);
		for(int k = 0; k<M*N; k++)
		{
			result.data[k] = data[k] - m.data[k];
		}	
		return result;	
	}

	void Mat::operator-=(const Mat& m)
	{
		assert(M == m.M && N == m.N);
		hasQR = false;
		for(int k = 0; k<M*N; k++)
		{
			data[k]-=m.data[k];
		}	
	}
	
	//Multiply by a scalar
	void Mat::operator*=(const Real& f)
	{
		hasQR = false;
		for(int k = 0; k<M*N; k++)
		{
			data[k]*=f;
		}
	}

	/** Solve the LU factorization of Ax = b with b = this */
	Mat Mat::operator/(Mat& A) const
	{
		return solveQR(A);

		/*
		assert(A.M == M);

		if(A.M == A.N)
		{
			return solveLU(A);
		}

		if(A.M > A.N)
		{
			return (A.transposeMult(*this))*
						solveLU(A.transposeMult(A));
		}
		else
		{
			return A.transposeMult(*this)*
						solveLU(A.multTranspose(A));
		}
		*/
	}

	Mat Mat::solveLU(const Mat& A)
	{
		assert(M == A.M && M == A.N);

		Mat LU(M,M);

		double sum = 0.0;

		for(int i = 0; i<M; i++)
		{
			for(int j = i; j<M; j++)
			{
				sum = 0.0;
				for(int k = 0;k<i;k++)
				{
					sum += LU(i,k)*LU(k,j);
				}
				LU(i,j)=A(i,j)-sum;
			}
			for(int j = i+1; j<M; j++)
			{
				sum = 0.0;
				for(int k = 0; k<i; k++)
				{
					sum += LU(j,k)*LU(k,i);
				}
				LU(j,i)=(A(j,i)-sum)/LU(i,i);
			}
		}

		Mat y(M,1);
		for(int i = 0; i<M; i++)
		{
			sum = 0.0;
			for(int k = 0;k<i;k++)
			{
				sum+=LU(i,k)*y(k);
			}
			y(i) = (*this)(i)-sum;
		}

		Mat x(M,1);
		for(int i = M-1;i>=0;i--)
		{
			sum=0.0;
			for(int k = i+1;k<M; k++)
			{
				sum+=LU(i,k)*x(k);
			}
			x(i) = (y(i)-sum)/LU(i,i);
		}
		return x;
	}
/*
	Mat Mat::solveGaussJordan(const Mat& A)
	{
		int n = A.M;

		Mat b = (*this);

		for(int k = 0; k<n-1; k++)
		{
	//		 Find the best row permutation 
			int ip = k;
			Real maxValue = abs(A(k,k));
			for(int j = k; j<n; j++)
			{
				if(abs(A(j,k)) > maxValue)
				{
					ip = j;
					maxValue = abs(A(j,k));
				}	
			}
			
	//		 Perform the permutation 
			for(int j = k; j<n; j++)
			{
				Real temp = A(k,j);
				A(k,j) = A(ip, j);
				A(ip,j) = temp;		
			}
			Real temp = b(k);
			b(k) = b(ip);
			b(ip) = temp;

	//		 Solve Gaussian elimination 
			for(int i = k+1; i<n; i++)
			{
				Real piv = A(i,k)/A(k,k);

				for(int j = k+1; j<n; j++)
				{
					A(i,j) -= piv * A(k,j);
				}
				b(i) -= piv * b(k);
			}
		}

		return backSolve(A, b);
	}*/

	
	void Mat::initQR()
	{
		Q = std::make_shared<Mat>(M,N);
		R = std::make_shared<Mat>(N,N);
		hasQR = true;

		for(int j = 0; j<N; j++)
		{
			Mat qj = getColumn(j);
			R->set(j,j,1.0);

			for(int i = 0; i<=j-1; i++)
			{
				Mat qi = Q->getColumn(i);

				Real rij = 0.0;

				for(int k = 0; k<M; k++)
					rij += qj(k) * qi(k);

				for(int k = 0; k<M; k++)
					qj(k) -= rij * qi(k);

				R->set(i,j,rij);
			}

			Real rjj = 0.0;
			for(int k = 0; k<M; k++)
				rjj += qj(k) * qj(k);

			rjj = std::sqrt(rjj);

			if(rjj > 1e-8)
			{
				for(int k = 0; k<M; k++)
					Q->set(k,j,qj(k)/rjj);

				R->set(j,j,rjj);
			}
		}
	}

	Mat Mat::backSolve(const Mat& A) const
	{
		Mat x(A.N,1);

		for(int i = A.N-1; i >= 0; i--)
		{
			Real t = (*this)(i);

			for(int j = i+1; j < A.N; j++)
			{
				t -= A(i,j) * x(j);
			}

			x(i) = t / A(i,i);
		}

		return x;
	}	

	Mat Mat::solveQR(Mat& A) const
	{
		if(!A.hasQR) A.initQR();

		return (A.Q->transposeMult(*this)).backSolve(*(A.R));
	}

	Mat Mat::getColumn(const int index) const
	{
		Mat c(M,1);

		for(int i = 0; i<M; i++)
			c(i,0) = (*this)(i,index);

		return c;
	}
			
	//Get the coefficient at row i and column j
	Real Mat::operator()(const int& i, const int& j) const
	{
		assert(i<M && j<N);
		return data[i*N+j];
	}
	
	//Set the coefficient at row i and column j
	Real& Mat::operator()(const int& i, const int& j)
	{
		assert(i<M && j<N);
		hasQR = false;
		return data[i*N+j];
	}

	void Mat::set(const int& i, const int& j, const Real& value)
	{
		assert(i<M && j<N);
		hasQR = false;
		data[i*N+j] = value;
	}

	std::size_t Mat::getNbRows() const
	{
		return M;
	}
	
	std::size_t Mat::getNbCols() const
	{
		return N;
	}

	Complex::Complex(const Real& angle):c(cos(angle)), s(sin(angle))
	{}
	
	Complex::Complex(const Real& real, const Real& imag):
					c(real), s(imag)
	{}

	Complex::Complex(const Mat& m):Complex(m(0,0), m(0,1))
	{}

	Real Complex::getAngle() const
	{
		return atan2(s,c);
	}
	
	Complex Complex::operator*(const Complex& q) const
	{
		Real a=c;Real b=s;
		Real c=q.c;Real d=q.s;
		return Complex(a*c-b*d, b*c+a*d);
	}
	
	Complex Complex::operator=(const Complex& q)
	{
		c=q.c;
		s=q.s;
		return *this;
	}
	
	Real Complex::length() const
	{
		return sqrt(c*c+s*s);
	}
	
	Complex Complex::normalize()
	{
		Real length = (*this).length();
		c/=length;
		s/=length;
		return *this;
	}
	
	Complex Complex::conjugate()
	{
		s*=-(Real)(1.0);
		return *this;
	}
	
	Real Complex::getReal() const
	{
		return c;
	}
	
	Real Complex::getImaginary() const
	{
		return s;
	}
	
	void Complex::addAngle(const Real& dtheta)
	{
		Complex q(dtheta);
		Real ctemp = c*q.c - s*q.s;
		s = c*q.s + s*q.c;
		c=ctemp;
	}


	Quaternion::Quaternion():w((Real)1.0), x(0.0), y(0.0), z(0.0)
	{}
	
	Quaternion::Quaternion(const Real& w, const Real& x, 
					  	   const Real& y, const Real& z):
							w(w), x(x), y(y), z(z)
	{}
	
	Quaternion::Quaternion(const Real& roll, const Real& pitch, 
						const Real& yaw)
	{
		Real cy = cos(yaw*0.5);
		Real sy = sin(yaw*0.5);
		Real cp = cos(pitch*0.5);
		Real sp = sin(pitch*0.5);
		Real cr = cos(roll*0.5);
		Real sr = sin(roll*0.5);
		
		w = cy*cp*cr+sy*sp*sr;
		x = cy*cp*sr-sy*sp*cr;
		y = sy*cp*sr+cy*sp*cr;
		z = sy*cp*cr-cy*sp*sr;
		
		//Source: Wikipedia, conversion between 
		// quaternion and Euler angles
	}

	Quaternion::Quaternion(const Mat& m)
	{
		Real trace = m(0, 0) + m(1, 1) + m(2, 2);
		
		if(trace > 0)
		{
			Real s = 0.5 / (Real)std::sqrt(trace+ 1.0);
			w = 0.25 / s;
			x = (m(1, 2) - m(2, 1)) * s;
			y = (m(2, 0) - m(0, 2)) * s;
			z = (m(0, 1) - m(1, 0)) * s;
		}
		else
		{
			if(m(0, 0) > m(1, 1) && m(0, 0) 
					> m(2, 2))
			{
				Real s = 2.0f * (Real)std::sqrt(1.0f + m(0, 0)
							- m(1, 1) - m(2, 2));
				w = (m(1, 2) - m(2, 1)) / s;
				x = 0.25f * s;
				y = (m(1, 0) + m(0, 1)) / s;
				z = (m(2, 0) + m(0, 2)) / s;
			}
			else if(m(1, 1) > m(2, 2))
			{
				Real s = 2.0f * (Real)std::sqrt(1.0f + m(1, 1)
						- m(0, 0) - m(2, 2));
				w = (m(2, 0) - m(0, 2)) / s;
				x = (m(1, 0) + m(0, 1)) / s;
				y = 0.25f * s;
				z = (m(2, 1) + m(1, 2)) / s;
			}
			else
			{
				Real s = 2.0f * (Real)std::sqrt(1.0f + m(2, 2)
						- m(0, 0) - m(1, 1));
				w = (m(0, 1) - m(1, 0) ) / s;
				x = (m(2, 0) + m(0, 2) ) / s;
				y = (m(1, 2) + m(2, 1) ) / s;
				z = 0.25f * s;
			}
		}

		Real length = -(Real)std::sqrt(x*x + y*y + z*z + w*w);
		x /= length;
		y /= length;
		z /= length;
		w /= length;
		
	}

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
			w*q.z + x*q.y - y*q.x + z*q.w  //z
		);
	}
	
	void Quaternion::operator*=(const Quaternion& q)
	{
		Real w0 = w*q.w - x*q.x - y*q.y - z*q.z;
		Real x0 = w*q.x + x*q.w + y*q.z - z*q.y;
		Real y0 = w*q.y - x*q.z + y*q.w + z*q.x;
		Real z0 = w*q.z + x*q.y - y*q.x + z*q.w;
		w = w0;
		x = x0;
		y = y0;
		z = z0;
	}
	
	void Quaternion::addAngle(const Real& roll, const Real& pitch, 
						 const Real& yaw)
	{
		
		Quaternion q = (*this)*Quaternion(0.0, 
										  -roll*0.5, 
										  -pitch*0.5, 
										  -yaw*0.5);
		w+=q.w;
		x+=q.x;
		y+=q.y;
		z+=q.z;
				
		normalize();
	}
	
	Real Quaternion::length() const
	{
		return sqrt(w*w+x*x+y*y+z*z);
	}
	
	Quaternion Quaternion::normalize()
	{
		Real length = (*this).length();
		w/=length;
		x/=length;
		y/=length;
		z/=length;

		return *this;
	}
	
	Quaternion Quaternion::conjugate()
	{
		x*=-1.0;
		y*=-1.0;
		z*=-1.0;

		return *this;
	}
	
	Quaternion Quaternion::getConjugate() const
	{
		return Quaternion(w, -x, -y, -z);
	}
	
	Vec<2> Quaternion::rotate(const Vec<2>& v2) const
	{
		Vec<3> r = rotate(Vec<3>(v2(0), v2(1), 0.0));
		return Vec<2>(r(0), r(1));
	}

	Vec<3> Quaternion::rotate(const Vec<3>& v2) const
	{
		if(v2.dot(v2) < 1e-10) return Vec<3>(v2(0),v2(1),v2(2));
		Quaternion v((Real)0.0, v2(0), v2(1), v2(2));
		Quaternion p(*this);
		p.normalize();
		Quaternion pConj=p.getConjugate();
		Quaternion result = p*v*pConj;
		Vec<3> v3;
		v3(0)=result.getX();
		v3(1)=result.getY();
		v3(2)=result.getZ();
		return v3;
	}
	
	Vec<2> Quaternion::rotateInv(const Vec<2>& v2) const
	{
		Vec<3> r = rotateInv(Vec<3>(v2(0), v2(1), 0.0));
		return Vec<2>(r(0), r(1));
	}

	Vec<3> Quaternion::rotateInv(const Vec<3>& v2) const
	{
		if(v2.dot(v2) < 1e-10) return Vec<3>(v2(0),v2(1),v2(2));
		Quaternion v((Real)0.0, v2(0), v2(1), v2(2));
		Quaternion p(*this);
		p.normalize();
		Quaternion pConj=p.getConjugate();
		Quaternion result = pConj*v*p;
		Vec<3> v3;
		v3(0)=result.getX();
		v3(1)=result.getY();
		v3(2)=result.getZ();
		return v3;
	}

	Vec<3> Quaternion::toEulerAngle() const
	{
		Real sinr_cosp = 2.0*(w*x+y*z);
		Real cosr_cosp = 1.0-2.0*(x*x+y*y);
		Real roll = (Real)(std::atan2(sinr_cosp, cosr_cosp));

		Real sinp = 2.0*(w*y-z*x);
		Real pitch = std::asin(sinp);
		if(std::abs(sinp) >= 1)
			pitch = std::copysign(3.14159/2.0, sinp);

		Real siny_cosp = 2.0*(w*z+x*y);
		Real cosy_cosp = 1.0-2.0*(y*y+z*z);
		Real yaw = (Real)(std::atan2(siny_cosp, cosy_cosp));

		return Vec<3>(roll, pitch, yaw);
	}

	Real Quaternion::angleWith(const Quaternion& other)
	{
		Vec<3> v1 = rotate(Vec<3>(1,0,0)); 
		Vec<3> v2 = other.rotate(Vec<3>(1,0,0)); 
		Vec<3> v3 = cross(v1,v2);
	
		if(dot(v3,Vec<3>(0,0,1)) > 0)
			return -std::acos(dot(v1,v2));
		
		return std::acos(dot(v1,v2));
	}

	void Quaternion::lookAt(const Vec<2>& lookAt)
	{
		Vec<2> fwd = Circe::normalize(lookAt);

		Mat m(3,3);

		m(0,0) = fwd(0);
		m(1,0) = fwd(1);
		m(0,1) = -fwd(1);
		m(1,1) = fwd(0);
		m(2,2) = 1.0;

		Quaternion q(m);

		w = q.w;
		x = q.x;
		y = q.y;
		z = q.z;
	}

	void Quaternion::lookAt(const Vec<3>& lookAt, const Vec<3>& up0)
	{
		Vec<3> fwd = Circe::normalize(lookAt);
		Vec<3> right = Circe::normalize(cross(up0, fwd));
		Vec<3> up = cross(fwd, right);

		Mat m(3,3);
		for(int i = 0; i<3; i++)
		{
			m(i,0) = fwd(i);
			m(i,1) = right(i);
			m(i,2) = up(i);
		}

		Quaternion q(m);

		w = q.w;
		x = q.x;
		y = q.y;
		z = q.z;
	}

	Real Quaternion::getW() const
	{
		return w;
	}
	
	Real Quaternion::getX() const
	{
		return x;
	}
	
	Real Quaternion::getY() const
	{
		return y;
	}
	
	Real Quaternion::getZ() const
	{
		return z;
	}
}
