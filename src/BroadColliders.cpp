#include "BroadColliders.h"

namespace Circe
{
	AABB::AABB()
	{}

	AABB::AABB(const Vec<DIMENSION>& center, const Vec<DIMENSION>& width)
			: center(center)
	{
		Real maxWidth = width(0)*1.5;
		for(int i = 1; i<DIMENSION; i++)
		{
			maxWidth = std::max(maxWidth, width(i)*1.5);
		}
		for(int i = 0; i<DIMENSION; i++)
		{
			halfWidth = maxWidth;
		}
	}

	void AABB::refit(const AABB& v1, const AABB& v2)
	{
		for(int i = 0; i<DIMENSION; i++)
		{
			Real minBound = std::min(
					v1.center(i)-v1.halfWidth(i)
					+std::min(v1.margin(i),0.0),
					v2.center(i)-v2.halfWidth(i)
					+std::min(v2.margin(i),0.0));
			Real maxBound = std::max(
					v1.center(i)+v1.halfWidth(i)
					+std::max(v1.margin(i),0.0),
					v2.center(i)+v2.halfWidth(i)
					+std::max(v2.margin(i),0.0));

			center(i)    = (maxBound + minBound) * 0.5;
			halfWidth(i) = (maxBound - minBound) * 0.5;
		}
	}

	bool AABB::isInside(const AABB& v) const
	{
		for(int i = 0; i<DIMENSION; i++)
		{
			if(std::abs(center(i)-v.center(i))+halfWidth(i) 
						> v.halfWidth(i)) return false;
		}
		return true;
	}

	/** AABB to AABB intersection **/ 
	bool AABB::intersects(const AABB& other) const
	{
		for(int i = 0; i<DIMENSION; i++)
			if(std::abs(center(i)-other.center(i)) > 
			   halfWidth(i) + other.halfWidth(i)) return false;
		
		return true;
	}

	Real AABB::getUnionArea(const AABB& s1) const
	{
		AABB aabb(center, halfWidth);
		aabb.refit(aabb, s1);
		return aabb.getArea();
	}

	Real AABB::getArea() const
	{
		Real A = 1.0;
		for(int i = 0; i<DIMENSION; i++)
			A *= 2.0*halfWidth(i);

		return A;
	}

	void AABB::setPosition(const Vec<DIMENSION>& position)
	{
		center = position;
		update();
	}

	void AABB::setSize(const Vec<DIMENSION>& width)
	{
		Real maxWidth = width(0)*1.5;
		for(int i = 1; i<DIMENSION; i++)
		{
			maxWidth = std::max(maxWidth, width(i)*1.5);
		}
		for(int i = 0; i<DIMENSION; i++)
		{
			halfWidth = maxWidth;
		}
	}

	void AABB::setMargin(const Vec<DIMENSION>& v)
	{
		margin = v;
	}


	Point::Point(const Vec<DIMENSION>& location) : position(location)
	{}

	bool Point::intersects(const AABB& box) const
	{
		for(int i = 0; i<DIMENSION; i++)
			if(std::abs(box.center(i)-position(i)) > box.halfWidth(i))
				return false;

		return true;
	}

	Ray::Ray(const Vec<DIMENSION>& origin, const Vec<DIMENSION>& direction)
		: origin(origin), dir(direction)
	{}

	bool Ray::intersects(const AABB& box) const
	{
		bool inside = true;
		int quadrant[DIMENSION];
		int whichPlane;
		Vec<DIMENSION> maxT;
		Vec<DIMENSION> candidatePlane;
		Vec<DIMENSION> minB = box.center-box.halfWidth;
		Vec<DIMENSION> maxB = box.center+box.halfWidth;
		Vec<DIMENSION> coord;

		for(int i = 0; i<DIMENSION; i++)
		{
			if(origin(i) < minB(i))
			{
				quadrant[i] = 1;
				candidatePlane(i) = minB(i);
				inside = false;
			}
			else if(origin(i) > maxB(i))
			{
				quadrant[i] = 0;
				candidatePlane(i) = maxB(i);
				inside = false;
			}
			else
			{
				quadrant[i] = 2;
			}	
		}

		if(inside)
		{
			return true;
		}

		for(int i = 0; i<DIMENSION; i++)
		{
			if(quadrant[i] != 2 && dir(i) != 0.0)
				maxT(i) = (candidatePlane(i) - origin(i)) / dir(i);
			else
				maxT(i) = -1.0;
		}

		whichPlane = 0;
		for(int i = 1; i<DIMENSION; i++)
			if(maxT(whichPlane) < maxT(i))
				whichPlane = i;

		if(maxT(whichPlane) < 0.0) return false;
		for(int i = 0; i<DIMENSION; i++)
		{
			if(whichPlane != i)
			{
				coord(i) = origin(i) + maxT(whichPlane) * dir(i);
				if(coord(i) < minB(i) || coord(i) > maxB(i))
					return false;
			}
			else
			{
				coord(i) = candidatePlane(i);
			}
		}

		return true;
	}

	Segment::Segment(const Vec<DIMENSION>& pStart, 
					 const Vec<DIMENSION>& pEnd)
		: p0(pStart), p1(pEnd)
	{}

	bool Segment::intersects(const AABB& box) const
	{
		Vec<DIMENSION> minB = box.center-box.halfWidth;
		Vec<DIMENSION> maxB = box.center+box.halfWidth;

		for(int i = 0; i<DIMENSION; i++)
		{
			if(p0(i) < minB(i) && p1(i) < minB(i)) return false;
			if(p0(i) > maxB(i) && p1(i) > maxB(i)) return false;
		}

		return Ray(p0, Circe::normalize(p1-p0)).intersects(box);
	}



}
