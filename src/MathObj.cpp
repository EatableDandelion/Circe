#include "MathObj.h"
#include <math.h>

namespace Circe{
	
	
	Position2::Position2(const REF_FRAME& frame, const float& x, const float& y):frame(frame), x(x), y(y)
	{}
	
	Position2::Position2(const REF_FRAME& frame, const Vec2 vector):frame(frame), x(vector(0)), y(vector(1))
	{}
	
	Position2::Position2(const Position2& other):frame(other.frame), x(other(0)), y(other(1))
	{}
	
	Position2& Position2::operator=(const Position2& other)
	{
		x = other.x;
		y = other.y;
		frame = other.frame;
		return *this;
	}
	
	Direction2 Position2::operator-(const Position2& p2)
	{
		return Direction2(frame, Vec2(x-p2.x, y-p2.y));
	}
	
	Position2 Position2::operator+(const Direction2& d)
	{
		return Position2(frame, x+d(0), y+d(1));
	}
	
	Position2 Position2::operator*(const float& f)
	{
		return Position2(frame, x*f, y*f);
	}
	
	void Position2::operator+=(const Direction2& d)
	{
		x+=d(0);
		y+=d(1);
	}
	
	Position2 Position2::operator-(const Direction2& d)
	{
		return Position2(frame, x-d(0), y-d(1));
	}
	
	void Position2::operator-=(const Direction2& d)
	{
		x-=d(0);
		y-=d(1);
	}
	
	float& Position2::operator()(const unsigned int& index)
	{
		if(index == 0)
		{
			return x;
		}
		else
		{
			return y;
		}
	}
	
	float Position2::operator()(const unsigned int& index) const
	{
		if(index == 0)
		{
			return x;
		}
		else
		{
			return y;
		}
	}
	
	float Position2::distance2(const Position2& other)
	{
		return x*other.x+y*other.y;
	}
	
	float Position2::distance(const Position2& other)
	{
		return sqrt(distance2(other));
	}
	
	REF_FRAME Position2::getFrame() const
	{
		return frame;
	}
	
	Vec2 Position2::getValue() const
	{
		return Vec2(x, y);
	}
	
	
	
	Position3::Position3(const REF_FRAME& frame, const float& x, const float& y, const float& z):frame(frame), x(x), y(y), z(z)
	{}
	
	Position3::Position3(const REF_FRAME& frame, const Vec3 vector):frame(frame), x(vector(0)), y(vector(1)), z(vector(2))
	{}
	
	Position3::Position3(const Position3& other):frame(other.frame), x(other.x), y(other.y), z(other.z)
	{}
	
	Position3& Position3::operator=(const Position3& other)
	{
		x = other.x;
		y = other.y;
		z = other.z;
		frame = other.frame;
		return *this;
	}
	
	Direction3 Position3::operator-(const Position3& p2)
	{
		return Direction3(frame, Vec3(x-p2.x, y-p2.y, z-p2.z));
	}
	
	Position3 Position3::operator+(const Direction3& d)
	{
		return Position3(frame, x+d(0), y+d(1), z+d(2));
	}
	
	Position3 Position3::operator*(const float& f)
	{
		return Position3(frame, x*f, y*f, z*f);
	}
	
	void Position3::operator+=(const Direction3& d)
	{
		x+=d(0);
		y+=d(1);
		z+=d(2);
	}
	
	Position3 Position3::operator-(const Direction3& d)
	{
		return Position3(frame, x-d(0), y-d(1), z-d(2));
	}
	
	void Position3::operator-=(const Direction3& d)
	{
		x-=d(0);
		y-=d(1);
		z-=d(2);
	}
	
	float& Position3::operator()(const unsigned int& index)
	{
		if(index == 0)
		{
			return x;
		}
		else if(index == 1)
		{
			return y;
		}
		else
		{
			return z;
		}
	}
	
	float Position3::operator()(const unsigned int& index) const
	{
		if(index == 0)
		{
			return x;
		}
		else if(index == 1)
		{
			return y;
		}
		else
		{
			return z;
		}
	}
	
	float Position3::distance2(const Position3& other)
	{
		return x*other.x+y*other.y;
	}
	
	float Position3::distance(const Position3& other)
	{
		return sqrt(distance2(other));
	}
	
	REF_FRAME Position3::getFrame() const
	{
		return frame;
	}
	
	Vec3 Position3::getValue() const
	{
		return Vec3(x, y, z);
	}
	
	
	Direction2::Direction2(const REF_FRAME& frame, const float& x, const float& y):frame(frame), x(x), y(y)
	{}
	
	Direction2::Direction2(const REF_FRAME& frame, const Vec2 vector):frame(frame), x(vector(0)), y(vector(1))
	{}
	
	Direction2::Direction2(const Direction2& other):frame(other.frame), x(other.x), y(other.y)
	{}
	
	Direction2& Direction2::operator=(const Direction2& other)
	{
		x = other.x;
		y = other.y;
		frame = other.frame;
		return *this;
	}
	
	Direction2 Direction2::operator+(const Direction2& d2)
	{
		return Direction2(frame, x+d2.x, y+d2.y);
	}
	
	Direction2 Direction2::operator-(const Direction2& d2)
	{
		return Direction2(frame, x-d2.x, y-d2.y);
	}
	
	void Direction2::operator*=(const float& f)
	{
		x*=f;
		y*=f;
	}
	
	Position2 Direction2::operator+(const Position2& d2)
	{
		return Position2(frame, x+d2(0), y+d2(1));
	}
	
	Position2 Direction2::operator-(const Position2& d2)
	{
		return Position2(frame, x-d2(0), y-d2(1));
	}
	
	float& Direction2::operator()(const unsigned int& index)
	{
		if(index == 0)
		{
			return x;
		}
		else
		{
			return y;
		}
	}
	
	float Direction2::operator()(const unsigned int& index) const
	{
		if(index == 0)
		{
			return x;
		}
		else
		{
			return y;
		}
	}
	
	REF_FRAME Direction2::getFrame() const
	{
		return frame;
	}
	
	Vec2 Direction2::getValue() const
	{
		return Vec2(x, y);
	}	


	Direction3::Direction3(const REF_FRAME& frame, const float& x, const float& y, const float& z):frame(frame), x(x), y(y), z(z)
	{}
	
	Direction3::Direction3(const REF_FRAME& frame, const Vec3 vector):frame(frame), x(vector(0)), y(vector(1)), z(vector(2))
	{}
	
	Direction3::Direction3(const Direction3& other):frame(other.frame), x(other.x), y(other.y), z(other.z)
	{}
	
	Direction3& Direction3::operator=(const Direction3& other)
	{
		x = other.x;
		y = other.y;
		z = other.z;
		frame = other.frame;
		return *this;
	}
	
	Direction3 Direction3::operator+(const Direction3& d2)
	{
		return Direction3(frame, x+d2.x, y+d2.y, z+d2.z);
	}
	
	Direction3 Direction3::operator-(const Direction3& d2)
	{
		return Direction3(frame, x-d2.x, y-d2.y, z-d2.z);
	}
	
	void Direction3::operator*=(const float& f)
	{
		x*=f;
		y*=f;
		z*=f;
	}
	
	Position3 Direction3::operator+(const Position3& d2)
	{
		return Position3(frame, x+d2(0), y+d2(1), z+d2(2));
	}
	
	Position3 Direction3::operator-(const Position3& d2)
	{
		return Position3(frame, x-d2(0), y-d2(1), z-d2(2));
	}
	
	float& Direction3::operator()(const unsigned int& index)
	{
		if(index == 0)
		{
			return x;
		}
		else if(index == 1)
		{
			return y;
		}
		else
		{
			return z;
		}
	}
	
	float Direction3::operator()(const unsigned int& index) const
	{
		if(index == 0)
		{
			return x;
		}
		else if(index == 1)
		{
			return y;
		}
		else
		{
			return z;
		}
	}
	
	REF_FRAME Direction3::getFrame() const
	{
		return frame;
	}
	
	Vec3 Direction3::getValue() const
	{
		return Vec3(x, y, z);
	}
	
	

	Transform3::Transform3(const Position3& position0, const Quaternion& rotation, const Direction3& scale):position(position0), rotation(rotation), scale(scale)
	{
		setFramePosition(position0);
	}	
	
	Transform3::Transform3(const Transform3& transform):position(transform.position), rotation(transform.rotation), scale(transform.scale)
	{}
	
	Transform3& Transform3::operator=(const Transform3& transform)
	{
		position=transform.position;
		rotation=transform.rotation;
		scale=transform.scale;
		return *this;
	}
	
	Position3 Transform3::getFramePosition(const REF_FRAME& reference) const
	{
		if(reference == REF_FRAME::LOCAL)
		{
			if(std::shared_ptr<Transform3>  parent = m_parent.lock())
			{
				return parent->toGlobalFrame(position);
			}
		}
		return position;
	}
	
	Quaternion Transform3::getFrameRotation() const
	{
		return rotation;
	}
	
	Direction3 Transform3::getFrameScale() const
	{
		return scale;
	}
	
	void Transform3::setFramePosition(const Position3& newPosition)
	{				
		if(newPosition.getFrame() == REF_FRAME::LOCAL)
		{
			if(std::shared_ptr<Transform3>  parent = m_parent.lock())
			{
				position = parent->toLocalFrame(newPosition);
				return;
			}
		}
		position = newPosition;			
	}
	
	void Transform3::setFramePosition(const float& x, const float& y, const float& z)
	{
		position(0) = x;
		position(1) = y;
		position(2) = z;
	}
	
	void Transform3::setFrameRotation(const Vec2& fwdAxis)
	{
		rotation = (Mat44::rotationMatrix(fwdAxis)).toQuaternion();
	}
	
	void Transform3::setFrameRotation(const Vec3& leftAxis, const Vec3& fwdAxis)
	{
		rotation = (Mat44::rotationMatrix(leftAxis, fwdAxis)).toQuaternion();
	}
	
	void Transform3::setFrameScale(const Direction3& newScale)
	{
		scale = newScale;
	}

	void Transform3::setFrameScale(const float& x, const float& y, const float& z)
	{
		scale(0) = x;
		scale(1) = y;
		scale(2) = z;
	}	
	
	void Transform3::rotate(const Direction3& eulerAngles)
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
	
	void Transform3::translate(const Direction3& v)
	{
		Direction3 v2 = v;
		if(v.getFrame() == REF_FRAME::LOCAL)
		{
			v2 = toLocalFrame(v);
		}
		position += v2;
	}
	
	void Transform3::resize(const float& scaleRatio)
	{
		scale*=scaleRatio;
	}
	
	Mat<4> Transform3::getTransformMatrix() const
	{
		Mat<4> res = (Mat44::positionMatrix(position)) * (Mat44::rotationMatrix(rotation)) * (Mat44::scaleMatrix(scale.getValue()));
		if(std::shared_ptr<Transform3> parent = m_parent.lock())
		{
			return parent->getTransformMatrix()*res;
		}
		else
		{
			return res;
		}
	}
	
	void Transform3::attachTo(const std::shared_ptr<Transform3>& parent)
	{
		position = parent->toLocalFrame(position);
		m_parent = parent;
	}
	
	void Transform3::detachFrom(const std::shared_ptr<Transform3>& parent)
	{
		position = parent->toGlobalFrame(position);
		m_parent.reset();
	}
	
	Direction3 Transform3::toLocalFrame(const Direction3& direction)
	{
		if(std::shared_ptr<Transform3>  parent = m_parent.lock())
		{
			return inThisFrame(parent->toLocalFrame(direction));
		}
		else
		{
			return inThisFrame(direction);
		}
	}
	
	Direction3 Transform3::toGlobalFrame(const Direction3 direction)
	{
		if(std::shared_ptr<Transform3>  parent = m_parent.lock())
		{
			return parent->toGlobalFrame(inParentFrame(direction));
		}
		else
		{
			return inParentFrame(direction);
		}
	}
	
	Position3 Transform3::toLocalFrame(const Position3& pos)
	{
		if(std::shared_ptr<Transform3>  parent = m_parent.lock())
		{
			return inThisFrame(parent->toLocalFrame(pos));
		}
		else
		{
			return inThisFrame(pos);
		}
	}
	
	Position3 Transform3::toGlobalFrame(const Position3 pos)
	{
		if(std::shared_ptr<Transform3>  parent = m_parent.lock())
		{
			return parent->toGlobalFrame(inParentFrame(pos));
		}
		else
		{
			return inParentFrame(pos);
		}
	}
	
	std::weak_ptr<Transform3> Transform3::getParent()
	{
		return m_parent;
	}

	Direction3 Transform3::inThisFrame(const Direction3& direction)
	{
		Direction3 newDirection(REF_FRAME::LOCAL, direction.getValue().rotateInv(rotation));
		return newDirection;
	}
	
	Direction3 Transform3::inParentFrame(const Direction3& direction)
	{
		Direction3 newDirection = Direction3(REF_FRAME::GLOBAL, direction.getValue().rotate(rotation));
		return newDirection;
	}
	
	Position3 Transform3::inThisFrame(const Position3& pos)
	{
		Position3 newPosition(REF_FRAME::LOCAL, (pos.getValue()-position.getValue()).rotate(rotation));
		return newPosition;
	}
	
	Position3 Transform3::inParentFrame(const Position3& pos)
	{
		
		Position3 newPosition = Position3(REF_FRAME::GLOBAL, pos.getValue().rotateInv(rotation) + position.getValue());
		return newPosition;
	}
	
	
	std::ostream& operator<<(std::ostream &strm, const Position2 &v)
	{
		std::string frame = "(Global position)";
		if(v.getFrame() == REF_FRAME::LOCAL)
		{
			frame = "(Local position)";
		}
		return strm << v.getValue() << frame;
	}
	
	std::ostream& operator<<(std::ostream &strm, const Direction2 &v)
	{
		std::string frame = "(Global direction)";
		if(v.getFrame() == REF_FRAME::LOCAL)
		{
			frame = "(Local direction)";
		}
		return strm << v.getValue() << frame;
	}
	
	std::ostream& operator<<(std::ostream &strm, const Position3 &v)
	{
		std::string frame = "(Global position)";
		if(v.getFrame() == REF_FRAME::LOCAL)
		{
			frame = "(Local position)";
		}
		return strm << v.getValue() << frame;
	}
	
	std::ostream& operator<<(std::ostream &strm, const Direction3 &v)
	{
		std::string frame = "(Global direction)";
		if(v.getFrame() == REF_FRAME::LOCAL)
		{
			frame = "(Local direction)";
		}
		return strm << v.getValue() << frame;
	}
	
}