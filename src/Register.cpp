#include "Register.h"

namespace Circe
{
	void Register::RegisterKey::addListener
						(const Listener<std::vector<Real>> listener)
	{
		if(!emitter)
		{
			emitter=std::make_shared<Emitter<std::vector<Real>>>();
		}

		emitter->addListener(listener);
	}

	void Register::RegisterKey::callListeners
								(const std::vector<Real>& values)
	{
		if(emitter)
		{
			std::vector<Real> vec(values.begin()+offset,
								  values.begin()+offset+size);
			emitter->broadcast(vec);
		}
	}

	void Register::setVariable(const std::string& name,
					const std::vector<Real> value)
	{
		std::size_t nameID = CIRCE_STRING_ID(name);

		if(hasVariable(nameID))
		{
			RegisterKey key = m_positions.at(nameID);
			unsigned int size = key.size;

			if(size != value.size()) 
				CIRCE_ERROR("The variable "
							+name
							+" has changed size.");

			unsigned int offset = key.offset;
			
			for(int i = 0; i < size; i++)
				m_values[offset + i] = value[i];

			key.callListeners(m_values);
		}
		else
		{
			RegisterKey position;

			position.offset = m_values.size();
			position.size 	= value.size();

			m_positions.insert(
				std::pair<std::size_t, RegisterKey>
						 (nameID, position));

			for(int i = 0; i < value.size(); i++)
				m_values.push_back(value[i]);
		}
	}

	void Register::setVariable(const std::string& name, const Real& value)
	{
		setVariable(name, std::vector<Real>({value}));	
	}

	void Register::setVariable(const std::string& name, 
					 const Circe::Vec2& value)
	{
		setVariable(name, std::vector<Real>({value(0), value(1)}));	
	}

	void Register::setVariable(const std::string& name, 
					 		   const Circe::Vec3& value)
	{
		setVariable(name, std::vector<Real>
					({value(0),value(1),value(2)}));
	}

	void Register::setVariable(const std::string& name, 
					 		   const Mat& value)
	{
		std::vector<Real> v(
				{value(0,0),value(0,1),value(0,2),value(0,3),
				 value(1,0),value(1,1),value(1,2),value(1,3),
				 value(2,0),value(2,1),value(2,2),value(2,3),
				 value(3,0),value(3,1),value(3,2),value(3,3)});

		setVariable(name, v);	
	}

	std::vector<Real> Register::operator()(const std::string& name) const
	{
		return getVariable(CIRCE_STRING_ID(name));
	}

	std::vector<Real> Register::getVariable(const std::size_t nameID) const
	{
		std::vector<Real> uniform;

		if(hasVariable(nameID))
		{	
			unsigned int size = m_positions.at(nameID).size;
			unsigned int offset = m_positions.at(nameID).offset;

			for(int i = 0; i < size; i++)
				uniform.push_back(m_values[offset + i]);
		}

		return uniform;
	}

	bool Register::hasVariable(const std::string& name) const
	{
		return hasVariable(CIRCE_STRING_ID(name));
	}

	bool Register::hasVariable(const std::size_t nameID) const
	{
		return m_positions.count(nameID);
	}

	void Register::addListener(const std::string& name,
					 		   const Listener<std::vector<Real>> lisnr)
	{
		if(hasVariable(name))
		{
			std::size_t nameID = CIRCE_STRING_ID(name);
			m_positions[nameID].addListener(lisnr);
			m_positions[nameID].callListeners(getVariable(nameID));
		}
	}
}
