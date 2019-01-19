#include "Utils.h"

namespace Circe
{
	size_t getId(const std::string& name)
	{
		return std::hash<std::string>{}(name);
	}
}