#include "FileUtils.h"

namespace Circe
{
	DataFileWriter::DataFileWriter(const std::string& name):fileName(name)
	{
		
	}
	
	DataFileWriter::~DataFileWriter()
	{
		if(file.is_open())
		{
			file.close();
		}
	}
}
