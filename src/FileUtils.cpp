#include "FileUtils.h"

namespace Circe
{
	FileReader::FileReader(const std::string& name):fileName(name), file(fileName)
	{}
	
	FileReader::~FileReader()
	{
		if(file.is_open())
		{
			file.close();
		}
	}
	
	void FileReader::read()
	{
		std::string line;
		if(!file.is_open())
		{
			file.open(fileName);	
			while (std::getline(file, line)) {
				lines.push_back(line);
			}
		}
	}
	
	bool FileReader::isEmpty()
	{
		if(!file.is_open())
		{
			file.open(fileName);	
		}
		return file.peek() == std::ifstream::traits_type::eof();
	}
	
	std::vector<std::string> FileReader::getLines() const
	{
		return lines;
	}
	
	void FileReader::clearFile()
	{
		std::fstream ofs;
		ofs.open(fileName, std::ios::out | std::ios::trunc);
		ofs.close();
	}
	
	
	KillSwitch::KillSwitch(const std::string& name):reader(name)
	{}
			
	bool KillSwitch::isActivated()
	{
		if(!reader.isEmpty())
		{
			reader.clearFile();
			return true;
		}
		return false;
	}
	
	
	
	DataFileWriter::DataFileWriter(const std::string& name, const char& delimiter):fileName(name), delimiter(delimiter)
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
