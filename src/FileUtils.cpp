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
	
	void FileReader::close()
	{
		if(file.is_open())
		{
			file.close();
		}
	}
	
	KillSwitch::KillSwitch():reader("killswitch.txt")
	{
		std::ofstream file("killswitch.txt");
		file.close();
		
	}
	
	KillSwitch::~KillSwitch()
	{
		reader.close();
		std::remove("killswitch.txt");
	}
			
	bool KillSwitch::isActivated()
	{
		if(!reader.isEmpty())
		{
			reader.close();
			std::remove("killswitch.txt");
			return true;
		}
		reader.close();
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
