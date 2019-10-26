#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

namespace Circe
{
	class DataFileWriter
	{
		public:
			DataFileWriter(const std::string& name);
			~DataFileWriter();
			
			template<typename... Args>
			void writeLine(const Args... variables)
			{
				if(!file.is_open())
				{
					file.open(fileName);
				}

				file << addLine(variables...) << std::endl;
			}
		
		private:
			const std::string fileName;
			std::ofstream file;
			
			template<typename T, typename... Args>
			std::string addLine(const T variable, const Args... variables)
			{
				return std::to_string(variable)+"\t"+addLine(variables...);
			}
			
			template<typename T>
			std::string addLine(const T variable)
			{
				return std::to_string(variable);
			}
	};
	
	
}
