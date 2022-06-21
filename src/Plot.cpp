#include "Plot.h"

namespace Circe
{
	Series2d::Series2d(const std::string& fileName, const std::string& seriesName, const int& xIndex, const int& yIndex, const bool& visibl, const LineType& type):fileName(fileName), seriesName(seriesName), xIndex(xIndex), yIndex(yIndex), visible(visible), type(type)
	{}
		
	std::string Series2d::getSeriesName() const
	{
		return seriesName;
	}
	
	std::string Series2d::getFileName() const
	{
		return fileName;
	}
	
	int Series2d::getXIndex() const
	{
		return xIndex;
	}
	
	int Series2d::getYIndex() const
	{
		return yIndex;
	}
	
	bool Series2d::isVisible() const
	{
		return visible;
	}
			
	void Series2d::setSeriesName(const std::string& newName)
	{
		seriesName = newName;
	}
	
	void Series2d::setVisible(const bool& isVisible)
	{
		visible = isVisible;
	}
	
	LineType Series2d::getType() const
	{
		return type;
	}
	
	
	Plot2d::Plot2d(const std::string& title):title(title)
	{
		pipe = popen("gnuplot -persist", "w");
		setTitle(title);
	}
	
	Plot2d::~Plot2d()
	{
		fprintf(pipe, "exit\n");
		pclose(pipe);
	}
	
	void Plot2d::setTitle(const std::string& title)
	{
		printCommand("set title \""+title+"\"");
	}
	
	void Plot2d::setXAxisLabel(const std::string& label)
	{
		printCommand("set xlabel \""+label+"\"");
	}
	
	void Plot2d::setYAxisLabel(const std::string& label)
	{
		printCommand("set ylabel \""+label+"\"");
	}
	
	Series2d Plot2d::addSeries(const std::string& fileName, const std::string& seriesName, const int& xIndex, const int& yIndex, const LineType& type)
	{
		Series2d newSeries(fileName, seriesName, xIndex+1, yIndex+1, true, type);
		allSeries.push_back(newSeries);
		return newSeries;
	}
	
	void Plot2d::display()
	{
		bool hasASeries = false;
		std::string command = "plot";
		for(Series2d series : allSeries)
		{
			if(series.isVisible())
			{
				if(hasASeries) // If true, this is the second series of the plot, then add a comma
				{
					command += ",";
				}
				hasASeries = true;
				command +=  " \""+series.getFileName()+"\" using ";
				command += std::to_string(series.getXIndex())+":"+std::to_string(series.getYIndex());
				
				if(series.getType() == LineType::LINE)
				{
					command += " with lines";
				}
				else if(series.getType() == LineType::STEPS)
				{
					command += " with steps";
				}
				command += " title \""+series.getSeriesName()+"\"";
			}
		}
		
		if(hasASeries) // Only display if there is a series
		{
			printCommand(command);
		}
	}
	
	void Plot2d::printCommand(const std::string& command)
	{
		fprintf(pipe, "%s\n", command.c_str());
	}
}