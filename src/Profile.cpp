#include "Profile.h"

namespace Circe
{
	
	unsigned int Debug::markIndex=0;
	high_resolution_clock::time_point Debug::timeStamp=high_resolution_clock::now();

	
	/** Log **/
	
	void Log::print(const std::string& type, const char* file, const int line, const std::string text) const
	{
		std::string fileName(file);
		std::cout << ("["+fileName+": "+std::to_string(line)+"] "+type+": "+text+"\n") << std::endl;
	}
	
	
	/** ProfileMarker **/
	
	ProfileMarker::ProfileMarker(const char* file, const char* functionName, const int line):name(functionName)
	{
		ProfilerService::get()->push(file, functionName, line);
	}
			
	ProfileMarker::~ProfileMarker()
	{
		ProfilerService::get()->pull();
	}
	
	
	/** Profile **/
	
	Profile::Profile(const std::string& fileName, const std::string& functionName, const int line, std::shared_ptr<Profile> parent):m_fileName(fileName), m_functionName(functionName), m_line(line), m_parent(parent), m_count(0), startTime(high_resolution_clock::now())
	{}
	
	void Profile::printHead(const double& scale)
	{
		for(int i =0; i< 20; i++)
			std::cout << std::endl;
		const std::size_t n(110);
		const std::size_t m(10);
		char data[n][m];
		for(int i =0; i<n; i++){
			for(int j =0; j<m; j++){
				data[i][j]=' ';
			}		
		}
		
		int x(0);
		int y(0);
		int ymax(0);
		m_time = 0.0;
		for(const auto& entry : m_children)
		{		
			m_time += entry.second->m_time;
		}
		 
		std::vector<int> limits;
		for(const auto& entry : m_children)
		{
			limits.push_back(x++);
			int y0 = y;
			std::shared_ptr<Profile> profile = entry.second;
			profile->print<n, m>(x, y0, scale, m_time, data);
			ymax = std::max(ymax, y0);
		}
		limits.push_back(x++);
		 
		for(int i : limits)
		{
			for(int j = 0; j < ymax+1; j++)
			{
				data[i][j] = '|';
			}
		}
		
		Circe::log.printArray<n,m>(data);
	}
	
	void Profile::printTreeHead()
	{
		m_time = 0.0;
		for(const auto& entry : m_children)
		{		
			m_time += entry.second->m_time;
		}
				
		printTree(0, m_time);
				
	}
	
	void Profile::printTree(const int depth, const double& totalTime)
	{		
		if(m_children.empty())
		{
			std::cout << "-----> ";
		}else{
			std::cout << "--+--> ";
		}
		double percentage(m_time*100.0/totalTime);
		std::cout << m_time << " ms - " << m_count << " calls - " << (m_time/m_count) << " ms/call (" << percentage << "%)";
		std::cout << "  ||  "<< m_fileName << " - " << m_functionName << " - Line " << m_line << std::endl;
		
		for(const auto& entry : m_children)
		{
			for(int i = 0; i<depth+1; i++)
			{
				std::cout << "  |";
			}std::cout << std::endl;
			
			for(int i = 0; i<depth; i++)
			{
				std::cout << "  |";
			}std::cout << "  +";

			entry.second->printTree(depth+1, totalTime);			
		}
		
	}
	
	
	/** Profiler **/
	
	Profiler::Profiler():m_head(std::make_shared<Profile>("","root", 0))
	{}
	
	Profiler::~Profiler()
	{
		m_head->printHead(100.0);
		m_head->printTreeHead();
	}

	void Profiler::push(const char* file, const char* functionName, const int line)
	{	
		
		std::string name = std::string(file)+" "+std::string(functionName)+" "+std::to_string(line);
		std::size_t id = Circe::getId(name);
		if(m_head->m_children.find(id) == m_head->m_children.end())
		{
			m_head->m_children[id] = std::make_shared<Profile>(std::string(file), std::string(functionName), line, m_head);
		}
		else
		{
			m_head->m_children[id]->startTime = high_resolution_clock::now();
		}
		m_head = m_head->m_children[id];
	}
	
	void Profiler::pull()
	{
		auto dt = duration_cast<microseconds>(high_resolution_clock::now() - m_head->startTime);
		m_head->m_count++;
		m_head->m_time += dt.count()/1000.0;
		m_head = m_head->m_parent;
	}
	
	std::shared_ptr<Profiler> ProfilerService::m_profiler;
	
}
