

#include "Circe.h"
#include <memory>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>



int main()
{

	using namespace Circe;
	
/*
	CIRCE_INITPROFILER;
	{
		CIRCE_PROFILEBLOCK;
		{
			for(int k = 0; k<5; k++){
				CIRCE_PROFILEBLOCK;
				for(int i = 0; i<1000; i++)
				{
					std::cout << "1";
				}
			}
		}
		{
			CIRCE_PROFILEBLOCK;
			
			{
				CIRCE_PROFILEBLOCK;
				for(int i = 0; i<600; i++)
				{
					std::cout << "2";
				}
			}
			
			{
				CIRCE_PROFILEBLOCK;
				for(int i = 0; i<800; i++)
				{
					std::cout << "2";
				}
			}
		}
	}*/
	
	Mat44 m = Mat44::rotationMatrix(Vec<3>(0.8,0.5,0), Vec<3>(-0.6,0.5,0));
	std::cout << m << std::endl;
	std::cout << Mat44::rotationMatrix(m.toQuaternion()) << std::endl;
	
	
	return 0;
}
