

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
	Quaternion q;
	Mat44 m;
	for(int i = 0; i<100; i++)
	{
		q.addAngle(0.01f, 0.1f, 0.1f, false);
		
		
	}
	float prev = 0;
	float theta = 0;
	for(int i = 0; i<100; i++)
	{
		q.addAngle(0.0f, 0.0f, 0.1f, true);
		m = Mat44::rotationMatrix(q);
		theta = atan2(m(1,0), m(0,0));
		std::cout << theta << std::endl;
		prev = theta;
	}
	//
	
	return 0;
}
