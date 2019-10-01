

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
	
	/*Mat44 m = Mat44::rotationMatrix(Vec<3>(0.8,0.5,0), Vec<3>(-0.6,0.5,0));
	std::cout << m << std::endl;
	std::cout << Mat44::rotationMatrix(m.toQuaternion()) << std::endl;*/
	

	//Direction<3> direction({0.0f,1.0f, 0.0f}, REF_FRAME::LOCAL);
	Direction<3> direction(REF_FRAME::LOCAL, 0.0f,1.0f,0.0f);
	Transform<3> transform;
	transform.setRotation(Vec<3>(1,1,0), Vec<3>(-1,1,0));
	transform.translate(direction);
	std::cout << transform.getPosition() << std::endl;
	return 0;
}
