

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
	
	
	/*Plot2d plot("1st plot");
	plot.setXAxisLabel("Time (s)");
	plot.addSeries("test.txt", "y(t)", 0, 1, LineType::LINE);
	plot.addSeries("test.txt", "exp(t)", 0, 2, LineType::LINE);
	plot.display();*/
	
	Polynomial p({0.0f, 0.0f, 1.5f});
	Polynomial p2({0.2f, 0.5f, 0.5f});
	Polynomial p3 = Circe::LegendrePolynomial(5);
	
	Vec<2> u({1.0f,1.0f});
	Mat<2> m;
	m(0,0) = 1.0f;
	m(0,1) = -0.5f;
	m(1,1) = 1.0f;
	m(1,0) = -0.2f;
	float t = 0.0f;
	for(int i = 0; i<100; i++)
	{
		RK4(u, m, 0.1f);
		t += 0.1f;
		std::cout << u << " " << exp(t) << std::endl;		
	}
	
	return 0;
}
