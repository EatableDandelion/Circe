#include "Circe.h"
#include <memory>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "BVH.h"
#include "BroadColliders.h"

int main()
{

	using namespace Circe;

/*	CIRCE_INIT_PROFILER();

	for(int i = 0; i<10; i++){
	CIRCE_START_PROFILE(TEST);
	for(int i = 0; i<1000; i++)
	{
		std::cout << "1";
	}
	CIRCE_STOP_PROFILE(TEST);}
*/
/*	
	Mat A(3,4);
	A(0,0) = 0.0; A(0,1) = 0.0; A(0,2) = 10.0; A(0,3) = -1.0;
	A(1,0) = 10.0; A(1,1) = -10.0; A(1,2) = 10.0; A(1,3) = -5.0;
	A(2,0) = 10.0; A(2,1) = 10.0; A(2,2) = 0.0; A(2,3) = 4.0;

	Mat b(3,1);
	b(0) = 1.0;
	b(1) = 1.0;
	b(2) = 1.0;

	std::cout << b/A << std::endl;*/

	BVH<2, AABB> tree;


	tree.insert(1, Vec2(0,0),  Vec2(1,1), Vec2(0,0), 1000.0);
	tree.insert(2, Vec2(10,0), Vec2(1,1), Vec2(10,0),1.0);
	tree.insert(3, Vec2(20,0), Vec2(1,1), Vec2(5,0), 10.0);
	tree.insert(4, Vec2(21,0), Vec2(1,1), Vec2(6,0), 0.1);

	std::set<GravityWell> wells = tree.getGravityWells(1);

	for(GravityWell well : wells)
		std::cout << well.mass << " " << well.mass/well.d2 << std::endl;


	//std::cout << reg("test")[0] << std::endl;

/*	CIRCE_INITPROFILER;
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
	}
*/	
	
	/*Plot2d plot("1st plot");
	plot.setXAxisLabel("Time (s)");
	plot.addSeries("test.txt", "y(t)", 0, 1, LineType::LINE);
	plot.addSeries("test.txt", "exp(t)", 0, 2, LineType::LINE);
	plot.display();*/
	
	/*Polynomial p({0.0f, 0.0f, 1.5f});
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
	}*/
	
	/*KillSwitch killswitch;
	int i = 0;
	while(i<30 && !killswitch.isActivated())
	{
		std::cin.get();
		std::cout << i << std::endl;
		i++;
	}*/
	//FileReader reader("FoilMesh.dat");
	
	return 0;
}
