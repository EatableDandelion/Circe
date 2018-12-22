#include <iostream>
#include "MathObj.h"
#include "Event.h"
#include <memory>

using namespace std;
using namespace Circe;


int main()
{

	Transform<3> transform;

	
	for(int i = 0; i<100; i++)
	{
		cout << transform.getRotation() << endl;
		transform.rotate(1.0f, 0.0f, 0.0f);
		
	}
	
	
	return 0;
}
