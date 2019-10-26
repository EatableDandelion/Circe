#pragma once

#include <iostream>
#include <MathObj.h>

namespace Circe
{
	template<std::size_t N>
	void RK4(Vec<N>& stateVector, const Mat<N>& systemMatrix, const float& dt)
	{
		Vec<N> k1 = systemMatrix*stateVector*dt;
		Vec<N> k2 = systemMatrix*(k1/2.0f+stateVector)*dt;
		Vec<N> k3 = systemMatrix*(k2/2.0f+stateVector)*dt;
		Vec<N> k4 = systemMatrix*(k3+stateVector)*dt;
		stateVector += (k1+2.0f*k2+2.0f*k3+k4)/6.0f;
	}
}