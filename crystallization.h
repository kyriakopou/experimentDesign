/****Crystallization Model*********/
/****R1: 2A -> B********/
/****R2: A + C -> D*****/


#ifndef __CRYSTALLIZATION_H__
#define __CRYSTALLIZATION_H__

#include "common.h"


namespace attraction
{
	class Crystallization
	{
	public:
		static const unsigned int modes = 1;		/*number of modes*/
		static const unsigned int dimension = 4;	/*dimension of the state space*/
		static const unsigned int maxFanout = 2;	/*number of reactions*/

		/*State: |A | B | C | D |*/
		
		static const uint S = 0;

		static const uint A = 0;
		static const uint B = 1;
		static const uint C = 2;
		static const uint D = 3;
		
	
		static inline std::string getName()
		{
			std::string name("Crystallization");
			return name;
		}

		static inline void exploreState(const State<dimension>& state, Row<dimension,maxFanout>& row, double* unknownParams, double* paramValues, int M)
		{
			/*give the default values of the parameters*/
			real params[maxFanout];
			params[0] = 4;	/*theta1*/
			params[1] = 0.1;	/*theta2*/
			
			/*get the values for the unknown params externally*/
			for (int i=0; i<M; i++)
				params[int(unknownParams[i])] = paramValues[i];
			
			
			row.fanout = 0;
			switch(state.mode){
				case S: {

					/* 2A -> B      (c1) */
					if (state.var.c[A] > 1) {
						row.nonZeros[row.fanout].number = 0;
						row.nonZeros[row.fanout].rate = state.var.c[A] * (state.var.c[A]-1)/2 * params[0];
						row.nonZeros[row.fanout].constant = params[0];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[A] = row.nonZeros[row.fanout].successor.var.c[A]-2;
						row.nonZeros[row.fanout].successor.var.c[B]++;
						row.fanout ++;
					}


					/* A + C -> D	(c2) */
					if (state.var.c[A] > 0 && state.var.c[C] > 0) {
						row.nonZeros[row.fanout].number = 1;
						row.nonZeros[row.fanout].rate = state.var.c[A] * state.var.c[C] * params[1];
						row.nonZeros[row.fanout].constant = params[1];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[A]--;
						row.nonZeros[row.fanout].successor.var.c[C]--;
						row.nonZeros[row.fanout].successor.var.c[D]++;
						row.fanout ++;
					}



				break;
				}
			}
		}

		/*returns the initial state of the specific model*/
		static State<dimension> getInitialState()
		{
			State<dimension> initial;
			
			initial.mode = S;
			initial.var.c[A] = 4;
			initial.var.c[B] = 3;
			initial.var.c[C] = 2;
			initial.var.c[D] = 1;


			return initial;
		}
	};
}

#endif
