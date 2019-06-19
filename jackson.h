/****Tandem Jackson Model*********/
/****R1: 0 -> A ********/
/****R2: A -> B *****/
/****R3: B -> 0 *****/

#ifndef __Jackson_H__
#define __Jackson_H__


#include "common.h"



namespace attraction
{
	class Jackson
	{
	public:
		static const unsigned int modes = 1;		/*number of modes*/
		static const unsigned int dimension = 2;	/*dimension of the state space*/
		static const unsigned int maxFanout = 3;	/*number of reactions*/

		/*State: |A | B | C | D |*/
		
		static const uint S = 0;

		static const uint A = 0;
		static const uint B = 1;

		static const real lambda = 0.04;
		static const real mu1 = 0.48;
		static const real mu2 = 0.48;
		
	
		static inline std::string getName()
		{
			std::string name("Jackson");
			return name;
		}

		static inline void exploreState(const State<dimension>& state, Row<dimension,maxFanout>& row, double* unknownParams, double* paramValues, int M)
		{
			
			real params[maxFanout];
			params[0] = 0.04;
			params[1] = 0.48;
			params[2] = 0.48;
			
			
			
			row.fanout = 0;
			switch(state.mode){
				case S: {

					/* 0 -> A   (lambda) */
					if (true) {
						row.nonZeros[row.fanout].number = 0;
						row.nonZeros[row.fanout].constant = params[0];
						row.nonZeros[row.fanout].rate = lambda;
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[A] = row.nonZeros[row.fanout].successor.var.c[A]+1;
						row.fanout++;
					}

					/* A -> B	(mu1) */
					if (state.var.c[A] > 0) {
						row.nonZeros[row.fanout].number = 1;
						row.nonZeros[row.fanout].constant = params[1];
						row.nonZeros[row.fanout].rate = mu1;
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[A]--;
						row.nonZeros[row.fanout].successor.var.c[B]++;	
						row.fanout++;
					}
					break;
					
					/* B -> 0	(mu2) */
					if (state.var.c[B] > 0) {
						row.nonZeros[row.fanout].number = 2;
						row.nonZeros[row.fanout].constant = params[2];
						row.nonZeros[row.fanout].rate = mu2;
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[B]--;
						row.fanout++;
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
			initial.var.c[A] = 0;
			initial.var.c[B] = 0;

			return initial;
		}
	};
}

#endif