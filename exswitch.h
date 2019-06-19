/*
 *  exswitch.h - Exclusive Switch Model
 *  Attraction
 *
 *  Created by David Spieler on 1/27/12.
 *  Copyright 2012 Saarland University. All rights reserved.
 *
 */

#ifndef __EXSWITCH_H__
#define __EXSWITCH_H__

#include "common.h"


namespace attraction
{
	class ExclusiveSwitch
	{
	public:
		static const unsigned int modes     = 3;	/*number of model modes*/
		static const unsigned int dimension = 2;	/*dimension of the state space*/
		static const unsigned int maxFanout = 8;    /*maximum number of reactions*/

		/* State: | DNA | P1 | P2 | */

		static const uint U   = 0;
		static const uint B1  = 1;
		static const uint B2  = 2;

		static const uint P1  = 1;
		static const uint P2  = 0;


		static inline std::string getName()
		{
			std::string name("Exclusive Switch");
			return name;
		}
		
		static inline void exploreState(const State<dimension>& state, Row<dimension,maxFanout>& row, double* unknownParams, double* paramValues, int M)
		{	
			
			/*the default values of the params*/
			real params[maxFanout];
			params[0] = 0.05;	//lambda1	
			params[1] = 0.05;	//lambda2
			params[2] = 0.001;	//beta1
			params[3] = 0.001;	//beta2
			params[4] = 0.008;	//nu1
			params[5] = 0.008;  //nu2
			params[6] = 0.0005;	//delta1	
			params[7] = 0.0005;	//delta2
			
			/*get the values for the unknown params externally*/
			for (int i=0; i<M; i++)
				params[int(unknownParams[i])] = paramValues[i];
			
			row.fanout = 0;
			
			switch(state.mode) {
				case U: {
					/* G -> G + P1 (lambda1) */
					row.nonZeros[row.fanout].number = 0;
					row.nonZeros[row.fanout].rate = params[0];
					row.nonZeros[row.fanout].constant = params[0];
					row.nonZeros[row.fanout].successor = state;
					row.nonZeros[row.fanout].successor.var.c[P1]++;
					row.fanout++;
					
					/* G -> G + P2  (lambda2) */
					row.nonZeros[row.fanout].number = 1;
					row.nonZeros[row.fanout].rate = params[1];
					row.nonZeros[row.fanout].constant = params[1];
					row.nonZeros[row.fanout].successor = state;
					row.nonZeros[row.fanout].successor.var.c[P2]++;
					row.fanout++;
					
					/* G + P1 -> G.P1 (beta1)*/
					if (state.var.c[P1] > 0) {
						row.nonZeros[row.fanout].number = 2;
						row.nonZeros[row.fanout].rate = state.var.c[P1] * params[2];
						row.nonZeros[row.fanout].constant = params[2];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[P1]--;
						row.nonZeros[row.fanout].successor.mode = B1;
						row.fanout++;
					}
					
					/* G + P2 -> G.P2  (beta2) */
					if (state.var.c[P2] > 0) {
						row.nonZeros[row.fanout].number = 3;
						row.nonZeros[row.fanout].rate = state.var.c[P2] * params[3];
						row.nonZeros[row.fanout].constant = params[3];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[P2]--;
						row.nonZeros[row.fanout].successor.mode = B2;
						row.fanout++;
					}
				} break;
					
				case B1: {
					/* G.P1 -> G.P1 + P1  (lambda1)*/
					row.nonZeros[row.fanout].number = 0;
					row.nonZeros[row.fanout].rate = params[0];
					row.nonZeros[row.fanout].constant = params[0];
					row.nonZeros[row.fanout].successor = state;
					row.nonZeros[row.fanout].successor.var.c[P1]++;
					row.fanout++;
					
					/* G.P1 -> G + P1 (nu1)*/
					row.nonZeros[row.fanout].number = 4;
					row.nonZeros[row.fanout].rate = params[4];
					row.nonZeros[row.fanout].constant = params[4];
					row.nonZeros[row.fanout].successor = state;
					row.nonZeros[row.fanout].successor.var.c[P1]++;
					row.nonZeros[row.fanout].successor.mode = U;
					row.fanout++;
				} break;
					
				case B2: {
					/* G.P2 -> G.P2 + P2 (lambda2)*/
					row.nonZeros[row.fanout].number = 1;
					row.nonZeros[row.fanout].rate = params[1];
					row.nonZeros[row.fanout].constant = params[1];
					row.nonZeros[row.fanout].successor = state;
					row.nonZeros[row.fanout].successor.var.c[P2]++;
					row.fanout++;
					
					/* G.P2 -> G + P2  (nu2)*/
					row.nonZeros[row.fanout].number = 5;
					row.nonZeros[row.fanout].rate = params[5];
					row.nonZeros[row.fanout].constant = params[5];
					row.nonZeros[row.fanout].successor = state;
					row.nonZeros[row.fanout].successor.var.c[P2]++;
					row.nonZeros[row.fanout].successor.mode = U;
					row.fanout++;
				} break;
					
				default: break;
			}
			
			/* P1 -> 0  (delta1)*/
			if (state.var.c[P1] > 0) {
				row.nonZeros[row.fanout].number = 6;
				row.nonZeros[row.fanout].rate = state.var.c[P1] * params[6];
				row.nonZeros[row.fanout].constant = params[6];
				row.nonZeros[row.fanout].successor = state;
				row.nonZeros[row.fanout].successor.var.c[P1]--;
				row.fanout++;
			}
			
			/* P2 -> 0  (delta2)*/
			if (state.var.c[P2] > 0) {
				row.nonZeros[row.fanout].number = 7;
				row.nonZeros[row.fanout].rate = state.var.c[P2] * params[6];
				row.nonZeros[row.fanout].constant = params[7];
				row.nonZeros[row.fanout].successor = state;
				row.nonZeros[row.fanout].successor.var.c[P2]--;
				row.fanout++;
			}
		}
		
		/*returns the initial state of the specific model*/
		static State<dimension> getInitialState()
		{
			State<dimension> initial;
			
			initial.mode = U;
			initial.var.c[P1] = 0;
			initial.var.c[P2] = 0;
			
			return initial;
		}
	};
}
#endif
