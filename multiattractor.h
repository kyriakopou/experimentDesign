/*
 *  exswitch.h - MUltiAttractor Model
 *  Attraction
 *
 *  Created by Charis Kyriakopoulos on 8/30/14.
 *  Copyright 2014 Saarland University. All rights reserved.
 *
 */

#ifndef __MULTIATTRACTOR_H__
#define __MULTIATTRACTOR_H__

#include "common.h"


namespace attraction
{
	class MultiAttractor
	{
	public:
		static const unsigned int modes     = 3;	/*number of model modes*/
		static const unsigned int dimension = 2;	/*dimension of the state space*/
		static const unsigned int maxFanout = 8;    /*maximum number of reactions*/

		/* State: | DNA | P1 | P2 | */

		static const uint U   = 0;
		static const uint B1  = 1;
		static const uint B2  = 2;

		static const uint MAFAProt  = 0;
		static const uint DeltaProt  = 1;
		static const uint PaxProt = 2;
		static const uint PaxDna = 3;
		static const uint MAFADna = 4;
		static const uint DeltaDna = 5;
		static const uint PaxDnaDeltaProt = 6;
		static const uint MAFADnaPaxProt = 7;
		static const uint MAFADnaMAFAProt = 8;
		static const uint MAFADnaDeltaProt = 9;
		static const uint DeltaDnaPaxProt = 10;
		static const uint DeltaDnaMAFAProt =11;
		static const uint DeltaDnaDeltaProt = 12;


		static inline std::string getName()
		{
			std::string name("Multi Attractor");
			return name;
		}
		
		static inline void exploreState(const State<dimension>& state, Row<dimension,maxFanout>& row, double* unknownParams, double* paramValues, int M)
		{	
			real params[maxFanout];
			params[0] = 5.0;	//p
			params[1] = 0.10;	//d
			params[2] = 1.0;	//b
			params[3] = 1.0;	//u

			
			row.fanout = 0;
			
			switch(state.mode) {
				case S: {
					/* PaxDna -> PaxDna + PaxProt  @ mass_action(p) */
					if (state.var.c[PaxDna] > 0){
						row.nonZeros[row.fanout].number = 0;
						row.nonZeros[row.fanout].rate = params[0] * state.var.c[PaxDna];
						row.nonZeros[row.fanout].constant = params[0];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[PaxProt]++;
						row.fanout++;
					}
					
					
					/* PaxProt -> 0  @ mass_action(d)*/
					if (state.var.c[PaxProt] > 0) {
						row.nonZeros[row.fanout].number = 1;
						row.nonZeros[row.fanout].rate = state.var.c[PaxProt] * params[1];
						row.nonZeros[row.fanout].constant = params[1];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[PaxProt]--;
						row.fanout++;
					}
					
					/*  PaxDna + DeltaProt -> PaxDnaDeltaProt  @ mass_action(b)*/
					if (state.var.c[PaxDna] > 0 && state.var.c[DeltaProt] > 0) {
						row.nonZeros[row.fanout].number = 2;
						row.nonZeros[row.fanout].rate = state.var.c[PaxDna] * state.var.c[DeltaProt] *  params[2];
						row.nonZeros[row.fanout].constant = params[2];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[PaxDna]--;
						row.nonZeros[row.fanout].successor.var.c[DeltaProt]--;
						row.nonZeros[row.fanout].successor.var.c[PaxDnaDeltaProt]++;
						row.fanout++;
					}
					
					/*  PaxDna + DeltaProt <- PaxDnaDeltaProt  @  mass_action(u)*/
					if (state.var.c[PaxDnaDeltaProt] > 0) {
						row.nonZeros[row.fanout].number = 3;
						row.nonZeros[row.fanout].rate = state.var.c[PaxDna] * state.var.c[DeltaProt] *  params[2];
						row.nonZeros[row.fanout].constant = params[2];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[PaxDna]++;
						row.nonZeros[row.fanout].successor.var.c[DeltaProt]++;
						row.nonZeros[row.fanout].successor.var.c[PaxDnaDeltaProt]--;
						row.fanout++;
					}

					/*  PaxDna + DeltaProt <- PaxDnaDeltaProt  @  mass_action(u)*/
					if (state.var.c[PaxDnaDeltaProt] > 0) {
						row.nonZeros[row.fanout].number = 4;
						row.nonZeros[row.fanout].rate = state.var.c[PaxDna] * state.var.c[DeltaProt] *  params[2];
						row.nonZeros[row.fanout].constant = params[2];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[PaxDna]++;
						row.nonZeros[row.fanout].successor.var.c[DeltaProt]++;
						row.nonZeros[row.fanout].successor.var.c[PaxDnaDeltaProt]--;
						row.fanout++;
					}


					/*  MAFADna -> MAFADna + MAFAProt  @ mass_action(p)*/
					if (state.var.c[MAFADna] > 0) {
						row.nonZeros[row.fanout].number = 5;
						row.nonZeros[row.fanout].rate = state.var.c[MAFADna] *  params[0];
						row.nonZeros[row.fanout].constant = params[0];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[MAFAProt]++;
						row.fanout++;
					}

					/*   MAFAProt -> 0  @ mass_action(d)*/
					if (state.var.c[MAFAProt] > 0) {
						row.nonZeros[row.fanout].number = 6;
						row.nonZeros[row.fanout].rate = state.var.c[MAFAProt] *  params[1];
						row.nonZeros[row.fanout].constant = params[1];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[MAFAProt]--;
						row.fanout++;
					}


					/* MAFADna + PaxProt -> MAFADnaPaxProt  @ mass_action(b)*/
					if (state.var.c[MAFADna] > 0 && state.var.c[PaxProt] > 0) {
						row.nonZeros[row.fanout].number = 7;
						row.nonZeros[row.fanout].rate = state.var.c[MAFAProt] * state.var.c[PaxProt] * params[2];
						row.nonZeros[row.fanout].constant = params[2];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[MAFADna]--;
						row.nonZeros[row.fanout].successor.var.c[PaxProt]--;
						row.nonZeros[row.fanout].successor.var.c[MAFADnaPaxProt]++;
						row.fanout++;
					}

					/* MAFADna + PaxProt <- MAFADnaPaxProt  @ mass_action(u)*/
					if (state.var.c[MAFADnaPaxProt] > 0) {
						row.nonZeros[row.fanout].number = 8;
						row.nonZeros[row.fanout].rate = state.var.c[MAFADnaPaxProt] * params[3];
						row.nonZeros[row.fanout].constant = params[3];
						row.nonZeros[row.fanout].successor = state;
						row.nonZeros[row.fanout].successor.var.c[MAFADna]++;
						row.nonZeros[row.fanout].successor.var.c[PaxProt]++;
						row.nonZeros[row.fanout].successor.var.c[MAFADnaPaxProt]--;
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
			
			initial.mode = U;
			initial.var.c[P1] = 0;
			initial.var.c[P2] = 0;
			
			return initial;
		}
	};
}
#endif
