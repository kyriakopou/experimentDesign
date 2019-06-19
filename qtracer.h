/*
 *  qtracer.h
 *  Attraction
 *
 *  Created by David Spieler on 4/2/12.
 *  Copyright 2012 Saarland University. All rights reserved.
 *
 */

#ifndef __QTRACER_H__
#define __QTRACER_H__
#define  MAX_REACTIONS 8

#include <iostream>
#include <fstream>
#include "tracer.h"
#include <iomanip>


extern int numOfUnParams;

namespace attraction
{	

	struct QStateInfo
	{
		real p;
		real p_dt;
		real pder[MAX_REACTIONS];
		real pder_dt[MAX_REACTIONS];

		real k[4];
		real kder[MAX_REACTIONS][4];

		QStateInfo()
		{
			p = 0;
			p_dt = 0;


			k[0]   = 0;
			k[1]   = 0;
			k[2]   = 0;
			k[3]   = 0;

			kder[MAX_REACTIONS][0] = 0;
			kder[MAX_REACTIONS][1] = 0;
			kder[MAX_REACTIONS][2] = 0;
			kder[MAX_REACTIONS][3] = 0;

		}
	};
	
	template<class Model>
	class QTracer : public Tracer<Model,QStateInfo>
	{		
	public:



		void trace() {};

		void trace(real max_h, int sample_num, int M, double t_prev_max, double t_cur_max, double* unknownParams, double* paramValues)
		{

			/******************************************************************/
			/* Exploration of new states and computation of maximal exit rate */
			/******************************************************************/
			typename Tracer<Model,QStateInfo>::StateSet layer[4];
			const typename Tracer<Model,QStateInfo>::StateSet* lastLayer[4] = {&(this->m_sets[this->m_active]), &(layer[0]), &(layer[1]), &(layer[2])};
			
			/*For each step of Runge Kutta*/
			for (uint i = 0; i < 4; i++)
			{
				typename Tracer<Model,QStateInfo>::StateSet& lastSet = (i > 0) ? layer[i-1] : this->m_sets[this->m_active];	/*if i>0 set lastSet to layer[i-1], otherwise set it to current set*/
				
				/* iterate through lastSet of states*/
				typename Tracer<Model,QStateInfo>::StateIterator it;
				for (it = lastSet.begin(); it != lastSet.end(); it++)
				{
					/* retrieve corresponding entry from state*/
					State<Model::dimension> state = *it;
					typename Tracer<Model,QStateInfo>::DynSpaceType::Entry* entry = this->m_store.lookup(state);	
					
					/*stores entry with the info about the neighbours, rates, exitrate, etc */
					this->ensureExplored(entry, unknownParams, paramValues, M);
					
					/* add successor states to current layer */
					for (uint t = 0; t < entry->fanout; t++)
					{
						bool newState = true;
						for (uint j = 0; j <= i; j++)
						{	
							/*if the neighbor is already in the lastSet do nothing*/
							if (lastLayer[j]->find(entry->s[t]->key) != lastLayer[j]->end())
							{
								newState = false;
								break;
							}
						}
						/*otherwise put the state of the new discovered neighbour to the current layer*/
						if (newState)		
							layer[i].insert(entry->s[t]->key);
					}
				}
			}
			
			/*decide the timestep h regarding the maxRate
			 without exceeding the timeHorizon  */
			real h = 1.0 / this->m_maxRate / 10;
			h = (h < max_h) ? h : max_h;



			/***********************/
			/* Runge Kutta 4 steps */
			/***********************/
			
			/*define the Set of the active States*/
			typename Tracer<Model,QStateInfo>::StateSet& activeSet = this->m_sets[this->m_active];

			/*compute the intermediate derivative estimations of RK method k0, ..., k3*/
			for (uint i = 0; i < 4; i++)
			{
				/*for every state s belonging in the active set*/
				typename Tracer<Model,QStateInfo>::StateIterator it;
				for (it = activeSet.begin(); it != activeSet.end(); it++)
				{

					/*retrieve the corresponding entry of the current state */
					State<Model::dimension> state = *it;
					typename Tracer<Model,QStateInfo>::DynSpaceType::Entry* entry = this->m_store.lookup(state);
	
					/*compute the probability and the derivatives
					outflow for every possible successor of the state s*/
					for (uint t = 0; t < entry->fanout; t++)
					{	
						/*pp is the argument of F function, i.e. the approximation of pi_s(t+h/2), pi_s(t+h) for k's computation*/
						real pp = entry->value.p;
						real ppder[MAX_REACTIONS];
						for (int j = 0; j < numOfUnParams; j++)
							ppder[j] = entry->value.pder[j];

						
						if ((i == 1) || (i == 2)){
							pp += 0.5 * h * entry->value.k[i-1];
							for (int j = 0; j < numOfUnParams; j++)
								ppder[j] += 0.5 * h * entry->value.kder[j][i-1];
						}
						else if (i == 3){
							pp += h * entry->value.k[i-1];
							for (int j = 0; j < numOfUnParams; j++)
								ppder[j] += h* entry->value.kder[j][i-1];
						}

						/*for all k[i]'s*/
						/*compute the outcoming probability flow to t-th successor*/
						real flow = pp * entry->tp[t];
						/*subtract the outcoming probability flow from the current state's k[i] value*/
						entry->value.k[i] -= flow;
						/*add the outcoming flow to k[i] value of successor t of current state*/
						entry->s[t]->value.k[i] += flow;
						
						/*for all k[i]'s*/
						/*compute the outcoming j-th derivative flow to t-th successor*/
						real flow_der[MAX_REACTIONS];
						for (int j = 0; j < numOfUnParams; j++){
							flow_der[j] = ppder[j] * entry->tp[t];
							if (unknownParams[j] == entry->tn[t])
								flow_der[j] += pp * entry->tp[t] / entry->tc[t];
							/*subtract the outcoming derivative flow from the current state's k[i] value*/
							entry->value.kder[j][i] -= flow_der[j];
							/*add the outcoming derivative flow to k[i] value of successor t of current state*/
							entry->s[t]->value.kder[j][i] += flow_der[j];
						}
					}
				}
				
				/*update the active set of states*/
				/*we add all the new neighbours discovered above in the activeSet*/
				activeSet.insert(layer[i].begin(), layer[i].end());
			}
			
			/********************************************************************/
			/* Computation of next step probabilities and compaction of state space */
			/********************************************************************/

			/*for every state s belonging in the active set*/
			/*first keep track of its previous probability value and then update this value*/
			/*compute pi_s(t+h) = pi_s(t) + 1/6*(k0 +k1 + k2 +k3) -k's here include time h */
			/*and delete from the active set all the states that have probability smaller than the threshold*/

			//update the values of p, p_dt, pder, pder_dt
			typename Tracer<Model,QStateInfo>::StateIterator it;
			for (it = activeSet.begin(); it != activeSet.end();)
			{
				State<Model::dimension> state = *it;
				typename Tracer<Model,QStateInfo>::DynSpaceType::Entry* entry = this->m_store.lookup(state);
				

				/*computation of the next probability value pi(t+h). The weighted sum of the k[i]'s is finally the prob_mass_change for each state (i.e. inflow - outflow) */
				entry->value.p_dt = (entry->value.k[0] + 2 * entry->value.k[1] + 2 * entry->value.k[2] + entry->value.k[3]) / 6;
				entry->value.p += entry->value.p_dt * h ;

				/*computation of the next value of d/dc(pi(t+h))*/
				for (int j = 0; j < numOfUnParams; j++){
					entry->value.pder_dt[j] = (entry->value.kder[j][0] + 2 * entry->value.kder[j][1] + 2 * entry->value.kder[j][2] + entry->value.kder[j][3]) / 6;
					entry->value.pder[j] += entry->value.pder_dt[j] * h;
				}

				/*reinitialization of k values for R.K.*/
				entry->value.k[0] = 0;
				entry->value.k[1] = 0;
				entry->value.k[2] = 0;
				entry->value.k[3] = 0;
				
				/*reinitialization of k values for derivatives for R.K.*/
				for (int j = 0; j < numOfUnParams; j++){
					entry->value.kder[j][0] = 0;
					entry->value.kder[j][1] = 0;
					entry->value.kder[j][2] = 0;
					entry->value.kder[j][3] = 0;
				}

				/*compaction of the state space*/
				/*take out the states that have absolute information (and der_information) smaller than m_delta*/
				//if (!this->checkStateInfo(state, M))
				if (entry->value.p < this->m_delta)
					activeSet.erase(it++);
				else
					it++;



				
			}


			//update time
			this->m_time += h;


			/******************************************/
			/*Append the History of the new time point*/
			/******************************************/
			//if we don't have the history append the new entries to the end
			//of the history file
			if(t_prev_max < this->m_time && this->m_time <= t_cur_max)
				this->keepHistory(sample_num, M);



		}




	};
}
	
#endif
