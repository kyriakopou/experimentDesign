/*
 *  tracer.h
 *  Attraction
 *
 *  Created by David Spieler on 5/23/12.
 *  Copyright 2012 Saarland University. All rights reserved.
 *
 */

#ifndef __TRACER_H__
#define __TRACER_H__
#define  MAX_REACTIONS 8

#include "common.h"
#include "dynspace.h"
#include <cmath>
#include <set>
#include <iomanip>
#include "misc.h"
#include <iostream>
#include <fstream>


extern int numOfUnParams;


namespace attraction 
{

	struct Information {
		real information;
		real derivative;
	};


	/*struct storedState{
		State<Model::dimension> state;
		real prob;
		real prob_der[MAX_REACTIONS];
	};*/

	template<class Model, typename StateInfo>
	class Tracer
	{
	protected:
		typedef DynSpace<State<Model::dimension>, StateInfo, StateEquals<Model::dimension>, StateHash<Model::dimension>, Model::maxFanout > DynSpaceType;
		typedef std::set<State<Model::dimension> > StateSet;
		typedef typename std::set<State<Model::dimension> >::iterator StateIterator;
		
		DynSpaceType    m_store;

		real			m_delta;
		real			m_time;
		real			m_maxRate;
		real			m_maxderRate;
		
		StateSet		m_sets[2];	/*Statesets*/
		uint	 		m_active;	/*number of active sets*/
		
		void ensureExplored(typename DynSpaceType::Entry* entry)
		{
			if (entry->fanout == INVALID)
			{
				
				//with exploreState  we store in struct row all the necessary information
				Row<Model::dimension,Model::maxFanout> row;
				Model::exploreState(entry->key, row);
				
				real exitRate = row.exitRate();
				/*find the maximum exit rate*/
				m_maxRate = (exitRate > m_maxRate) ? exitRate : m_maxRate;	
				
				/*copy all the info for the state from row struct to entry struct*/
				entry->exitRate = exitRate;
				entry->iexitRate = 1.0 / exitRate;
				entry->fanout = row.fanout;
				for (uint i = 0; i < entry->fanout; i++)
				{
					/*obtain transition rates*/
					entry->tp[i] = row.nonZeros[i].rate;
					/*obtain transition number*/
					entry->tn[i] = row.nonZeros[i].number;
					/*obtain transition constant*/
					entry->tc[i] = row.nonZeros[i].constant;
					/*obtain the successor entries*/
					entry->s[i] = m_store.ensureExistence(row.nonZeros[i].successor);
				}

				/*consider also the maxRate w.r.t. derivatives*/
				real derRate[MAX_REACTIONS];

				for (int j = 0; j < numOfUnParams; j++){
					derRate[j] = entry->tp[j] / entry->tc[j];
					m_maxderRate = (derRate[j] > m_maxderRate) ? derRate[j] : m_maxderRate;
				}

				m_maxRate = (m_maxderRate > m_maxRate) ? m_maxderRate : m_maxRate;


			}
		}
		
	public:
		Tracer()
		{
			//set up the hash tables that give the ptrs to the entry struct
			m_store.setup(1024*1024, 1024*1024);
			
			m_delta = 1e-20;
			m_time = 0;
			m_maxRate = 0;
			m_maxderRate = 0;
			m_active = 0;

		}
		
		void setDelta(real delta)
		{
			m_delta = delta;
		}
		
		/*returns the total probability mass of the truncated state space */ 
		real pMass()
		{
			real mass = 0;
			
			StateIterator it;
			for (it = m_sets[m_active].begin(); it != m_sets[m_active].end(); it++)
			{
				typename DynSpaceType::Entry* entry = m_store.lookup(*it);
				mass += entry->value.p;
			}
			
			return mass;
		}
		
		/*returns the size of the truncated state space*/
		uint activeStates()
		{
			return m_sets[m_active].size();
		}
		
		/*returns a pointer to the beginning of the truncated state space table*/
		StateIterator begin()
		{
			return m_sets[m_active].begin();
		}
		
		/*returns a pointer to the end of the truncated state space table*/
		StateIterator end()
		{
			return m_sets[m_active].end();
		}
		
		/*get the probability of a state*/
		real getProbability(const State<Model::dimension>& state)
		{
			return m_store.lookup(state)->value.p;
		}
		
		/*get the total time elapsed*/
		real getTime()
		{
			return m_time;
		}
		
		/*put initial state into the active set*/
		void putInitialState(const State<Model::dimension>& state, real p, double* unknownParams, double* paramValues, int M)
		{
			typename DynSpaceType::Entry* entry = m_store.lookup(state);
			
			/* maybe the state was added when exploring another one, add it to the active set */
			if (entry)
			{
				entry->value.p += p;

				m_sets[m_active].insert(state); 
			}
			else
			{
				entry = m_store.add(state);
				entry->value.p = p;

				ensureExplored(entry, unknownParams, paramValues, M);
				m_sets[m_active].insert(state);
			}



		}
		

		/*print the whole active state space*/
		/*and the probability of each state*/
		void printResults()
		{

			std::cout << "START" << std::setw(15) << "Prob";
			for (int i = 0; i < numOfUnParams; i++)
				std::cout << std::setw(15) << "der" << i;
			std::cout << std::endl;
			StateIterator it;
			for (it = m_sets[m_active].begin(); it != m_sets[m_active].end(); it++)
			{
				State<Model::dimension> state = *it;
				typename DynSpaceType::Entry* entry = m_store.lookup(state);
				
				//print the transient probability values
				std::cout << entry->key.mode << " | ";
				for (uint i = 0; i < Model::dimension; i++)
				{
					std::cout << entry->key.var.c[i] << " ";
				}
				std::cout << ":" <<  std::setw(15) << entry->value.p;

				//print the derivative values
				for (int i = 0; i < numOfUnParams; i++)
					std::cout << std::setw(15) << entry->value.pder[i];
				std::cout << std::endl;
			}

			std::cout << "END" << std::endl << std::endl;
		}
		

		struct Moments
		{
			real firstMoment[Model::dimension];
			real secondMoment[Model::dimension];
			real thirdMoment[Model::dimension];
			real fourthMoment[Model::dimension];
			real secondCentralMoment[Model::dimension];
			real thirdCentralMoment[Model::dimension];
			real fourthCentralMoment[Model::dimension];
			real firstMomentDer[Model::dimension];
			real secondMomentDer[Model::dimension];
			real secondCentralMomentDer[Model::dimension];
		};



		/*print some statistics*/
		Moments computeStatistics()
		{

			Moments stat;

			//initialize all the moments
			for (uint i = 0; i < Model::dimension; i++)
			{
				stat.firstMoment[i] = 0;
				stat.secondMoment[i] = 0;
				stat.thirdMoment[i] = 0;
				stat.fourthMoment[i] = 0;
				stat.secondCentralMoment[i] = 0;
				stat.thirdCentralMoment[i] = 0;
				stat.fourthCentralMoment[i] = 0;
				stat.firstMomentDer[i] = 0;
				stat.secondMomentDer[i] = 0;
				stat.secondCentralMomentDer[i] = 0;
			}
			
			StateIterator it;
			for (it = m_sets[m_active].begin(); it != m_sets[m_active].end(); it++)
			{
				State<Model::dimension> state = *it;
				typename DynSpaceType::Entry* entry = m_store.lookup(state);
				
				//compute the moments of the distribution for each species
				//and the derivatives of the first and the second moment wrt to unknown parameter
				for (uint i = 0; i < Model::dimension; i++)
				{
					stat.firstMoment[i] += entry->key.var.c[i] * entry->value.p;
					stat.secondMoment[i] += pow(entry->key.var.c[i],2) * entry->value.p;
					stat.thirdMoment[i] += pow(entry->key.var.c[i], 3) * entry->value.p;
					stat.fourthMoment[i] += pow(entry->key.var.c[i], 4) * entry->value.p;

					//derivatives of the first two moments
					stat.firstMomentDer[i] += entry->key.var.c[i] * entry->value.pder;
					stat.secondMomentDer[i] += pow(entry->key.var.c[i], 2) * entry->value.pder;
				}
			}
			
			//compute the central moments of the distribution for each species
			for (uint i = 0; i < Model::dimension; i++)
			{
				real mean = stat.firstMoment[i];
				stat.secondCentralMoment[i] = stat.secondMoment[i] - pow(mean,2);
				stat.thirdCentralMoment[i] = stat.thirdMoment[i] - 3 * mean * stat.secondMoment[i] + 2 * pow(mean,3);
				stat.fourthCentralMoment[i] = stat.fourthMoment[i] - 4 * mean * stat.thirdMoment[i] + 6 * pow(mean,2) * stat.secondMoment[i] - 3 * pow(mean,4);
				//derivative of the second central moment
				stat.secondCentralMomentDer[i] = stat.secondMomentDer[i] - 2 * stat.firstMomentDer[i] * mean;
			}


			/*std::cout << "STATISTICS" << std::endl << std::endl;
			for (uint i = 0; i < Model::dimension; i++)
			{
				real mean = fmoment[i];
				real variance = smoment[i] - (mean * mean);
				real stddev = sqrt(variance);
				std::cout << "var " << i << " : " << mean << " +- " << stddev << std::endl;
			}*/

			return stat;
		}



		void computeInformation(double *Fisher, double *derFisher, int M)
		{


			/*Initialize the Fisher Information and its derivative */
			for (int m = 0; m < M; m++){
				for (int n = 0; n< M; n++){
					Fisher[m*M + n] = 0;
					derFisher[m*M + n] = 0;
				}
			}

			StateIterator it;
			for (it = m_sets[m_active].begin(); it != m_sets[m_active].end(); it++)
			{
				State<Model::dimension> state = *it;
				typename DynSpaceType::Entry* entry = m_store.lookup(state);

				//Fill each entry of the Fisher Information Matrix and its derivative
				for (int m = 0; m < M; m++){
					for (int n = 0; n < M; n++){
						real state_info = (entry->value.pder[m] * entry->value.pder[n]) / entry->value.p;
						real state_der_info = ((entry->value.pder_dt[m] * entry->value.pder[n] + entry->value.pder[m] * entry->value.pder_dt[n]) * entry->value.p - entry->value.pder[m] * entry->value.pder[n] * entry->value.p_dt ) / pow(entry->value.p, 2);
						if ( !std::isnan(state_info) )
							Fisher[m*M + n] += state_info;
						if ( !std::isnan(state_der_info) )
							derFisher[m*M + n] += state_der_info;
					}
				}

			}

			return;
		}
		
		bool checkStateInfo(State<Model::dimension> state, int M)
		{
			typename DynSpaceType::Entry* entry = m_store.lookup(state);
			bool cond = 0;
			int m = 0;
			int n = 0;
			while ( m < M  && cond == 0){
				while (n < M && cond == 0){
					real state_info = (entry->value.pder[m] * entry->value.pder[n]) / entry->value.p;
					real state_der_info = ((entry->value.pder_dt[m] * entry->value.pder[n] + entry->value.pder[m] * entry->value.pder_dt[n]) * entry->value.p - entry->value.pder[m] * entry->value.pder[n] * entry->value.p_dt ) / pow(entry->value.p, 2);
					if ( (state_info > this->m_delta || state_der_info > this->m_delta || state_der_info < -this->m_delta) && !std::isnan(state_info) && !std::isnan(state_der_info) ) 
						cond = 1;
					n++;
				}	
				m++;
		
			}
			
			return cond;
			
		}

		Information computeInformationLygeros(uint species)
		{
			Information info;
			Moments stat = computeStatistics();

			info.information = 0;
			/*keep the information and the derivative of it
			 in the table inform*/
			real sample_mean_info = pow(stat.firstMomentDer[species], 2) / stat.secondCentralMoment[species];
			real nominator = pow(stat.secondCentralMoment[species] * stat.secondCentralMomentDer[species] - stat.firstMomentDer[species] * stat.thirdCentralMoment[species], 2);
			real denominator = pow(stat.secondCentralMoment[species], 2) * (stat.fourthCentralMoment[species] - pow(stat.secondCentralMoment[species], 2)) - stat.secondCentralMoment[species] * pow(stat.thirdCentralMoment[species], 2);
			if ( !std::isnan(nominator / denominator) && !std::isnan(sample_mean_info) )
				info.information =  sample_mean_info + nominator / denominator;

			return info;
		}




		void keepHistory(int sample_num, int M){


			//Write the new times and the new information into files
			//in order to take the plots in MATLAB
			//Comment if you run the optimization procedure
			std::ofstream times, fim, fim_dt;
			std::string name1 = "data_files/times", name2 = "data_files/fim", name3 = "data_files/fim_dt";

			char numstr[21];
			sprintf(numstr, "%d", sample_num);
			//define the file names
			name1 = name1 + numstr + ".dat";
			name2 = name2 + numstr + ".dat";
			name3 = name3 + numstr + ".dat";

			//open the files
			times.open(name1.c_str(), std::ios::app | std::ios::binary);
			fim.open(name2.c_str(), std::ios::app | std::ios::binary);
			fim_dt.open(name3.c_str(), std::ios::app | std::ios::binary);


			//define the arrays and the pointers to return the history
			double FisherMatrix[M][M], derFisherMatrix[M][M];
			double *p1, *p2;
			p1 = &FisherMatrix[0][0];
			p2 = &derFisherMatrix[0][0];

			//compute FIM and FIM_der
			this->computeInformation(p1, p2, M);

			//write FIM and FIM_der in the corresponding files
			if (fim.is_open() && fim_dt.is_open() ){

				fim.write((char*) &FisherMatrix[0], sizeof(double) * M * M );
				fim_dt.write((char*) &derFisherMatrix[0], sizeof(double) * M * M );

			}

			//write in times
			if (times.is_open())
				times.write((char*) &this->m_time, sizeof(double) );

			//close all the files
			times.close();
			fim.close();
			fim_dt.close();

		}



		void revokeHistory(double t, double *Fisher, double *derFisher, int sample_num, int M){

			//open files times and
			std::ifstream times, fim, fim_dt;
			std::string name1 = "data_files/times", name2 = "data_files/fim", name3 = "data_files/fim_dt";

			char numstr[21];
			sprintf(numstr, "%d", sample_num);
			//define the file names
			name1 = name1 + numstr + ".dat";
			name2 = name2 + numstr + ".dat";
			name3 = name3 + numstr + ".dat";
			//open the files
			times.open(name1.c_str(), std::ios::app | std::ios::binary);
			fim.open(name2.c_str(), std::ios::app |std::ios::binary);
			fim_dt.open(name3.c_str(), std::ios::app | std::ios::binary);

			//take the number of lines of times
			std::streampos sizeInBytes;
			times.seekg(0, std::ios::end);
			sizeInBytes = times.tellg();
			//std::cout << size << std::endl;
			times.seekg(0, std::ios::beg);
			int num_lines = sizeInBytes  / sizeof(double);
			//std::cout << num_lines << std::endl;

			//temporally map file to an array of size = lines
			double timePoints[num_lines];
			//store the timePoints from the file into a double array
			times.read((char*) &timePoints[0], sizeInBytes );

			//do a binary search into timePoints to find
			//the line_number of the time points closest to t_new
			int next = binarySearch(timePoints, 0, num_lines, t);
			int prev = next - 1;
			//std::cout << prev << std::endl;

			//find the corresponding rows in Fim, Fim_dt file
			//t_next, t_prev and place there the position of the file
			int size = M * M;
			fim.seekg(prev * sizeof(double) * size);
			fim_dt.seekg(prev * sizeof(double) * size);

			//copy the entries of the files into 4 temporal arrays
			double fim_next[size], fim_dt_next[size], fim_prev[size], fim_dt_prev[size];
			fim.read((char*) &fim_prev[0], size * sizeof(double));
			fim_dt.read((char*) &fim_dt_prev[0], size * sizeof(double));

			fim.read((char*) &fim_next[0], size * sizeof(double));
			fim_dt.read((char*) &fim_dt_next[0], size * sizeof(double));


			times.close();
			fim.close();
			fim_dt.close();

			//Take the information matrices for the intermediate time point t by doing interpolation
			//between the matrices of the previous and next time points
			for (int i = 0; i < size; i++){
					Fisher[i] = linearInterpolation(timePoints[prev], timePoints[next] , fim_prev[i], fim_next[i], t);
					derFisher[i] = linearInterpolation(timePoints[prev], timePoints[next] , fim_dt_prev[i], fim_dt_next[i], t);

			}



		}

		void clear(){

			m_sets[0].clear();
			m_sets[1].clear();
		}


		//THIS HAS TO BE DONE AFTER THE PAPER----------------------------
		//write two structs into two files: one containing the states and the other containing the prob, and prob_der[i]
		/*void keepProbHistory(int sample_num, int M){

			int size = M * M;

			std::ofstream storedInfo;
			std::string name = "storedInfo";
			char numstr[21];
			sprintf(numstr, "%d", sample_num);
			//define the file names
			name = name + numstr + ".dat";

			storedInfo.open(name.c_str(), std::ios::binary);


			StateIterator it;
			for (it = m_sets[m_active].begin(); it != m_sets[m_active].end(); it++)
			{
				State<Model::dimension> state = *it;
				typename DynSpaceType::Entry* entry = m_store.lookup(state);

				//store the transient probability values into the file
				storedInfo.write((char*) entry->key.mode, sizeof(uint));
				for (uint i = 0; i < Model::dimension; i++)
				{
					storedInfo.write((char*) entry->key.var.c[i], sizeof(int));
				}
				storedInfo.write((char*) entry->value.p, sizeof(double));

				//store the derivative values into the file
				for (int i = 0; i < numOfUnParams; i++)
					storedInfo.write((char*) entry->value.pder[i], sizeof(double));


			}

			storedInfo.close();


		}



		storedState revokeState(int sample_num, int M){

			storedState stored;

			std::ifstream storedInfo;
			std::string name = "storedInfo";
			char numstr[21];
			sprintf(numstr, "%d", sample_num);
			//define the file names
			name = name + numstr + ".dat";

			storedInfo.open(name.c_str(), std::ios::binary);

			State<Model::dimension> state;

			storedInfo.read((char*) state.mode, sizeof(uint));
			for (uint i = 0; i < Model::dimension; i++)
			{
				storedInfo.read((char*) &state.var.c[i], sizeof(int));
			}
			stored.state = state;

			//read the prob and the prob_der of each state
			storedInfo.read((char*) &stored.prob, sizeof(real));
			for (int i = 0; i < numOfUnParams; i++){

				storedInfo.read((char*) &stored.prob_der[i], sizeof(real));
			}



			return stored;


		}*/












	};




}

#endif
