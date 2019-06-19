/*
 *  common.h
 *  Attraction
 *
 *  Created by David Spieler on 2/1/12.
 *  Copyright 2012 Saarland University. All rights reserved.
 *
 */
/* Here the state is defined as a struct and operations btw states are also defined*/ 


#ifndef __COMMON_H__
#define __COMMON_H__

#include <sys/time.h>
#include <cmath>
#include "vector.h"



namespace attraction
{
	typedef double		real;
	typedef unsigned int	uint;
	
	const uint INVALID = 0xFFFFFFFF;
	
	/*                                    */
	/* State Descriptor                   */
	/*                                    */
	/* | mode | var[0] | ... | var[d-1] | */
	/*                                    */
	template<uint d>
	struct State
	{
		uint		mode;
		Vector<uint,d>	var;
		
		inline bool operator < (const State& o) const
		{
			if (mode < o.mode)
				return true;
			else if (mode > o.mode)
				return false;
			
			for (uint i = 0; i < d; i++)
			{
				if (var.c[i] < o.var.c[i])
					return true;
				else if (var.c[i] > o.var.c[i])
					return false;
			}
			
			return false;
		}
		
		/*return the state as a stream of characters*/
		/*friend inline std::ostream& operator << (std::ostream& output, const State<d>& s)
		{
			output << s.mode << " " << s.var;
			return output;
		}*/
	};
	
	/* State Hash Functor */
	template<uint d>
	struct StateHash
	{
		inline size_t operator () (const State<d>& state) const
		{
			size_t hash = state.mode;
			
			for (uint i = 0; i < d; i++)
				hash = hash * 7919 + state.var.c[i];
			
			return hash;
		}
	};
	
	/* State Equality Functor */
	template<uint d>
	struct StateEquals
	{
		inline bool operator () (const State<d>& s1, const State<d>& s2) const
		{
			if (s1.mode != s2.mode)
				return false;
			
			for (uint i = 0; i < d; i++)
			{
				if (s1.var.c[i] != s2.var.c[i])
					return false;
			}
			
			return true;
		}
	};
	
	/*            */
	/* Transition */
	/*            */
	
	template<uint d>
	struct Transition
	{
		int 	number;
		real	constant;
		real    rate;
		State<d> successor;

		/*
		friend inline std::ostream& operator << (std::ostream &output, const Transition<d> &t)
		{
			output << "--" << t.rate << "--" << t.successor;
			return output;
		}*/
	};
	
	
	/*     */
	/* Row */
	/*     */
	template<uint d, uint f>
	struct Row
	{
		uint          fanout;
		Transition<d> nonZeros[f];
		
		/*returns the exit rate of a state*/
		inline real exitRate()
		{
			real er = 0;
			for (uint i = 0; i < fanout; i++) 
				er += nonZeros[i].rate;
			return er;
		}



		/*
		friend inline std::ostream& operator << (std::ostream& output, const Row<d,f>& row)
		{
			for (uint i = 0; i < row.fanout; i++)
				output << row.nonZeros[i] << "\t\t";
			return output;
		}*/
	};	
	
	
}

#endif
