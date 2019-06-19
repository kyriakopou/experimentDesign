/*
 *  vector.h
 *  Attraction
 *
 *  Created by David Spieler on 1/27/12.
 *  Copyright 2012 Saarland University. All rights reserved.
 *
 */
/*Vector is defined as a struct and operations between vectors are also defined*/

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <iostream>

namespace attraction {
	
	template<typename T, unsigned int n>
	struct Vector
	{		
		T c[n];
		
		/*CONSTRUCTORS*/
		inline Vector<T,n>()
		{
			for (unsigned int i=0; i < n; i++)
				c[i] = 0;
		}
		
		inline Vector<T,n>(T s)
		{
			for (unsigned int i=0; i < n; i++)
				c[i] = s;
		}
		
		/*OPERATORS*/
		/*multiplication*/
		inline Vector<T,n> operator * (T s)
		{
			Vector<T,n> r;
			for (unsigned int i=0; i < n; i++)
				r.c[i] = c[i] * s;
			return r;
		}
		/*division*/
		inline Vector<T,n> operator / (T s)
		{
			Vector<T,n> r;
			for (unsigned int i=0; i < n; i++)
				r.c[i] = c[i] / s;
			return r;
		}
		/*subtraction*/
		inline void operator -= (Vector<T,n> o)
		{
			for (unsigned int i=0; i < n; i++)
				c[i] -= o.c[i];
		}
		/*addition*/
		inline void operator += (Vector<T,n> o)
		{
			for (unsigned int i=0; i < n; i++)
				c[i] += o.c[i];
		}
		/*substraction*/
		inline Vector<T,n> operator - (Vector<T,n> o)
		{
			Vector<T,n> r;
			for (unsigned int i=0; i < n; i++)
				r.c[i] = c[i] - o.c[i];
			return r;
		}
		/*addition*/
		inline Vector<T,n> operator + (Vector<T,n> o)
		{
			Vector<T,n> r;
			for (unsigned int i=0; i < n; i++)
				r.c[i] = c[i] + o.c[i];
			return r;
		}
		/*increase*/
		inline void increase(T by)
		{
			for (unsigned int i=0; i < n; i++)
				c[i] += by;
		}
		/*sum of the entries*/
		inline T sum()
		{
			T s = 0;
			for (unsigned int i=0; i < n; i++)
				s += c[i];
			return s;
		}
		
		friend inline std::ostream& operator << (std::ostream& output, const Vector<T,n>& v)
		{
			for (unsigned int i=0; i < n; i++)
				output << v.c[i] << " ";
			return output;
		}
	};
}

#endif
