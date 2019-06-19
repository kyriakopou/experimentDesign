/*
 * misc.h
 *
 *  Created on: Mar 18, 2014
 *      Author: harry
 */

#ifndef MISC_H_
#define MISC_H_


#endif /* MISC_H_ */

int binarySearch(double array[], int first, int last, double search_key){

	int index, mid;
	if (first > last)
		index = last + 1;

	else
	{
		mid = (first + last)/2;
		if (search_key == array[mid])
			index = mid;
		else
			//divide and conquer
			if (search_key < array[mid])
				index = binarySearch(array, first, mid - 1, search_key);
			else
				index = binarySearch(array, mid + 1, last, search_key);

	}


	return index;

}


double linearInterpolation(double x0, double x1, double y0, double y1, double x){

	double y;
	return y = y0 + (y1 - y0) * (x- x0)/ (x1 - x0);


}
