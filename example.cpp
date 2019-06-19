#include <iostream>
#include <fstream>
#include "qtracer.h"
#include "tracer.h"
#include "crystallization.h"



using namespace attraction;

int numOfUnParams = 2;


int main(int argc, char* argv[]) {
	//create the object of tracer of class Qtracer
	QTracer<Crystallization> tracer;

	/*set the truncation threshold*/
	tracer.setDelta(0);


	/*define the arrays that we will use for the output*/
	int M = numOfUnParams;
	double unknownParams[M];
	double paramValues[M];
	
	unknownParams[0] = 0;
	unknownParams[1] = 1;
	paramValues[0] = 4;
	paramValues[1] = 0.1;
	
	double FisherMatrix[M][M], derFisherMatrix[M][M];
	double *p1, *p2;
	p1 = &FisherMatrix[0][0];
	p2 = &derFisherMatrix[0][0];


	/*take the next time point(taken from matlab)*/
	double t = 10;
	int sample_num = 1;



	//define the file to store tmax so far
	double tmax = 0;

	



	tracer.putInitialState(Crystallization::getInitialState(), 1.0, &unknownParams[0], &paramValues[0], M);

	/*compute the distribution at t time units*/
	while(tracer.getTime() < t) {
		double maxdt = t - tracer.getTime();
		tracer.trace(maxdt, sample_num, M, tmax, t, &unknownParams[0], &paramValues[0]);

	}

	//Print the results
	std::cout << "time=" << tracer.getTime() << std::endl;
	tracer.printResults();


	tracer.computeInformation(p1, p2, M);

	//write the new tmax in the file
	std::ofstream t_max;
	t_max.open("data_files/t_max.dat");
	if (t_max.is_open() ){
		t_max << t << std::endl;
		t_max.close();
	}




	std::cout << "Fisher Information Matrix at time point " << t << " is: " << std::endl << std::endl;
	for (int i = 0; i < M; i++){
			for (int j = 0; j < M; j++)
				std::cout << std::setw(10) << FisherMatrix[i][j] <<  std::setw(10)<<'\t';
			std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << "Time Derivative of the Fisher Information Matrix at time point " << t << " is: " << std::endl <<std::endl;
	for (int i = 0; i < M; i++){
		for (int j = 0; j < M; j++)
			std::cout << std::setw(10) << derFisherMatrix[i][j] <<  std::setw(10)<<'\t';
		std::cout << std::endl;
	}

	//clear all the memory
	tracer.clear();




	return 0;
}
