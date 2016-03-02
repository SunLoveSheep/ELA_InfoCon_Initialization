//Header file for sampling functions.

#ifndef _Sampling_H
#define _Sampling_H

class Sampling
{
public:
	Sampling();
	virtual~Sampling();
	//void LHD_IndexSamples(int *OutIndex);//This is the function that generates samples in terms of index. Each index represents which section of that dimension is sampled. Output is saved in *OutIndex
	//void LHD_RealSampling(int *InIndex, double *OutSamples);//This is the function that takes the index samples and each dimension range, and generate real samples. Results are saved in *OutSamples
	//void RNGSampling(double *OutSamples);//This is the function that randomly generate samples

	void SamplingProcess();//This is the sampling process that can be called by other .cpp files. The data.SamplingMethod decides which sampling method to use
};

#endif