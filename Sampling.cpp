//This is the functions that generate samples

#include "Sampling.h"
#include "Data.h"
#include <iostream>
#include <list>

extern Data data;

using namespace std;

Sampling::Sampling()
{
}

Sampling::~Sampling()
{
}

//this function is used to return a random double variable within the range of min and max
double RandomGen(double min, double max)
{
	double Min = (min*1000000);
    double Max = (max*1000000);
    double Rand = rand()*rand();
    double Result;

	if (min!=max)
	{
		Result = ((int)(Rand)%(int)(Max-Min))+Min;
	}
	else
	{
		Result=min;
	}
    return Result/1000000.0;
}

//this function returns a randomly swapped integer vector containing values from Min to Max. Min should be 0 and Max is data.NumSample-1
void RandomIntVec(int Min, int Max, int *Vec)
{
	for (int i=0;i<data.NumSample;i++){
		Vec[i]=i;
	}

	//swap two random numbers in the vec for NumSample times
	int swap_time = 0;
	int swap_temp = 0;
	while (swap_time<data.NumSample){
		int rNum1 = (int)(rand()*rand())%(Max+1-Min)+Min-1;//random generate two integers
		int rNum2 = (int)(rand()*rand())%(Max+1-Min)+Min-1;//random generate two integers
		//swap the numbers on that integer position of Vec
		swap_temp = Vec[rNum1];
		Vec[rNum1] = Vec[rNum2];
		Vec[rNum2] = swap_temp;

		swap_time++;//counter update
	}
}

//this function generate index samples for LHD
void LHD_IndexSamples(int **OutIndex)
{
	int *temp_vec = new int[data.NumSample];
	
	for (int d=0;d<data.D;d++)
	{
		RandomIntVec(1,data.NumSample,temp_vec);
		for (int n=0;n<data.NumSample;n++)
		{
			OutIndex[n][d] = temp_vec[n];
		}
	}

	delete []temp_vec;
}

void LHD_RealSampling(int **InIndex, double **OutSamples)//this function takes index samples and generate actual samples by means of LHD
{
	for (int d=0;d<data.D;d++)
	{
		double Sec_Length = (double)((data.Max[d]-data.Min[d]))/(double)(data.NumSample);//to calculate the length of each section, dividing according to the number of samples
		
		for (int n=0;n<data.NumSample;n++)
		{
			//randomly generate real sample values according to index.
			double randnum = RandomGen(0.0,Sec_Length);
			OutSamples[n][d] = data.Min[d] + (InIndex[n][d]) * Sec_Length +randnum;
		}
	}
}

void RNGSampling(double **OutSamples)//this function randomly generate samples
{
	for (int n=0;n<data.NumSample;n++)
	{
		for (int d=0;d<data.D;d++)
		{
			OutSamples[n][d] = RandomGen(data.Min[d],data.Max[d]);
		}
	}
}

//LHD only generates samples. This function sorts the samples in the sense of nearest neighbors put together.
void NearestSorting(double **InSamples, double **OutSequences)
{
	//thought: Randomly pick the 1st point, then calculate Euclidean distances from all points, pick the shortest one as the 2nd.
	//Eliminate the selected points and repeat this procedure to select all points into the sequence.

	//Randomly select the first point:
	int Index_FirstSample = (int)(RandomGen(0,data.NumSample-1));
	int Current_Sample = Index_FirstSample; // record the current sample being considered
	int NumRemainSample = data.NumSample;//to record the number of unsorted samples
	while (NumRemainSample>0)
	{
		//add the selected sample to the current position of the sequence
		//and eliminating the current evaluated sample
		for (int d=0;d<data.D;d++)
		{
			//add the selected sample to the current position of the sequence
			OutSequences[data.NumSample - NumRemainSample][d] = InSamples[Current_Sample][d];
			
			//and eliminating the current evaluated sample
			for (int n=Current_Sample;n<NumRemainSample-1;n++)
			{
				InSamples[n][d] = InSamples[n+1][d];
			}
			InSamples[NumRemainSample-1][d] = 0;
		}

		//reduce number of unsorted samples by 1
		NumRemainSample -- ;
		//calculate Euclidean Distance from all other remained samples from the current one and find the shorted distanced one
		double *EuclideanSum = new double[NumRemainSample];//here store the Euclidean distance from current sample to all remaining samples respectively
		double Min_Distance = 100000000000000.0;
		for (int n=0;n<NumRemainSample;n++)
		{
			EuclideanSum[n]=0;
			for (int d=0;d<data.D;d++)
			{
				EuclideanSum[n]+=(OutSequences[data.NumSample - NumRemainSample-1][d] - InSamples[n][d])
					*(OutSequences[data.NumSample - NumRemainSample-1][d] - InSamples[n][d]);
			}

			//compare and find the min distance
			if (EuclideanSum[n]<Min_Distance)
			{
				Min_Distance = EuclideanSum[n];
				//update the current sample being considered
				Current_Sample = n;
			}
		}

		delete []EuclideanSum;
	}
}

void Sampling::SamplingProcess()
{	
	if (data.SamplingMethod == "LHD")//sampling using LHD
	{
		//initiate variables
		int **LHD_index = new int*[data.NumSample];
		double **LHD_samples = new double*[data.NumSample];
		for (int i=0;i<data.NumSample;i++)
		{
			LHD_index[i] = new int[data.D];
			LHD_samples[i] = new double[data.D];
			for (int d=0;d<data.D;d++)
			{
				LHD_index[i][d]=0;
				LHD_samples[i][d]=0;
			}
		}

		//LHD indexing and sampling
		LHD_IndexSamples(LHD_index);
		LHD_RealSampling(LHD_index,LHD_samples);

		/*cout<<"Unsorted: "<<endl;
		for (int i=0;i<data.NumSample;i++)
		{
			for (int d=0;d<data.D;d++)
			{
				cout<<LHD_samples[i][d]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;*/

		if (data.SequencingMethod == "NN")
		{
			//LHD samples are randomly generated, use nearest neighbor to sort and sequencing samples
			NearestSorting(LHD_samples, data.SampleSequence);
		}
		else if (data.SequencingMethod == "RNG")
		{
			//LHD samples are already randomly generated. Directly copy to output sequence
			for (int i=0;i<data.NumSample;i++)
			{
				for (int d=0;d<data.D;d++)
				{
					data.SampleSequence[i][d] = LHD_samples[i][d];
				}
			}
		}

		/*cout<<endl<<"Sorted: "<<endl;
		for (int i=0;i<data.NumSample;i++)
		{
			for (int d=0;d<data.D;d++)
			{
				cout<<data.SampleSequence[i][d]<<" ";
			}
			cout<<endl;
		}*/

		//release RAM
		for (int i=0;i<data.NumSample;i++)
		{
			delete []LHD_index[i];
			delete []LHD_samples[i];
		}
		delete []LHD_index;
		delete []LHD_samples;
	}
	else if (data.SamplingMethod == "RNG")//sampling using RNG
	{
		RNGSampling(data.SampleSequence);
	}
}