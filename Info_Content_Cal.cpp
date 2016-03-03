#include "Info_Content_Cal.h"
#include "Data.h"
#include <iostream>
#include <math.h>

using namespace std;
extern Data data;

Info_Content_Cal::Info_Content_Cal()
{}

Info_Content_Cal::~Info_Content_Cal()
{}

//This function calculates the 1-,0,1 vector given the current epsilon value, sample sequence and objective function values
void CalPsiVector(int *psivector, double epsilon)
{
	for (int n=0;n<data.NumSample-1;n++)
	{
		//calculation of difference of control variable distances and objective function values between conjective samples in sample sequence
		double DeltaY = data.ObjectiveSequence[n+1]-data.ObjectiveSequence[n];
		double DeltaX = 0;
		for (int d=0;d<data.D;d++)
		{
			DeltaX += sqrt((data.SampleSequence[n+1][d]-data.SampleSequence[n][d])*(data.SampleSequence[n+1][d]-data.SampleSequence[n][d]));
		}
		
		//judge if the current bit of psivector should be 1-, 0, or 1.
		if (DeltaY/DeltaX < -epsilon)
			psivector[n] = -1;
		else if ((DeltaY/DeltaX < epsilon) && (DeltaY/DeltaX > -epsilon))
			psivector[n] = 0;
		else if (DeltaY/DeltaX > epsilon)
			psivector[n] = 1;
	}
}

//These two functions calculate corresponding H and M values given the psi vector
double CalH(int *psivector)
{
	double H = 0;

	//define 6 counters to count all 6 combinations of 1-,0,1. M for -1, Z for 0, P for 1. Situation that Either value appears in two consecutive positions is not considered.
	int CountMZ = 0, CountMP = 0, CountZP = 0, CountZM = 0, CountPZ = 0, CountPM = 0;

	for (int i=0;i<data.NumSample-2;i++)//need to minus 2. Since the length of 1-,0,1 sequence is already data.NumSample-1
	{
		if (psivector[i] == -1)
		{
			if (psivector[i+1] == 0)
				CountMZ ++;
			else if (psivector[i+1] == 1)
				CountMP ++;
		}
		else if (psivector[i] == 0)
		{
			if (psivector[i+1] == -1)
				CountZM ++;
			else if (psivector[i+1] == 1)
				CountZP ++;
		}
		else if (psivector[i] == 1)
		{
			if (psivector[i+1] == -1)
				CountPM ++;
			else if (psivector[i+1] == 0)
				CountPZ ++;
		}
	}
	
	int CountSum = CountMZ + CountMP + CountZP + CountZM + CountPM + CountPZ;
	double LOGCountMZ = (CountMZ==0) ? 0:(log(CountMZ));//if Count = 0, then set its log to 0, otherwise use its log value
	double LOGCountMP = (CountMP==0) ? 0:(log(CountMP));//if Count = 0, then set its log to 0, otherwise use its log value
	double LOGCountZP = (CountZP==0) ? 0:(log(CountZP));//if Count = 0, then set its log to 0, otherwise use its log value
	double LOGCountZM = (CountZM==0) ? 0:(log(CountZM));//if Count = 0, then set its log to 0, otherwise use its log value
	double LOGCountPM = (CountPM==0) ? 0:(log(CountPM));//if Count = 0, then set its log to 0, otherwise use its log value
	double LOGCountPZ = (CountPZ==0) ? 0:(log(CountPZ));//if Count = 0, then set its log to 0, otherwise use its log value
	//check if CountSum is 0. If so, epsilon is big enough, stop and return
	//if not 0, calculate H value
	H = (CountSum == 0) ? 0 : -((CountMZ/CountSum) * (LOGCountMZ/log(6)) - ((double)(CountMP)/(double)CountSum) * (LOGCountMP/log(6)) - ((double)CountZP/(double)CountSum) * (LOGCountZP/log(6))
		- ((double)CountZM/(double)CountSum) * (LOGCountZM/log(6)) - ((double)CountPM/(double)CountSum) * (LOGCountPM/log(6)) - ((double)CountPZ/(double)CountSum) * (LOGCountPZ/log(6)));
	
	/*cout<<CountMZ <<" "<< LOGCountMZ<< " "<<log(6) <<" "<< LOGCountMZ/log(6)<<endl;
	cout<<CountMP <<" "<< LOGCountMP<< " "<<log(6) <<" "<< LOGCountMP/log(6)<<endl;
	cout<<CountZP <<" "<< LOGCountZP<< " "<<log(6) <<" "<< LOGCountZP/log(6)<<endl;
	cout<<CountZM <<" "<< LOGCountZM<< " "<<log(6) <<" "<< LOGCountZM/log(6)<<endl;
	cout<<CountPM <<" "<< LOGCountPM<< " "<<log(6) <<" "<< LOGCountPM/log(6)<<endl;
	cout<<CountPZ <<" "<< LOGCountPZ<< " "<<log(6) <<" "<< LOGCountPZ/log(6)<<endl;
	cout<<CountSum<<" "<<H<<endl;
	getchar();*/

	return H;
}
double CalM(int *psivector)
{
	int *nonzero_psi_vector = new int[data.NumSample-1];
	int *norepeat_psi_vector = new int[data.NumSample-1];
	int nonzero_vec_length = 0;
	//First, remove all 0s
	for (int i=0;i<data.NumSample-1;i++)
	{
		if (psivector[i] != 0)
		{
			nonzero_psi_vector[nonzero_vec_length] = psivector[i];
			nonzero_vec_length++;//copy the non-zero value into temp_psi_vector and increase its length by 1
		}
	}

	//Then, remove all repeated symbols
	bool Repeat = false;
	int repeat_check = 1;
	int no_repeat_counter = 0;
	while (repeat_check < nonzero_vec_length)
	{
		//if the two consecutive position have equal bits
		if (nonzero_psi_vector[repeat_check] - nonzero_psi_vector[repeat_check-1] == 0)
		{
			if (repeat_check == 1)//first bit
			{
				norepeat_psi_vector[no_repeat_counter] = nonzero_psi_vector[repeat_check-1];
				no_repeat_counter++;
			}
			else//not the first bit
			{
				//do nothing
			}
			repeat_check++;
		}
		else //if the two consecutive bits are different:
		{
			if (repeat_check == 1)//first bit
			{
				norepeat_psi_vector[no_repeat_counter] = nonzero_psi_vector[repeat_check-1];
				no_repeat_counter++;
				norepeat_psi_vector[no_repeat_counter] = nonzero_psi_vector[repeat_check];
				no_repeat_counter++;
			}
			else//not the first bit
			{
				norepeat_psi_vector[no_repeat_counter] = nonzero_psi_vector[repeat_check];
				no_repeat_counter++;
			}
			repeat_check++;
		}
	}

	//Now, the norepeat_psi_vector should contains the psi without zero and repeated bits.
	
	//calculate M:
	double M=0;
	M = (double)no_repeat_counter/(double)(data.NumSample-1);

	delete []nonzero_psi_vector;
	delete []norepeat_psi_vector;

	return M;
}

//Functions to find the 4 featured values: (note that FindM0 has no content, since M0 = MVector[0])
double FindHMax(double *HVector)
{
	double HMax = -100000000000000;
	for (int i=0;i<data.NumEpsilon+1;i++)
	{
		if (HMax < HVector[i])
			HMax = HVector[i];
	}

	return HMax;
}
double FindEpsilon_s(double *HVector, double epsilon_step)
{
	//loop through the HVector to find the first position that has H(epsilon) < 0.05
	double Epsilon_s = 0;
	int PositionEpsilon_s = 0;

	for (int i=0;i<data.NumEpsilon+1;i++)
	{
		if (HVector[i]<0.05)
		{
			PositionEpsilon_s = i;
			break;
		}
	}

	//calculate the correspondent epsilon value
	if (PositionEpsilon_s == 0)
		Epsilon_s = 0;
	else
		Epsilon_s = data.EpsilonMin + epsilon_step * (PositionEpsilon_s-1);

	Epsilon_s = log10(Epsilon_s);

	return Epsilon_s;
}
double FindM0()
{}
double FindEpsilon_05(double *MVector, double epsilon_step)
{
	int PositionEpsilon_05 = 0;
	double Epsilon_05 = 0;
	
	for (int i=data.NumEpsilon+1;i>0;i--)
	{
		if (MVector[i] > 0.5 * MVector[0])
		{
			PositionEpsilon_05 = i;
			break;
		}
	}

	//calculate the correspondent epsilon value
	if (PositionEpsilon_05 == 0)
		Epsilon_05 = 0;
	else
	{
		Epsilon_05 = data.EpsilonMin + epsilon_step * (PositionEpsilon_05-1);
		Epsilon_05 = log10(Epsilon_05);
	}

	//Epsilon_05 = log10(Epsilon_05);

	return Epsilon_05;
}

//Given a function and a sample sequence, this function calculate the 4 low level features based on information contents
void Info_Content_Cal::Info_Con_Process(int funcnum)
{
	//Divided the Epsilon range into data.NumEpsilon pieces
	double LogEpsilonMin = log10(data.EpsilonMin); //log base of 10 of Epsilon max
	double LogEpsilonMax = log10(data.EpsilonMax); //log base of 10 of Epsilon min
	double EpsilonStep = (LogEpsilonMax-LogEpsilonMin)/data.NumEpsilon; //dividing steps
	double LogEpsilon = LogEpsilonMin; //current epsilon value (log 10 based)
	double *HVector = new double[data.NumEpsilon+1];//H vector to record all H values, +1 is because need to count Epsilon = 0 as well.
	double *MVector = new double[data.NumEpsilon+1];//M vector to record all M values, +1 is because need to count Epsilon = 0 as well.

	//the 1-,0,1 vector. -1 for 1-, 0 for 0, and 1 for 1
	int *PsiVector = new int[data.NumSample-1];
	for (int i=0;i<data.NumSample-1;i++)
	{
		PsiVector[i]=0;
	}

	int EpsilonCount = 0;
	//First deal with the case Epsilon = 0:
	double Epsilon = 0;
	CalPsiVector(PsiVector,Epsilon);

	double H = CalH(PsiVector);
	if (H == 0)
	{
		cout<<"Epsilon Big Enough! All 0 in Psi Vector."<<endl;
		cout<<"Current epsilon: "<<Epsilon<<endl;
		double BreakEpsilon = 0;
	}
	double M = CalM(PsiVector);
	HVector[EpsilonCount] = H;
	MVector[EpsilonCount] = M;
	EpsilonCount++;

	//Then for epsilon values among all Epsilon values:
	while (LogEpsilon<LogEpsilonMax)
	{
	//{
		//read sample sequence and calculate 1-,0,1 sequence
		Epsilon = pow(10, LogEpsilon);
		CalPsiVector(PsiVector,Epsilon);

		//calculate H and M
		H = CalH(PsiVector);
		if (H == 0)
		{
			cout<<"Epsilon Big Enough! All 0 in Psi Vector."<<endl;
			cout<<"Current epsilon: "<<LogEpsilon<<" "<<Epsilon<<endl;
			double BreakEpsilon = Epsilon;
			break;
		}
		M = CalM(PsiVector);

		//assign H,M into the H_vector and M_vector
		HVector[EpsilonCount] = H;
		MVector[EpsilonCount] = M;

		LogEpsilon += EpsilonStep;
		EpsilonCount++;
	//}
	}

	//From the H_vector and M_vector, conclude the values for the 4 low level features.
	double Hmax = FindHMax(HVector);
	double Epsilon_s = FindEpsilon_s(HVector, EpsilonStep);
	double M0 = MVector[0];
	double Epsilon_05 = FindEpsilon_05(MVector, EpsilonStep);

	/*cout<<"Hmax: "<<Hmax<<endl;
	cout<<"Epsilon_s: "<<Epsilon_s<<endl;
	cout<<"M0 : "<<M0<<endl;
	cout<<"Epsilon_05: "<<Epsilon_05<<endl;*/

	//record the four features.
	data.HMaxVector[funcnum] = Hmax;
	data.Epsilon_sVector[funcnum] = Epsilon_s;
	data.M0Vector[funcnum] = M0;
	data.Epsilon_05Vector[funcnum] = Epsilon_05;

	delete []PsiVector;
	delete []HVector;
	delete []MVector;
}
