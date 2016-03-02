#ifndef _Info_Content_Cal_H_
#define _Info_Content_Cal_H_

class Info_Content_Cal
{
public:
	Info_Content_Cal();
	virtual~Info_Content_Cal();

	void Info_Con_BinarySequence();//this function will read in an Epsilon value and sample sequence, and output the 1-,0,1 sequence using information content
	void Info_Con_H_M_Vector();//this function will read the whole range of Epsilon and the 1-,0,1 sequence and output the H and M vector for all Epsilon values

	void Info_Con_Process();//this is the general function that returns the four low level features of a given function, range of Epsilon, and sample sequences.
};

#endif