/*
This .cpp file contains the headers to build the result output functions
*/

#ifndef _OUTPUT_H_
#define _OUTPUT_H_

class Output
{
public:
	Output();
	virtual~Output();

	void ResultOutput(int funcnum);//to write the result to .csv files
};

#endif