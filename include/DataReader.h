#ifndef _DATAREADER_H_
#define _DATAREADER_H_
#include <stdexcept>
#include <iostream>
#include <istream>
#include <fstream>
#include <string>

class DataReader final
{
	private:
		std::string pathToFile;
		std::filebuf file;
		std::istream *input;
	public:
		DataReader();
		DataReader(std::string path);
		~DataReader();

		void setPath(std::string path);
		std::string getLine();
};


#endif
