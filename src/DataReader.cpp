#include "../include/DataReader.h"
#include "../include/Debug.h"

DataReader::DataReader()
{
#ifdef DEBUG
	std::cerr << "DEFAULT CONSTRUCTOR - DataReader" << std::endl;
#endif
}

DataReader::DataReader(std::string path)
{
#ifdef DEBUG
	std::cerr << "INIT CONSTRUCTOR - DataReader(path)" << std::endl;
#endif
	pathToFile = path;
	if (!(file.open(pathToFile.c_str(), std::ios::in)))
	{
#ifdef DEBUG
		std::cerr << "ERROR DataReader - Constructor: couldnt open file " << path.c_str() << std::endl;
#endif
		std::string err = "DataReader::DataReader(string): The path specified does not point to a valid file location: "+path;
		throw std::invalid_argument(err);
		return;
	}
	input = new std::istream(&file);
}

DataReader::~DataReader()
{
#ifdef DEBUG
	std::cerr << "DESTRUCTOR - DataReader" << std::endl;
#endif
	if (file.is_open())
	{
		delete input;
		file.close();
	}
}

void DataReader::setPath(std::string path)
{
#ifdef DEBUG
	std::cerr << "setPath(path) - DataReader" << std::endl;
#endif
	pathToFile = path;
	if (!(file.open(pathToFile.c_str(), std::ios::in)))
	{
		std::cerr << "ERROR DataReader - setPath: couldnt open file " << path.c_str() << std::endl;
		return;
	}
	input = new std::istream(&file);
}

std::string DataReader::getLine()
{
#ifdef DEBUG
	std::cerr << "getLine() - DataReader" << std::endl;
#endif
	std::string ret;
	if (pathToFile.empty())
	{
		std::cerr << "ERROR DataReader - getLine: No Path assigned" << std::endl; 
		return ret;
	}
	if (file.is_open())
	{
		if (!(input->good()))
		{
			if (!(input->eof()))
				std::cerr << "ERROR DataReader - getLine: Inputstream not good anymore" << std::endl; 
			return ret;
		} else 
		{
			std::getline(*input, ret);
		}
	}
	return ret;
}
