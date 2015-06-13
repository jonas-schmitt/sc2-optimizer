#include "../include/DataReader.h"

DataReader::DataReader()
{

}

DataReader::DataReader(std::string path)
{

	pathToFile = path;
	if (!(file.open(pathToFile.c_str(), std::ios::in)))
	{

		std::string err = "DataReader::DataReader(string): The path specified does not point to a valid file location: "+path;
		throw std::invalid_argument(err);
		return;
	}
	input = new std::istream(&file);
}

DataReader::~DataReader()
{

	if (file.is_open())
	{
		delete input;
		file.close();
	}
}

void DataReader::setPath(std::string path)
{

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
