#include <iostream>
#include <fstream>
#include <set>
#include <string>

using namespace std;

int main(int argc, const char* argv[])
{
	if(argc == 1)
	{
		cout
			<< "Usage: " << *argv << " file1 [file2 [file3 [...]]]" << endl
			<< "Function:  Sorts and merges all lines in the specified input " << endl
			<< "files and merges them to stdout" << endl
			<< endl
			<< "Note:  This function must read each of the files completely into memory " << endl
			<< "before it can resolve duplicates.  Ensure that you have enough memory " << endl
			<< "available before proceeding on very large files." << endl;
		return -1;
	}

	typedef set<string> t_stType;
	t_stType rs;
	for(int i = 1; i < argc; i++)
	{
		ifstream in(argv[i]);
		if(!in)
		{
			cerr << "Failed to open input file " << argv[i] << endl;
			continue;
		}

		string line;
		while(!in.eof())
		{
			getline(in, line);
			rs.insert(line);
		}
	}
	
	// Write the results to stdout:
	for(t_stType::const_iterator q = rs.begin(); q != rs.end(); q++)
		cout << *q << endl;
}
