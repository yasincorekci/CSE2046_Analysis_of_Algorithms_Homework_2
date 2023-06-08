
#include <iostream>

#include "tsp.h"
#include "util.h"

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " <input file>" << std::endl;
		return -1;
	}

	auto points = Util::readFile(argv[1]);

	TSP tsp(points);
	tsp.run();

	return 0;
}
