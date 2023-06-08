
#ifndef __UTIL_H
#define __UTIL_H

#include <vector>
#include <fstream>

class Util {
public:
	static std::vector<std::pair<int, int>> readFile(const std::string& filename);
	static void writeFile(const std::string& filename, std::vector<int>& path, uint64_t distance);
	static int randi();
	static int randi(int start, int end);
	static double rand();

private:
	static bool seeded;
};

#endif
