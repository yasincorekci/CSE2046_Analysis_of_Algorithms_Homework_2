
#ifndef __TSP_H
#define __TSP_H

#include <map>
#include <vector>
#include <algorithm>
#include <string>

struct Path {
public:
	Path(const std::vector<int>& vertexes, double distance, double rank) : vertexes(vertexes), distance(distance), rank(rank) {}
		
public:
	double rank;
	double distance;
	std::vector<int> vertexes;
};

class TSP {
public:
	TSP(std::vector<std::pair<int, int>>& points);
	void run(double mutationRate = 0.01, int poolSize = 100, double selectRate = 0.5);

private:
	double rank(std::vector<int>& vertexes);
	double distance(std::vector<int>& vertexes);
	std::vector<int> crossover(std::vector<int>& vertexes1, std::vector<int>& vertexes2);
	void mutate(std::vector<int>& vertexes, double mutationRate);
	std::vector<Path> selectPopulation(std::vector<Path>& paths, int selectSize);
	void selectPopulation(std::vector<Path>& paths, std::vector<Path>& newPaths, int selectSize);
	void shuffle(std::vector<int>& vertexes);
	std::vector<int> extractBestHalf();
	std::vector<int> calculateBestInitials();
	std::vector<int> extractHalf(int index);

private:
	std::map<int, std::map<int, double>> adjPointToDistance;
	std::map<int, std::multimap<double, int>> adjDistanceToPoint;
	std::vector<std::pair<int, int>> points;
	int numOfPoints;
	bool isLargeData;
};

#endif
