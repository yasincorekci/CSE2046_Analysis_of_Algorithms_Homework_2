
#include "tsp.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cfloat>

#include "util.h"

TSP::TSP(std::vector<std::pair<int, int>>& points) : points(points), isLargeData(false)
{
	double d, dx, dy;

	numOfPoints = (int) points.size();

    if (numOfPoints > 3000) {
        isLargeData = true;
        return;
    }

	for (int i = 0; i < (int) points.size(); ++i) {
		adjPointToDistance[i] = std::map<int, double>();
		adjDistanceToPoint[i] = std::multimap<double, int>();
	}

	for (int i = 0; i < (int) points.size(); ++i) {
        std::cout << "Creating " << (i + 1) << std::endl;

		for (int j = i + 1; j < (int) points.size(); ++j) {
			dx = points[i].first - points[j].first;
			dy = points[i].second - points[j].second;

			d = std::round(std::sqrt(dx * dx + dy * dy));
			
			adjPointToDistance[i][j] = d;
			adjPointToDistance[j][i] = d;

			adjDistanceToPoint[i].insert({d, j});
			adjDistanceToPoint[j].insert({d, i});
		}
	}
}


void TSP::run(double mutationRate, int poolSize, double selectRate)
{
    Path bestPath(std::vector<int>(), DBL_MAX, 0);

    int loopCount = 2000;

    if (numOfPoints == 0)
        return;

    if (numOfPoints > 1000)
        loopCount = 1000;

    std::vector<int> bestInitials;
    int limit = 1;

    if (!isLargeData) {
        bestInitials = calculateBestInitials();
    } 
    else {
        bestInitials.resize(limit, 0);
        poolSize = 10;
    }

    for (int t = 0; t < limit; ++t) {
        std::vector<Path> pathPool;

        std::vector<int> bestHalf = extractHalf(bestInitials[t]);

        for (int i = 0; i < poolSize; ++i) {
            std::vector<int> vertexes = bestHalf;
            shuffle(vertexes);

            double dist = distance(vertexes);
            double rank = 1.0 / dist;

            pathPool.push_back(Path(vertexes, dist, rank));
        }

        int selectSize = (int)(poolSize * selectRate);

        for (int i = 0; i < loopCount; i++) {
            std::cout << "Loop " << (i + 1) << std::endl;

            std::sort(pathPool.begin(), pathPool.end(),
                [](const Path& p1, const Path& p2) -> bool
                {
                    return p1.rank > p2.rank;
                });

            if (pathPool[0].distance < bestPath.distance) {
                bestPath = pathPool[0];
            }

            pathPool = selectPopulation(pathPool, selectSize);

            for (int k = 0; k < poolSize - selectSize; ++k) {
                int r1 = Util::randi(0, selectSize - 1);
                int r2 = Util::randi(0, selectSize - 1);

                auto vertexes = crossover(pathPool[r1].vertexes, pathPool[r2].vertexes);
                double dist = distance(vertexes);
                double rank = 1.0 / dist;

                pathPool.push_back(Path(vertexes, dist, rank));
            }

            for (int k = 0; k < selectSize; ++k) {
                mutate(pathPool[k].vertexes, mutationRate);

                pathPool[k].distance = distance(pathPool[k].vertexes);
                pathPool[k].rank = 1.0 / pathPool[k].distance;
            }
        }
    }

    Util::writeFile("output.txt", bestPath.vertexes, (uint64_t) std::round(bestPath.distance));

    for (auto i : bestPath.vertexes)
        std::cout << i << " ";

    std::cout << std::endl;
    std::cout << "Distance: " << distance(bestPath.vertexes) << std::endl;
}

double TSP::rank(std::vector<int>& vertexes)
{
    return 1.0 / distance(vertexes);
}

double TSP::distance(std::vector<int>& vertexes)
{
    double total = 0.0;

    if (isLargeData) {
        double dx, dy, d;

        for (int i = 0; i < (int)vertexes.size() - 1; ++i) {
            dx = points[vertexes[i]].first - points[vertexes[i+1]].first;
            dy = points[vertexes[i]].second - points[vertexes[i+1]].second;

            d = std::round(std::sqrt(dx * dx + dy * dy));

            total += d;
        }

        dx = points[vertexes[(int)vertexes.size() - 1]].first - points[vertexes[0]].first;
        dy = points[vertexes[(int)vertexes.size() - 1]].second - points[vertexes[0]].second;

        d = std::round(std::sqrt(dx * dx + dy * dy));

        total += d;
    }
    else {
        total += adjPointToDistance[vertexes[(int)vertexes.size() - 1]][vertexes[0]];

        for (int i = 0; i < (int)vertexes.size() - 1; ++i) {
            total += adjPointToDistance[vertexes[i]][vertexes[i + 1]];
        }
    }

    return total;
}

std::vector<int> TSP::crossover(std::vector<int>& vertexes1, std::vector<int>& vertexes2)
{
    std::vector<int> newVertexes(vertexes1.size());

    int start = Util::randi(1, (int) vertexes1.size() - 1);
    int end = Util::randi(1, (int) vertexes1.size() - 1);

    if (start > end)
        std::swap(start, end);

    std::vector<int> subVertexes;

    for (int i = start; i <= end; ++i)
        subVertexes.push_back(vertexes1[i]);

    int i = 1;
    int current = 1;

    newVertexes[0] = vertexes1[0];

    while (i < start) {
        int item = vertexes2[current];

        if (std::find(subVertexes.begin(), subVertexes.end(), item) == subVertexes.end()) {
            newVertexes[i] = item;
            ++i;
        }

        ++current;
    }

    for (int i = start; i <= end; ++i)
        newVertexes[i] = subVertexes[i-start];

    i = end + 1;

    while (i < (int) vertexes1.size()) {
        int item = vertexes2[current];

        if (std::find(subVertexes.begin(), subVertexes.end(), item) == subVertexes.end()) {
            newVertexes[i] = item;
            ++i;
        }

        ++current;
    }

    return newVertexes;
}

void TSP::mutate(std::vector<int>& vertexes, double mutationRate)
{
    for (int i = 1; i < (int) vertexes.size(); ++i) {
        if (Util::rand() < mutationRate) {
            int j = Util::randi(1, (int) vertexes.size() - 1);

            int item = vertexes[i];
            vertexes[i] = vertexes[j];
            vertexes[j] = item;
        }
    }
}

void TSP::selectPopulation(std::vector<Path>& paths, std::vector<Path>& newPaths, int selectSize)
{
    if (selectSize == 1) {
        newPaths.push_back(paths[0]);
        paths.erase(paths.begin());
        return;
    }
    
    double totalRank = 0.0;

    for (int i = 0; i < (int) paths.size(); ++i)
        totalRank += paths[i].rank;

    int selected = 0;
    double total = 0.0;

    double selectRandom = Util::rand() * totalRank;

    for (auto& path : paths) {
        total += path.rank;

        if (total > selectRandom)
            break;

        ++selected;
    }

    if (selected >= (int) paths.size())
        selected = (int) paths.size() - 1;

    newPaths.push_back(paths[selected]);
    paths.erase(paths.begin() + selected);

    selectPopulation(paths, newPaths, selectSize - 1);

    return;
}

void TSP::shuffle(std::vector<int>& vertexes)
{
    for (int i = 0; i < (int) vertexes.size() * 16; i++)
    {
        int first = Util::randi(1, (int) vertexes.size() - 1);
        int second = Util::randi(1, (int)vertexes.size() - 1);

        int t = vertexes[first];
        vertexes[first] = vertexes[second];
        vertexes[second] = t;
    }
}

std::vector<int> TSP::extractBestHalf()
{
    std::vector<int> vertexes;

    int minIndex = 0;
    double minDistance = DBL_MAX;

    int i;
    int limit = (int) adjDistanceToPoint.size() / 2;

    for (auto& kv : adjDistanceToPoint) {
        double total = 0.0;

        i = 0;

        for (auto& e : kv.second) {
            if (i++ == limit)
                break;

            total += e.first;
        }

        if (total < minDistance) {
            minDistance = total;
            minIndex = kv.first;
        }
    }

    vertexes.push_back(minIndex);

    i = 1;
    
    for (auto& e : adjDistanceToPoint[minIndex]) {
        if (i++ == limit)
            break;

        vertexes.push_back(e.second);
    }

    return vertexes;
}

std::vector<int> TSP::calculateBestInitials()
{
    std::multimap<double, int> distToPoints;
    std::vector<int> bestInitials;

    int i;
    int limit = (int) adjDistanceToPoint.size() / 2;

    for (auto& kv : adjDistanceToPoint) {
        double total = 0.0;

        i = 1;
        limit = (int) kv.second.size() / 2;

        for (auto& e : kv.second) {
            if (i++ == limit)
                break;

            total += e.first;
        }

        distToPoints.insert({ total, kv.first });
    }

    for (auto& e : distToPoints) {
        bestInitials.push_back(e.second);
    }

    return bestInitials;
}

std::vector<int> TSP::extractHalf(int index)
{
    std::vector<int> vertexes;

    if (isLargeData) {
        std::vector<int> v;

        for (int i = 0; i < (int) points.size(); i++) {
            v.push_back(i);
        }

        for (int i = 0; i < (int) v.size() * 16; i++)
        {
            int first = Util::randi(0, (int) v.size() - 1);
            int second = Util::randi(0, (int) v.size() - 1);

            int t = v[first];
            v[first] = v[second];
            v[second] = t;
        }

        for (int i = 0; i < (int) std::ceil((double) v.size() / 2); i++) {
            vertexes.push_back(v[i]);
        }
    }
    else {
        vertexes.push_back(index);

        int i = 1;
        int limit = (int)adjDistanceToPoint.size() / 2;

        for (auto& e : adjDistanceToPoint[index]) {
            if (i++ == limit)
                break;

            vertexes.push_back(e.second);
        }
    }

    return vertexes;
}

std::vector<Path> TSP::selectPopulation(std::vector<Path>& paths, int selectSize)
{
    std::vector<Path> newPaths;

    if ((int) paths.size() <= selectSize)
        return paths;

    selectPopulation(paths, newPaths, selectSize);

    return newPaths;
}
