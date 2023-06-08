
#include "util.h"

#include <cstdlib>
#include <ctime>

bool Util::seeded = false;

std::vector<std::pair<int, int>> Util::readFile(const std::string& filename)
{
	std::vector<std::pair<int, int>> points;

    std::ifstream ifs(filename);

    if (ifs.is_open()) {
        int id, x, y;

        while (ifs >> id >> x >> y)
        {
            points.push_back({x, y});
        }

        ifs.close();
    }

    return points;
}

void Util::writeFile(const std::string& filename, std::vector<int>& path, uint64_t distance)
{
    std::ofstream ofs(filename);

    if (ofs.is_open()) {
        ofs << distance << std::endl;

        for (auto i : path)
            ofs << i << std::endl;

        ofs.close();
    }
}

int Util::randi()
{
    if (!seeded) {
        std::srand((unsigned) time(NULL));
        seeded = true;
    }

    return std::rand();
}

int Util::randi(int start, int end)
{
    if (start > end)
        std::swap(start, end);

    return (randi() % (end - start + 1)) + start;
}

double Util::rand()
{
    return ((double) randi() / (RAND_MAX));
}
