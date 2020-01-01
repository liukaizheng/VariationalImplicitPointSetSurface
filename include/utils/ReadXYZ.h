#ifndef READXYZ_H
#define READXYZ_H

#include <string>
#include <fstream>
#include <sstream>
#include <cctype>
#include <vector>

template <typename Poin3,typename Alloc>
inline bool ReadXYZ(const std::string& name,
        std::vector<Poin3, Alloc>& points)
{
    std::ifstream in(name);
    if(!in) return false;

    points.clear();

    std::string line;
    std::stringstream ss;
    Poin3 p;

    unsigned i = 0;
    while(!in.eof())
    {
        std::getline(in, line);
        if(in.bad()) return false;

        ss.str(line);
        ss.clear();
        ss >> std::ws;
        if(ss.eof()) continue;
        auto c = ss.peek();
        if(!std::isdigit(c) && c != '-' && c != '+') continue;

        for(i = 0; i < 3; i++)
        {
            if(ss.eof()) break;
            ss >> p[i];
        }
        if(i != 3) continue;

        points.emplace_back(p);
    }
    return true;
}
#endif
