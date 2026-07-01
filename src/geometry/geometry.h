#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

enum dir { pos, neg };

class Geometry {
public:
    int L;
    int V;

private:
    std::vector<int> neighbors;
    std::vector<std::pair<int, int>> links_staples;

public:
    Geometry(int L_);

    int index_site(int x, int t) const;
    int index_neigbors(int site, int mu, dir d) const;
    int index_staples(int site, int mu, int i_staple, int i_link) const;
    int get_neighbor(int site, int mu, dir d) const;
    std::pair<int, int> get_link_staple(int site, int mu, int i_staple, int i_link) const;
};

#endif
