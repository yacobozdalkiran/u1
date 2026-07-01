#include "geometry.h"

Geometry::Geometry(int L_) : L(L_), V(L * L), neighbors(4 * V) {
    for (int x = 0; x < L; x++) {
        for (int t = 0; t < L; t++) {
            int site = index_site(x, t);
            neighbors[index_neigbors(site, 0, pos)] = index_site((x + 1) % L, t);
            neighbors[index_neigbors(site, 0, neg)] = index_site((x - 1 + L) % L, t);
            neighbors[index_neigbors(site, 1, pos)] = index_site(x, (t + 1) % L);
            neighbors[index_neigbors(site, 1, neg)] = index_site(x, (t - 1 + L) % L);
        }
    }

    links_staples.resize(V * 2 * 2 * 3, std::make_pair(-1, -1));
    for (int x = 0; x < L; x++) {
        for (int t = 0; t < L; t++) {
            int site = index_site(x, t);
            for (int mu = 0; mu < 2; mu++) {
                int nu = 1-mu;
                int xmu = get_neighbor(site, mu, pos);     // x+mu
                int xnu = get_neighbor(site, nu, pos);     // x+nu
                int xmunu = get_neighbor(xmu, nu, neg);  // x+mu-nu
                int xmnu = get_neighbor(site, nu, neg);  // x-nu

                links_staples[index_staples(site, mu, 0, 0)] = {xmu, nu};
                links_staples[index_staples(site, mu, 0, 1)] = {xnu, mu};
                links_staples[index_staples(site, mu, 0, 2)] = {site, nu};
                links_staples[index_staples(site, mu, 1, 0)] = {xmunu, nu};
                links_staples[index_staples(site, mu, 1, 1)] = {xmnu, mu};
                links_staples[index_staples(site, mu, 1, 2)] = {xmnu, nu};
            }
        }
    }
};

int Geometry::index_site(int x, int t) const { return t * L + x; };
int Geometry::index_neigbors(int site, int mu, dir d) const { return site * 4 + mu * 2 + d; };
int Geometry::get_neighbor(int site, int mu, dir d) const {
    return neighbors[index_neigbors(site, mu, d)];
};
int Geometry::index_staples(int site, int mu, int i_staple, int i_link) const {
    return site * 3 * 2 * 2 + mu * 3 * 2 + i_staple * 3 + i_link;
}
std::pair<int, int> Geometry::get_link_staple(int site, int mu, int i_staple, int i_link) const {
    return links_staples[index_staples(site, mu, i_staple, i_link)];
}
