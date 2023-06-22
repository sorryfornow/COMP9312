#include <iostream>
#include <vector>
#include <unordered_map>
#include <utility>
#include <set>
#include <queue>
#include <sstream>

template <typename Vertex>
class disjoint_set {
private:
    std::unordered_map<Vertex, Vertex> parent{};
    std::unordered_map<Vertex, int> size_{};
public:
    explicit disjoint_set(const std::unordered_map<Vertex, std::set<Vertex>>& neighbors) {
        for (const auto &[v,_]: neighbors) {
            parent[v] = v;
            size_[v] = 1;
        }

        std::size_t union_cnt = 1;
        for (const auto &[u, nbrs]: neighbors) {
            for (const auto &v: nbrs) {
                std::cout << "\nunion_cnt: " << union_cnt++ << std::endl;
                std::cout << "u: " << u << ", v: " << v << std::endl;
                union_set(u, v);
                std::cout << showUp() << std::endl;
            }
        }
    }

    explicit disjoint_set(const std::vector<std::pair<Vertex,Vertex>>& edges){
        for (const auto &[u,v]: edges) {
            parent[u] = u;
            parent[v] = v;
            size_[u] = 1;
            size_[v] = 1;
        }
        std::size_t union_cnt = 1;
        for (const auto &[u,v]: edges) {
            std::cout << "\nunion_cnt: " << union_cnt++ << std::endl;
            std::cout << "u: " << u << ", v: " << v << std::endl;
            union_set(u, v);
            std::cout << showUp() << std::endl;
        }
    }

    void union_set(Vertex u, Vertex v) {    // union by rank
        auto u_parent = find_set(u);
        auto v_parent = find_set(v);
        if (u_parent == v_parent) return;
        if (size_[u_parent] <= size_[v_parent]) {
            parent[u_parent] = v_parent;
            size_[v_parent] += size_[u_parent];
        } else if (size_[u_parent] > size_[v_parent]) {
            parent[v_parent] = u_parent;
            size_[u_parent] += size_[v_parent];
        } else {
            throw std::runtime_error("union_set() error");
        }

    }

    Vertex find_set(Vertex u) {
        if (parent[u] == u) return u;
        parent[u] = find_set(parent[u]);
        return parent[u];
    }

    // current tree print
    std::string showUp() {
        std::stringstream ss;
        ss << "son -> parent" << std::endl;
        for (const auto &[u,_]: parent) {
            ss << u << " -> " << parent[u] << std::endl;
        }
        return ss.str();
    }
};

template <typename Vertex>
class undirectedUnweightedGraph {
    using Edge = std::pair<Vertex, Vertex>;
private:
    std::unordered_map<Vertex, std::set<Vertex>> neighbors{};
public:

    explicit undirectedUnweightedGraph(std::vector<Edge> edges){
        for (const auto &[u,v]: edges) {
            neighbors[u].insert(v);
            neighbors[v].insert(u);
        }
    }

    std::unordered_map<Vertex, std::set<Vertex>> get_neighbors() {
        return neighbors;
    }

    std::vector<Vertex> bfs(Vertex start) {
        std::vector<Vertex> result;
        std::queue<Vertex> q;
        std::unordered_map<Vertex, bool> visited;
        q.push(start);
        visited[start] = true;
        while (!q.empty()) {
            auto u = q.front();
            q.pop();
            result.push_back(u);
            for (const auto &v: neighbors[u]) {
                if (!visited[v]) {
                    q.push(v);
                    visited[v] = true;
                }
            }
        }
        return result;
    }

    void dfs_helper(Vertex u, std::unordered_map<Vertex, bool> &visited, std::vector<Vertex> &result) {
        visited[u] = true;
        result.push_back(u);
        for (const auto &v: neighbors[u]) {
            if (!visited[v]) {
                dfs_helper(v, visited, result);
            }
        }
    }
    std::vector<Vertex> dfs(Vertex start) {
        std::vector<Vertex> result;
        std::unordered_map<Vertex, bool> visited;
        dfs_helper(start, visited, result);
        return result;
    }

    std::vector<Vertex> topological_sort() {
        std::vector<Vertex> result;
        std::unordered_map<Vertex, bool> visited;
        for (const auto &[u,_]: neighbors) {
            if (!visited[u]) {
                dfs_helper(u, visited, result);
            }
        }
        std::reverse(result.begin(), result.end());
        return result;
    }
};


int main() {
    std::cout << "Hello, World!" << std::endl;
    std::vector<std::pair<char, char>> Fig1 = { {'M','P'},{'P','N'},{'M','N'},

                                                {'I','A'},{'I','C'}, {'I','G'},
                                                {'A','C'},{'A','B'}, {'A','E'},
                                                {'C','B'}, {'C','F'}, {'B','G'},
                                                {'G','F'}, {'E','F'}, {'F','B'},

                                                {'D','H'}, {'D','K'}, {'H','J'},
                                                {'H','K'}, {'J','K'}
    };
    std::vector<std::pair<char, char>> Fig2 = { {}

    };

    undirectedUnweightedGraph<char> g1(Fig1);

    disjoint_set<char> ds1(Fig1);
//    disjoint_set<char> ds1(g1.get_neighbors());


    return 0;
}
