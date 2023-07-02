#include <iostream>
#include <vector>
#include <unordered_map>
#include <utility>
#include <set>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <sstream>

namespace ASMT1 {
    template<typename V, typename T>
    using triplet = std::tuple<V, V, T>;

    template<typename Vertex>
    class disjoint_set {
    private:
        std::unordered_map<Vertex, Vertex> parent{};
        std::unordered_map<Vertex, int> size_{};
    public:
        explicit disjoint_set(const std::unordered_map<Vertex, std::set<Vertex>> &neighbors) {
            for (const auto &[v, _]: neighbors) {
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

        explicit disjoint_set(const std::vector<std::pair<Vertex, Vertex>> &edges) {
            for (const auto &[u, v]: edges) {
                parent[u] = u;
                parent[v] = v;
                size_[u] = 1;
                size_[v] = 1;
            }
            std::size_t union_cnt = 1;
            for (const auto &[u, v]: edges) {
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
            } else [[unlikely]] {
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
            for (const auto &[u, _]: parent) {
                ss << u << " -> " << parent[u] << std::endl;
            }
            return ss.str();
        }
    };

    template<typename Vertex>
    class undirectedUnweightedGraph {
        using Edge = std::pair<Vertex, Vertex>;
    private:
        std::unordered_map<Vertex, std::set<Vertex>> neighbors{};
    public:

        explicit undirectedUnweightedGraph(std::vector<Edge> edges) {
            for (const auto &[u, v]: edges) {
                neighbors[u].insert(v);
                neighbors[v].insert(u);
            }
        }

        /** The function below is used to calculate the distance between two vertices sets.
         * The distance is defined as the number of edges in the shortest path between two vertices sets.
         * The function is implemented by rebuild the adjacency list and using BFS to probe the distance.
         * **/
        std::size_t setDistance(std::unordered_set<Vertex> x, std::unordered_set<Vertex> y) {

            // intersection
            std::unordered_set<Vertex> intersection = std::set_intersection(x.begin(), x.end(), y.begin(), y.end());
            if (!intersection.empty()) {
                return 0;
            }

            // rebuild graph
            std::unordered_map<Vertex, std::unordered_set<Vertex>> neighbors_new{};

            for (auto [u,nbrs]: neighbors) {
                if (x.count(u)) u = 'x';
                if (y.count(u)) u = 'y';
                for (auto v: nbrs) {
                    if (x.count(v)) v = 'x';
                    if (y.count(v)) v = 'y';
                    if (u == v) continue;
                    neighbors_new[u].insert(v);
                    neighbors_new[v].insert(u);
                }
            }

            // BFS
            Vertex start = 'x';
            std::size_t distance = 0;
            std::queue<Vertex> q;
            std::unordered_set<Vertex> visited;

            q.push(start);
            visited.insert(start);

            std::size_t level = 0;

            while (!q.empty()) {
                ++level;
                std::size_t Q_size = q.size();
                for (std::size_t i = 0; i < Q_size; ++i) {
                    auto u = q.front();
                    q.pop();
                    for (const auto &v: neighbors_new[u]) {
                        if (visited.count(v)) continue;
                        q.push(v);
                        visited.insert(v);
                        if (v == 'y') return level;
                    }
                }

            }

            // cannot reach
            return std::numeric_limits<std::size_t>::max();
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


    };


    template<typename Vertex>
    class directedUnweightedGraph {
        using Edge = std::pair<Vertex, Vertex>;
    private:
        std::unordered_map<Vertex, std::set<Vertex>> next_{};
    public:
        explicit directedUnweightedGraph(std::vector<Edge> edges) {
            for (const auto &[u, v]: edges) next_[u].insert(v);
        }

        std::unordered_map<Vertex, std::set<Vertex>> get_next_() {
            return next_;
        }


        std::vector<Vertex> topological_sort() {
            std::vector<Vertex> result;
            std::unordered_set<Vertex> visited;
            std::queue<Vertex> Q;
            std::unordered_map<Vertex, int> in_degree;

            for (const auto &[u, adj]: next_) in_degree[u] = 0;
            for (const auto &[u, adj]: next_) {
                for (const auto &v: adj) in_degree[v]++;
            }
            for (const auto &[u, _]: in_degree) {
                if (in_degree[u] == 0) Q.push(u);
            }

            while (!Q.empty()) {
                auto u = Q.front();
                visited.insert(u);
                result.push_back(u);
                Q.pop();
                for (const auto &v: next_[u]) {
                    if (visited.count(v)) continue;
                    if (--in_degree[v] == 0) Q.push(v);
                }
            }
            return result;
        }

    };

    template<typename Vertex, typename T>
    class undirectedTemporalGraph {
        using Edge = triplet<Vertex, T>;
    private:
        std::unordered_map<Vertex, std::set<Vertex>> neighbors{};
        std::queue<Edge> edges_queue{}; // edges awaiting to be inserted
        std::queue<Edge> edges_cur{};   // cur edges in the sliding window
        T L;  // sliding window size
        T cur_time;
    public:
        explicit undirectedTemporalGraph(std::vector<Edge> edges, T L_) {
            L = L_;
            cur_time = L_;
            if (L == 0 || edges.empty()) throw std::runtime_error("L should be greater than 0");
            for (const auto &e: edges) edges_queue.push(e);

            // first sliding window
            while (!edges_queue.empty() && static_cast<T>(get<2>(edges_queue.front())) <= cur_time) {
                insert_edge();
                showUp();
            }

            // maintain the window
            while (!edges_queue.empty() && static_cast<T>(get<2>(edges_queue.front())) <= ++cur_time){
                insert_edge();
                delete_edge();
                showUp();
            }

        }

        void showUp() {
            // visualize the graph
            std::cout << "cur_time: " << cur_time << std::endl;
            for (const auto &[u, adj]: neighbors) {
                std::cout << u << ": ";
                for (const auto &v: adj) std::cout << v << " ";
                std::cout << std::endl;
            }
        }

        // operations:

        // 1) scanning all neighbors of the vertex in the sliding window
        auto get_neighbors() {
            return neighbors;
        }

        // 2) inserting a new neighbor when a new edge arrives
        void insert_edge() {
            while (!edges_queue.empty() && static_cast<T>(get<2>(edges_queue.front())) <= cur_time) insert_one_edge();
        }

        void insert_one_edge(){
            auto [u, v, t] = edges_queue.front();
            edges_cur.push(edges_queue.front());

            edges_queue.pop();
            neighbors[u].insert(v);
            neighbors[v].insert(u);
        }

        // 3) deleting the expired neighbor
        void delete_edge() {
            while (!edges_cur.empty() && static_cast<T>(get<2>(edges_cur.front())) <= cur_time - L) {
                auto [u, v, t] = edges_cur.front();
                edges_cur.pop();
                neighbors[u].erase(v);
                neighbors[v].erase(u);
            }
        }

    };
}

int main() {
    // undirected graph
    std::vector<std::pair<char, char>> Fig1 = { {'M','P'},{'P','N'},{'M','N'},

                                                {'I','A'},{'I','C'}, {'I','G'},
                                                {'A','C'},{'A','B'}, {'A','E'},
                                                {'C','B'}, {'C','F'}, {'B','G'},
                                                {'G','F'}, {'E','F'}, {'F','B'},

                                                {'D','H'}, {'D','K'}, {'H','J'},
                                                {'H','K'}, {'J','K'}
    };

    // directed graph with edges from u1 to u2
    std::vector<std::pair<char, char>> Fig2 = { {'A', 'D'}, {'A', 'F'},
                                                {'B', 'A'}, {'B', 'C'},
                                                {'C', 'D'},
                                                {'D', 'E'},
                                                {'F', 'C'}, {'F', 'D'},
                                                {'G', 'E'}, {'G', 'F'}, {'G', 'I'},
                                                {'H', 'A'}, {'H', 'G'}, {'H', 'E'}, {'H', 'I'},{'H', 'J'},
                                                {'I', 'J'},
                                                {'J', 'F'},
    };

    // Undirected Temporal Graph
    // triplet (u1, u2, t) means that there is an edge between u1 and u2 at time t
    // edges are naturally sorted by time
    std::vector<std::tuple<char,char,std::size_t>> Fig3 = { {'A', 'B', 1u}, {'D', 'C', 2u}, {'A', 'C', 3u}, {'B', 'C', 4u},};

    using namespace ASMT1;
    // Q1
    undirectedUnweightedGraph<char> g1(Fig1);

    // Q2
    disjoint_set<char> ds1(Fig1);

    // Q3
    directedUnweightedGraph<char> g2(Fig2);
    auto res = g2.topological_sort();
    std::cout<< "Topological sort: ";
    for (auto& u: res) {
        std::cout<< u << " ";
    }
    std::cout<<std::endl;

    // Q4 line 105
    /** std::size_t directedUnweightedGraph.setDistance(std::unordered_set<Vertex> x, std::unordered_set<Vertex> y) **/

    // Q5
    undirectedTemporalGraph<char, std::size_t> g3(Fig3, 3u);



    return 0;
}
