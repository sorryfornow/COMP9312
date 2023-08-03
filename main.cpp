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

void asmt1(){
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

}

namespace ASMT2 {

    template<typename Vertex>
    class disjoint_set {
    private:
        std::unordered_map<Vertex, Vertex> parent{};
        std::unordered_map<Vertex, std::size_t> rank{}; //  union by rank
    public:
        explicit disjoint_set(std::vector<Vertex> vertices) {
            for (const auto &v: vertices) {
                parent[v] = v;
                rank[v] = 0u;
            }
        }

        explicit disjoint_set() = default;

        void make_set(Vertex v) {
            if (parent.count(v)) return;
            parent[v] = v;
            rank[v] = 0u;
        }

        Vertex find(Vertex u) {
            if (parent[u] != u) parent[u] = find(parent[u]);
            return parent[u];
        }

        void merge(Vertex u, Vertex v) {
            auto u_root = find(u);
            auto v_root = find(v);
            if (u_root == v_root) return;
            if (rank[u_root] < rank[v_root]) {
                parent[u_root] = v_root;
            } else if (rank[u_root] > rank[v_root]) {
                parent[v_root] = u_root;
            } else {
                parent[v_root] = u_root;
                rank[u_root]++;
            }
        }

        std::unordered_map<Vertex, std::vector<Vertex>> get_components() {
            std::unordered_map<Vertex, std::vector<Vertex>> components{};
            for (const auto &[u, _]: parent) {
                components[find(u)].emplace_back(u);
            }
            return components;
        }

        std::unordered_map<Vertex, Vertex> get_parent() {
            return parent;
        }

    };

    template<typename Vertex>
    class undirectedUnweightedGraph {
    private:
        // type definitions
        using edge_t = std::pair<Vertex, Vertex>;   // edge (u, v) type
        using truss_tri_t = std::tuple<Vertex, Vertex, std::size_t>;    // edge with its truss number (u, v, k) type
        using adjacency_t = std::unordered_map<Vertex, std::vector<Vertex>>;    // adjacency list type
        // using path_t = std::vector<Vertex>
        // using Paths = std::vector<path>

        // member variables
        adjacency_t neighbors{};
        std::unordered_map<Vertex, std::vector<std::vector<Vertex>>> shortest_paths{};

        // member functions
        void dfs_(Vertex start,
                 Vertex current,
                 std::unordered_map<Vertex, int>& distance,
                 std::unordered_map<Vertex, std::vector<std::vector<Vertex>>>& allPaths,
                 std::vector<Vertex>& currentPath) {
            currentPath.push_back(current);

            if (current == start) {
                std::reverse(currentPath.begin(), currentPath.end());
                allPaths[start].push_back(currentPath);
                std::reverse(currentPath.begin(), currentPath.end());
            } else {
                for (Vertex prev : neighbors.at(current)) {
                    if (distance[current] == distance[prev] + 1) {
                        dfs_(start, prev, distance, allPaths, currentPath);
                    }
                }
            }

            currentPath.pop_back();
        }

        std::unordered_map<Vertex, std::vector<std::vector<Vertex>>> shortestPathsFrom( Vertex start) {
            std::unordered_map<Vertex, std::vector<std::vector<Vertex>>> allPaths;
            std::unordered_map<Vertex, int> distance;
            std::queue<Vertex> q;

            for (const auto& entry : neighbors) {
                Vertex vertex = entry.first;
                distance[vertex] = -1; // Initialize distance as -1 (unreachable)
            }

            distance[start] = 0;
            q.push(start);

            while (!q.empty()) {
                Vertex current = q.front();
                q.pop();

                for (Vertex neighbor : neighbors.at(current)) {
                    if (distance[neighbor] == -1) {
                        distance[neighbor] = distance[current] + 1;
                        q.push(neighbor);
                    }
                }
            }

            for (const auto& entry : distance) {
                Vertex vertex = entry.first;
                if (vertex != start) {
                    std::vector<Vertex> currentPath;
                    dfs_(start, vertex, distance, allPaths, currentPath);
                }
            }

            return allPaths;
        }

        std::unordered_map<Vertex, std::vector<std::vector<Vertex>>> shortestPaths_all(){
            // yield all shortest paths in the graph
            std::unordered_map<Vertex, std::vector<std::vector<Vertex>>> all_paths{};
            // initialize the paths dictionary
            for (const auto &[u, adj]: neighbors) all_paths[u] = {};
            // BFS probe all shortest paths from each vertex
            for (const auto &[u, adj]: neighbors) {
                std::unordered_map<Vertex, std::vector<std::vector<Vertex>>> paths = shortestPathsFrom(u);
                for (const auto &[_, path]: paths)
                    for (const auto & p: path)
                        all_paths[u].emplace_back(p);
            }

            for (auto &[_, paths]: all_paths) {
                // sort paths regards the last element in the path
                std::sort(paths.begin(), paths.end(), [](const auto &lhs, const auto &rhs) {
                    return lhs.back() < rhs.back();
                });
            }

            return all_paths;
        }

    public:
        undirectedUnweightedGraph() = delete;
        // constructor with initialization
        explicit undirectedUnweightedGraph(std::vector<edge_t> edges) {
            std::unordered_map<Vertex, std::set<Vertex>> neighbors_{};
            for (const auto &[u, v]: edges) {
                neighbors_[u].insert(v);
                neighbors_[v].insert(u);
            }
            for (const auto &[u, adj]: neighbors_) {
                neighbors[u] = std::vector<Vertex>(adj.begin(), adj.end());
            }
            shortest_paths = {};
        }

        // return the adjacency list
        std::unordered_map<Vertex, std::vector<Vertex>> get_neighbors() {
            return neighbors;
        }

        // return the shortest paths dictionary
        std::unordered_map<Vertex, std::vector<std::vector<Vertex>>> get_shortest_paths() {
            if (shortest_paths.empty()) {
                shortest_paths = shortestPaths_all();
            }
            return shortest_paths;
        }

        // solution for question 2
        std::vector<std::vector<Vertex>>
        get_truss_set([[maybe_unused]] adjacency_t G, std::vector<truss_tri_t> KTruss_arr, int k_ ){
            auto k = static_cast<std::size_t>(k_);
            // initializing an empty disjoint set
            disjoint_set<Vertex> DS;
            // merge edges with truss number >= k where the edges are sorted regards truss number in descending order
            for (const auto &[u, v, n]: KTruss_arr) {
                if (n<k) break;
                DS.make_set(u);
                DS.make_set(v);
                DS.merge(u, v);
            }
            // get the disjoint set
            std::unordered_map<Vertex, Vertex> parent = DS.get_parent();

            // yield the components
            std::unordered_map<Vertex, std::vector<Vertex>> components{};
            for (const auto &[u, _]: parent) components[DS.find(u)].emplace_back(u);

            // transforms the components dictionary into the truss_set vector.
            std::vector<std::vector<Vertex>> truss_set{};
            for (const auto &[u, comp]: components) truss_set.emplace_back(comp);

            // sort the vertices in each component if needed
            // for (auto &comp: truss_set) std::sort(comp.begin(), comp.end());

            return truss_set;
        }

        // solution for question 4
        double betweenness(Vertex v) {

            double betweenness = 0.0;
            auto all_paths_dict = get_shortest_paths();
            std::unordered_set<Vertex> S{}; // store the visited vertices
            S.insert(v);

            for(const auto &[s, _]: neighbors) {
                if (S.count(s)) continue; // if the vertex is visited, skip
                S.insert(s);    // marked as visited

                auto paths_list = all_paths_dict.at(s);
                std::unordered_map<Vertex, std::vector<std::vector<Vertex>>> paths_dict{};
                for (const auto &path: paths_list) {
                    if (S.count(path.back())) continue; // if the last vertex in the path is visited, skip
                    auto& PD = paths_dict[path.back()];
                    PD.emplace_back(path);
                }

                for (const auto &[t, paths]: paths_dict) {
                    double numerator = 0.0;
                    double denominator = paths.size();
                    for (const auto& path: paths){
                        if (std::find(path.begin(), path.end(), v) != path.end()) numerator += 1.0;
                    }
                    betweenness += numerator / denominator;
                }
            }

            return betweenness;
        }

        double closeness(Vertex v) {

            double closeness_reciprocal = 0.0;
            auto all_paths_dict = get_shortest_paths();
            std::unordered_set<Vertex> S{}; // store the visited vertices
            S.insert(v);

            auto paths_list = all_paths_dict.at(v);
            for (const auto &path: paths_list) {
                if (S.count(path.back())) continue; // if the last vertex in the path is visited, skip
//                printf("t: %c, length: %lu\n", path.back(), path.size()-1);
                closeness_reciprocal += path.size()-1;
                S.insert(path.back());
            }

//            printf("closeness_reciprocal : %f\n", closeness_reciprocal);
            return 1.0 / closeness_reciprocal;
        }

    };

}


void asmt2(){

    using namespace ASMT2;

    // graph of Figure 4
    std::vector<std::pair<std::size_t, std::size_t>> Fig4{
            {1u, 2u}, {1u, 3u}, {1u, 4u}, {1u, 5u},
            {2u, 3u}, {2u, 4u}, {2u, 5u},
            {3u, 4u}, {4u, 5u},

            {4u, 11u}, {5u, 11u}, {6u, 7u}, {6u, 8u},
            {6u, 10u}, {6u, 12u}, {7u, 12u},
            {8u, 9u}, {8u, 10u}, {9u, 10u},

            {2u, 6u}, {3u, 12u}, {7u, 13u},
    };

    // precomputed truss number for each edge
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> truss_tri{
            {1u, 2u, 4u}, {1u, 3u, 4u}, {1u, 4u, 4u}, {1u, 5u, 4u},
            {2u, 3u, 4u}, {2u, 4u, 4u}, {2u, 5u, 4u},
            {3u, 4u, 4u}, {4u, 5u, 4u},

            {4u, 11u, 3u}, {5u, 11u, 3u}, {6u, 7u, 3u}, {6u, 8u, 3u},
            {6u, 10u, 3u}, {6u, 12u, 3u}, {7u, 12u, 3u},
            {8u, 9u, 3u}, {8u, 10u, 3u}, {9u, 10u, 3u},

            {2u, 6u, 2u}, {3u, 12u, 2u}, {7u, 13u, 2u},
    };

    // Q2
    undirectedUnweightedGraph<std::size_t> g1(Fig4);
    auto res = g1.get_truss_set(g1.get_neighbors(), truss_tri, 3);

    // print the result
    for (auto& comp: res) {
        std::cout<< "[";
        for (auto& u: comp) {
            std::cout<< u << ", ";
        }
        std::cout<<"]"<<std::endl;
    }

    // Q4
    std::vector<std::pair<char, char>> Fig5{
            {'A', 'B'}, {'A', 'C'}, {'A', 'D'}, {'A', 'H'},
            {'C', 'H'}, {'C', 'D'}, {'C', 'F'}, {'C', 'E'},
            {'D', 'F'}, {'D', 'E'}, {'E', 'F'}, {'E', 'G'},
            {'F', 'G'},
    };

    undirectedUnweightedGraph<char> g2(Fig5);
    // Compute the betweenness centrality and closeness centrality of nodes C and D.

    auto betweenness_centrality_C = g2.betweenness('C');
    printf("betweenness centrality of node C: %f\n", betweenness_centrality_C);
    auto betweenness_centrality_D = g2.betweenness('D');
    printf("betweenness centrality of node D: %f\n", betweenness_centrality_D);
    auto closeness_centrality_C = g2.closeness('C');
    printf("closeness centrality of node C: %f\n", closeness_centrality_C);
    auto closeness_centrality_D = g2.closeness('D');
    printf("closeness centrality of node D: %f\n", closeness_centrality_D);

    auto allpath = g2.get_shortest_paths();

    // E's graphlet type 1 paths
    auto e_paths = allpath.at('E');
    std::vector<std::vector<char>> res_e = {};
    for (const auto& path: e_paths) {
        if (path.back() == 'B') res_e.push_back(path);
    }
    std::cout<< "all paths from E to B: ";
    for (const auto& path: res_e) {
        std::cout<< "[";
        for (const auto& v: path) {
            std::cout<< v << ", ";
        }
        std::cout<<"]"<<std::endl;
    }

    // A's graphlet type 1 paths
    auto a_paths = allpath.at('A');
    std::vector<std::vector<char>> res_a = {};
    for (const auto& path: a_paths) {
        if (path.back() == 'G') res_a.push_back(path);
    }
    std::cout<< "all paths from A to G: ";
    for (const auto& path: res_a) {
        std::cout<< "[";
        for (const auto& v: path) {
            std::cout<< v << ", ";
        }
        std::cout<<"]"<<std::endl;
    }


//    //   print all paths
//    auto allpath = g2.get_shortest_paths();
//    for (const auto &[u, paths]: allpath) {
//        std::cout<< u << ": ";
//        for (const auto &path: paths) {
//            std::cout<< "[";
//            for (const auto &v: path) {
//                std::cout<< v << ", ";
//            }
//            std::cout<<"]";
//        }
//        std::cout<<std::endl;
//    }


}

int main() {

    asmt2();

    return 0;
}
