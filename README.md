# COMP9312
code for asmt1

According to the lecture, performing an id-to-name mapping before graph processing is necessary since the time comlexity can be reduced by doing so given the overhead and the real performance of hashmap is unsatisfactory somewhat. 
Furthermore, after mapping any type of vertex id into integer, the vertex is now naturally std::size_t type, which can be directly used as the subscript of array.

Instead of doing bijection by 2 std::unordered_map between vertex id and vertex name, hereby, std::unordered_map <Vertex, std::set> is used for coding simplification.
Likewise, in Disjoint set and other else class, the size array and parents array (std::vector or std::array)are replaced by std::unordered_map. Since there will only be one more std::unordered_map used and amortized time complexity can only have a negligible effect.

The purpose of this is just to simplify the tediousness of coding and explanation above illustrated that the expected time complexity is unaltered.
