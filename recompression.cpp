/**
 * Copyright 2018 Christopher Osthues
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithm>
#include <limits>
#include <fstream>
#include <climits>
#include <iostream>
#include <list>
#include <string>
#include <tuple>
#include <vector>
#include <utility>
#include <sstream>
#include <string>

#include <chrono>
#include <map>

#include "rlslp.hpp"
//#include <sdsl/int_vector.hpp>
//#include <sdsl/rank_support_v5.hpp>

//#include "radix_sort.hpp"
//#include "../rlslp.hpp"

namespace lce {

typedef std::int32_t variable_t;
typedef std::uint32_t terminal_count_t;
typedef std::vector<variable_t> text_t;
//typedef std::list<variable_t>::iterator iterator_t;
//typedef std::tuple<variable_t, variable_t, iterator_t> triple_t;
typedef std::tuple<variable_t, variable_t, bool> multiset_t;

//template<typename T, typename V>
long radix_sort(std::vector<multiset_t> &vec) {
    const auto startTime = std::chrono::system_clock::now();
    size_t mask = 1;
    for (std::uint32_t j = 0; j < sizeof(variable_t) * CHAR_BIT; ++j) {
        std::vector<multiset_t> zeroes;
        std::vector<multiset_t> ones;
        for (const auto &v : vec) {
            const auto &sec = std::get<1>(v);
            size_t bucket = mask & sec;
            if (bucket) {
                ones.emplace_back(v);
            } else {
                zeroes.emplace_back(v);
            }
        }
        std::uint32_t k = 0;
        for (; k < zeroes.size(); ++k) {
            vec[k] = zeroes[k];
        }
        for (std::uint32_t l = 0; l < ones.size(); ++l, ++k) {
            vec[k] = ones[l];
        }
        mask <<= 1;
    }
    mask = 1;
    for (std::uint32_t j = 0; j < sizeof(variable_t) * CHAR_BIT; ++j) {
        std::vector<multiset_t> zeroes;
        std::vector<multiset_t> ones;
        for (const auto &v : vec) {
            const auto &sec = std::get<0>(v);
            size_t bucket = mask & sec;
            if (bucket) {
                ones.emplace_back(v);
            } else {
                zeroes.emplace_back(v);
            }
        }
        std::uint32_t k = 0;
        for (; k < zeroes.size(); ++k) {
            vec[k] = zeroes[k];
        }
        for (std::uint32_t l = 0; l < ones.size(); ++l, ++k) {
            vec[k] = ones[l];
        }
        mask <<= 1;
    }
    const auto endTime = std::chrono::system_clock::now();
        const auto timeSpan = endTime - startTime;
        std::cout << "Time for radix sort triple: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
	return std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count();
}

/**
 * @brief Implementation of radix sort for tuples according to the first two parameters.
 *
 * This implementation sorts the numbers in the vector at bit level.
 *
 * @param vec[in,out] The vector to sort
 */
/*template<typename t_size>
long radix_sort(std::vector<t_size> &vec) {
const auto startTime = std::chrono::system_clock::now();
    size_t mask = 1;
    for (std::uint32_t j = 0; j < sizeof(t_size) * CHAR_BIT; ++j) {
        std::vector<t_size> zeroes;
        std::vector<t_size> ones;
        for (const auto &v : vec) {
            const auto &sec = std::get<0>(v);
            size_t bucket = mask & sec;
            if (bucket) {
                ones.emplace_back(v);
            } else {
                zeroes.emplace_back(v);
            }
        }
        std::uint32_t k = 0;
        for (; k < zeroes.size(); ++k) {
            vec[k] = zeroes[k];
        }
        for (std::uint32_t l = 0; l < ones.size(); ++l, ++k) {
            vec[k] = ones[l];
        }
        mask <<= 1;
    }
const auto endTime = std::chrono::system_clock::now();
        const auto timeSpan = endTime - startTime;
        std::cout << "Time for tuple: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
	return std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count();
}*/

//struct pair_hash {
//    template<typename P1, typename P2>
//    size_t operator()(const std::pair<P1, P2>& pair) const {
//        std::hash<P1> h1;//(pair.first);
//        std::hash<P2> h2;//(pair.second);
//        return h1(pair.first) ^ (h2(pair.second) << 1);
//    }
//};
//
//struct tuple_hash {
//    template<typename T1, typename T2, typename T3>
//    size_t operator()(const std::tuple<T1, T2, T3>& tup) const {
//        std::hash<T1> h1;//(pair.first);
//        std::hash<T2> h2;//(pair.second);
//        std::hash<T3> h3;
//        return h1(std::get<0>(tup)) ^ (h2(std::get<1>(tup)) << 1) ^ (h3(std::get<2>(tup)) << 2);
//    }
//};

    

struct lceq {
    variable_t nt;
    size_t pos;
    size_t traverse;

    lceq(variable_t nt, size_t pos, size_t traverse) {
        this->nt = nt;
        this->pos = pos;
        this->traverse = traverse;
    }
};

/**
 * @brief Replaces all letters in the text with new non-terminals.
 *
 * @param text[in] The text to compress
 * @param t[out] The new text where all letters are replaced with non-terminals
 * @param alphabet_size[out] The size of the alphabet used in the text
 * @param mapping[out] The mapping of the symbols in the text to the non-terminal
 */
void replace_letters(text_t &t, rlslp& rlslp, variable_t& alphabet_size, std::vector<variable_t> &mapping) {
    const auto startTime = std::chrono::system_clock::now();
    std::map<variable_t, variable_t> alpha;
    for (size_t i = 0; i < t.size(); ++i) {
        alpha[t[i]] = 0;
    }

    size_t i = 0;
    for (auto& a : alpha) {
        mapping.emplace_back(i);
        a.second = i++;
        rlslp.non_terminals.emplace_back(a.first, 1);
    }

    /*rlslp.blocks.resize(alpha.size());
    for (size_t j = 0; j < rlslp.blocks.size(); ++j) {
        rlslp.blocks[j] = false;
    }*/
    alphabet_size = alpha.size();
    rlslp.terminals = alphabet_size;

    for (size_t i = 0; i < t.size(); ++i) {
        t[i] = alpha[t[i]];
    }
    const auto endTime = std::chrono::system_clock::now();
    const auto timeSpan = endTime - startTime;
    std::cout << "Time for replace letters: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
}

/**
 * @brief Detects all blocks in the text and replaces them with new non-terminals.
 *
 * @param t[in,out] The text
 * @param alphabet_size[in,out] The size of the alphabet used in the text
 * @param mapping[in,out] The mapping of the symbols in the text to the non-terminal
 */
void bcomp(text_t& t, rlslp& rlslp, size_t& text_size, variable_t& alphabet_size,std::vector<variable_t>& mapping) {
    std::cout << "Text size (Input BComp): " << text_size << std::endl;
    const auto startTime = std::chrono::system_clock::now();
    variable_t next_nt = rlslp.non_terminals.size();
    std::vector<bool> block_vec(alphabet_size);
    variable_t block_len = 1;
    size_t block_count = 0;
    size_t new_text_size = text_size;
    std::vector<std::map<variable_t, variable_t>> blocks(alphabet_size);
    std::vector<std::pair<variable_t, size_t>> positions;
    size_t copy_i = 0;
    bool copy = false;
    const auto startTimeBlock = std::chrono::system_clock::now();
    for (size_t i = 1; i < text_size; ++i, ++copy_i) {
        while (i < text_size && t[i-1] == t[i]) {
            block_len++;
            i++;
            new_text_size--;
        }
        if (block_len > 1) {
            if (!copy) {
                t[copy_i] = t[i-1];
                copy_i++;
                copy = true;
            }
            positions.emplace_back(block_len, copy_i-1);
            block_vec[t[i - 1]] = true;
            blocks[t[i-1]][block_len] = 1;
            block_count++;
            block_len = 1;
        }
        if (copy) {
            t[copy_i] = t[i];
        }
    }
    if (copy_i < new_text_size) {
        t[copy_i] = t[text_size-1];
    }
    const auto endTimeBlock = std::chrono::system_clock::now();
    const auto timeSpanBlock = endTimeBlock - startTimeBlock;
    std::cout << "Time for finding blocks: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpanBlock).count() << "[ms]" << std::endl;

    block_count = 0;

    const auto startTimeAss = std::chrono::system_clock::now();
    for (size_t i = 0; i < block_vec.size(); ++i) {
        if (block_vec[i]) {
            for (auto& block : blocks[i]) {
                block.second = alphabet_size++;
                block_count++;
                rlslp.non_terminals.emplace_back(mapping[i], -block.first, rlslp.non_terminals[mapping[i]].len * block.first);
                mapping.emplace_back(next_nt++);
            }
        }
    }
    const auto endTimeAss = std::chrono::system_clock::now();
    const auto timeSpanAss = endTimeAss - startTimeAss;
    std::cout << "Time for block nts: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpanAss).count() << "[ms]" << std::endl;


    const auto startTimeRep = std::chrono::system_clock::now();
    text_size = new_text_size;
    for (const auto& pos : positions) {
        t[pos.second] = blocks[t[pos.second]][pos.first];
    }
    const auto endTimeRep = std::chrono::system_clock::now();
    const auto timeSpanRep = endTimeRep - startTimeRep;
    std::cout << "Time for replacing blocks: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpanRep).count() << "[ms]" << std::endl;
    t.resize(new_text_size);
    t.shrink_to_fit();

    /*if (block_count > 0) {
        rlslp.blocks.resize(rlslp.blocks.size() + block_count, true);
    }*/

    const auto endTime = std::chrono::system_clock::now();
    const auto timeSpan = endTime - startTime;
    std::cout << "Time for bcomp: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
    std::cout << "Text size (Ouput BComp): " << text_size << std::endl;
}


/**
 * @brief Computes the multiset representing the adjacency graph of the symbols in the text.
 *
 * @param t[in] The text
 * @return The multiset
 */
inline void compute_multiset(const text_t& t, size_t& text_size, std::map<std::pair<variable_t, variable_t>, std::pair<size_t, size_t>>& multiset) {
    const auto startTime = std::chrono::system_clock::now();
    // Compute adjacency graph of the symbols in the current text
    for (size_t i = 1; i < text_size; ++i) {
        std::pair<variable_t, variable_t> adj;
        if (t[i-1] > t[i]) {
            adj.first = t[i-1];
            adj.second = t[i];
        } else {
            adj.first = t[i];
            adj.second = t[i-1];
        }
        auto found = multiset.find(adj);
        if (found == multiset.end()) {
            if (t[i-1] > t[i]) {
                multiset[adj].first = 1;
                multiset[adj].second = 0;
            } else {
                multiset[adj].first = 0;
                multiset[adj].second = 1;
            }
        } else {
            if (t[i-1] > t[i]) {
                (*found).second.first += 1;
            } else {
                (*found).second.second += 1;
            }
        }
    }

    const auto endTime = std::chrono::system_clock::now();
    const auto timeSpan = endTime - startTime;
    std::cout << "Time for multiset: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
}

/**
 * @brief Computes the partition based on the adjacency graph of the text given by the multiset.
 *
 * @param multiset[i] Multiset representing the adjacency graph of the text
 * @param alphabet_size[in] The size of the alphabet used in the text
 *
 * @return The partition of the symbols in the alphabet to maximize the number of pairs to be compressed
 */
inline void compute_partition(const std::map<std::pair<variable_t, variable_t>, std::pair<size_t, size_t>>& multiset, std::vector<bool>& partition, const variable_t& alphabet_size) {
    const auto startTime = std::chrono::system_clock::now();
    //int lr_count = 0;
    //int rl_count = 0;
    variable_t c = 0;
    int left = 0, right = 0;
    for (const auto& adj : multiset) {
        auto first = adj.first;
        while (c < alphabet_size && first.first > c) {
            partition[c] = left > right;
            c++;
            left = 0;
            right = 0;
        }
        if (partition[first.second]) {
            right += adj.second.first + adj.second.second;
        } else {
            left += adj.second.first + adj.second.second;
        }
    }

    while (c < alphabet_size) {
        partition[c] = left > right;
        c++;
        left = 0;
        right = 0;
    }
    
    /*if (rl_count > lr_count) {
        partition.flip();
    }*/

    const auto endTime = std::chrono::system_clock::now();
    const auto timeSpan = endTime - startTime;
    std::cout << "Time for computing partition: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
}

/**
 * @brief Counts the number of pairs.
 *
 * Counts the number of pairs build by symbols from the left set with symbols from the right set of the
 * partition and vice versa. If the number of pairs build by symbols from the right set with symbols of
 * the left set is greater than the number of symbol pairs from left to right then the sets of the
 * partition are swapped.
 *
 * @param multiset[in] Multiset representing the adjacency graph of the text
 * @param partition[in,out] The partition of the letters in the current alphabet represented by a bitvector
 */
inline void count_pairs(const std::map<std::pair<variable_t, variable_t>, std::pair<size_t, size_t>>& multiset, std::vector<bool>& partition) {
    const auto startTime = std::chrono::system_clock::now();
    int lr_count = 0;
    int rl_count = 0;
    for (const auto &adj : multiset) {
        if (!partition[adj.first.first] && partition[adj.first.second]) {  // (c,b,0) -> cb in text
            lr_count += adj.second.first;
            rl_count += adj.second.second;
        } else {
            rl_count += adj.second.first;
            lr_count += adj.second.second;
        }
    }

    // If there are more pairs in the current text from right set to left set swap partition sets
    if (rl_count > lr_count) {
        partition.flip();
    }
    const auto endTime = std::chrono::system_clock::now();
    const auto timeSpan = endTime - startTime;
    std::cout << "Time for count pairs: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
}

/**
 * @brief Determines a partition of the letters in the current alphabet.
 *
 * This is an implementation of a greedy 1/4 approximation algorithm.
 *
 * At first the adjacency graph of the text will be computed by creating the multiset
 * of the tuples (t[i], t[i+1], 0) for t[i] > t[i+1] and (t[i+1], t[i], 1) for
 * t[i+1] > t[i]. This multiset will be sorted with radix sort. Then for all increasingly
 * sorted letters c in the current alphabet the sums l and r of freq(c, b, cdot)
 * for all b in the left set of the partition and freq(c, b, cdot) for all b
 * in the right set of the partition will be computed. freq{c, bar{c}, cdot} counts all local
 * pairs of cb and bc in the text. If l is greater then r the letter will be added
 * to the right set of the partition else to the left set. At last the total pairs of the
 * partition from the left set to the right set and vice versa occurring in the text will be
 * counted. If there are more pairs from the right set to the left set occurring in the text
 * the sets will be swapped.
 *
 * @param t[in] The text
 * @param alphabet_size[in] The size of the alphabet used in the text
 *
 * @return The partition of the letters in the current alphabet represented by a bitvector
 */
inline void partition(const text_t& t, size_t& text_size, const variable_t& alphabet_size, std::vector<bool>& partition) {
    const auto startTime = std::chrono::system_clock::now();
    std::map<std::pair<variable_t, variable_t>, std::pair<size_t, size_t>> multiset;
    compute_multiset(t, text_size, multiset);

    compute_partition(multiset, partition, alphabet_size);

    // Count pairs in the current text based on the pairs build by the partition
    // from left set to right set and vice versa
    count_pairs(multiset, partition);

    const auto endTime = std::chrono::system_clock::now();
    const auto timeSpan = endTime - startTime;
    std::cout << "Time for partition: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
}

/**
 * @brief Detects all pairs in the text and replaces them with new non-terminals.
 *
 * @param t[in,out] The text
 * @param alphabet_size[in,out] The size of the alphabet used in the text
 * @param mapping[in,out] The mapping of the symbols in the text to the non-terminal
 */
void pcomp(text_t& t, rlslp& rlslp, size_t& text_size, variable_t& alphabet_size, std::vector<variable_t>& mapping) {
    std::cout << "Text size (Input PComp): " << text_size << std::endl;
    const auto startTime = std::chrono::system_clock::now();
    std::vector<bool> part(alphabet_size, false);
    partition(t, text_size, alphabet_size, part);


    variable_t next_nt = rlslp.non_terminals.size();
    size_t pair_count = 0;
    size_t new_text_size = text_size;
    std::vector<std::map<variable_t, std::vector<variable_t>>> pairs(alphabet_size);
    size_t copy_i = 0;
    bool copy = false;
    std::cout << "Alphabet: " << alphabet_size << std::endl;
    if (text_size == 2) {
        for (const auto& b : part) {
            std::cout << b;
        }
        std::cout << std::endl;
        for (size_t i = 0; i < text_size; ++i) {
            std::cout << t[i] << " ";
        }
        std::cout << std::endl;
    }
    for (size_t i = 1; i < text_size; ++i, ++copy_i) {
        if (!part[t[i-1]] && part[t[i]]) {
            copy = true;
            auto found = pairs[t[i-1]].find(t[i]);
            if (found == pairs[t[i-1]].end()) {
                mapping.emplace_back(next_nt++);
                pair_count++;
            }
            pairs[t[i-1]][t[i]].push_back(copy_i);
            t[copy_i] = t[i-1];
            i++;
            new_text_size--;
        } else if (copy) {
            t[copy_i] = t[i-1];
        }
    }
    if (copy_i < new_text_size) {
        t[copy_i] = t[text_size-1];
    }
    std::cout << "All pairs computed" << std::endl;

    text_size = new_text_size;

    size_t alpha_size = alphabet_size;
    for (size_t i = 0; i < alpha_size; ++i) {
        if (!part[i]) {
            for (const auto& pair : pairs[i]) {
                for (const auto& pos : pair.second) {
                    t[pos] = alphabet_size;
                }
                rlslp.non_terminals.emplace_back(mapping[i], mapping[pair.first], rlslp.non_terminals[mapping[i]].len + rlslp.non_terminals[mapping[pair.first]].len);
                alphabet_size++;
            }
        }
    }
    t.resize(new_text_size);
    t.shrink_to_fit();

    std::cout << "Replaced all pairs" << std::endl;

    /*if (pair_count > 0) {
        rlslp.blocks.resize(rlslp.blocks.size() + pair_count, false);
    }*/

    // Introduce new non-terminals for all found pairs
    const auto endTime = std::chrono::system_clock::now();
    const auto timeSpan = endTime - startTime;
    std::cout << "Time for pcomp: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
    
    std::cout << "Text size (Ouput PComp): " << text_size << std::endl;
}


/**
 * @brief Updates the alphabet used in the current text.
 *
 * @param t[in,out] The current text
 * @param alphabet_size[in,out] The size of the alphabet used in the text
 * @param mapping[in,out] The mapping of the symbols in the text to the non-terminal
 */
void compute_alphabet(text_t& t, size_t& text_size, variable_t &alphabet_size, std::vector<variable_t> &mapping) {
    const auto startTime = std::chrono::system_clock::now();
    // Create bitvector for all letters to determine which are used in the current text
    std::vector<bool> used(alphabet_size, 0);
    std::vector<size_t> ranks(alphabet_size, 0);
    for (size_t i = 0; i < text_size; ++i) {
        used[t[i]] = true;
    }
    // Build rank structure over the bitvector
    ranks[0] = 0;//(used[0]? 1 : 0);
    for (size_t i = 1; i < ranks.size(); ++i) {
        if (used[i-1]) {
            ranks[i] = ranks[i-1] + 1;
        } else {
            ranks[i] = ranks[i-1];
        }
    }

    alphabet_size = 0;
    auto m = 0;
    auto i = 0;
    for (const auto& u : used) {
        // If letter is used in the current text add it to the new alphabet and set the rank of the letter
        if (u) {
            alphabet_size++;
            mapping[i++] = mapping[m];
        }
        ++m;
    }
    mapping.resize(i);
    mapping.shrink_to_fit();

    // overwrite text with ranks of symbols
    for (size_t i = 0; i < text_size; ++i) {
        t[i] -= (t[i] - ranks[t[i]]);
    }
    const auto endTime = std::chrono::system_clock::now();
    const auto timeSpan = endTime - startTime;
    std::cout << "Time for alphabet: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
}

/**
 * @brief Determines the length of the remaining substring to be extracted.
 *
 * @param i The start index within the current subtree
 * @param len The length to extract
 * @param derived_string_len The length of the substring derived by the current non-terminal
 *
 * @return The length of the remaining substring to be extracted
 */
inline size_t rest_len(const size_t i, const size_t len, const size_t derived_string_len) {
    if (i + len > derived_string_len) {
        return len - derived_string_len + i;
    }
    return 0;
}

/**
 * @brief Extracts a substring of length @code{len} beginning at position @code{i}.
 *
 * @param i The start position in the subtree
 * @param len The length to be extracted
 * @param nt The current node
 * @param s[out] The extracted (sub)string
 */
void extract(const size_t i, const size_t len, const variable_t nt, text_t& s, rlslp& rlslp) {
    if (len > 0) {
        if (nt < rlslp.terminals && i == 0) {
            s.push_back(rlslp.non_terminals[nt].production[0]);
        } else {
            auto first = rlslp.non_terminals[nt].production[0];
            auto second = rlslp.non_terminals[nt].production[1];
            auto first_len = rlslp.non_terminals[first].len;

            if (second < 0/*rlslp.blocks[nt]*/) {
                auto block_len = first_len;
                size_t idx = 0;
                // Compute correct position inside the block
                while (i >= first_len) {
                    first_len += block_len;
                    idx++;
                }
                size_t r_len = rest_len(i, len, first_len);

                // Extract the terminals of the current non-terminal inside the block
                extract(i - idx * block_len, len - r_len, first, s, rlslp);

                // Extract remaining terminals of the block if they are part of the substring to be extracted
                while (r_len > 0) {
                    if (r_len > block_len) {
                        r_len = block_len;
                    }
                    extract(0, r_len, first, s, rlslp);
                    first_len += block_len;
                    r_len = rest_len(i, len, first_len);
                }
            } else {
                size_t r_len = rest_len(i, len, first_len);

                // Check if position is in left subtree
                if (i < first_len) {
                    extract(i, len - r_len, first, s, rlslp);

                    // Compute also right subtree if i + len >= first_len
                    if (r_len > 0) {
                        extract(0, r_len, second, s, rlslp);
                    }
                } else {  // Check if position is in left subtree
                    extract(i - first_len, len, second, s, rlslp);
                }
            }
        }
    }
}

size_t lcequery(rlslp& rlslp,
                const size_t i,
                const size_t j,
                const variable_t i_nt,
                const variable_t j_nt,
                size_t traverse,
                std::vector<size_t>& s_len,
                std::vector<lceq>& i_nts,
                std::vector<lceq>& j_nts);

/**
     * @brief Determines the next non-terminal und position inside the subtree to match.
     *
     * @param pos[in,out] The position inside the next non-terminal to match
     * @param nt[in,out] The next non-terminal to match
     */
inline void subtree(rlslp& rlslp, size_t& pos, variable_t& nt) {
    if (nt >= rlslp.terminals) {
        auto left = rlslp.non_terminals[nt].production[0];
        auto right = rlslp.non_terminals[nt].production[1];
        auto len = rlslp.non_terminals[left].len;

        if (right < 0/*rlslp.blocks[nt]*/) {
            auto block_len = 0;
            while (block_len + len <= pos) {
                block_len += len;
            }
            pos = pos - block_len;
            nt = left;
        } else {
            if (pos < len) {
                nt = left;
            } else {
                pos = pos - len;
                nt = right;
            }
        }
    }
}

/**
 * @brief Traverses the derivation tree upwards to find the next subtrees containing the next leaves to match.
 *
 * @param i The position of the first suffix to compare next
 * @param j The position of the second suffix to compare next
 * @param i_nt The current non-terminal to compare for the first suffix
 * @param j_nt The current non-terminal to compare for the second suffix
 * @param traverse The number of the next traverse
 * @param i_nts Seen non-terminal for the first suffix
 * @param j_nts Seen non-terminal for the second suffix
 *
 * @return The length of the longest common extension
 */
size_t find_next_nt(rlslp& rlslp,
                    size_t i,
                    size_t j,
                    const variable_t i_nt,
                    const variable_t j_nt,
                    std::vector<size_t>& s_len,
                    size_t traverse,
                    std::vector<lceq>& i_nts,
                    std::vector<lceq>& j_nts) {
    variable_t p_i_nt = i_nt;

    // Length of the last matched non-terminals
    size_t len = s_len[traverse - 1] - s_len[traverse - 2];
    size_t comp_i_len = len;
    // If the parents traverse is older than the traverse of the current non-terminal we need to take the
    // difference of the actual lce length and the lce length of the traverse of the parent.
    if (i_nts[p_i_nt].traverse + 1 != traverse) {
        comp_i_len = s_len[traverse - 1] - s_len[i_nts[p_i_nt].traverse - 1];
    }
//        std::cout << "s_len:" << s_len[traverse-1] << ", i:" << i << ", len: " << rlslp.non_terminals[i_nt].len << std::endl;
//        std::cout << "comp_i_len: " << comp_i_len << std::endl;

    // While the position inside the non-terminal plus the comp_i_len is greater or equal the length of the string
    // derived by this non-terminal we have to traverse upward.
    // We reached the root if p_i_nt is equal to the parent of p_i_nt
    while (i + comp_i_len >= rlslp.non_terminals[p_i_nt].len && p_i_nt != i_nts[p_i_nt].nt) {
//            std::cout << "up for i" << std::endl;
//            std::cout << p_i_nt << ", pos: " << i_nts[p_i_nt].pos << ", nt: " << i_nts[p_i_nt].nt << ", tr: " << i_nts[p_i_nt].traverse << std::endl;
        p_i_nt = i_nts[p_i_nt].nt;
        i = i_nts[p_i_nt].pos;

        // Compute next comp_i_len
        if (i_nts[p_i_nt].traverse + 1 == traverse) {
            comp_i_len = len;
        } else {
            comp_i_len = s_len[traverse - 1] - s_len[i_nts[p_i_nt].traverse - 1];
        }
//            std::cout << "comp_i_len: " << comp_i_len << std::endl;
    }
    i_nts[p_i_nt].traverse = traverse;
    i_nts[p_i_nt].pos = i + comp_i_len;
//        std::cout << "point to: " << p_i_nt << ", pos: " << i_nts[p_i_nt].pos << ", nt: " << i_nts[p_i_nt].nt << ", tr: " << i_nts[p_i_nt].traverse << std::endl;

    variable_t p_j_nt = j_nt;
    size_t comp_j_len = len;
    if (j_nts[p_j_nt].traverse + 1 != traverse) {
        comp_j_len = s_len[traverse - 1] - s_len[j_nts[p_j_nt].traverse - 1];
    }
//        std::cout << "s_len:" << s_len[traverse-1] << ", j:" << j << ", len: " << rlslp.non_terminals[j_nt].len << std::endl;
//        std::cout << "comp_j_len: " << comp_j_len << std::endl;
    while (j + comp_j_len >= rlslp.non_terminals[p_j_nt].len && p_j_nt != j_nts[p_j_nt].nt) {
//            std::cout << "up for j" << std::endl;
//            std::cout << p_j_nt << ", pos: " << j_nts[p_j_nt].pos << ", nt: " << j_nts[p_j_nt].nt << ", tr: " << j_nts[p_j_nt].traverse << std::endl;
        p_j_nt = j_nts[p_j_nt].nt;
        j = j_nts[p_j_nt].pos;

        if (j_nts[p_j_nt].traverse + 1 == traverse) {
            comp_j_len = len;
        } else {
            comp_j_len = s_len[traverse - 1] - s_len[j_nts[p_j_nt].traverse - 1];
        }
//            std::cout << "comp_j_len: " << comp_j_len << std::endl;
    }
    j_nts[p_j_nt].traverse = traverse;
    j_nts[p_j_nt].pos = j + comp_j_len;
//        std::cout << "point to: " << p_j_nt << ", pos: " << j_nts[p_j_nt].pos << ", nt: " << j_nts[p_j_nt].nt << ", tr: " << j_nts[p_j_nt].traverse << std::endl;

    // If we reached the root and are looking for a position greater than the maximum length of the string, we are finished
    if ((p_i_nt == i_nts[p_i_nt].nt && i + comp_i_len >= rlslp.non_terminals[p_i_nt].len) ||
        (p_j_nt == j_nts[p_j_nt].nt && j + comp_j_len >= rlslp.non_terminals[p_j_nt].len)) {
        return 0;
    } else {  // Traverse down to the i + comp_i_lenth and j + comp_j_lenth leaves
//            std::cout << "call" << std::endl;
//            std::cout << "i: " << i << ", s_len: " << comp_i_len << ", nt: " << p_i_nt << std::endl;
//            std::cout << "j: " << j << ", s_len: " << comp_j_len << ", nt: " << p_j_nt << std::endl;
        return lcequery(rlslp, i + comp_i_len, j + comp_j_len, p_i_nt, p_j_nt, traverse, s_len, i_nts, j_nts);
    }
}


/**
 * @brief Returns the length of the longest common extension.
 *
 * @param i Begin position of the first suffix
 * @param j Begin position of the second suffix
 * @param i_nt Current non-terminal of the first suffix
 * @param j_nt Current non-terminal of the second suffix
 * @param traverse The number of the downward traverse
 * @param s_len The cumulated sum of the lengths of the matching non-terminals
 * @param i_nts Visited non-terminals for the first suffix
 * @param j_nts Visited non-terminals for the second suffix
 *
 * @return The length of the longest common extension
 */
size_t lcequery(rlslp&/*<variable_t, terminal_count_t>&*/ rlslp,
                const size_t i,
                const size_t j,
                const variable_t i_nt,
                const variable_t j_nt,
                size_t traverse,
                std::vector<size_t>& s_len,
                std::vector<lceq>& i_nts,
                std::vector<lceq>& j_nts) {
    size_t pos_i = i, pos_j = j;
    variable_t nt_i = i_nt, nt_j = j_nt;
    subtree(rlslp, pos_i, nt_i);
    subtree(rlslp, pos_j, nt_j);

    if (i_nt != nt_i) {
        i_nts[nt_i].nt = i_nt;
        i_nts[nt_i].pos = pos_i;
        i_nts[nt_i].traverse = traverse;
    }
    if (j_nt != nt_j) {
        j_nts[nt_j].nt = j_nt;
        j_nts[nt_j].pos = pos_j;
        j_nts[nt_j].traverse = traverse;
    }

//        std::cout << pos_i << ": " << nt_i << ", " << pos_j << ": " << nt_j << std::endl;

    if (nt_i == nt_j && pos_i == pos_j) {
        // if both non-terminals are the same and starting at the same position
//            std::cout << "Begin same level" << std::endl;
//            print(i_nts, j_nts);

        auto len = rlslp.non_terminals[nt_i].len - pos_i;
        if (s_len.empty()) {
            s_len.push_back(len);
        } else {
            s_len.push_back(len + s_len[s_len.size()-1]);
        }
        return len + find_next_nt(rlslp, i, j, i_nt, j_nt, s_len, traverse + 1, i_nts, j_nts);
    } else if (rlslp.non_terminals[nt_i].production[1] < 0 && rlslp.non_terminals[nt_j].production[1] < 0 && /*rlslp.blocks[nt_i] && rlslp.blocks[nt_j] &&*/
               rlslp.non_terminals[nt_i].production[0] == rlslp.non_terminals[nt_j].production[0] &&
               pos_i % rlslp.non_terminals[rlslp.non_terminals[nt_i].production[0]].len ==
               pos_j % rlslp.non_terminals[rlslp.non_terminals[nt_i].production[0]].len) {
        // if two blocks with same (non-)terminal but different length start at the same position within one of
        // the block symbols
//            std::cout << "block same symbol" << std::endl;
//            std::cout << nt_i << "," << nt_j << std::endl;
//            print(i_nts, j_nts);

        auto len = std::min(rlslp.non_terminals[nt_i].len - pos_i, rlslp.non_terminals[nt_j].len - pos_j);
        if (s_len.empty()) {
            s_len.push_back(len);
        } else {
            s_len.push_back(len + s_len[s_len.size()-1]);
        }
        return len + find_next_nt(rlslp, pos_i, pos_j, nt_i, nt_j, s_len, traverse + 1, i_nts, j_nts);
    } else if (j_nts[nt_i].nt != rlslp.non_terminals.size() && j_nts[nt_i].pos == pos_i && j_nts[nt_i].traverse == traverse) {
        // if we found a non-terminal we have seen for the other suffix starting at the same position
//            std::cout << "Begin i in j" << std::endl;
//            print(i_nts, j_nts);
        auto len = rlslp.non_terminals[nt_i].len - pos_i;
        if (s_len.empty()) {
            s_len.push_back(len);
        } else {
            s_len.push_back(len + s_len[s_len.size()-1]);
        }
//            std::cout << nt_i << ", " << pos_i << ": " << rlslp.non_terminals[nt_i].len << std::endl;
        return len + find_next_nt(rlslp, pos_i, pos_j, nt_i, nt_j, s_len, traverse + 1, i_nts, j_nts);
    } else if (i_nts[nt_j].nt != rlslp.non_terminals.size() && i_nts[nt_j].pos == pos_j && i_nts[nt_j].traverse == traverse) {
        // if we found a non-terminal we have seen for the other suffix starting at the same position
//            std::cout << "Begin j in i" << std::endl;
//            print(i_nts, j_nts);
        auto len = rlslp.non_terminals[nt_j].len - pos_j;
        if (s_len.empty()) {
            s_len.push_back(len);
        } else {
            s_len.push_back(len + s_len[s_len.size()-1]);
        }
//            std::cout << nt_j << ", " << pos_j << ": " << rlslp.non_terminals[nt_j].len << std::endl;
        return len + find_next_nt(rlslp, pos_i, pos_j, nt_i, nt_j, s_len, traverse + 1, i_nts, j_nts);
    } else if (nt_i < rlslp.terminals && nt_j < rlslp.terminals) {
        // if we traversed down to the leaves for both suffices
//            std::cout << "Begin terminals" << std::endl;
//            print(i_nts, j_nts);
        if (nt_i == nt_j) {
            // if both leaves contain the same terminal
            auto len = rlslp.non_terminals[nt_i].len;
            if (s_len.empty()) {
                s_len.push_back(len);
            } else {
                s_len.push_back(len + s_len[s_len.size()-1]);
            }
            return len + find_next_nt(rlslp, i, j, nt_i, nt_j, s_len, traverse + 1, i_nts, j_nts);
        } else {
            // The terminals of the leaves don't match
            return 0;
        }
    } else {
        return lcequery(rlslp, pos_i, pos_j, nt_i, nt_j, traverse, s_len, i_nts, j_nts);
    }
}



void print(const std::vector<lceq>& i_nts,
           const std::vector<lceq>& j_nts) {
    for (size_t i = 0; i < i_nts.capacity(); ++i) {
        std::cout << i << ": " << i_nts[i].nt << ", " << i_nts[i].pos << ", " << i_nts[i].traverse << ": " << j_nts[i].nt << ", " << j_nts[i].pos << ", " << j_nts[i].traverse << std::endl;
    }
    std::cout << "End" << std::endl;
}


/**
 * @brief Builds the straight-line program generating the given text using the recompression technique.
 *
 * A straight-line program (SLP) is in particular a context free grammar in Chomsky normal form.
 * First all letters in the text are replaced by non-terminals which derive the letters. Then
 * the block compression bcomp and the pair compression pcomp alternately compress blocks and local
 * pairs in the texts resulting of the previous compression. This will be repeated until the text
 * contains only one letter.
 *
 * @param text The text
 */
void recomp(text_t& t, rlslp& rlslp) {
    variable_t alphabet_size = 0;
    size_t text_size = t.size();
    std::vector<variable_t> mapping;

    replace_letters(t, rlslp, alphabet_size, mapping);

    bool not_finished = text_size > 1;

    while (not_finished) {
//            std::cout << "BComp" << std::endl;
        bcomp(t, rlslp, text_size, alphabet_size, mapping);
//            std::cout << "Alpha" << std::endl;
        compute_alphabet(t, text_size, alphabet_size, mapping);
        not_finished = text_size > 1;
        if (not_finished) {
//                std::cout << "PComp" << std::endl;
            pcomp(t, rlslp, text_size, alphabet_size, mapping);
//                std::cout << "Alpha" << std::endl;
            compute_alphabet(t, text_size, alphabet_size, mapping);
            not_finished = text_size > 1;
        }
    }
}

/**
 * @brief Extracts the substring beginning at position @code{i} and length @code{len}.
 *
 * Note that if i + length is greater than the length of the text represented by
 * the SLP the substring from i to the end of the text will be returned.
 *
 * If i is greater than the length of the text represented by the SLP an empty
 * string will be returned.
 *
 * @param i The start index
 * @param len The length of the substring to extract
 *
 * @return The substring at position @code{i} and length @code{len}
 */
text_t extract(rlslp& rlslp, const size_t i, const size_t len) {
    text_t extracted_string;
    if (!rlslp.non_terminals.empty() && i < rlslp.non_terminals[rlslp.non_terminals.size()-1].len && len > 0) {
        extract(i, len, rlslp.non_terminals.size()-1, extracted_string, rlslp);
    }
    return extracted_string;
}

/**
 * @brief Returns the length of the longest common extension.
 *
 * @param i The start position of the first suffix
 * @param j The start position of the second suffix
 * @return The length of the longest common extension
 */
size_t lcequery(rlslp& rlslp, const size_t i, const size_t j) {
    if (rlslp.non_terminals.empty()) {
        return 0;
    }

    variable_t nt_root = rlslp.non_terminals.size() - 1;
    auto root = rlslp.non_terminals[nt_root];
    if (i >= root.len || j >= root.len) {
        return 0;
    }
    if (i == j) {
        return root.len - i;
    }

    variable_t size = rlslp.non_terminals.size();
    std::vector<lceq> i_non_terms(size, lceq{size, (size_t)size, 0});
    std::vector<lceq> j_non_terms(size, lceq{size, (size_t)size, 0});
    i_non_terms[nt_root].nt = nt_root;
    j_non_terms[nt_root].nt = nt_root;
    i_non_terms[nt_root].pos = i;
    j_non_terms[nt_root].pos = j;
    i_non_terms[nt_root].traverse = 1;
    j_non_terms[nt_root].traverse = 1;

    std::vector<size_t> len;
    len.push_back(0);

    return lcequery(rlslp, i, j, nt_root, nt_root, 1, len, i_non_terms, j_non_terms);
}

}  // namespace lce


void read_file(const std::string& file_name, std::string& text, bool term) {
    std::cout << "Reading file: " << file_name << std::endl;
    std::ifstream ifs(file_name.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
    auto file_size = ifs.tellg();
    if (file_size < 0) {
        std::cerr << "Cannot read file" << std::endl;
        return;
    }
    ifs.seekg(0, std::ios::beg);

    if (term) {
        text.resize(((size_t)file_size)+1, '\0');
    } else {
        text.resize(file_size, '\0');
    }
    ifs.read((char*)text.data(), file_size);
    if (term) {
        text[file_size] = '$';
    }
    ifs.close();
    std::cout << "Read " << file_size << " bytes" << std::endl;
    std::cout << "Finished reading file" << std::endl;
}

void write_to_file(const std::string& file_name, rlslp& rlslp) {
    std::cout << "Wrinting to file " << file_name << ".rcmp" << std::endl;
    std::ofstream ofs(file_name + ".rcmp", std::ios::out | std::ios::binary);
    lce::terminal_count_t terms = rlslp.terminals;
    size_t size = rlslp.non_terminals.size();
    ofs.write((char*)&terms, sizeof(lce::terminal_count_t));
    ofs.write((char*)&size, sizeof(size_t));
    size_t i = 0;
    char term;
    while (i < terms) {
        term = (char)rlslp.non_terminals[i++].production[0];
        ofs.write(&term, sizeof(char));
    }
    for (; i < rlslp.non_terminals.size(); ++i) {
        ofs.write((char*)&rlslp.non_terminals[i].production[0], sizeof(lce::variable_t));
        ofs.write((char*)&rlslp.non_terminals[i].production[1], sizeof(lce::variable_t));
//        ofs << factors[i].first << factors[i].second;
    }
    ofs.close();
    std::cout << "Finished writing file" << std::endl;
}

void load_from_file(const std::string& file_name, rlslp& rlslp) {
    std::cout << "Reading file " << file_name << ".rcmp" << std::endl;
    std::ifstream istream(file_name + ".rcmp", std::ios::in | std::ios::binary);
    lce::terminal_count_t terms;
    size_t size;
    istream.read((char*)&terms, sizeof(lce::terminal_count_t));
    istream.read((char*)&size, sizeof(size_t));
    rlslp.terminals = terms;
    rlslp.non_terminals.resize(size);
    size_t i = 0;
    lce::variable_t var;
    char term;
    while (i < terms) {
        istream.read(&term, sizeof(char));
        rlslp.non_terminals[i++].production.emplace_back((lce::variable_t)((unsigned char)term));
    }
    for (; i < size; ++i) {
        istream.read((char*)&var, sizeof(lce::variable_t));
        rlslp.non_terminals[i].production.emplace_back(var);
        istream.read((char*)&var, sizeof(lce::variable_t));
        rlslp.non_terminals[i].production.emplace_back(var);
    }
    istream.close();
    std::cout << "Finished reading file" << std::endl;
}

bool equal(const rlslp& exp_rlslp, const rlslp& chck_rlslp) {
    if (exp_rlslp.terminals != chck_rlslp.terminals) {
        std::cerr << "Terminals not equal" << std::endl;
        return false;
    }
    if (exp_rlslp.non_terminals.size() != chck_rlslp.non_terminals.size()) {
        std::cerr << "Number of non-terminals not equal" << std::endl;
        return false;
    }
    size_t i = 0;
    while (i < exp_rlslp.terminals) {
        if (exp_rlslp.non_terminals[i].production[0] != chck_rlslp.non_terminals[i].production[0]) {
            std::cerr << "Wrong terminal" << std::endl;
            return false;
        }
        i++;
    }

    for (; i < exp_rlslp.non_terminals.size(); ++i) {
        if (exp_rlslp.non_terminals[i].production[0] != chck_rlslp.non_terminals[i].production[0] ||
                exp_rlslp.non_terminals[i].production[1] != chck_rlslp.non_terminals[i].production[1]) {
            std::cerr << "Wrong non-terminal" << std::endl;
            return false;
        }
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "arg error" << std::endl;
        return -1;
    }
    std::string filename = argv[1];

    std::string str;
    read_file(filename, str, false);

    std::vector<lce::variable_t> istr(str.length());
    for (uint i = 0; i < str.length(); ++i) {
        istr[i] = (lce::variable_t)((unsigned char)str[i]);
    }

    std::cout << "vector copied" << std::endl;

    rlslp slp;
    std::cout << "Text size: " << istr.size() << std::endl;
    const auto startTime = std::chrono::system_clock::now();
    lce::recomp(istr, slp);
    const auto endTime = std::chrono::system_clock::now();
    const auto timeSpan = endTime - startTime;
    std::cerr << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
    std::cerr << "rlslp size: " << slp.non_terminals.size() << std::endl;

    { //// Checking correctness.
        std::cerr << "Checking if rlslp can be decompressed correctly..." << std::endl;
        std::string decomp_str;
        std::vector<lce::variable_t> dec = lce::extract(slp, 0, str.size());
        for (auto i : dec) {
            decomp_str.push_back(((char)i));
        }
        if (decomp_str == str) {
            std::cerr << "Decompressed correctly." << std::endl;
        } else {
            std::cerr << "Decompression failed." << std::endl;
        }
    }

    /*write_to_file(filename, slp);

    rlslp load_slp;
    load_from_file(filename, load_slp);

    if (equal(slp, load_slp)) {
        std::cout << "RLSLP correctly stored and loaded." << std::endl;
    } else {
        std::cerr << "Storing and loading failed." << std::endl;
    }*/

    return 0;
}
