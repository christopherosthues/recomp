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
 
#pragma once

#include <cstdint>
#include <string>
#include <vector>

//template<typename V = size_t, typename T = std::uint32_t>
struct rlslp {
    typedef std::int32_t variable_t;
    typedef std::uint32_t terminal_count_t;
//    typedef V variable_t;
//    typedef T terminal_count_t;

    struct non_terminal {
        // Length of the (sub)string derived by the non-terminal
        size_t len;
        
        // Derivation rule of the non-terminal
        std::vector<variable_t> production;
        
        non_terminal();

        non_terminal(variable_t term, size_t len);

        non_terminal(variable_t first, variable_t second, size_t len);

        bool operator==(const non_terminal &nt) const;

        std::string to_string() const;
    };

    // All non-terminals are implicitly given by their index
    std::vector<non_terminal> non_terminals;
    
    // Bitvector to indicate whether a derivation rule of the non-terminal given by its index is a block
    //std::vector<bool> blocks;
    //sdsl::bit_vector blocks;
    
    // The number of terminals (terminals - 1 is the last production which derives a terminal)
    terminal_count_t terminals = 0;
    
    size_t block_count = 0;

    std::string to_string() const;
};

