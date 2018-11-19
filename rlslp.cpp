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
 
#include "rlslp.hpp"

#include <iostream>
#include <sstream>

rlslp::non_terminal::non_terminal() = default;

rlslp::non_terminal::non_terminal(variable_t term, size_t len) : production(1) {
    this->production[0] = term;//.emplace_back(term);
    //this->production.resize(1);
    //this->production.shrink_to_fit();
    this->len = len;
}

rlslp::non_terminal::non_terminal(variable_t first, variable_t second, size_t len) : production(2) {
    this->production[0] = first;//.emplace_back(first);
    this->production[1] = second;//.emplace_back(second);
    //this->production.resize(2);
    //this->production.shrink_to_fit();
    this->len = len;
}

bool rlslp::non_terminal::operator==(const non_terminal &nt) const {
    return this->production == nt.production && this->len == nt.len;
}

std::string rlslp::non_terminal::to_string() const {
    std::string sec = (production.size() > 1)? std::to_string(production[1]) : "";
    return "(" + std::to_string(production[0]) + "," + sec + ")," +
           std::to_string(len);
}

std::string rlslp::to_string() const {
    std::stringstream sstream;
    sstream << "Number of terminals: " << terminals << std::endl;
    /*sstream << "Blocks: ";
    for (const auto& b : blocks) {
        std::cout << b << std::endl;
    }*/
    sstream << "Non-terminals: " << std::endl;
    size_t i = 0;
    for (const auto& nt : non_terminals) {
        sstream << i << ": " << nt.to_string() << std::endl;
        i++;
    }
    return sstream.str();
}
