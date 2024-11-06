#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <cassert>
#include "common.h"


namespace {
    // Hash used to place genome in the hashtable.
    hash_t hasher(chrom_t chrom, position_t pos, char ref, char alt) {
        hash_t hash = 0;
        hash_t xorshift_const = 0x9e3779b9;
        hash ^= std::hash<chrom_t>()(chrom) + xorshift_const + (hash << 6) + (hash >> 2);
        hash ^= std::hash<position_t>()(pos) + xorshift_const + (hash << 6) + (hash >> 2);
        hash ^= std::hash<char>()(ref) + xorshift_const + (hash << 6) + (hash >> 2);
        hash ^= std::hash<char>()(alt) + xorshift_const + (hash << 6) + (hash >> 2);
        hash %= Consts::MAX_HASH;
        return hash;
    }

    // One to one hash. Encrypting (Chrom, Pos, Ref, Alt) on 64 bits.
    hash_t unique_representation(chrom_t chrom, position_t pos, char ref, char alt) {
        hash_t hash = 0;
        // 8 bits per chrom / alt / ref.
        hash = chrom;
        hash ^= hash_t(ref) << 8;
        hash ^= hash_t(alt) << 16;
        // Rest of the bits for pos.
        hash ^= (pos << 24);
        return hash;
    }
}


void check_open(std::ifstream &file) {
    if (!file.is_open()) {
        std::cerr << "ERROR: File is not open." << std::endl;
        exit(2);
    }
}

void check_open(std::ofstream &file) {
    if (!file.is_open()) {
        std::cerr << "ERROR: File is not open." << std::endl;
        exit(2);
    }
}


GenomeInfo::GenomeInfo(const std::string &line, bool is_variant) {
    // Used to store chromosoms as integers.
    const static std::map<char, int> chroms = {
        {'M', 23},
        {'X', 24},
        {'Y', 25},
    };
    
    std::istringstream iss(line);
    std::string id, chrom_str;
    chrom_t chrom;
    position_t pos;
    char ref, alt;

    if (!is_variant)
        iss >> chrom_str >> pos >> ref >> alt;
    else
        iss >> chrom_str >> pos >> id >> ref >> alt;

    // Parse chromosome to integer.
    if (chrom_str.size() == 1) {
        auto iter = chroms.find(chrom_str[0]);
        if (iter != chroms.end())
            chrom = iter->second;
        else
            chrom = chrom_str[0] - '0';
    }
    else {
        chrom = (chrom_str[0] - '0') * 10 + chrom_str[1] - '0';
    }

    this->hash = hasher(chrom, pos, ref, alt);
    this->representation = unique_representation(chrom, pos, ref, alt);
}

hash_t GenomeInfo::get_hash() {
    return this->hash;
}

hash_t GenomeInfo::get_representation() {
    return this->representation;
}


InputData::InputData(int argc, char *argv[]) {
    const std::string OPTION_INDEX = "-i";

    if (argc != 4) {
        std::cerr << "ERROR: Wrong number of arguments!" << std::endl;
        exit(1);
    }

    if (std::string(argv[1]) == OPTION_INDEX)
        *this = {argv[2], argv[3]};
    else
        *this = {argv[1], argv[2], argv[3]};
}

bool InputData::is_indexing() {
    return this->indexing;
}

std::string InputData::get_database() const {
    return this->database;
}

std::string InputData::get_variant() const {
    return this->variant;
}

std::string InputData::get_index() const {
    return this->index;
}

std::string InputData::get_output() const {
    return this->output;
}


DatabaseParser::DatabaseParser(const InputData &input_data) {
    *this = {input_data.get_database(), input_data.get_index()};
}

void DatabaseParser::parse_database() {
    std::ifstream database(this->database);
    check_open(database);
    std::ofstream index(this->index);
    check_open(index);

    std::string line;
    int max_collisions = 0, processed = 0;
    std::vector<int> collisions(Consts::MAX_HASH);
    // Skip the first line with columns names.
    std::getline(database, line);
    while (std::getline(database, line)) {
        auto gen_info = GenomeInfo(line, false);
        auto hash = gen_info.get_hash();
        auto repr = gen_info.get_representation();
        index << hash << " " << repr << " " << line << std::endl;
        collisions[hash]++;
        max_collisions = std::max(max_collisions, collisions[hash]);

        processed++;
        if (processed == Consts::CHUNK_SIZE) {
            processed = 0;
            std::fill(collisions.begin(), collisions.end(), 0);
        }
    }

    // Make sure that the hashing function is good enough.
    assert((size_t) max_collisions <= Consts::MAX_COLLISIONS);

    database.close();
    index.close();
}
