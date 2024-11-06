// Common parts of GPU and CPU code.

#ifndef __COMMON_H__
#define __COMMON_H__

#include <fstream>
#include <string>

using hash_t = std::size_t;
using chrom_t = uint32_t;
using position_t = uint64_t;


void check_open(std::ifstream &file);
void check_open(std::ofstream &file);


class Consts {
public:
    static const std::size_t CHUNK_SIZE = 512 * 1024;
    static const std::size_t MAX_HASH = 1024 * 1024;
    static const std::size_t MAX_COLLISIONS = 10;
    static const std::size_t MAX_HASHTABLE_SIZE = MAX_HASH * MAX_COLLISIONS;
    static const std::size_t INF = MAX_HASHTABLE_SIZE + 1;
    static const int BLOCKS = 32;
    static const int THREADS_PER_BLOCK = 256;
};


class InputData {
private:
    bool indexing;
    std::string database, variant, index, output;

public:
    InputData(std::string database, std::string index):
        indexing(true),
        database(database),
        index(index) {}

    InputData(std::string variant, std::string index, std::string output):
        indexing(false),
        variant(variant),
        index(index),
        output(output) {}

    InputData(int argc, char *argv[]);

    bool is_indexing();

    std::string get_database() const;

    std::string get_variant() const;

    std::string get_index() const;

    std::string get_output() const;
};


class DatabaseParser {
private:
    std::string database, index;

public:
    DatabaseParser(std::string database, std::string index):
        database(database), index(index) {}

    DatabaseParser(const InputData &input_data);

    void parse_database();
};


// Information which will be matched from variant to database.
class GenomeInfo {
private:
    hash_t hash, representation;

public:
    GenomeInfo(const std::string &line, bool is_variant);

    hash_t get_hash();

    hash_t get_representation();

    bool operator==(const GenomeInfo&) const = default;
};


#endif /* __COMMON_H__ */
