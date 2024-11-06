#include <sstream>
#include <vector>
#include "common.h"


void iterate_over_variant(const InputData &input_data,
    const std::vector<std::string> &hashtable, std::vector<std::string> &results) {
    std::ifstream variant(input_data.get_variant());
    check_open(variant);

    std::string line;
    // Skip the first two lines with non-important data.
    std::getline(variant, line);
    std::getline(variant, line);
    int curr_line = 1;
    while (std::getline(variant, line)) {
        curr_line++;
        if (results[curr_line] != "")
            continue;

        auto gen_info = GenomeInfo(line, true);
        auto hash = gen_info.get_hash();
        hash *= Consts::MAX_COLLISIONS;
        
        // Iterate over collisions as long as needed.
        for (std::size_t i = 0; i < Consts::MAX_COLLISIONS; i++) {
            if (hashtable[hash] == "")
                break;
            if (gen_info == GenomeInfo(hashtable[hash], false)) {
                results[curr_line] = hashtable[hash];
                break;
            }
            hash++;
        }
    }

    variant.close();
}


void process_data(const InputData &input_data) {
    std::ifstream index(input_data.get_index());
    check_open(index);
    std::ifstream variant(input_data.get_variant());
    check_open(variant);

    int variant_size = 0;
    std::string line;
    while (std::getline(variant, line))
        variant_size++;
    variant.close();

    hash_t hash, repr;
    std::string row;
    int processed = 0;
    std::vector<std::string> hashtable(Consts::MAX_HASHTABLE_SIZE);
    // Array used for saving the results which later will be written to output file.
    std::vector<std::string> results(variant_size);

    while (std::getline(index, line)) {
        std::istringstream iss(line);
        iss >> hash >> repr;
        std::getline(iss >> std::ws, row);
        // Scale the hash.
        hash *= Consts::MAX_COLLISIONS;
        while (hashtable[hash] != "") {
            hash++;
        }
        hashtable[hash] = row;

        processed++;
        if (processed == Consts::CHUNK_SIZE) {
            iterate_over_variant(input_data, hashtable, results);
            // Reset the hashmap.
            processed = 0;
            std::fill(hashtable.begin(), hashtable.end(), "");
        }
    }

    if (processed != 0)
        iterate_over_variant(input_data, hashtable, results);

    std::ofstream output(input_data.get_output());
    check_open(output);

    for (int i = 0; i < variant_size; i++) {
        if (results[i] != "")
            output << results[i] << std::endl;
    }

    output.close();
    index.close();
}


int main(int argc, char *argv[]) {
    InputData input_data(argc, argv);

    if (input_data.is_indexing()) {
        DatabaseParser parser(input_data);
        parser.parse_database();
    }
    else {
        process_data(input_data);
    }

    return 0;
}
