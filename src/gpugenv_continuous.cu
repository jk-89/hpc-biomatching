#include <fstream>
#include <sstream>
#include <vector>
#include <functional>
#include "common.h"


__global__ void iterate_over_variant(const hash_t *hashtable, size_t variant_size,
    const hash_t *variant_hash, const hash_t *variant_repr, hash_t *results) {
    size_t total_blocks = gridDim.x;
    size_t block_id = blockIdx.x;
    size_t elements_per_block = (variant_size + total_blocks - 1) / total_blocks;

    size_t start_idx = elements_per_block * block_id;
    size_t end_idx = start_idx + elements_per_block;
    if (variant_size < end_idx)
        end_idx = variant_size;

    // In this version we make the for loop iterate over the continuous fragment in memory.
    size_t total_threads = blockDim.x;
    size_t thread_id = threadIdx.x;
    size_t elements_per_thread = (end_idx - start_idx + total_threads - 1) / total_threads;
    
    size_t start_delta = start_idx + elements_per_thread * thread_id;
    size_t end_delta = start_delta + elements_per_thread;
    if (variant_size < end_delta)
        end_delta = variant_size;

    for (size_t i = start_delta; i < end_delta; i++) {
        auto hash = variant_hash[i];
        hash *= Consts::MAX_COLLISIONS;

        // Iterate over collisions as long as needed.
        for (std::size_t j = 0; j < Consts::MAX_COLLISIONS; j++) {
            if (hashtable[hash] == Consts::INF)
                break;
            if (hashtable[hash] == variant_repr[i]) {
                results[i] = hash;
                break;
            }
            hash++;
        }
    }
}


void save_results(const std::vector<std::string> &hashtable,
    const std::vector<hash_t> &results, std::vector<std::string> &results_str) {
    int variant_size = results.size();
    for (int i = 0; i < variant_size; ++i) {
        if (results[i] != Consts::INF)
            results_str[i] = std::move(hashtable[results[i]]);
    }
}


// The main idea is to pass to GPU only hashes and avoid passing strings.
// Unfortunately this causes the code to be much more complicated.
void process_data(const InputData &input_data) {
    std::ifstream index(input_data.get_index());
    check_open(index);
    std::ifstream variant(input_data.get_variant());
    check_open(variant);

    int variant_size = 0;
    std::string line;
    while (std::getline(variant, line))
        variant_size++;
    // Skip the first two lines with non-important data.
    variant_size -= 2;

    // Save the content of the variant file in program memory.
    std::vector<hash_t> variant_hash(variant_size);
    std::vector<hash_t> variant_repr(variant_size);
    variant.clear();
    variant.seekg(0);
    // Skip the first two lines with non-important data.
    std::getline(variant, line);
    std::getline(variant, line);
    for (int i = 0; i < variant_size; i++) {
        std::getline(variant, line);
        auto gen_info = GenomeInfo(line, true);
        variant_hash[i] = gen_info.get_hash();
        variant_repr[i] = gen_info.get_representation();
    }
    variant.close();

    // Initialize cuda arrays.
    size_t variant_bytes = variant_size * sizeof(hash_t);
    hash_t *cuda_variant_hash, *cuda_variant_repr;
    hash_t *cuda_hashtable;
    hash_t *cuda_results;
    cudaMalloc((void**)&cuda_variant_hash, variant_bytes);
    cudaMalloc((void**)&cuda_variant_repr, variant_bytes);
    cudaMalloc((void**)&cuda_hashtable, Consts::MAX_HASHTABLE_SIZE * sizeof(hash_t));
    cudaMalloc((void**)&cuda_results, variant_bytes);
    std::vector<hash_t> results(variant_size, Consts::INF);
    std::vector<std::string> results_str(variant_size);

    cudaMemcpy(cuda_variant_hash, variant_hash.data(), variant_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_variant_repr, variant_repr.data(), variant_bytes, cudaMemcpyHostToDevice);

    hash_t hash, repr;
    std::string row;
    int processed = 0;
    // Local hashtable contains pairs (hash, line).
    // Other hashtable contains pairs (hash, representation).
    std::vector<std::string> hashtable(Consts::MAX_HASHTABLE_SIZE);
    std::vector<hash_t> hashtable_repr(Consts::MAX_HASHTABLE_SIZE, Consts::INF);

    std::function<void()> execute_kernel = [&]() {
        // Copy the hashtable to GPU.
        cudaMemcpy(
            cuda_hashtable, hashtable_repr.data(),
            Consts::MAX_HASHTABLE_SIZE * sizeof(hash_t), cudaMemcpyHostToDevice
        );
        cudaMemcpy(cuda_results, results.data(), variant_bytes, cudaMemcpyHostToDevice);
        iterate_over_variant<<<Consts::BLOCKS, Consts::THREADS_PER_BLOCK>>>(
            cuda_hashtable, variant_size,
            cuda_variant_hash, cuda_variant_repr, cuda_results
        );
        cudaMemcpy(results.data(), cuda_results, variant_bytes, cudaMemcpyDeviceToHost);
        save_results(hashtable, results, results_str);
    };

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
        hashtable_repr[hash] = repr;

        processed++;
        if (processed == Consts::CHUNK_SIZE) {
            execute_kernel();

            // Reset the arrays.
            processed = 0;
            std::fill(hashtable.begin(), hashtable.end(), "");
            std::fill(hashtable_repr.begin(), hashtable_repr.end(), Consts::INF);
            std::fill(results.begin(), results.end(), Consts::INF);
        }
    }

    if (processed != 0)
        execute_kernel();

    std::ofstream output(input_data.get_output());
    check_open(output);
    for (int i = 0; i < variant_size; ++i) {
        if (results_str[i] != "")
            output << results_str[i] << std::endl;
    }

    cudaFree(cuda_variant_hash);
    cudaFree(cuda_variant_repr);
    cudaFree(cuda_hashtable);
    cudaFree(cuda_results);
    index.close();
    output.close();
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
