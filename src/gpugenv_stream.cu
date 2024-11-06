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

    for (size_t i = start_idx + threadIdx.x; i < end_idx; i += blockDim.x) {
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
    hash_t* &results, std::vector<std::string> &results_str) {
    int variant_size = results_str.size();
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
    // Use page-locked memory.
    size_t variant_bytes = variant_size * sizeof(hash_t);
    size_t hashtable_bytes = Consts::MAX_HASHTABLE_SIZE * sizeof(hash_t);
    hash_t *variant_hash, *variant_repr;
    cudaHostAlloc((void**)&variant_hash, variant_bytes, cudaHostAllocDefault);
    cudaHostAlloc((void**)&variant_repr, variant_bytes, cudaHostAllocDefault);
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
    hash_t *cuda_variant_hash, *cuda_variant_repr;
    hash_t *cuda_hashtable;
    hash_t *cuda_results;
    cudaMalloc((void**)&cuda_variant_hash, variant_bytes);
    cudaMalloc((void**)&cuda_variant_repr, variant_bytes);
    cudaMalloc((void**)&cuda_hashtable, hashtable_bytes);
    cudaMalloc((void**)&cuda_results, variant_bytes);
    hash_t *results;
    cudaHostAlloc((void**)&results, variant_bytes, cudaHostAllocDefault);
    std::vector<std::string> results_str(variant_size);

    cudaStream_t stream;
    cudaStreamCreate(&stream);
    cudaMemcpyAsync(
        cuda_variant_hash, variant_hash, variant_bytes, cudaMemcpyHostToDevice, stream
    );
    cudaMemcpyAsync(
        cuda_variant_repr, variant_repr, variant_bytes, cudaMemcpyHostToDevice, stream
    );
    cudaStreamSynchronize(stream);

    hash_t hash, repr;
    std::string row;
    int processed = 0;
    // Local hashtable contains pairs (hash, line).
    // Other hashtable contains pairs (hash, representation).
    // We keep two copies of hashtable in order to concurrently read from file and process data.
    std::vector<std::vector<std::string>> hashtable(
        2, std::vector<std::string>(Consts::MAX_HASHTABLE_SIZE)
    );
    hash_t *hashtable_repr;
    cudaHostAlloc((void**)&hashtable_repr, hashtable_bytes, cudaHostAllocDefault);

    // `should_save` determines if there are some results from previous iteration to save.
    std::function<void(bool)> execute_kernel = [&](bool should_save) {
        // Synchronize stream, copy the hashtable.
        cudaStreamSynchronize(stream);
        cudaMemcpyAsync(
            cuda_hashtable, hashtable_repr, hashtable_bytes, cudaMemcpyHostToDevice, stream
        );
        // Save previous results.
        if (should_save)
            save_results(hashtable[1], results, results_str);
        for (int i = 0; i < variant_size; i++)
            results[i] = Consts::INF;
        cudaMemcpyAsync(
            cuda_results, results, variant_bytes, cudaMemcpyHostToDevice, stream
        );

        // Execute kernel, save the results.
        iterate_over_variant<<<Consts::BLOCKS, Consts::THREADS_PER_BLOCK, 0, stream>>>(
            cuda_hashtable, variant_size,
            cuda_variant_hash, cuda_variant_repr, cuda_results
        );
        cudaMemcpyAsync(results, cuda_results, variant_bytes, cudaMemcpyDeviceToHost, stream);
    };

    bool should_save = false;
    while (std::getline(index, line)) {
        std::istringstream iss(line);
        iss >> hash >> repr;
        std::getline(iss >> std::ws, row);
        // Scale the hash.
        hash *= Consts::MAX_COLLISIONS;
        while (hashtable[0][hash] != "") {
            hash++;
        }
        hashtable[0][hash] = row;
        hashtable_repr[hash] = repr;

        processed++;
        if (processed == Consts::CHUNK_SIZE) {
            execute_kernel(should_save);
            should_save = true;

            // Reset the arrays.
            processed = 0;
            std::swap(hashtable[0], hashtable[1]);
            std::fill(hashtable[0].begin(), hashtable[0].end(), "");
            for (std::size_t i = 0; i < Consts::MAX_HASHTABLE_SIZE; i++)
                hashtable_repr[i] = Consts::INF;
        }
    }

    if (processed != 0) {
        execute_kernel(should_save);
        cudaStreamSynchronize(stream);
        save_results(hashtable[0], results, results_str);
    }

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
    cudaFreeHost(variant_hash);
    cudaFreeHost(variant_repr);
    cudaFreeHost(results);
    cudaFreeHost(hashtable_repr);
    cudaStreamDestroy(stream);
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
