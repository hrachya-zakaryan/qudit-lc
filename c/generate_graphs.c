#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <string.h>
#include <stdbool.h>

typedef uint64_t PackedGraph;

void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void reverse(int *arr, int start, int end) {
    while (start < end) {
        swap(&arr[start], &arr[end]);
        start++;
        end--;
    }
}

bool nextPermutation(int *arr, int n) {
    int pivot = -1;
    for (int i = n - 2; i >= 0; i--) {
        if (arr[i] < arr[i + 1]) {
            pivot = i;
            break;
        }
    }

    if (pivot == -1) {
        reverse(arr, 0, n - 1);
        return false;
    }

    for (int i = n - 1; i > pivot; i--) {
        if (arr[i] > arr[pivot]) {
            swap(&arr[i], &arr[pivot]);
            break;
        }
    }

    reverse(arr, pivot + 1, n - 1);
    return true;
}

PackedGraph bitpack_encode(int n, int d, int *matrix) {
    PackedGraph packed = 0;
    int bit_length = 1;
    while ((1 << bit_length) < d) bit_length++;
    int shift = 0;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            packed |= ((PackedGraph)matrix[i * n + j] << shift);
            shift += bit_length;
        }
    }
    return packed;
}

void bitpack_decode(PackedGraph packed, int n, int d, int *matrix) {
    int bit_length = 1;
    while ((1 << bit_length) < d) bit_length++;
    int shift = 0;

    memset(matrix, 0, n * n * sizeof(int));

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            int weight = (packed >> shift) & ((1 << bit_length) - 1);
            matrix[i * n + j] = matrix[j * n + i] = weight;
            shift += bit_length;
        }
    }
}

void local_scaling(int n, int v, int k, int d, int *matrix) {
    for (int u = 0; u < n; u++) {
        if (u != v) {
            matrix[v * n + u] = matrix[u * n + v] = (k * matrix[v * n + u]) % d;
        }
    }
}

void local_complementation(int n, int v, int k, int d, int *matrix) {
    for (int u = 0; u < n; u++) {
        for (int w = u; w < n; w++) {
            if (u != w && u != v && w != v && matrix[v * n + u] && matrix[v * n + w]) {
                matrix[u * n + w] = matrix[w * n + u] = (matrix[u * n + w] + k * matrix[v * n + u] * matrix[v * n + w]) % d;
            }
        }
    }
}

PackedGraph find_minimum_encoding(int n, int d, int *matrix) {
    PackedGraph min_encoded = bitpack_encode(n, d, matrix);
    int *perm = malloc(n * sizeof(int));
    if (!perm) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++) perm[i] = i;

    do {
        PackedGraph encoded = 0;
        int bit_length = 1;
        while ((1 << bit_length) < d) bit_length++;
        int shift = 0;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                int weight = matrix[perm[i] * n + perm[j]];
                encoded |= ((PackedGraph)weight << shift);
                shift += bit_length;
            }
        }

        if (encoded < min_encoded) {
            min_encoded = encoded;
        }
    } while (nextPermutation(perm, n));

    free(perm);
    return min_encoded;
}

void generate_weighted(int n, int d, PackedGraph encoded) {
    int num_edges = 0;
    int *adj_matrix = calloc(n * n, sizeof(int));
    if (!adj_matrix) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    bitpack_decode(encoded, n, d, adj_matrix);

    int max_edges = n * (n - 1) / 2;
    int (*upper_indices)[2] = malloc(max_edges * sizeof(int[2]));
    if (!upper_indices) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (adj_matrix[i * n + j] == 1) {
                if (num_edges < max_edges) {
                    upper_indices[num_edges][0] = i;
                    upper_indices[num_edges][1] = j;
                    num_edges++;
                }
            }
        }
    }
    

    long long total_matrices = 1;
    for (int i = 0; i < num_edges; i++) {
        total_matrices *= (d - 1);
    }

    int seen_capacity = 25000;
    int seen_size = 0;
    PackedGraph *seen = malloc(seen_capacity * sizeof(PackedGraph));
    if (!seen) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (long long mask = 0; mask < total_matrices; mask++) {
        int *new_matrix = malloc(n * n * sizeof(int));
        if (!new_matrix) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }
        memcpy(new_matrix, adj_matrix, n * n * sizeof(int));

        long long temp_mask = mask;
        for (int k = 0; k < num_edges; k++) {
            int weight = (temp_mask % (d - 1)) + 1;
            temp_mask /= (d - 1);
            int i = upper_indices[k][0], j = upper_indices[k][1];
            new_matrix[i * n + j] = new_matrix[j * n + i] = weight;
        }

        PackedGraph min_encoded_value = find_minimum_encoding(n, d, new_matrix);

        bool is_unique = true;
        for (int i = 0; i < seen_size; i++) {
            if (seen[i] == min_encoded_value) {
                is_unique = false;
                break;
            }
        }

        if (is_unique) {
            printf("%llu\n", (unsigned long long)min_encoded_value);
            if (seen_size == seen_capacity) {
                seen_capacity *= 2;
                seen = realloc(seen, seen_capacity * sizeof(PackedGraph));
                if (!seen) {
                    fprintf(stderr, "Memory allocation failed\n");
                    exit(EXIT_FAILURE);
                }
            }
            seen[seen_size++] = min_encoded_value;
        }

        free(new_matrix);
    }

    free(upper_indices);
    free(seen);
    free(adj_matrix);
}

int main(int argc, char **argv) {
    PackedGraph encoded;
    int n = atoi(argv[1]);
    int d = atoi(argv[2]);
    encoded=strtoull(argv[3],NULL,10);
    generate_weighted(n, d, encoded);
    return 0;
}
