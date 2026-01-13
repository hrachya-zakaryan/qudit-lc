#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <string.h>
#include <stdbool.h>


typedef uint64_t PackedGraph;  // Use 64-bit integer for storage



void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void reverse(int arr[], int start, int end) {
    while (start < end) {
        swap(&arr[start], &arr[end]);
        start++;
        end--;
    }
}

bool nextPermutation(int *arr, int n) {
  
    // Find the pivot index
    int pivot = -1;
    for (int i = n - 2; i >= 0; i--) {
        if (arr[i] < arr[i + 1]) {
            pivot = i;
            break;
        }
    }

    // If pivot point does not exist, 
    // reverse the whole array
    if (pivot == -1) {
        reverse(arr, 0, n - 1);
        return false;
    }

    // Find the element from the right that
    // is greater than pivot
    for (int i = n - 1; i > pivot; i--) {
        if (arr[i] > arr[pivot]) {
            swap(&arr[i], &arr[pivot]);
            break;
        }
    }

    // Reverse the elements from pivot + 1 to the end
    reverse(arr, pivot + 1, n - 1);
    return true;
}

// Bitpack encoding of an adjacency matrix
PackedGraph bitpack_encode(int n, int d,int matrix[n][n]) {
    PackedGraph packed = 0;
    int bit_length = 1;
    while ((1 << bit_length) < d) bit_length++;
    int shift = 0;
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            packed |= ((PackedGraph)matrix[i][j] << shift);
            shift += bit_length;
        }
    }
    return packed;
}

// Bitpack decoding
void bitpack_decode(PackedGraph packed, int n, int d, int matrix[n][n]) {
    int bit_length = 1;
    while ((1 << bit_length) < d) bit_length++;
    int shift = 0;
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            int weight = (packed >> shift) & ((1 << bit_length) - 1);
            matrix[i][j] = matrix[j][i] = weight;
            shift += bit_length;
        }
    }
}

// Apply local scaling to vertex v
void local_scaling( int n, int v, int k, int d,int matrix[n][n]) {
    for (int u = 0; u < n; u++) {
        if (u != v) {
            matrix[v][u] = matrix[u][v] = (k * matrix[v][u]) % d;
        }
    }
}

// Apply local complementation to vertex v
void local_complementation(int n, int v, int k, int d, int matrix[n][n]) {
    for (int u = 0; u < n; u++) {
        for (int w = u; w < n; w++) {
            if (u != w && u != v && w != v && matrix[v][u] && matrix[v][w]) {
                matrix[u][w] = matrix[w][u] = (matrix[u][w] + k * matrix[v][u] * matrix[v][w]) % d;
            }
        }
    }
}

// Generate all permutations and find the minimal encoding
PackedGraph find_minimum_encoding(int n, int d, int matrix[n][n]) {
    PackedGraph min_encoded = bitpack_encode(n, d, matrix);
    int perm[n];
    for (int i = 0; i < n; i++) perm[i] = i;
    
    do {
        PackedGraph encoded = 0;
        int bit_length = 1;
        while ((1 << bit_length) < d) bit_length++;
        int shift = 0;

        // Encode only the upper triangular portion with permuted indices
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                int weight = matrix[perm[i]][perm[j]];
                encoded |= ((PackedGraph)weight << shift);
                shift += bit_length;
            }
        }

        if (encoded < min_encoded) {
            min_encoded = encoded;
        }
    } while (nextPermutation(perm, n));

    return min_encoded;
}

// Create a complementation layer and return the resultant graphs
void create_complementation_layer(PackedGraph current, int n, int d,PackedGraph results[]) {
    int matrix[n][n];
    int result_count = 0;
    bitpack_decode(current, n, d, matrix);
    
    for (int v = 0; v < n; v++) {
        for (int k = 2; k < d; k++) {  // Scaling
            int scaled[n][n];
            memcpy(scaled, matrix, sizeof(scaled));
            local_scaling(n, v, k, d,scaled);
            PackedGraph encoded_scaled = find_minimum_encoding(n, d, scaled);
            results[result_count] = encoded_scaled;
            result_count++;
        }
        for (int k = 1; k < d; k++) {  // Complementation
            int complemented[n][n];
            memcpy(complemented, matrix, sizeof(complemented));
            local_complementation(n, v, k, d, complemented);
            PackedGraph encoded_complemented = find_minimum_encoding(n, d,complemented);
            results[result_count] = encoded_complemented;
            result_count++;
        }
    }
}

int main(int argc, char **argv) {
    int n;
    int d;
    PackedGraph encoded;
    n = atoi(argv[1]);
    d = atoi(argv[2]);
    encoded=strtoull(argv[3],NULL,10);
    int complementations = n*(2*d-3);
    PackedGraph results[complementations];
    create_complementation_layer(encoded, n, d, results);
    for (int i = 0; i < complementations; i++) {
        printf("%llu\n", (unsigned long long)results[i]);
    }
    
    return 0;
}
