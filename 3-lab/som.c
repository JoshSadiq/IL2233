#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "data.h"

typedef struct {
    int x;
    int y;
} pos_t;

double euclidean_distance(double *v1, double *v2, int n) {
    double dist = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = v1[i] - v2[i];
        dist += diff * diff;
    }
    return sqrt(dist);
}

void update_weights(double *weights, double *input, int m, int x, int y, float lr, float sigma, double dist_factor) {
    for (int i = 0; i < m; i++) {
        double delta = lr * dist_factor * (input[i] - weights[y * m * x + x * m + i]);
        weights[y * m * x + x * m + i] += delta;
    }
}

void SOM(pos_t *assignment, double *data, int n_samples, int m_features, int height, int width, int max_iter, float lr, float sigma) {
    double *weights = malloc(height * width * m_features * sizeof(double));
    for (int i = 0; i < height * width * m_features; i++) {
        weights[i] = ((double) rand() / (RAND_MAX));
    }

    for (int iter = 0; iter < max_iter; iter++) {
        float current_sigma = sigma * exp(-(float)iter / (float)max_iter);

        float current_lr = lr * exp(-(float)iter / (float)max_iter);

        for (int s = 0; s < n_samples; s++) {
            double *sample = &data[s * m_features];
            double min_distance = INFINITY;
            pos_t min_pos;

            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    double *weight = &weights[y * m_features * width + x * m_features];
                    double distance = euclidean_distance(sample, weight, m_features);

                    if (distance < min_distance) {
                        min_distance = distance;
                        min_pos.x = x;
                        min_pos.y = y;
                    }
                }
            }

            double sigma_squared = current_sigma * current_sigma;
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    double dist_x = (double) (x - min_pos.x);
                    double dist_y = (double) (y - min_pos.y);
                    double distance = dist_x * dist_x + dist_y * dist_y;
                    double dist_factor = exp(-(distance / (2.0 * sigma_squared)));

                    update_weights(weights, sample, m_features, x, y, current_lr, current_sigma, dist_factor);
                }
            }
            assignment[s] = min_pos;
        }
    }

    free(weights);
}



int main(int argc, char *argv[]) {
    int width, height;
    int n_samples, m_features;
    int max_iter = 100;
    float lr = 0.5;
    float sigma = 5.0;

    if (argc != 2) {
        printf("Usage: %s rand/iris/bme\n", argv[0]);
        return 1;
    }

    char *dataset = argv[1];

    // seed the PRNG
    srand(time(NULL));

    // first level array: samples, 
    // second level array observations in each sample
    double *data;
    if (strcmp(dataset, "rand") == 0) {
        n_samples = 100;
        m_features = 5;
        height = width = 2;
        data = malloc(n_samples * m_features * sizeof(double));
        for (int i = 0; i < n_samples * m_features; i++)
            data[i] = ((double)rand() / (RAND_MAX)) * 10.0;
    } else if (strcmp(dataset, "iris") == 0) {
        data = &iris_data[0][0]; // data is 150 x 4
        m_features = 4;
        height = 1;
        width = 3;
        n_samples = 150;
    } else if (strcmp(dataset, "bme") == 0) {
        n_samples = 180;
        m_features = 128;
        height = 1;
        width = 3;
        data = &bme_data[0][0]; // data is 180 x 128
    } else {
        printf("invalid dataset! (rand/iris/bme)\n");
        exit(1);
    }

    pos_t *assignments = malloc(n_samples * sizeof(pos_t));
    clock_t start = clock();
    SOM(assignments, data, n_samples, m_features, height, width, max_iter, lr, sigma);
    clock_t end = clock();

    for (int s = 0; s < n_samples; s++) {
        printf("(%d, %d) ", assignments[s].x, assignments[s].y);
        if (s % 10 == 0 && s != 0) {
            printf("\n");
        }
    }

    printf("\n");
    if (strcmp(dataset, "rand") == 0)
        free(data);
    free(assignments);
    printf("\nTime elapsed: %.7f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    return 0;

}


