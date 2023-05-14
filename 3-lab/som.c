#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int x;
    int y;
} t_pos;

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

void SOM(t_pos *assignment, double *data, int n_samples, int m_features, int height, int width, int max_iter, float lr, float sigma) {
    double *weights = malloc(height * width * m_features * sizeof(double));
    for (int i = 0; i < height * width * m_features; i++) {
        weights[i] = ((double) rand() / (RAND_MAX));
    }

    for (int iter = 0; iter < max_iter; iter++) {
        float current_sigma = sigma * exp(-(float) iter / (float) max_iter);

        float current_lr = lr * exp(-(float) iter / (float) max_iter);

        for (int s = 0; s < n_samples; s++) {
            double *sample = &data[s * m_features];
            double min_distance = INFINITY;
            t_pos min_pos;

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


void compute_umatrix(double *map, double *umatrix, int height, int width) {
    for (int i = 0; i < height * width; i++) {
        double sum = 0;
        int count = 0;
        if (i - 1 >= 0) { // left neighbor
            sum += pow(map[i] - map[i-1], 2);
            count++;
        }
        if (i + 1 < height * width) { // right neighbor
            sum += pow(map[i] - map[i+1], 2);
            count++;
        }
        if (i - width >= 0) { // top neighbor
            sum += pow(map[i] - map[i-width], 2);
            count++;
        }
        if (i + width < height * width) { // bottom neighbor
            sum += pow(map[i] - map[i+width], 2);
            count++;
        }
        umatrix[i] = sqrt(sum / count);
    }
}



int main() {
    // Example usage
    int n_samples = 5;
    int m_features = 3;
    int height = 4;
    int width = 4;
    int max_iter = 100;
    float lr = 0.1;
    float sigma = 5.0;
    double *data = malloc(n_samples * m_features * sizeof(double));
    for (int i = 0; i < n_samples * m_features; i++) {
        data[i] = ((double) rand() /    (RAND_MAX)) * 10.0;
}
t_pos *assignment = malloc(n_samples * sizeof(t_pos));
SOM(assignment, data, n_samples, m_features, height, width, max_iter, lr, sigma);


for (int s = 0; s < n_samples; s++) {
    printf("(%d, %d) ", assignment[s].x, assignment[s].y);
}
    printf("\n");
    free(data);
    free(assignment);
    return 0;
}


