
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <float.h>
#include "data.h"

/* *********************************************************************
*********************** HELPER FUNCTIONS *******************************
********************************************************************* */

void displayTimeSeries(double *data, int samples, int observations) {
    for (int s = 0; s < samples; s++) {
        printf("Sample %d: ", s+1);
        for (int o = 0; o < observations; o++) {
            if (o % 8 == 0 && o != 0)
                printf("\n");
            printf("%.2f ", data[s * observations + o]);
        }
        printf("\n");
    }
}

void displayCentroids(int clusters, int observations, double centroids[][observations]) {
    for (int c = 0; c < clusters; c++) {
        printf("Centroid %d: ", c+1);
        for (int o = 0; o < observations; o++) {
            if (o % 8 == 0 && o != 0)
                printf("\n");
            printf("%.2f ", centroids[c][o]);
        }
        printf("\n");
    }
}


// helper function to avoid duplicate initial centroids
bool isAlreadyCentroid(int sampleNo, int *sampleNoInCluster, int clusters) {
    for (int i = 0; i < clusters; i++)
        if (sampleNoInCluster[i] == sampleNo)
            return true;
    return false;
}

// function to generate random samples
void generateSamples(double *data, int samples, int observations) {
    for (int s = 0; s < samples; s++)
        for (int o = 0; o < observations; o++)
            data[s * observations + o] = (double)rand() / RAND_MAX * 10.0 - 5.0;
}

/* *********************************************************************
******************************* K-MEANS *******************************
********************************************************************* */

// "loss func": euclidean diff between observations in series and centroid
double distance(double *series, double *centroid, int observations) {
    double dist = 0.0;
    for (int o = 0; o < observations; o++)
        dist += fabs(series[o] - centroid[o]);
    return dist;
}

void initCentroids(double *centroids, double *data, int clusters, int samples, int observations) {
    // to avoid duplicate initial centroids
    int randomSampleNo;
    int sampleNoInCentroids[clusters];
    memset(sampleNoInCentroids, -1, clusters * sizeof(int));

    // random samples become centroids
    for (int c = 0; c < clusters; c++) {
        while (isAlreadyCentroid(randomSampleNo = rand() % samples, sampleNoInCentroids, clusters));
        sampleNoInCentroids[c] = randomSampleNo;
        memcpy(
            &centroids[c * observations], 
            &data[randomSampleNo * observations], 
            observations * sizeof(double));
    }
}

void updateCentroids(double *centroids, double *data, int *clusterAssignments, int clusters, int samples, int observations) {
    for (int c = 0; c < clusters; c++) {
        int samplesInCluster[samples]; // at most all samples
        int clusterSize = 0;
        
        // find which samples assigned to the cluster
        for (int s = 0; s < samples; s++)
            if (clusterAssignments[s] == c)
                samplesInCluster[clusterSize++] = s;

        // for (int i = 0; i < clusterSize; i++)
        //     printf("Found mapping: TS %d -> C %d\n", samplesInCluster[i]+1, c+1);

        // average of observations at each time step
        if (clusterSize > 0)
            for (int o = 0; o < observations; o++) {
                double sum = 0.0;
                for (int s = 0; s < clusterSize; s++)
                    sum += data[samplesInCluster[s] * observations + o];
                centroids[c * observations + o] = sum / clusterSize;
            }
    }
}

double updateAssignments(double *centroids, double *data, int *clusterAssignments, int clusters, int samples, int observations) {
    double cost = 0.0;
    // re-assign samples to the cluster whose centroid is closest
    for (int s = 0; s < samples; s++) {
        double minDistance = DBL_MAX;
        // printf("Sample %d closest to %d\n", s+1,  clusterAssignments[s]);
        
        for (int c = 0; c < clusters; c++) {
            double currentDistance = 
                distance(&data[s*observations], &centroids[c*observations], observations);
            if (currentDistance < minDistance) {
                minDistance = currentDistance;
                clusterAssignments[s] = c;
            }
        }
        // printf("Sample %d closest to %d\n", s+1, closestCluster);
        cost += minDistance;
    }

    return cost;
}

// ORIGINAL: void kmeans(int *assignment, int K, int max_iter, int n_samples, int m_features, double *data)
// assignment -> clusterAssignments (symbol table)
// K -> clusters
// n_samples -> samples
// m_features -> observations
void kmeans(double *data, int *clusterAssignments, int max_iter, int clusters, int samples, int observations) {
    double *centroids = malloc(clusters * observations * sizeof(double));
    double newCost, previousCost = 0.0;

    initCentroids(centroids, data, clusters, samples, observations);
    displayCentroids(clusters, observations, centroids);
    
    // assign samples to random clusters
    for (int s = 0; s < samples; s++)
        clusterAssignments[s] = rand() % clusters;

    // display initial cluster assignments
    // for (int s = 0; s < samples; s++)
    //     printf("Initial mapping: TS %d -> C %d\n", s+1, clusterAssignments[s]+1);

    for (int iter = 0; iter < max_iter; iter++) {
        newCost = 
            updateAssignments(centroids, data, clusterAssignments, clusters, samples, observations);

        // printf("Cost: %f\n", newCost);
        if (fabs(newCost - previousCost) < 1.0) {
            printf("Converged after %d iterations\n", iter+1);
            break;
        } else {
            previousCost = newCost;
        }

        updateCentroids(centroids, data, clusterAssignments, clusters, samples, observations);
        // displayCentroids(clusters, observations, centroids);
    }

    displayCentroids(clusters, observations, centroids);

    free(centroids);
}

bool testImpl() {
    double data[3][6] = {
        {1.0, 2.0, 1.0, 2.0, 1.0, 2.0},
        {8.0, 9.0, 8.0, 9.0, 8.0, 9.0},
        {15.0, 14.0, 15.0, 14.0, 15.0, 14.0}};
    int samples = 3;
    int observations = 6;
    int clusters = 2;
    double centroids[clusters][observations];
    int clusterAssignments[samples];

    memcpy(centroids[0], data[0], observations * sizeof(double));
    memcpy(centroids[1], data[2], observations * sizeof(double));

    printf("Initial data:\n");
    displayTimeSeries(&data[0][0], samples, observations);

    printf("\nInitial centroids:\n");
    displayCentroids(clusters, observations, centroids);

    // "random assignments"
    clusterAssignments[0] = 1; // sample 1 -> cluster 2 (sample 3)
    clusterAssignments[1] = 1; // sample 2 -> cluster 2 (sample 3)
    clusterAssignments[2] = 0; // sample 3 -> cluster 1 (sample 1)

    printf("\nInitial 'random' mapping:\n");
    for (int i = 0; i < samples; i++)
        printf("time series %d -> cluster %d\n", i+1, clusterAssignments[i]+1);

    // test cost function
    double expectedCost = 7 + 7 + 7 + 7 + 7 + 7;
    double actualCost = distance(data[1], centroids[0], observations);
    if (expectedCost != actualCost) {
        printf("Cost(distance) between sample 2 and cluster 1 is %f, expected %f\n", actualCost, expectedCost);
        return false;
    } else {
        printf("\n- Cost is working correctly\n");
        printf("Cost of sample 2 to cluster 1 is %.2f\n", actualCost);
    }

    updateAssignments(&centroids[0][0], &data[0][0], clusterAssignments, clusters, samples, observations);
    if (clusterAssignments[0] != 0 || clusterAssignments[1] != 1 || clusterAssignments[2] != 1) {
        return false;
    } else {
        printf("\n- Cluster assignments are working correctly, updated:\n");
        for (int i = 0; i < samples; i++)
            printf("time series %d -> cluster %d\n", i+1, clusterAssignments[i]+1);
    }


    updateCentroids(&centroids[0][0], &data[0][0], clusterAssignments, clusters, samples, observations);
    double expectedCentroid[6] = {11.5, 11.5, 11.5, 11.5, 11.5, 11.5};
    for (int i = 0; i < observations; i++)
        if (expectedCentroid[i] != centroids[1][i]) {
            printf("centroid 1, feature %d: %.2f, expected %.2f\n", i+1, centroids[0][i], expectedCentroid[i]);
            return false;
        }
    printf("\n- Centroids are updated correctly, now they are:\n");
    displayCentroids(clusters, observations, centroids);

    printf("\n");

    return true;
}

int main(int argc, char *argv[]) {
    // if (!testImpl())
    //     return 1;

    if (argc != 3 && argc != 5) {
        printf("Usage: %s rand/iris/bme clusters [samples] [observations]\n", argv[0]);
        return 1;
    }

    srand(time(NULL));

    double *data;
    int samples, observations;
    char *dataset = argv[1];
    int clusters = atoi(argv[2]);

    // first level array: samples, 
    // second level array observations in each sample
    if (strcmp(dataset, "rand") == 0) {
        samples = atoi(argv[3]);
        observations = atoi(argv[4]);
        data = malloc(samples * observations * sizeof(double));
        generateSamples(data, samples, observations);
    } else if (strcmp(dataset, "iris") == 0) {
        samples = sizeof(iris_data) / sizeof(iris_data[0]);
        observations = sizeof(iris_data[0]) / sizeof(iris_data[0][0]);
        data = &iris_data[0][0];
    } else if (strcmp(dataset, "bme") == 0) {
        samples = sizeof(bme_data) / sizeof(bme_data[0]);
        observations = sizeof(bme_data[0]) / sizeof(bme_data[0][0]);
        data = &bme_data[0][0];
    } else {
        printf("invalid dataset selection '%s'! (rand/iris/bme)\n", dataset);
        exit(1);
    }

    // display initial data
    // displayTimeSeries(data, samples, observations);

    int *clusterAssignments = malloc(samples * sizeof(int));
    int max_iter = 100;

    clock_t start = clock();
    kmeans(data, clusterAssignments, max_iter, clusters, samples, observations);
    clock_t end = clock();

    if (strcmp(dataset, "rand") == 0)
        free(data);

    // print final mapping as S:# -> C:# with four entries per line and use padding for better alignment
    printf("\nFinal mapping:\n");
    for (int i = 0; i < samples; i++) {
        printf("S:%-3d -> C:%d\t", i+1, clusterAssignments[i]+1);
        if ((i+1) % 4 == 0) printf("\n");
    }

    free(clusterAssignments);
    
    printf("\nTime elapsed: %.7f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    return 0;
}