#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <float.h>

#define THRESHOLD 0.0001

void genTimeSeriesSamples(double *data, int samples, int observations) {
    srand(time(NULL));
    for (int sample = 0; sample < samples; sample++)
        for (int observation = 0; observation < observations; observation++)
            data[sample*samples+observation] = (double)rand() / RAND_MAX * 10.0 - 5.0;
}

void displayTimeSeries(double *data, int samples, int observations) {
    for (int sample = 0; sample < samples; sample++) {
        printf("Sample %d: ", sample+1);
        for (int observation = 0; observation < observations; observation++) {
            if (observation % 8 == 0 && observation != 0)
                printf("\n");
            printf("%f ", data[sample*samples+observation]);
        }
        printf("\n");
    }
}

void displayCentroids(int clusters, int observations, double centroids[][observations]) {
    for (int cluster = 0; cluster < clusters; cluster++) {
        printf("Centroid %d: ", cluster+1);
        for (int observation = 0; observation < observations; observation++) {
            if (observation % 8 == 0 && observation != 0)
                printf("\n");
            printf("%f ", centroids[cluster][observation]);
        }
        printf("\n");
    }
}

// loss func: euclidean diff between observations in series and centroid
double distance(double *series, double *centroid, int observations) {
    double dist = 0.0;
    for (int o = 0; o < observations; o++)
        dist += fabs(series[o] - centroid[o]);
    return dist;
}

bool isAlreadyCentroid(int sampleNo, int *sampleNoInCluster, int clusters) {
    for (int i = 0; i < clusters; i++)
        if (sampleNoInCluster[i] == sampleNo)
            return true;
    return false;
}

// ORIGINAL: void kmeans(int *assignment, int K, int max_iter, int n_samples, int m_features, double *data)
// assignment -> clusterAssignments (symbol table)
// K -> clusters
// n_samples -> samples
// m_features -> observations
void kmeans(double *data, int *clusterAssignments, int max_iter, int clusters, int samples, int observations) {
    double centroids[clusters][observations];
    double previousCost = 0.0;

    // to keep track of which samples are already centroids
    int sampleNoInCentroids[clusters];
    for (int i = 0; i < clusters; i++)
        sampleNoInCentroids[i] = -1;

    // randomly select initial centroids
    for (int i = 0; i < clusters; i++) {
        int randomSample = -1;
        while (isAlreadyCentroid(randomSample = rand() % samples, sampleNoInCentroids, clusters));
        sampleNoInCentroids[i] = randomSample;
        memcpy(centroids[i], &data[randomSample], observations * sizeof(double));
    }

    // displayCentroids(clusters, observations, centroids);
    
    // assign samples to random clusters
    for (int s = 0; s < samples; s++)
        clusterAssignments[s] = rand() % clusters;

    // display initial cluster assignments
    // for (int s = 0; s < samples; s++) {
    //     printf("Initial mapping: TS %d -> C %d\n", s+1, clusterAssignments[s]+1);
    // }
    
    for (int iter = 0; iter < max_iter; iter++) {
        // update each cluster's centroid
        for (int c = 0; c < clusters; c++) {
            int samplesInCluster[samples];
            int clusterSize = 0;
            
            // find the samples belonging to the cluster
            for (int s = 0; s < samples; s++)
                if (clusterAssignments[s] == c)
                    samplesInCluster[clusterSize++] = s;

            // display found samples in cluster
            // for (int i = 0; i < clusterSize; i++)
            //     printf("Found mapping: TS %d -> C %d\n", samplesInCluster[i]+1, c+1);

            // calculate new centroids by summing observation values for each
            // sample at each time step then dividing by the number of samples
            for (int o = 0; o < observations; o++) {
                double sum = 0.0;
                for (int s = 0; s < clusterSize; s++)
                    sum += data[samplesInCluster[s] * samples + o];
                centroids[c][o] = sum / clusterSize;
            }
        }
        
        // displayCentroids(clusters, observations, centroids);

        double newCost = 0.0;
        // re-assign samples to the cluster whose centroid is closest
        for (int s = 0; s < samples; s++) {
            double minDistance = DBL_MAX;
            int closestCluster = clusterAssignments[s];
            // printf("Sample %d closest to %d\n", s+1, closestCluster);
            
            for (int c = 0; c < clusters; c++) {
                double currentDistance = 
                    distance(&data[s * samples], centroids[c], observations);
                if (currentDistance < minDistance) {
                    minDistance = currentDistance;
                    closestCluster = c;
                }
            }
            // printf("Sample %d closest to %d\n", s+1, closestCluster);
            newCost += minDistance;
            clusterAssignments[s] = closestCluster;
        }

        printf("Cost: %f\n", newCost);
        if (fabs(newCost - previousCost) == 0) {
            printf("Converged after %d iterations\n", iter+1);
            break;
        } else {
            previousCost = newCost;
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <samples> <observations> <clusters>\n", argv[0]);
        return 1;
    }

    int samples = atoi(argv[1]);
    int observations = atoi(argv[2]);
    int clusters = atoi(argv[3]);

    // first level array: samples, 
    // second level array observations in each sample
    double *data = malloc(samples * observations * sizeof(double*));
    genTimeSeriesSamples(data, samples, observations);

    // display initial data
    // displayTimeSeries(data, samples, observations);

    int *clusterAssignments = malloc(samples * sizeof(int));

    int max_iter = 100;
    kmeans(data, clusterAssignments, max_iter, clusters, samples, observations);
    
    free(data);

    printf("\nSymbol Table:\n");
    for (int i = 0; i < samples; i++) {
        printf("Time series %d -> Cluster %d\n", i+1, clusterAssignments[i]+1);
    }

    free(clusterAssignments);
    
    return 0;
}
