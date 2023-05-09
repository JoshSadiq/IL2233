#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <float.h>

void genTimeSeriesSamples(double **data, int observations, int samples) {
    srand(time(NULL));
    for (int i = 0; i < samples; i++) {
        data[i] = malloc(observations * sizeof(double));
        for (int j = 0; j < observations; j++) {
            // random value between -5 and 5
            data[i][j] = (double)rand() / RAND_MAX * 10.0 - 5.0;
        }
    }
}

double distance(double *p1, double *p2) {
    return fabs(*p1 - *p2);
}

double lossFunc(double** dataPoints, int numPoints) {
    // must be implemented
    return 0.0;
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

    // initialize centroids to random time series to
    // keep track of which time series that have been used
    int sampleNoInCentroids[clusters];
    for (int i = 0; i < clusters; i++)
        sampleNoInCentroids[i] = -1;
    for (int i = 0; i < clusters; i++) {
        int randomSample = -1;
        while (isAlreadyCentroid(randomSample = rand() % samples, sampleNoInCentroids, clusters));
        sampleNoInCentroids[i] = randomSample;
        for (int j = 0; j < observations; j++) // copy whole time series
            centroids[i][j] = data[randomSample][j];
    }
    
    // random cluster assignments
    for (int i = 0; i < samples; i++)
        clusterAssignments[i] = rand() % clusters;
    
    double lastLossResult = DBL_MAX;
    while (1) {
        // re-do
        for (int c = 0; c < clusters; c++) {
            int numPointsInCluster = 0;
            double sumX = 0.0, sumY = 0.0;
            
            for (int p = 0; p < observations; p++) {
                if (symbolTable[p] == c) {
                    // Assuming 2D space for centroids
                    sumX += centroids[p][0];
                    sumY += centroids[p][1];
                    numPointsInCluster++;
                }
            }
            
            if (numPointsInCluster > 0) {
                centroids[c][0] = sumX / numPointsInCluster;
                centroids[c][1] = sumY / numPointsInCluster;
            }
        }
        
        // Update symbolTable with the best cluster for each data point
        for (int p = 0; p < samples; p++) {
            int bestCluster = clusterAssignments[p];
            double currentlyBestDist = distance(centroids[bestCluster][0], centroids[bestCluster][1],
                                                centroids[p][0], centroids[p][1]);
            
            for (int c = 0; c < numClusters; c++) {
                double distToCluster = distance(centroids[c][0], centroids[c][1],
                                                centroids[p][0], centroids[p][1]);
                if (distToCluster < currentlyBestDist) {
                    bestCluster = c;
                    currentlyBestDist = distToCluster;
                }
            }
            
            symbolTable[p] = bestCluster;
        }
        
        for (int c = 0; c < numClusters; c++) {
            int* pointsInCluster = malloc(numDataPoints * sizeof(int));
            int numPointsInCluster = 0;
            
            for (int p = 0; p < numDataPoints; p++) {
                if (symbolTable[p] == c) {
                    pointsInCluster[numPointsInCluster] = p;
                    numPointsInCluster++;
                }
            }
            
            double newLossResult = lossFunc(pointsInCluster, numPointsInCluster);
            free(pointsInCluster);
            
            if (fabs(newLossResult - lastLossResult) < 0.0001) {
                break;
            }
            
            lastLossResult = newLossResult;
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <numSamples> <numDataPoints> <numClusters>\n", argv[0]);
        return 1;
    }

    int samples = atoi(argv[1]);
    int observations = atoi(argv[2]);
    int clusters = atoi(argv[3]);

    // first level: samples, second dim: number of data points
    double **data = malloc(observations * sizeof(double*));
    genTimeSeriesSamples(dataPoints, numDataPoints);
    
    // Print the generated data points
    printf("Data Points:\n");
    for (int i = 0; i < numDataPoints; i++) {
        printf("Data Point %d: ", i);
        for (int j = 0; j < observations; j++) {
            printf("%f ", dataPoints[i][j]);
        }
        printf("\n");
    }

    exit(0);
    int max_iter = 100;
    // void kmeans(double *data, int *clusterAssignments, int clusters, int max_iter, int samples, int observations) {
    kmeans(data, clusterMapping, max_iter, clusters, observations, );
    
    printf("Symbol Table:\n");
    for (int i = 0; i < ; i++) {
        printf("Time series %d -> Cluster %d\n", i, clusterMapping[i]);
    }
    
    for (int i = 0; i < numClusters; i++)
        free(centroids[i]);
    for (int i = 0; i < numDataPoints; i++)
        free(dataPoints[i]);
    free(dataPoints);
    free(centroids);
    free(clusterMapping);
    
    return 0;
}
