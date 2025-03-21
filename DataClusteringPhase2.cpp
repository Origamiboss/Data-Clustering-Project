//Coding Practices: https://www.geeksforgeeks.org/7-common-programming-principles-that-every-developer-must-follow/

/*
Random Partition Method Vs Random Selection Method
In Theory, the random partition method would be more consistent than its random selection counterpart. 
Since the random partition method assigns each data point to a random cluster, it gives a lot more opportunities for the clusters to not be outliers.
With the random selection method, each cluster is merely a randomly selected data point. This greatly increases the opportunity for the clusters to become outliers if a bad data point is selected.

The results of this project revealed this to be the case. The random selection method resulted in high variance results, but within 100 runs it would eventually create optimal results.
The random partition method would have less outlier cases, but its best results are often worse than the random selection's best results.

As such, if you are going to run through an entire data set a few times then the random partition method would be a good method to choose,
but if you are going to run through an entire data set many times (ex: 100) then random selection method would be a good method due to its ability to generate optimal outliers.
*/

//Phase 1 libraries
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <set>


//Phase 2 libraries
#include <random>

//Phase 4 library
#include <array>


using namespace std;



//Phase 1 (Gather the data)
static int checkTheArguments(string fileName, int maxIterations, double convergenceThreshold, int numOfRuns, int typeOfClustering);
double** readData(string fileName, int& numOfInstances, int& sizeOfInstances);
double** setClusters(double** data, int numOfClusters, int numOfInstances, int sizeOfInstance);


//Phase 2 (Run the K-Means)
double runIterations(int& maxIterations, double& convergenceThreshold, double** clusters, double** data, double& initialSSE, int& iterationsRan, int* labels, double* distancesToClosestCluster, int* clusterSizes, int numOfInstances, int sizeOfInstance, int numOfClusters);
double calculateSquaredDistance(const double* x, const double* c, int sizeOfInstance);

//Phase 3 (Normalization and Initialization)
void minMaxNorm(double** data, int numOfInstances, int sizeOfInstance);
void zScoreNorm(double** data, int numOfInstances, int sizeOfInstance);
double** randomParitionClusters(double** data, int numOfClusters, int numOfInstances, int sizeOfInstance);
double** maximumMethodClusters(double** data, int numOfClusters, int numOfInstances, int sizeOfInstance);

//Phase 4
double calinski_validity(double** data, double** clusters, double& finalSSE, int* clusterSizes, int numClusters, int numOfInstances, int sizeOfInstance);
double silhouette_width(double** data, double** clusters, int* labels, int* clusterSizes, int numClusters, int numOfInstances, int sizeOfInstance);
double dunn(double** data, double** clusters, int* labels, int numClusters, int numOfInstances, int sizeOfInstance);
double davies_bouldin(double** data, double** clusters, int* labels, int numClusters, double* distancesToClosestCluster, int* clustersSizes, int numOfInstances, int sizeOfInstance);

void runClusterForRunAndNumClusters(int runs, int numOfClusters, double** data, double** clusters,
    int maxIterations, double convergenceThreshold,
    int typeOfClustering, int numOfInstances, int sizeOfInstances);


//overAll Mean
double* overallMean;

double* bestIndexCH;
double* bestIndexSW;
double* bestIndexDU;
double* bestIndexDB;
double* minDistCluster;
double* maxDistCluster;

int main(int argc, char* argv[])
{
    string fileName;
    int minNumClusters = -1;
    int maxIterations = -1;
    double convergenceThreshold = -1;
    int numOfRuns = -1;
    int typeOfClustering = -1;
    

    // make sure there are 5 arguments
    if (argc != 6)
    {
        cout << "Usage: <F> <I> <T> <R> <V> <K>" << endl;
        cout << "F: Name of the data file" << endl;
        cout << "I: Maximum number of iterations (positive integer)" << endl;
        cout << "T: Convergence threshold (non-negative real number)" << endl;
        cout << "R: Number of runs (positive integer)" << endl;
        cout << "V: Cluster Generation (random | random partition | maximum method) (0 | 1 | 2)" << endl;
        return 1;
    }

    // collect the data from arguments
    fileName = argv[1];
    minNumClusters = 2;
    maxIterations = stoi(argv[2]);
    convergenceThreshold = stod(argv[3]);
    numOfRuns = stoi(argv[4]);
    typeOfClustering = stoi(argv[5]);
    
    //check the data
    int result = checkTheArguments(fileName, maxIterations, convergenceThreshold, numOfRuns, typeOfClustering);

    if (result == 1) {
        //error
        return 1;
    }
    // Play with the data here!
    int sizeOfInstances, numOfInstances;
    //Read the data from the file
    double** data = readData(fileName, numOfInstances, sizeOfInstances);
    //normalize the data
    minMaxNorm(data, numOfInstances, sizeOfInstances);
    
    
    // Compute overall mean
    overallMean = new double [sizeOfInstances];
    for (int i = 0; i < sizeOfInstances; i++) {
        overallMean[i] = 0.0;
    }
    for (int j = 0; j < numOfInstances; j++) {
        for (int i = 0; i < sizeOfInstances; i++) {
            overallMean[i] += data[j][i];
        }
    }
    for (int i = 0; i < sizeOfInstances; i++) {
        overallMean[i] /= numOfInstances;
    }



    int maxNumClusters = sqrt(numOfInstances / 2);
    bestIndexCH = new double [maxNumClusters];
    bestIndexSW = new double[maxNumClusters];
    bestIndexDU = new double[maxNumClusters];
    bestIndexDB = new double[maxNumClusters];
    for (int i = 0; i < maxNumClusters; i++) {
        bestIndexCH[i] = -INFINITY;
        bestIndexSW[i] = -INFINITY;
        bestIndexDU[i] = -INFINITY;
        bestIndexDB[i] = INFINITY;
    }
    double** clusters = nullptr;

    maxDistCluster = new double[maxNumClusters];
    minDistCluster = new double[maxNumClusters];
    
    for (int runs = 1; runs < numOfRuns + 1; runs++) {
        for (int numOfClusters = minNumClusters; numOfClusters < maxNumClusters + 1; numOfClusters++) {
            //reset maxDist and minDist
            for (int i = 0; i < numOfClusters; i++) {
                maxDistCluster[i] = -INFINITY;
                minDistCluster[i] = INFINITY;
            }
            //deallocate cluster memory
            if (clusters != nullptr) {
                for (int i = 0; i < numOfClusters; i++) {
                    delete[] clusters[i];  // Free each cluster's data
                }
                delete[] clusters;  // Free the array of pointers
            }
            runClusterForRunAndNumClusters(runs, numOfClusters, data, clusters,
                maxIterations, convergenceThreshold, typeOfClustering, numOfInstances, sizeOfInstances);
        }
    }
    
    
    //Header
    //cout << "File Name\tCluster Num\tCH\tSW\tDU\tDB" << endl;
    for (int i = minNumClusters; i < maxNumClusters + 1; i++) {
        cout << fileName << "\t" << i << "\t"  << bestIndexCH[i] << "\t" << bestIndexSW[i] << "\t" << bestIndexDU[i] <<"\t" <<  bestIndexDB[i] << endl;
    }
    //deallocate memory
    for (int i = 0; i < numOfInstances; i++) {
        delete[] data[i];
    }
    delete[] data;
    //deallocate cluster memory
    if (clusters != nullptr) {
        for (int i = 0; i < maxNumClusters; i++) {
            delete[] clusters[i];  // Free each cluster's data
        }
        delete[] clusters;  // Free the array of pointers
    }
    delete[] bestIndexCH;
    delete[] bestIndexSW;
    delete[] bestIndexDU;
    delete[] bestIndexDB;
    delete[] overallMean;
    delete[] minDistCluster;
    delete[] maxDistCluster;

    return 0;
}
// The function you want to run in each thread if you so choose
void runClusterForRunAndNumClusters(int runs, int numOfClusters, double** data, double** clusters,
    int maxIterations, double convergenceThreshold,
    int typeOfClustering, int numOfInstances, int sizeOfInstance){
    
    int* labels = new int[numOfInstances];
    double* distancesToClosestCluster = new double [numOfInstances];
    int* clusterSizes = new int[numOfClusters];
    
    // Obtain clusters based on the chosen method
    switch (typeOfClustering) {
    case 0:
        clusters = setClusters(data, numOfClusters, numOfInstances, sizeOfInstance);
        break;
    case 1:
        clusters = randomParitionClusters(data, numOfClusters, numOfInstances, sizeOfInstance);
        break;
    case 2:
        clusters = maximumMethodClusters(data, numOfClusters, numOfInstances, sizeOfInstance);
        break;
    }
    
    // Run the iterations
    double initialSSE;
    int iterations;
    double finalSSE = runIterations(maxIterations, convergenceThreshold, clusters, data, initialSSE, iterations, labels, distancesToClosestCluster, clusterSizes, numOfInstances, sizeOfInstance, numOfClusters);
    //gather indexes
    double index = 0;
    //Calinski
    index = calinski_validity(data, clusters, finalSSE, clusterSizes, numOfClusters, numOfInstances, sizeOfInstance);
    bestIndexCH[numOfClusters] = max(bestIndexCH[numOfClusters], index);
    //silhouette width
    index = silhouette_width(data, clusters, labels, clusterSizes, numOfClusters, numOfInstances, sizeOfInstance);
    bestIndexSW[numOfClusters] = max(bestIndexSW[numOfClusters], index);
    //Dunn
    index = dunn(data, clusters, labels, numOfClusters, numOfInstances, sizeOfInstance);
    bestIndexDU[numOfClusters] = max(bestIndexDU[numOfClusters], index);

    //Davies
    index = davies_bouldin(data, clusters, labels, numOfClusters, distancesToClosestCluster, clusterSizes, numOfInstances, sizeOfInstance);
    bestIndexDB[numOfClusters] = min(bestIndexDB[numOfClusters], index);
    //deallocate memory
    delete[] labels;
    delete[] distancesToClosestCluster;
    delete[] clusterSizes;
}

//makes sure the arguments are valid
static int checkTheArguments(string fileName, int maxIterations, double convergenceThreshold, int numOfRuns, int typeOfClustering) {
    try
    {
        // Validate the arguments
        ifstream file(fileName);
        if (!file)
        {
            throw invalid_argument(("File: " + fileName + " does not exist.").c_str());
        }
        file.close();
        if (maxIterations <= 0)
        {
            throw invalid_argument("Maximum number of iterations (I) must be positive.");
        }
        if (convergenceThreshold < 0)
        {
            throw invalid_argument("Convergence threshold (T) cannot be negative.");
        }
        if (numOfRuns <= 0)
        {
            throw invalid_argument("Number of runs (R) must be positive.");
        }
        if (typeOfClustering != 0 && typeOfClustering != 1 && typeOfClustering != 2) {
            throw invalid_argument("Type cluster generation (V) must be 0, 1, or 2");
        }
        return 0;
    }
    catch (const invalid_argument& e) {
        cout << "Error: " << e.what() << endl;
        return 1;
    }
    catch (const exception& e) {
        cout << "Unexpected error: " << e.what() << endl;
        return 1;
    }
}
//reads the data from the file
double** readData(string fileName, int& numOfInstances, int& sizeOfInstances) {
    ifstream file(fileName);
    string line;

    // Read the first line to get the number of lines and the size of each line
    getline(file, line);
    stringstream ss(line);
    ss >> numOfInstances >> sizeOfInstances;

    // Create a 2D array with numOfInstances rows and sizeOfInstances columns
    double** data = new double*[numOfInstances];
    //create the rows inside of data
    for (int i = 0; i < numOfInstances; i++) {
        data[i] = new double[sizeOfInstances];
    }
    int dataIndex = 0;

    // Go through the rest of the file and store the data
    while (getline(file, line) && dataIndex < numOfInstances) {
        stringstream point(line);
        double newDouble;
        int index = 0;

        while (point >> newDouble && index < sizeOfInstances) {
            data[dataIndex][index] = newDouble;
            index++;
        }
        dataIndex++;
    }

    file.close();
    return data;
}
//sets up the clusters
double** setClusters(double** data, int numOfClusters, int numOfInstances, int sizeOfInstance) {
    //Randomly pick data for the clusters
    double** clusters = new double* [numOfClusters];

    // Use random_device for a better random seed
    random_device rd;
    mt19937 gen(rd());  // Initialize random number generator with seed
    uniform_int_distribution<int> dist(0, numOfInstances - 1);

    // Make a set for the data positions
    set<double*> clusterData;

    for (int i = 0; i < numOfClusters; i++) {
        int randomNum;
        do {
            randomNum = dist(gen);
        } while (clusterData.find(data[randomNum]) != clusterData.end());  // Ensure no duplicates
        clusterData.insert(data[randomNum]);

        // Create a new cluster and copy the data from the chosen instance
        double* newCluster = new double[sizeOfInstance];
        for (int j = 0; j < sizeOfInstance; j++) {
            newCluster[j] = data[randomNum][j];
        }
        clusters[i] = newCluster;
    }

    return clusters;
}

//Run through the Iterations and returns Initial, Final SSE and iterations | & are to save on memory managment and speed
double runIterations(int& maxIterations, double& convergenceThreshold, double** clusters, double** data, double& initialSSE, int& iterationsRan, int* labels, double* distancesToClosestCluster, int* clusterSizes, int numOfInstances, int sizeOfInstance, int numOfClusters) {
    double oldSSE = 0;
    //data for the loop | Initialize here to save on memory allocation
    //make an array to hold squared distances
    double** sqauredDistances = new double* [numOfClusters];
    double** newClusters = new double* [numOfClusters];
    for (int i = 0; i < numOfClusters; i++) {
        sqauredDistances[i] = new double[numOfInstances];
    }

    //i is the iteration we are on
    for (int i = 1; i <= maxIterations; i++) {
        double SSE = 0.0;
        //generate new clusters by taking the average of the seperated data points
        //reset the newClusters list
        //Reset an array to hold the calculated squared distances
        for (int k = 0; k < numOfClusters; k++) {
            newClusters[k] = new double[sizeOfInstance];
            // Initialize all elements to 0.0
            for (int j = 0; j < sizeOfInstance; j++) {
                newClusters[k][j] = 0.0;
                sqauredDistances[k][j] = 0.0;
            }
            //reset cluster sizes
            clusterSizes[k] = 0;
        }


        //An Iteration
        for (int j = 0; j < numOfInstances; j++) {
            //Find which cluster is closer, initializing with the first cluster distance
            int closest = 0;
            // -1 means that there is no closest distance yet
            double closestDist = calculateSquaredDistance(data[j], clusters[0], sizeOfInstance);

            //store the distances
            sqauredDistances[0][j] = closestDist;

            // // Loop through all clusters to find the closest
            for (int h = 1; h < numOfClusters; h++) {  // Start from 1 since 0 is already checked
                double dist = calculateSquaredDistance(data[j], clusters[h], sizeOfInstance);
                //store the distances
                sqauredDistances[h][j] = dist;
                if (dist < closestDist) {
                    closestDist = dist;  // Update closest distance
                    closest = h;         // Update closest cluster index
                }

            }

            // Calculate the SSE
            SSE += closestDist;

            //generate new clusters
            //add the data to the cluster
            for (int k = 0; k < sizeOfInstance; k++) {
                newClusters[closest][k] += data[j][k];
                //update who is the closest globally
                labels[j] = closest;
                distancesToClosestCluster[j] = closestDist;
                maxDistCluster[closest] = max(maxDistCluster[closest], sqrt(closestDist));
                minDistCluster[closest] = min(minDistCluster[closest], sqrt(closestDist));
            }

            //update the number of data points for this cluster
            clusterSizes[closest] += 1;
        }




        //check if the convergenceThreshold is reached and kill the run if so
        if (oldSSE != 0 && (oldSSE - SSE) / oldSSE < convergenceThreshold) {
            //save the iterations and Final SSE
            iterationsRan = i - 1;
            break;
        }
        //save initial
        if (i == 1) {
            initialSSE = SSE;
        }
        oldSSE = SSE;


        //once all of the data has been added together find the average for each cluster
        //Also handle singleton clusters
        for (int i = 0; i < numOfClusters; i++) {
            if (clusterSizes[i] != 1) {
                //is not a singleton cluster
                for (int j = 0; j < sizeOfInstance; j++) {
                    if (clusterSizes[i] > 0)
                        newClusters[i][j] /= clusterSizes[i];
                }
                //remove the old cluster
                delete[] clusters[i];
                //save the new clusters as the old
                clusters[i] = newClusters[i];
            }
            else {
                //Is a singleton cluster
                // Find the point contributing most to the SSE
                int maxErrorPointIndex = -1;
                double maxError = -1;
                for (int j = 0; j < numOfInstances; j++) {
                    double dist = sqauredDistances[i][j];
                    if (dist > maxError) {
                        maxError = dist;
                        maxErrorPointIndex = j;
                    }
                }
                // Reassign the cluster center to this point
                if (maxErrorPointIndex != -1) {
                    //create a new instance of this location
                    double* maxErrorPoint = new double[sizeOfInstance];
                    for (int u = 0; u < sizeOfInstance; u++) {
                        maxErrorPoint[u] = data[maxErrorPointIndex][u];
                    }
                    //unallocate old memory
                    delete[] clusters[i];
                    //assign clusters to its new point
                    clusters[i] = maxErrorPoint;
                }
            }
        }
    }
    //unallocate memory
    for (int i = 0; i < numOfClusters; i++) {
        delete[] sqauredDistances[i];
        delete[] newClusters[i];
    }
    delete[] sqauredDistances;
    delete[] newClusters;
    return oldSSE;
}
//returns a squared distance and quits if the closest Dist is exceeded to save on time (& avoids making copies to make it faster)
double calculateSquaredDistance(const double* x, const double* c, int sizeOfInstance) {
    double distance = 0;
    for (int i = 0; i < sizeOfInstance; i++) {
        double diff = x[i] - c[i];
        //add the squared difference
        distance += diff * diff;
    }
    return distance;  // Return the squared distance
}

//min-max normalization
void minMaxNorm(double** data, int numOfInstances, int sizeOfInstance) {
    // Initialize the min and max arrays for each column
    double* min = new double[sizeOfInstance];
    double* max = new double[sizeOfInstance];

    // Initialize min and max with the values from the first row
    for (int j = 0; j < sizeOfInstance; j++) {
        min[j] = data[0][j];
        max[j] = data[0][j];
    }

    // Find the min and max for each column
    for (int i = 1; i < numOfInstances; i++) {  
        for (int j = 0; j < sizeOfInstance; j++) {
            if (data[i][j] > max[j]) {
                max[j] = data[i][j];
            }
            if (data[i][j] < min[j]) {
                min[j] = data[i][j];
            }
        }
    }

    // Normalize the data
    for (int i = 0; i < numOfInstances; i++) {
        for (int j = 0; j < sizeOfInstance; j++) {
            if (max[j] != min[j]) {  // Avoid division by zero
                data[i][j] = (data[i][j] - min[j]) / (max[j] - min[j]);
            }
        }
    }

    // Free memory used for min and max arrays
    delete[] min;
    delete[] max;
}
//z-score normalization function
void zScoreNorm(double** data, int numOfInstances, int sizeOfInstance) {

    double* average = new double[sizeOfInstance];
    double* sd = new double[sizeOfInstance];
    for (int i = 0; i < sizeOfInstance; i++) {
        average[i] = 0.0;
        sd[i] = 0.0;
    }

    //find mean
    for (int i = 0; i < numOfInstances; i++) {
        for (int j = 0; j < sizeOfInstance; j++) {
            average[j] += data[i][j];
        }
    }
    //average the mean
    for (int j = 0; j < sizeOfInstance; j++) {
        average[j] /= numOfInstances;
    }
    
    //find the standard deviation
    for (int i = 0; i < numOfInstances; i++) {
        for (int j = 0; j < sizeOfInstance; j++) {
            sd[j] += (data[i][j] - average[j]) * (data[i][j] - average[j]);
        }
    }
    //divide the sd by the number of data points and take the sqrt
    for (int j = 0; j < sizeOfInstance; j++) {
        //Bessel's Correction
        sd[j] /= (numOfInstances -1);
        sd[j] = sqrt(sd[j]);
    }
   
    //now apply our newly found mean and standard deviation to all the data points to find their z score
    for (int i = 0; i < numOfInstances; i++) {
        for (int j = 0; j < sizeOfInstance; j++) {
            if (sd[j] != 0)
                data[i][j] = (data[i][j] - average[j]) / sd[j];
            else
                data[i][j] = 0;
        }
    }
    delete[] average;
    delete[] sd;
    //data will be saved and changed
}
//calculate the clusters with the random partition method
double** randomParitionClusters(double** data, int numOfClusters, int numOfInstances, int sizeOfInstance) {
    //initialize Data
    double** clusters = new double*[numOfClusters];
    int* clusterSize = new int[numOfClusters];

    for (int i = 0; i < numOfClusters; i++) {
        clusters[i] = new double[sizeOfInstance];
        for (int j = 0; j < sizeOfInstance; j++) {
            clusters[i][j] = 0.0;
        }
        clusterSize[i] = 0;
    }
    //Assign all data points to random clusters. Then compute the mean of them. These means will be the clusters that are returned
    int clusterNum = 0;
    //set up random generator
    // Use random_device for a better random seed
    random_device rd;
    mt19937 gen(rd());  // Initialize random number generator with seed
    uniform_int_distribution<int> dist(0, numOfClusters - 1);

    for (int i = 0; i < numOfInstances; i++) {
        //randomly pick a cluster
        clusterNum = dist(gen);
        //add the datapoint to the cluster slot and add one to the count
        for (int j = 0; j < sizeOfInstance; j++) {
            clusters[clusterNum][j] += data[i][j];
        }
        clusterSize[clusterNum]++;
    }
    //now find the means of all of the added up data
    for (int i = 0; i < numOfClusters; i++) {
        for (int j = 0; j < sizeOfInstance; j++) {
            clusters[i][j] /= clusterSize[i];
        }
    }
    //unallocate Data
    delete[] clusterSize;
    //return the newly made clusters
    return clusters;
}
//find the clusters using the maximum method
double** maximumMethodClusters(double** data, int numOfClusters, int numOfInstances, int sizeOfInstance) {
    double** clusters = new double* [numOfClusters];

    // Initialize clusters with 0.0
    for (int i = 0; i < numOfClusters; i++) {
        clusters[i] = new double[sizeOfInstance]();
    }

    // Random generator for choosing initial centers
    random_device rd;
    mt19937 gen(rd());  // Initialize random number generator with seed
    uniform_int_distribution<int> dist(0, numOfInstances - 1);

    // Choose the first cluster center (deep copy the data)
    int firstCenterIndex = dist(gen);
    for (int j = 0; j < sizeOfInstance; j++) {
        clusters[0][j] = data[firstCenterIndex][j];
    }

    // Select remaining clusters
    for (int i = 1; i < numOfClusters; i++) {
        double maxMinDist = -1;
        int selectedIndex = -1;

        // For each data point, find the minimum distance to any of the existing cluster centers
        for (int j = 0; j < numOfInstances; j++) {
            double minDist = numeric_limits<double>::max();

            // Find the closest center to data[j]
            for (int h = 0; h < i; h++) {
                double dist = calculateSquaredDistance(data[j], clusters[h], sizeOfInstance);
                // Keep track of the smallest distance
                if (minDist > dist) {
                    minDist = dist;
                }
            }

            // Track the point with the maximum minimum distance
            if (minDist > maxMinDist) {
                maxMinDist = minDist;
                selectedIndex = j;
            }
        }

        // Deep copy the selected point into the next cluster
        for (int j = 0; j < sizeOfInstance; j++) {
            clusters[i][j] = data[selectedIndex][j];
        }
    }

    return clusters;
}

// Function to compute the Calinski-Harabasz Index
double calinski_validity(double** data, double** clusters, double& finalSSE, int* clusterSizes, int numClusters, int numOfInstances, int sizeOfInstance) {
    //Check to see if it should run
    if (numClusters < 2 || numOfInstances <= numClusters) {
        return 0;
    }
    // Compute Between-cluster scatter (BSS)
    double BSS = 0.0;
    for (int i = 0; i < numClusters; i++) {
        if (clusterSizes[i] == 0) continue;
        
        double dist = 0.0;
        for (int j = 0; j < sizeOfInstance; j++) {
            double diff = clusters[i][j] - overallMean[j];
            dist += diff * diff;
        }
        BSS += clusterSizes[i] * dist;
    }
    // Compute Within-cluster scatter (WSS)
    double WSS = finalSSE;

    // Avoid division by zero
    if (WSS == 0) {
        return 0;
    }
    // Compute CH index
    double CH = (BSS / (numClusters - 1)) / (WSS / (numOfInstances - numClusters));
    return CH;
}
double silhouette_width(double** data, double** clusters, int* labels, int* clusterSizes, int numClusters, int numOfInstances, int sizeOfInstance) {
    // Step 1: Find the distance inside each cluster a(i)
    double* a_i = new double[numOfInstances];
    double total_silhouette = 0.0;

    // Step 2: Find the distance to the closest next cluster b(i)
    double* b_i = new double[numOfInstances];

    for (int i = 0; i < numOfInstances; i++) {
        a_i[i] = 0.0;
        b_i[i] = INFINITY;
    }
    // Step 1 - Calculate a(i) - average distance within the same cluster
    for (int i = 0; i < numOfInstances; i++) {
        int myCluster = labels[i];
        double totalDistToCluster = 0.0;
        int count = 0;

        // Loop through all data points to calculate the average distance to points in the same cluster
        for (int j = 0; j < numOfInstances; j++) {
            if (labels[j] == myCluster && i != j) {  // Do not include the point itself
                totalDistToCluster += calculateSquaredDistance(data[i], data[j], sizeOfInstance);
                count++;
            }
        }
        if (count > 0) {
            a_i[i] = totalDistToCluster / count;  // Calculate average distance within the cluster
        }
        // Step 2 - Calculate b(i) - minimum distance to the closest different cluster

        for (int j = 0; j < numClusters; j++) {
            if (j != myCluster) {  // Only consider other clusters
                double totalDistToOtherCluster = 0.0;
                int otherClusterCount = 0;

                // Loop through all data points to find those in cluster j
                for (int k = 0; k < numOfInstances; k++) {
                    if (labels[k] == j) {  // Only consider points in cluster j
                        totalDistToOtherCluster += calculateSquaredDistance(data[i], data[k], sizeOfInstance);
                        otherClusterCount++;
                    }
                }

                if (otherClusterCount > 0) {
                    double avgDistToOtherCluster = totalDistToOtherCluster / otherClusterCount;
                    b_i[i] = min(b_i[i], avgDistToOtherCluster);  // Keep the minimum distance to a different cluster
                }
            }
        }
    }
    


    // Step 3: Calculate the silhouette score for each point and sum it
    for (int i = 0; i < numOfInstances; i++) {
        double silhouette_i = (b_i[i] - a_i[i]) / max(a_i[i], b_i[i]);
        total_silhouette += silhouette_i;
    }
    //deallocate memory
    delete[] a_i;
    delete[] b_i;
    // Step 4: Return the average silhouette width
    return total_silhouette / numOfInstances;
}


double dunn(double** data, double** clusters, int* labels, int numClusters, int numOfInstances, int sizeOfInstance) {

    // Step 1: Calculate the maximum intra-cluster distance for each cluster
    double maxIntraClusterDist = 0.0;
    for (int i = 0; i < numClusters; i++) {
        // maxDistCluster[i] holds the maximum distance from a point in the cluster to the centroid
        maxIntraClusterDist = max(maxIntraClusterDist, maxDistCluster[i]);
    }

    // Step 2: Calculate the minimum inter-cluster distance
    double minInterClusterDist = numeric_limits<double>::infinity();
    for (int i = 0; i < numClusters; i++) {
        for (int j = i + 1; j < numClusters; j++) {
            // Calculate the Euclidean distance between centroids of clusters i and j
            double dist = sqrt(calculateSquaredDistance(clusters[i], clusters[j], sizeOfInstance));
            minInterClusterDist = min(minInterClusterDist, dist);
        }
    }

    // Step 3: Compute the Dunn Index
    double dunnIndex = minInterClusterDist / maxIntraClusterDist;
    return dunnIndex;
}

double davies_bouldin(double** data, double** clusters, int* labels, int numClusters, double* distancesToClosestCluster, int* clustersSizes, int numOfInstances, int sizeOfInstance) {
    double dbIndex = 0.0;

    // Step 1: Compute the compactness of each cluster
    double* avgCompact = new double[numClusters];
    for (int i = 0; i < numClusters; i++) {
        avgCompact[i] = 0.0;
    }
    for (int i = 0; i < numOfInstances; i++) {
        int clusterIndex = labels[i];
        avgCompact[clusterIndex] += distancesToClosestCluster[i];
    }

    // take average
    for (int i = 0; i < numClusters; i++) {
        avgCompact[i] /= clustersSizes[i];
    }

    // Step 2: Compute the Davies-Bouldin index
    for (int i = 0; i < numClusters; i++) {
        double maxRatio = -1.0;

        for (int j = 0; j < numClusters; j++) {
            if (i != j) {
                double dist = calculateSquaredDistance(clusters[i], clusters[j], sizeOfInstance);

                // no identicals
                if (dist > 0.0) {
                    double ratio = (avgCompact[i] + avgCompact[j]) / dist;
                    maxRatio = std::max(maxRatio, ratio);
                }
            }
        }

        //get the sum
        dbIndex += maxRatio;
    }

    // deallocate memory
    delete[] avgCompact;

    // Return the average Davies-Bouldin index
    return dbIndex / numClusters;
}



