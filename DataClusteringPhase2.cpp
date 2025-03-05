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
#include <vector>

using namespace std;

//Phase 1 (Gather the data)
int checkTheArguments(string fileName, int numOfClusters, int maxIterations, double convergenceThreshold, int numOfRuns, int typeOfClustering);
vector<vector<double>> readData(string fileName, int& numOfInstances, int& sizeOfInstances);
vector<vector<double>> setClusters(vector<vector<double>>& data, int numOfClusters);


//Phase 2 (Run the K-Means)
double runIterations(int& maxIterations, double& convergenceThreshold, vector<vector<double>>& clusters, vector<vector<double>>& data, double& initialSSE, int& iterationsRan);
double calculateSquaredDistance(const vector<double>& x, const vector<double>& c);

//Phase 3 (Normalization and Initialization)
vector<vector<double>> minMaxNorm(vector<vector<double>>& data);
vector<vector<double>> zScoreNorm(vector<vector<double>>& data);
double sqrtNewton(double num, double precision);
vector<vector<double>> randomParitionClusters(vector<vector<double>>& data, int numOfClusters);
vector<vector<double>> maximumMethodClusters(vector<vector<double>>& data, int numOfClusters);

int main(int argc, char* argv[])
{
    string fileName;
    int numOfClusters = -1;
    int maxIterations = -1;
    double convergenceThreshold = -1;
    int numOfRuns = -1;
    int typeOfClustering = -1;


    // make sure there are 6 arguments
    if (argc != 7)
    {
        cout << "Usage: <F> <K> <I> <T> <R> <V>" << endl;
        cout << "F: Name of the data file" << endl;
        cout << "K: Number of clusters (positive integer greater than 1)" << endl;
        cout << "I: Maximum number of iterations (positive integer)" << endl;
        cout << "T: Convergence threshold (non-negative real number)" << endl;
        cout << "R: Number of runs (positive integer)" << endl;
        cout << "V: Cluster Generation (random | random partition | maximum method) (0 | 1 | 2)" << endl;
        return 1;
    }

    // collect the data from arguments
    fileName = argv[1];
    numOfClusters = stoi(argv[2]);
    maxIterations = stoi(argv[3]);
    convergenceThreshold = stod(argv[4]);
    numOfRuns = stoi(argv[5]);
    typeOfClustering = stoi(argv[6]);

    //check the data
    int result = checkTheArguments(fileName, numOfClusters, maxIterations, convergenceThreshold, numOfRuns, typeOfClustering);

    if (result == 1) {
        //error
        return 1;
    }
    // Play with the data here!
    int sizeOfInstances, numOfInstances;
    //Read the data from the file
    vector<vector<double>> data = readData(fileName, numOfInstances, sizeOfInstances);
    data = zScoreNorm(data);

    vector<vector<double>> clusters;
    double bestFinalSSE = 0;
    double bestInitialSSE = 0;
    int bestRun = 0;
    int bestIteration = 0;
    //run for the number of runs
    for (int i = 1; i <= numOfRuns; i++) {
        //obtain clusters | check which type of clustering to do
        switch(typeOfClustering){
        case 0:
            clusters = setClusters(data, numOfClusters);
            break;
        case 1:
            clusters = randomParitionClusters(data, numOfClusters);
            break;
        case 2:
            clusters = maximumMethodClusters(data, numOfClusters);
            break;
        }
        //display which run we are on
        cout << "Run " << i << endl << "__________" << endl;
        //phase 2
        //run the iterations and display the SSE
        double initialSSE;
        int iterations;
        double finalSSE = runIterations(maxIterations, convergenceThreshold, clusters, data,initialSSE, iterations);
        if (bestFinalSSE == 0 || finalSSE < bestFinalSSE) {
            bestFinalSSE = finalSSE;
            bestRun = i;
        }
        if (bestInitialSSE == 0 || initialSSE < bestInitialSSE) {
            bestInitialSSE = initialSSE;
        }
        if (bestIteration == 0 || iterations < bestIteration) {
            bestIteration = iterations;
        }
    }
    //Print out what the best run was
    cout << "Best Run: " << bestRun << endl;
    cout << "Best Iterations: " << bestIteration << endl;
    cout << "Best Initial SSE: " << bestInitialSSE<< endl;
    cout << "Best Final SSE: " << bestFinalSSE << endl;

    
    return 0;
}

//makes sure the arguments are valid
int checkTheArguments(string fileName, int numOfClusters, int maxIterations, double convergenceThreshold, int numOfRuns, int typeOfClustering) {
    try
    {
        // Validate the arguments
        ifstream file(fileName);
        if (!file)
        {
            throw invalid_argument(("File: " + fileName + " does not exist.").c_str());
        }
        file.close();

        if (numOfClusters <= 1)
        {
            throw invalid_argument("Number of clusters (K) must be greater than 1.");
        }
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
vector<vector<double>> readData(string fileName, int& numOfInstances, int& sizeOfInstances) {
    ifstream file(fileName);
    string line;

    // Read the first line to get the number of lines and the size of each line
    getline(file, line);
    stringstream ss(line);
    ss >> numOfInstances >> sizeOfInstances;

    // Create a 2D vector with numOfInstances rows and sizeOfInstances columns
    vector<vector<double>> data(numOfInstances, vector<double>(sizeOfInstances));

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
    return std::move(data);
}
//sets up the clusters
vector<vector<double>> setClusters(vector<vector<double>>& data, int numOfClusters) {
    //Randomly pick data for the clusters
    vector<vector<double>> clusters(numOfClusters, vector<double>(data[0].size()));
    int numOfInstances = data.size();

    // Use random_device for a better random seed
    random_device rd;
    mt19937 gen(rd());  // Initialize random number generator with seed
    uniform_int_distribution<int> dist(0, numOfInstances - 1);

    // Make a set for the data psoition
    set<vector<double>> clusterData;

    
    for (int i = 0; i < numOfClusters; i++) {
        int randomNum;
        do {
            randomNum = dist(gen);
        } while (clusterData.find(data[randomNum]) != clusterData.end());
        clusterData.insert(data[randomNum]);


        
        //save cluster
        clusters[i] = data[randomNum];
        
    }

    return clusters;
}

//Run through the Iterations and returns Initial, Final SSE and iterations | & are to save on memory managment and speed
double runIterations(int& maxIterations, double& convergenceThreshold, vector<vector<double>>& clusters, vector<vector<double>>& data, double& initialSSE, int& iterationsRan) {
    
    
    double oldSSE = 0;
    int numOfInstances = data.size();
    int numOfClusters = clusters.size();
    int sizeOfInstances = data[0].size();

    //data for the loop | Initialize here to save on memory allocation
    vector<vector<double>> newClusters(numOfClusters, vector<double>(sizeOfInstances));
    vector<int> clusterSize(numOfClusters);


    //i is the iteration we are on
    for (int i = 1; i <= maxIterations; i++) {
        double SSE = 0.0;
        //generate new clusters by taking the average of the seperated data points
        //reset the newClusters list
        for (auto& cluster : newClusters) {
            fill(cluster.begin(), cluster.end(), 0.0);
        }
        //make a vector to keep track of the number of data points in each cluster
        fill(clusterSize.begin(), clusterSize.end(), 0);

        //Make an array to hold the calculated squared distances
        vector<vector<double>> sqauredDistances(numOfClusters, vector<double>(numOfInstances, 0.0));
        //An Iteration
        for (int j = 0; j < numOfInstances; j++) {
            //Find which cluster is closer, initializing with the first cluster distance
            int closest = 0;
            // -1 means that there is no closest distance yet
            double closestDist = calculateSquaredDistance(data[j], clusters[0]);
            //store the distances
            sqauredDistances[0][j] = closestDist;

            // // Loop through all clusters to find the closest
            for (int h = 1; h < numOfClusters; h++) {  // Start from 1 since 0 is already checked
                double dist = calculateSquaredDistance(data[j], clusters[h]);
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
            for (int k = 0; k < sizeOfInstances; k++) {
                newClusters[closest][k] += data[j][k];
            }

            //update the number of data points for this cluster
            clusterSize[closest] += 1;
        }

            
        
        
        //check if the convergenceThreshold is reached and kill the run if so
        if ((oldSSE - SSE) / oldSSE < convergenceThreshold && oldSSE != 0) {
            //save the iterations and Final SSE
            iterationsRan = i - 1;
            break;
        }
        //save initial
        if (i == 1) {
            initialSSE = SSE;
        }
        cout << "Iteration: " << i << "| SSE: " << SSE << endl;
        oldSSE = SSE;


        //once all of the data has been added together find the average for each cluster
        //Also handle singleton clusters
        for (int i = 0; i < numOfClusters; i++) {
            if (clusterSize[i] != 1) {
                //is not a singleton cluster
                for (int j = 0; j < sizeOfInstances; j++) {
                    if (clusterSize[i] > 0)
                        newClusters[i][j] /= clusterSize[i];
                }
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
                    clusters[i] = data[maxErrorPointIndex];
                }
            }
        }
        
    }
    return oldSSE;
}
//returns a squared distance and quits if the closest Dist is exceeded to save on time (& avoids making copies to make it faster)
double calculateSquaredDistance(const vector<double>& x, const vector<double>& c) {
    double distance = 0;
    for (int i = 0; i < x.size(); i++) {
        double diff = x[i] - c[i];
        //add the squared difference
        distance += diff * diff;
    }
    return distance;  // Return the squared distance
}

//min-max normalization
vector<vector<double>> minMaxNorm(vector<vector<double>>& data) {
    
    //initialize a starting data value
    vector<double> max = data[0];
    vector<double> min = data[0];
    int lineSize = data[0].size();

    vector<vector<double>> normData(data.size(), vector<double>(lineSize, 0.0));
    //find the max and min for  each column and store them
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < lineSize; j++) {
            double value = data[i][j];
            if (value > max[j]) {
                max[j] = value;
            }
            if (value < min[j]) {
                min[j] = value;
            }
        }
    }
    
    //Now normalize everything
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < lineSize; j++) {
            if(max[j] != min[j])
                normData[i][j] = (data[i][j] - min[j]) / (max[j] - min[j]);
        }
    }
    //return the normalized data
    return normData;
}
//z-score normalization function
vector<vector<double>> zScoreNorm(vector<vector<double>>& data) {
    int lineSize = data[0].size();

    vector<vector<double>> normData(data.size(), vector<double>(lineSize, 0.0));
    vector<double> average(lineSize, 0.0);
    //find mean
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < lineSize; j++) {
            average[j] += data[i][j];
        }
    }
    //average the mean
    for (int j = 0; j < lineSize; j++) {
        average[j] /= data.size();
    }
    vector<double> sd(lineSize, 0.0);
    //find the standard deviation
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < lineSize; j++) {
            sd[j] += (data[i][j] - average[j]) * (data[i][j] - average[j]);
        }
    }
    //divide the sd by the number of data points and take the sqrt
    for (int j = 0; j < lineSize; j++) {
        //Bessel's Correction
        sd[j] /= (data.size() -1);
        sd[j] = sqrtNewton(sd[j], 1e-6);
    }
   
    //now apply our newly found mean and standard deviation to all the data points to find their z score
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < lineSize; j++) {
            if (sd[j] != 0)
                normData[i][j] = (data[i][j] - average[j]) / sd[j];
            else
                normData[i][j] = 0;
        }
    }
    return normData;
}

//finding the sqrt without using sqrt(), Newton is a fast method
double sqrtNewton(double num, double precision) {
    if (num == 0) return 0;
    // Initial guess
    double x = num;  
    double root;

    while (true) {
        root = 0.5 * (x + (num / x));
        // Stop when close enough
        if (abs(root - x) < precision)  
            break;
        x = root;
    }

    return root;
}
//calculate the clusters with the random partition method
vector<vector<double>> randomParitionClusters(vector<vector<double>>& data, int numOfClusters) {
    vector<vector<double>> clusters(numOfClusters, vector<double>(data[0].size(), 0.0));
    //Assign all data points to random clusters. Then compute the mean of them. These means will be the clusters that are returned
    int clusterNum = 0;
    vector<int> clusterSize(numOfClusters, 0);

    //set up random generator
    // Use random_device for a better random seed
    random_device rd;
    mt19937 gen(rd());  // Initialize random number generator with seed
    uniform_int_distribution<int> dist(0, numOfClusters - 1);

    for (int i = 0; i < data.size(); i++) {
        //randomly pick a cluster
        clusterNum = dist(gen);
        //add the datapoint to the cluster slot and add one to the count
        for (int j = 0; j < data[i].size(); j++) {
            clusters[clusterNum][j] += data[i][j];
        }
        clusterSize[clusterNum]++;
    }
    //now find the means of all of the added up data
    for (int i = 0; i < numOfClusters; i++) {
        for (int j = 0; j < clusters[i].size(); j++) {
            clusters[i][j] /= clusterSize[i];
        }
    }
    
    //return the newly made clusters
    return clusters;
}
//find the clusters using the maximum method
vector<vector<double>> maximumMethodClusters(vector<vector<double>>& data, int numOfClusters) {
    vector<vector<double>> clusters(numOfClusters, vector<double>(data[0].size(), 0.0));

    int numOfInstances = data.size();
    //choose the first center
    //set up the random selector
    // Use random_device for a better random seed
    random_device rd;
    mt19937 gen(rd());  // Initialize random number generator with seed
    uniform_int_distribution<int> dist(0, numOfInstances - 1);

    clusters[0] = data[dist(gen)];
    for (int i = 1; i < numOfClusters; i++) {
        double maxMinDist = -1;
        int selectedIndex = -1;

        for (int j = 0; j < numOfInstances; j++) {
            double minDist = numeric_limits<double>::max(); //makes it the biggest value

            // Find the closest center to data[j]
            for (int h = 0; h < i; h++) {
                double dist = calculateSquaredDistance(data[j], clusters[h]);
                //picks the smallest distance and keeps it
                if (minDist > dist) {
                    minDist = dist;
                }
            }

            // Track the point with the maximum of these minimum distances
            if (minDist > maxMinDist) {
                maxMinDist = minDist;
                selectedIndex = j;
            }
        }

        // Assign the farthest valid point as the next center
        clusters[i] = data[selectedIndex];
    }
    return clusters;
}
