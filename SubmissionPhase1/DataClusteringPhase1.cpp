#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <set>

using namespace std;

int checkTheArguments(string fileName, int numOfClusters, int maxIterations, double convergenceThreshold, int numOfRuns);
double* readData(string fileName, int& numOfLines, int& sizeOfLine);
double* setClusters(double* data, int numOfLines, int sizeOfLine, int numOfClusters);
void outputClusters(string fileName, double* clusters, int numOfClusters, int sizeOfLine);

int main(int argc, char* argv[])
{
    string fileName;
    int numOfClusters = -1;
    int maxIterations = -1;
    double convergenceThreshold = -1;
    int numOfRuns = -1;



    // make sure there are 5 arguments
    if (argc != 6)
    {
        cout << "Usage: <F> <K> <I> <T> <R>" << endl;
        cout << "F: Name of the data file" << endl;
        cout << "K: Number of clusters (positive integer greater than 1)" << endl;
        cout << "I: Maximum number of iterations (positive integer)" << endl;
        cout << "T: Convergence threshold (non-negative real number)" << endl;
        cout << "R: Number of runs (positive integer)" << endl;
        return 1;
    }

    // collect the data from arguments
    fileName = argv[1];
    numOfClusters = stoi(argv[2]);
    maxIterations = stoi(argv[3]);
    convergenceThreshold = stod(argv[4]);
    numOfRuns = stoi(argv[5]);

    //check the data
    int result = checkTheArguments(fileName, numOfClusters, maxIterations, convergenceThreshold, numOfRuns);

    if (result == 1) {
        //error
        return 1;
    }

    // Play with the data here!
    int sizeOfLine, numOfLines;
    //Read the data from the file
    double* data = readData(fileName, numOfLines, sizeOfLine);

    //obtain clusters
    double* clusters = setClusters(data, numOfLines, sizeOfLine, numOfClusters);

    //print the clusters
    for (int i = 0; i < numOfClusters * sizeOfLine; i++) {
        
        if (i % sizeOfLine == 0) {
            cout << endl;
        }
        cout << clusters[i] << " ";
    }
    cout << endl;
    //Send the clusters to an output file
    outputClusters(fileName, clusters, numOfClusters, sizeOfLine);

    return 0;
}

//makes sure the arguments are valid
int checkTheArguments(string fileName, int numOfClusters, int maxIterations, double convergenceThreshold, int numOfRuns) {
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
double* readData(string fileName, int& numOfLines, int& sizeOfLine) {
    ifstream file(fileName);
    string line;
    //read the first line to see how big the data is
    getline(file, line);
    stringstream ss(line);
    ss >> numOfLines >> sizeOfLine;

    double* data = new double[numOfLines * sizeOfLine];
    int index = 0;
    //go through, read the data, and save it to a double array
    while (getline(file, line)) {
        stringstream point(line);
        double newDouble;
        while (point >> newDouble) {
            data[index] = newDouble;
            index++;
        }
    }
    file.close();
    return data;
}
//sets up the clusters
double* setClusters(double* data, int numOfLines, int sizeOfLine, int numOfClusters) {
    //randomly pick lines for the clusters
    double* clusters = new double[sizeOfLine * numOfClusters];

    // Seed the random number generator
    srand(static_cast<unsigned int>(time(nullptr)));

    //make a set for the line number
    set<int> clusterLines;

    for (int i = 0; i < numOfClusters; i++) {
        int randomNum;
        do {
            randomNum = rand() % numOfLines;
        } while (clusterLines.find(randomNum) != clusterLines.end());
        clusterLines.insert(randomNum);

        //start from that line and put all the data from that line into the cluster

        for (int j = 0; j < sizeOfLine; j++) {
            // Use sizeOfLine as the correct multiplier to place data in the correct spot
            clusters[i * sizeOfLine + j] = data[randomNum * sizeOfLine + j];
        }
    }

    return clusters;
}
//send the clusters to an output file
void outputClusters(string fileName, double* clusters, int numOfClusters, int sizeOfLine) {

    //make sure the output folder exists
    struct stat info;
    if (stat("Output", &info) != 0) {
        system("mkdir Output");
    }

    //get the last part of the path
    char lastSlash = fileName.find_last_of('/\\');
    string newName = (lastSlash == string::npos) ? fileName : fileName.substr(lastSlash + 1);

    //Make a new file
    ofstream outfile("Output/output_" + newName);

    //check if there is output file
    if (!outfile) {
        cout << "Error: Unable to create the output file." << endl;
        return;
    }
    //put the num of clusters and lineNum
    outfile << numOfClusters<< " " << sizeOfLine;
    for (int i = 0; i < numOfClusters * sizeOfLine; i++) {

        if (i % sizeOfLine == 0) {
            outfile << endl;
        }
        outfile << clusters[i] << " ";
    }
    cout << "Output File Created" << endl;
    outfile.close();
}
