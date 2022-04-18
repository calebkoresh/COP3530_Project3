#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <limits>
#include <queue>
#include <chrono>
using namespace std;

int BinarySearch(vector<string> names, string key) {
    int low = 0;
    int high = 199999;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (names[mid].substr(0, key.length()) < key) {
            low = mid + 1;
        }
        else if (names[mid].substr(0, key.length()) > key) {
            high = mid - 1;
        }
        else if (names[mid].substr(0, key.length()) == key) {
            return mid;
        }
    }
    return -1;
}

void ReformatKey(string &key) {
    size_t pos = key.find(" ");
    if (key[0] - 'a' >= 0) {
        key[0] = key[0] - 32;
    }
    if (key[pos + 1] - 'a' >= 0) {
        key[pos + 1] = key[pos + 1] - 32;
    }
}

int StringToInteger(string str) {
    for (int i = 0; i < str.length(); i++) {
        if (str[i] - '0' < 0 || str[i] - '0' > 9) {
            return -1;
        }
    }
    return stoi(str);
}

int Dijkstra(int source, int destination, map<int, vector<int>>& graph) {
    //Implementation without priority queue
    /*
    map<int, int> distances;
    set<int> completed;
    int count = graph.size();
    vector<int> vec;
    int minVertex = source;
    distances[minVertex] = 0;
    for (int i = 0; i < count; i++) {
        int min = 9999;
        vec = graph[minVertex];
        for (int i = 0; i < vec.size(); i++) {
            if (distances.count(vec[i]) == 0) {
                distances[vec[i]] = distances[minVertex] + 1;
            }
            else if (distances[minVertex] + 1 < distances[vec[i]]) {
                distances[vec[i]] = distances[minVertex] + 1;
            }
        }
        for (auto iter = distances.begin(); iter != distances.end(); iter++) {
            if (iter->second < min && completed.count(iter->first) == 0) {
                min = iter->second;
                minVertex = iter->first;
                break;
            }
        }
        completed.insert(minVertex);
    }
    return distances[destination];
    */

    vector<int> distances(graph.size(), INT32_MAX);
    distances[source] = 0;
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push({0,source});
    while (!pq.empty())
    {
        int vertex = pq.top().second;
        int degree = pq.top().first;
        pq.pop();
        if (distances[vertex] < degree)
            continue;
        vector<int> &adj = graph[vertex];
        for (int i = 0; i < adj.size(); i++)
        {
            if (distances[adj[i]] > degree + 1)
            {
                distances[adj[i]] = degree + 1;
                pq.push({distances[adj[i]], adj[i]});
            }
        }
    }
    return distances[destination];
}

int BreadthFirstSearch(int source, int destination, map<int, vector<int>>& graph) {
    set<int> completed;
    queue<int> q;
    vector<int> distances(graph.size(), INT32_MAX);
    distances[source] = 0;

    completed.insert(source);
    q.push(source);

    while (!q.empty())
    {
        int temp = q.front();
        q.pop();

        vector<int> v = graph[temp];
        for (int i = 0; i < v.size(); i++)
        {
            if (completed.count(v[i]) == 0)
            {
                completed.insert(v[i]);
                q.push(v[i]);
                distances[v[i]] = distances[temp] + 1;
                if (v[i] == destination)
                    return distances[v[i]];
            }
        }
    }
    return -1;
}

int main() {

    //=========================FILE INPUT=========================//

    ifstream firstNamesFile("NationalNames.csv");
    ifstream lastNamesFile("facebook-lastnames-withcount.txt");

    // Variables

    vector<string> firstNames;
    vector<int> firstPopularity;
    vector<string> lastNames;
    vector<int> lastPopularity;
    vector<string> names;
    vector<int> following;

    // Input for first names (from the 2000 most popular U.S. baby names in 1880)

    string first;
    firstNamesFile >> first;
    vector<int> firstABCOrder;
    for (int i = 0; i < 2000; i++) {
        firstNamesFile >> first;
        size_t pos1 = first.find(',');
        size_t pos2 = first.find(',', pos1+1);string temp = first.substr(pos1+1, pos2-pos1-1);
        int count = 0;
        for (int j = 0; j < firstNames.size(); j++) {
            if (temp < firstNames[j]) {
                firstABCOrder[j]++;
            }
            else {
                count++;
            }
        }
        firstABCOrder.push_back(count);
        firstNames.push_back(temp);
        size_t pos3 = first.find(',', pos2+1);
        size_t pos4 = first.find(',', pos3+1);
        firstPopularity.push_back(stoi(first.substr(pos4+1)));
    }

    // First names in alphabetical order to help for searching

    map<int, pair<string, int>> firstMap;
    for (int i = 0; i < 2000; i++) {
        firstMap[firstABCOrder[i]] = make_pair(firstNames[i], firstPopularity[i]);
    }

    // Input for last names (from the top 1000 Facebook last names)

    string last;
    vector<int> lastABCOrder;
    for (int i = 0; i < 100; i++) {
        lastNamesFile >> last;
        lastPopularity.push_back(stoi(last));
        lastNamesFile >> last;
        last[0] = last[0] - 32;
        int count = 0;
        for (int j = 0; j < lastNames.size(); j++) {
            if (last < lastNames[j]) {
                lastABCOrder[j]++;
            }
            else {
                count++;
            }
        }
        lastABCOrder.push_back(count);
        lastNames.push_back(last);
    }

    // Last names in alphabetical order to help for searching

    map<int, pair<string, int>> lastMap;
    for (int i = 0; i < 100; i++) {
        lastMap[lastABCOrder[i]] = make_pair(lastNames[i], lastPopularity[i]);
    }

    // Generating all combinations of first and last names
    // Generating a metric for popularity based on the counts for first and last names
    for (int i = 0; i < 2000; i++) {
        for (int j = 0; j < 100; j++) {
            string name = firstMap[i].first;
            name.append(" ");
            name.append(lastMap[j].first);
            names.push_back(name);

            // Formula for calculating how many people each name follows
            // Popularity of the first and last names directly correlates to the number of people they follow
            int pop = (int)(firstMap[i].second/10.0 + lastMap[j].second/1000.0);
            following.push_back(pop);
        }
    }

    // If you want to print the names

    /*
    for (int i = 0; i < 200000; i++) {
        cout << names[i] << endl;
        cout << following[i] << endl;
    }
    */

    //=========================CREATING GRAPH=========================//


    map<int, vector<int>> graph;
    for (int i = 0; i < 200000; i++) {
        vector<int> &vec = graph[i];
        for (int j = 0; j < following[i]; j++) {
            // RNG
            int random = rand() % 200000;
            // Use set data structure to ensure there are no repeats
            // Models a real social network because you can't follow the same person twice
            // Also not allowed to follow yourself
            set<int> randoms;
            while (randoms.count(random) != 0 || random == i) {
                random = (random + 1) % 200000;
            }
            randoms.insert(random);
            // Reverse of normal adjacency list to save time
            vec.push_back(random);
        }
    }
    cout << "Social Network Generated!" << endl << endl;

    // Printing first 10 names and who they follow
    /*
    for (int i = 0; i < 10; i++) {
        vector<int> vec = graph[i];
        cout << names[i] << ": ";
        for (int j = 0; j < vec.size(); j++) {
            cout << names[vec[j]] << ", ";
        }
        cout << endl;
    }
    */
    //=========================MENU=========================//
    cout << "1) Check degrees of seperation between two people in the social network" << endl;
    cout << "2) Speed test breadth first search versus Dijkstra's algorithm" << endl;
    int choice;
    cin >> choice;
    while (choice != 1 && choice != 2) {
        cout << "Invalid choice please pick 1 or 2" << endl;
        cin >> choice;
    }
    if (choice == 1) {

        //=========================SEARCHING FOR NAMES=========================//

        cout << "Choose two people from the social network to find their degree of separation" << endl;
        cout << "Search for the name of the first person below" << endl;

        string name1, name2;
        int foo = 0;
        string key;
        getline(cin, key);
        while (true) {
            getline(cin, key);
            ReformatKey(key);
            int index = BinarySearch(names, key);
            if (index == -1) {
                cout << "Name not found. Please enter another name to search for" << endl;
                continue;
            }

            int lower = index;
            while (names[lower].substr(0, key.length()) == key) {
                lower--;
            }
            lower++;

            int upper = index;
            while (names[upper].substr(0, key.length()) == key) {
                upper++;
            }
            for (int i = lower; i < upper; i++) {
                cout << i - lower + 1 << ": " << names[i] << endl;
            }

            cout << "Enter the number of the name you'd like to choose OR enter 0 to search again" << endl;
            int bar = 0;
            while (true) {
                string input;
                getline(cin, input);
                int choice = StringToInteger(input);
                if (choice == 0) {
                    bar = 1;
                    cout << "Search again" << endl;
                } else if (choice >= 1 && choice <= upper - lower) {
                    if (foo == 0) {
                        name1 = names[choice + lower - 1];
                        cout << "name 1: " << name1 << endl;
                        foo++;
                        cout << "Now search for the name of the second person" << endl;
                        bar = 1;
                    } else if (foo == 1) {
                        name2 = names[choice + lower - 1];
                        cout << "name 2: " << name2 << endl;
                    }
                } else {
                    cout << "Invalid choice, please pick another number" << endl;
                    continue;
                }
                break;
            }
            if (bar == 1) {
                continue;
            }
            break;

        }

        //=========================DIJKSTRA=========================//

        /* important for caleb
         * graph is called graph. its a map<int, vector<int>> and stores the indices of the names, not the actual names
         * names is a vector with all 200000 names in ABC order
         * following is a vector with the number of people each person follows, if u wanna use it
         * source name for dijsktra is name1, destination name is name2
         */
        int source = BinarySearch(names, name1);
        int destination = BinarySearch(names, name2);

        cout << endl << "Degrees of Seperation: ";

        //Calculate the duration of Dijkstra's Algorithm
        //Source: https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
        auto startDijkstra = std::chrono::high_resolution_clock::now();
        int degreesDijkstra = Dijkstra(source, destination, graph);
        auto stopDijkstra = std::chrono::high_resolution_clock::now();
        auto durationDijkstra = std::chrono::duration_cast<std::chrono::microseconds>(stopDijkstra - startDijkstra);
        double dijkstraTime = durationDijkstra.count() / 1000000.0;
        cout << degreesDijkstra << endl;
        cout << "Djikstra's Algorithm Execution Time: " << dijkstraTime << " seconds" << endl;

        //====================BREADTH FIRST SEARCH====================//

        //Calculate the duration of BFS
        //Source: https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
        auto startBfs = std::chrono::high_resolution_clock::now();
        int degreesBfs = BreadthFirstSearch(source, destination, graph);
        auto stopBfs = std::chrono::high_resolution_clock::now();
        auto durationBfs = std::chrono::duration_cast<std::chrono::microseconds>(stopBfs - startBfs);
        double bfsTime = durationBfs.count() / 1000000.0;
        cout << "Breadth First Search Execution Time: " << bfsTime << " seconds" << endl;
    }


    else if (choice == 2) {

        //==========================Speed Test========================//
        //Choose random source to avoid bias
        int source = rand() % 200000;

        //Ask user how many iterations to execute
        cout << "How many iterations should be tested?" << endl;
        int iterations;
        cin >> iterations;

        //Find average time it takes to reach a destination for Dijkstra's Algorithm
        double dijkstraSum = 0.0;
        for (int i = 0; i < iterations; i++)
        {
            //Random destinations to avoid bias
            int destination = rand() % 200000;

            //Test speed
            auto startDijkstra = std::chrono::high_resolution_clock::now();
            int degreesDijkstra = Dijkstra(source, destination, graph);
            auto stopDijkstra = std::chrono::high_resolution_clock::now();
            auto durationDijkstra = std::chrono::duration_cast<std::chrono::microseconds>(stopDijkstra - startDijkstra);
            double dijkstraTime = durationDijkstra.count() / 1000000.0;
            dijkstraSum += dijkstraTime;
        }
        double dijkstraAvg = dijkstraSum / double(iterations);
        cout << "Djikstra's Algorithm Average Execution Time: " << dijkstraAvg << " seconds" << endl;

        //Find average time it takes to reach a destination for Breadth First Search
        double bfsSum = 0.0;
        for (int i = 0; i < iterations; i++)
        {
            //Random destination to avoid bias
            int destination = rand() % 200000;

            //Test speed
            auto startBfs = std::chrono::high_resolution_clock::now();
            int degreesBfs = BreadthFirstSearch(source, destination, graph);
            auto stopBfs = std::chrono::high_resolution_clock::now();
            auto durationBfs = std::chrono::duration_cast<std::chrono::microseconds>(stopBfs - startBfs);
            double bfsTime = durationBfs.count() / 1000000.0;
            bfsSum += bfsTime;
        }
        double bfsAvg = bfsSum / double(iterations);
        cout << "Breadth First Search Average Execution Time: " << bfsAvg << " seconds" << endl;
    }
    return 0;
}