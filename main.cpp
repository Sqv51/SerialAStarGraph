#include <vector>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <queue>
#include <iostream>
#include <algorithm>
#include <set>
#include <random>
#include <chrono>
struct Node {
    [[maybe_unused]] int id{};
    double x{}, y{};
    std::vector<std::pair<int, double>> neighbors; // Pair of neighbor node id and edge weight
};

struct Graph {
    std::unordered_map<int, Node> nodes; //unordered_map to store nodes with id as key

    void addNode(int id, [[maybe_unused]] double x, [[maybe_unused]] double y) {
        nodes[id] = {id, x, y, {}};
    }

    void addEdge(int id1, int id2, double weight) {
        nodes[id1].neighbors.emplace_back(id2, weight);
        nodes[id2].neighbors.emplace_back(id1, weight); // Assuming undirected graph
    }
};

double euclideanDistance(const Node& a, const Node& b) {
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
    //Euclidean distance between two nodes will use as heuristic
}


struct AStarNode {
    int id;
    [[maybe_unused]] double f, g, h;

    bool operator>(const AStarNode& other) const {
        return f > other.f;
    }
};

std::vector<int> aStar(const Graph& graph, int startId, int goalId) {
    std::unordered_map<int, double> gScore;
    std::unordered_map<int, int> cameFrom;
    std::priority_queue<AStarNode, std::vector<AStarNode>, std::greater<>> openSet;
    std::set<int> closedSet;

    for (const auto& node : graph.nodes) {
        gScore[node.first] = std::numeric_limits<double>::infinity();
    }
    gScore[startId] = 0;

    Node startNode = graph.nodes.at(startId);
    Node goalNode = graph.nodes.at(goalId);

    openSet.push({startId, euclideanDistance(startNode, goalNode), 0, euclideanDistance(startNode, goalNode)});

    while (!openSet.empty()) {
        AStarNode current = openSet.top();
        openSet.pop();

        if (current.id == goalId) {
            std::vector<int> path;
            while (cameFrom.find(current.id) != cameFrom.end()) {
                path.push_back(current.id);
                current.id = cameFrom[current.id];
            }
            path.push_back(startId);
            std::reverse(path.begin(), path.end());
            return path;
        }

        closedSet.insert(current.id);

        for (const auto& neighbor : graph.nodes.at(current.id).neighbors) {
            int neighborId = neighbor.first;
            double tentative_gScore = gScore[current.id] + neighbor.second;

            if (closedSet.find(neighborId) != closedSet.end()) {
                continue;
            }

            if (tentative_gScore < gScore[neighborId]) {
                cameFrom[neighborId] = current.id;
                gScore[neighborId] = tentative_gScore;
                double hScore = euclideanDistance(graph.nodes.at(neighborId), goalNode);
                openSet.push({neighborId, gScore[neighborId] + hScore, gScore[neighborId], hScore});
            }
        }
    }

    return {}; // Return empty path if no path is found
}

std::vector<int> dijkstra(const Graph& graph, int startId, int goalId) {
    std::unordered_map<int, double> dist;
    std::unordered_map<int, int> prev;
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> queue;

    for (const auto& node : graph.nodes) {
        dist[node.first] = std::numeric_limits<double>::infinity();
    }
    dist[startId] = 0;

    queue.push({0, startId});

    while (!queue.empty()) {
        int currentId = queue.top().second;
        queue.pop();

        for (const auto& neighbor : graph.nodes.at(currentId).neighbors) {
            int neighborId = neighbor.first;
            double altDist = dist[currentId] + neighbor.second;

            if (altDist < dist[neighborId]) {
                dist[neighborId] = altDist;
                prev[neighborId] = currentId;
                queue.push({altDist, neighborId});
            }
        }
    }

    std::vector<int> path;
    for (int nodeId = goalId; nodeId != startId; nodeId = prev[nodeId]) {
        path.push_back(nodeId);
    }
    path.push_back(startId);
    std::reverse(path.begin(), path.end());

    return path;
}

int main() {
    Graph graph;
   //graph with a million nodes and edges all generated randomly
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0, 100000);

    //start the timer
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100000; i++) {
        graph.addNode(i, dis(gen), dis(gen));
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to add nodes: " << elapsed.count() << "s" << std::endl;


    //add edges and make the graph connected

    //start the timer
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100000; i++) {
        for (int j = 0; j < 100; j++) {
            int neighborId = dis(gen);
            graph.addEdge(i, neighborId, euclideanDistance(graph.nodes.at(i), graph.nodes.at(neighborId)));
        }
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time taken to add edges: " << elapsed.count() << "s" << std::endl;

    //check if the graph is connected
    //start the timer
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100000; i++) {
        if (graph.nodes.at(i).neighbors.empty()) {
            std::cout << "Graph is not connected" << std::endl;
            return 1;
        }
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time taken to check if graph is connected: " << elapsed.count() << "s" << std::endl;




    //select random start and goal nodes
    int startId = dis(gen);
    int goalId = dis(gen);



    std::cout << "Start node: " << startId << std::endl;
    std::cout << "Goal node: " << goalId << std::endl;


    start = std::chrono::high_resolution_clock::now();
    std::cout << "Running A* algorithm..." << std::endl;
    std::vector<int> path = aStar(graph, startId, goalId);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time taken to run A*: " << elapsed.count() << "s" << std::endl;

    std::cout << "Path: ";
    for (int nodeId : path) {
        std::cout << nodeId << " ";
    }
    std::cout << std::endl;

    start = std::chrono::high_resolution_clock::now();
    std::cout << "Running Dijkstra's algorithm..." << std::endl;
    path = dijkstra(graph, startId, goalId);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Time taken to run Dijkstra's: " << elapsed.count() << "s" << std::endl;

    std::cout << "Path: ";
    for (int nodeId : path) {
        std::cout << nodeId << " ";
    }
    std::cout << std::endl;



    return 0;
}