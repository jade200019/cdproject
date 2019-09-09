//
// Created by ubuntu on 4/8/19.
//

#ifndef TRYENV_GRAPHHANDLER_H
#define TRYENV_GRAPHHANDLER_H

#include "Config.h"

#include <cstdlib>

using namespace boost;
//typedef adjacency_list < vecS, vecS, undirectedS,
//        no_property, property < edge_weight_t, float >, disallow_parallel_edge_tag > Graph;
struct VertexData {
    float x_pixel;
    float y_pixel;
};
typedef adjacency_list < vecS, vecS, undirectedS,
        VertexData, property < edge_weight_t, float >, disallow_parallel_edge_tag > Graph;
typedef graph_traits < Graph >::edge_descriptor Edge;
typedef graph_traits < Graph >::vertex_descriptor Vertex;
typedef std::pair<int, int> E;

typedef property_map<Graph, vertex_index_t>::type IndexMap;



class GraphHandler {

    std::vector < Edge > spanning_tree; // vector of edges of MST result

    std::string printEdge(Edge e, const Graph & gMST);

public:
    Graph g; // Graph g(MAX_NUM_NODE);
    GraphHandler() = default;
    ~GraphHandler() = default;

    void test_image(){
        // all output: 1
        cout << addEdge(0, 2, 1) << endl;
        cout << addEdge(1, 3, 1) << endl;
        cout << addEdge(1, 4, 2) << endl;
        cout << addEdge(2, 1, 7.5) << endl;
        cout << addEdge(2, 3, 3) << endl;
        cout << addEdge(3, 4, 1) << endl;
        cout << addEdge(4, 0, 1) << endl;
        cout << addEdge(4, 1, 1) << endl;
        cout << addEdge(4, 5, 3) << endl;

        cout << addEdge(6, 7, 1) << endl;
        cout << addEdge(0, 6, 5) << endl;
        cout << addEdge(6, 8, 1) << endl;

        this->minimumSpanningTree();

        property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
        std::cout << "Print the edges in the MST:" << std::endl;
        for (std::vector < Edge >::iterator ei = spanning_tree.begin();
             ei != spanning_tree.end(); ++ei) {
            std::cout << source(*ei, g) << " <--> " << target(*ei, g)
                      << " with weight of " << weight[*ei]
                      << std::endl;
        }

        this->saveDotImage(this->g, "bin/kruskal-eg.dot");
    }

    void test_random_image(){
        srand(0);
        for (int i = 0; i < 30; i++){
            unsigned v1 = (rand()%30);
            unsigned v2 = (rand()%30);
            float w = (rand()%30);
            addEdge(v1, v2, w);
        }

        this->minimumSpanningTree();

        this->saveDotImage(g, "bin/kruskal-eg.dot");
    }

    // REQUIRES: none
    // EFFECTS: return true if adding succeed
    bool addEdge(unsigned v1, unsigned v2, float weight){
        auto epair = add_edge(v1, v2, weight, g);
        return epair.second;
    }

    // REQUIRES: pts[4] of rRects in pixels
    // EFFECTS: return true if adding succeed
    bool addEdge(unsigned v1, unsigned v2, float weight,
            float v1_x_pixel, float v1_y_pixel, float v2_x_pixel, float v2_y_pixel ){
        auto epair = add_edge(v1, v2, weight, g);
        g[v1].x_pixel = v1_x_pixel;
        g[v1].y_pixel = v1_y_pixel;
        g[v2].x_pixel = v2_x_pixel;
        g[v2].y_pixel = v2_y_pixel;
        return epair.second;
    }

    // REQUIRES: edges are added
    // EFFECTS: store MST edges in the vector spanning_tree
    void minimumSpanningTree(){
        kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));
    }

    // REQUIRES: MST is calculated
    // EFFECTS: cut MST so that degree <= 2
    void cutBranch(std::string pathMST, std::string pathCut);

    // REQUIRES: none
    // EFFECTS: save to a dot file
    // dot -Tpng bin/kruskal-eg.dot > bin/kruskal-eg.png
    void saveDotImage(const Graph & G, std::string path);

    // REQUIRES: a set of index to rRect: s, x, y, r
    // EFFECTS: a corresponding vertex index: v
    int sxyrIndex2vIndex(int s, int x, int y, int r);

    // REQUIRES: a vertex index: v
    // EFFECTS: a corresponding set of index to rRect: s, x, y, r, store in a vector
    vector<int> vIndex2sxyrIndex(int v);

    void visualizeGraph(const Graph & G, std::string path);
};


#endif //TRYENV_GRAPHHANDLER_H
