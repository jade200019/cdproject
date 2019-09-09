//
// Created by ubuntu on 4/8/19.
//

#include "GraphHandler.h"

void GraphHandler::cutBranch(std::string pathMST, std::string pathCut) {
    std::set<Edge> vecEdgeToRemove;                                           // store the edges to remove

    // get the property map for vertex indices
    IndexMap index = get(vertex_index, g);

    // weight map
    property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);

    // store MST in a new graph
    // Graph g(edge_array, edge_array + num_edges, weights, num_nodes);
    // but there's no weight according to spanning_tree (because I am lazy and do not want to create one)
    // so have to add edges one-by-one
    Graph gMST;
    for (std::vector < Edge >::iterator ei = spanning_tree.begin();
         ei != spanning_tree.end(); ++ei) {
//        std::cout << source(*ei, g) << " <--> " << target(*ei, g)
//                  << " with weight of " << weight[*ei]
//                  << std::endl;
        add_edge(source(*ei, g), target(*ei, g), weight[*ei], gMST);        // copy the weight to the new graph for MST
        gMST[source(*ei, g)].x_pixel = g[source(*ei, g)].x_pixel;
        gMST[source(*ei, g)].y_pixel = g[source(*ei, g)].y_pixel;           // also copy the pixel coordinate
        gMST[target(*ei, g)].x_pixel = g[target(*ei, g)].x_pixel;
        gMST[target(*ei, g)].y_pixel = g[target(*ei, g)].y_pixel;
    }
    //this->saveDotImage(gMST, "bin/gMST.dot");
    this->visualizeGraph(gMST, pathMST);

#ifdef DEBUG_PRINT_VERTICES_AND_OUT_DEGREE
    std::cout << "vertices(g) = \n";
#endif
    typename graph_traits<Graph>::vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(gMST); vi != vi_end; ++vi) {            // for each vertex in vertices(g)
#ifdef DEBUG_PRINT_VERTICES_AND_OUT_DEGREE
        std::cout << "vertex " << index[*vi] << std::endl;

        cout << "out-degree: " << out_degree(*vi, gMST) << endl;
#endif
        unsigned nodeDegree = out_degree(*vi, gMST);                        // skip node with degree <= 2
        if (nodeDegree <= 2){
            continue;
        }
#ifdef DEBUG_PRINT_OUTEDGE_OF_THE_VERTEX
        std::cout << "out-edges: ";
                cout << "e2 =" << printEdge(e2, gMST) << std::endl;
#endif
        typename graph_traits<Graph>::out_edge_iterator out_i, out_end;     // look for its out edges
        typename graph_traits<Graph>::edge_descriptor e;
        typename graph_traits<Graph>::edge_descriptor e1 = *(out_edges(*vi, gMST).first);
        typename graph_traits<Graph>::edge_descriptor e2 = *++(out_edges(*vi, gMST).first);
#ifdef DEBUG_PRINT_OUTEDGE_OF_THE_VERTEX
        cout << "e1 =" << printEdge(e1, gMST) << std::endl;
        cout << "e2 =" << printEdge(e2, gMST) << std::endl;
#endif
        for (tie(out_i, out_end) = out_edges(*vi, gMST);                    // for each edge e in out-edges of *vi
             out_i != out_end; ++out_i) {
            e = *out_i;
            Vertex src = source(e, gMST), targ = target(e, gMST);
#ifdef DEBUG_PRINT_OUTEDGE_OF_THE_VERTEX
            std::cout << "(" << index[src] << ","                           // index of source is always *vi
                      << index[targ] << ";" << weight[e] << ") ";           // the other node of the edge, and the weight
#endif
            if (weight[e] < weight[e1]){                                    // e1, e2 hold the two edges with lowest weight
                e1 = e;
            }
            else if (weight[e] < weight[e2]){
                e2 = e;
            }
        }
#ifdef DEBUG_PRINT_OUTEDGE_OF_THE_VERTEX
        std::cout << std::endl;
#endif
        for (tie(out_i, out_end) = out_edges(*vi, gMST);                    // for each edge e in out-edges of *vi
             out_i != out_end; ++out_i) {
            e = *out_i;
            if (e != e1 && e != e2){
                vecEdgeToRemove.insert(e);                               // push any edges other than e1 and e2
#ifdef DEBUG_PRINT_PUSH_EDGE_TO_REMOVE
                cout << "push edge: " << printEdge(e, gMST) << std::endl;
#endif
            }
        }

    }

    for (auto e : vecEdgeToRemove){                                         // remove the edges
#ifdef DEBUG_PRINT_REMOVED_EDGES
        cout << "remove edge: " << printEdge(e, gMST) << std::endl;
#endif
        remove_edge(e, gMST);                                               // cause segmentation fault if remove non-existing edge

    }

    //this->saveDotImage(gMST, "bin/gCut.dot");
    this->visualizeGraph(gMST, pathCut);

    // https://www.boost.org/doc/libs/1_65_0/libs/graph/example/connected_components.cpp
    std::vector<int> component(num_vertices(gMST));
    int num = connected_components(gMST, &component[0]);

#ifdef DEBUG_PRINT_CONNECTED_COMPONENTS
    std::vector<int>::size_type i;
    cout << "Total number of components: " << num << endl;
    for (i = 0; i != component.size(); ++i)
        cout << "Vertex " << i <<" is in component " << component[i] << endl;
    cout << endl;
#endif

}

void GraphHandler::saveDotImage(const Graph & G, std::string path) {
    // dot -Tpng bin/kruskal-eg.dot > bin/kruskal-eg.png
    std::ofstream fout(path);
    fout << "graph A {\n"
         << " rankdir=LR\n"
         << " size=\"10,10\"\n"
         << " ratio=\"filled\"\n"
         << " edge[style=\"bold\"]\n" << " node[shape=\"circle\"]\n";
    graph_traits<Graph>::edge_iterator eiter, eiter_end;
    for (boost::tie(eiter, eiter_end) = edges(G); eiter != eiter_end; ++eiter) {
        fout << source(*eiter, G) << " -- " << target(*eiter, G);
        if (std::find(spanning_tree.begin(), spanning_tree.end(), *eiter)
            != spanning_tree.end())
            fout << "[color=\"black\", label=\"" << get(edge_weight, G, *eiter)
                 << "\"];\n";
        else
            fout << "[color=\"gray\", label=\"" << get(edge_weight, G, *eiter)
                 << "\"];\n";
    }
    fout << "}\n";
}

std::string GraphHandler::printEdge(Edge e, const Graph & gMST) {
    std::stringstream ss;
    IndexMap index = get(vertex_index, this->g);
    property_map < Graph, edge_weight_t >::type weight = get(edge_weight, this->g); // compilation error if feed the second parameter with gMST
    Vertex src = source(e, gMST), targ = target(e, gMST);
    ss << "(" << index[src] << ","                           // index of source is always *vi
              << index[targ] << ";" << weight[e] << ") ";           // the other node of the edge, and the weight
    return ss.str();

}

int GraphHandler::sxyrIndex2vIndex(int s, int x, int y, int r) {
    return s * 1200 + x * 60 + y * 3 + r;
}

vector<int> GraphHandler::vIndex2sxyrIndex(int v) {
    int r = v % 3;
    v /= 3;
    int y = v % 20;
    v /= 20;
    int x = v % 20;
    int s = v / 20;
    vector<int> ret;
    ret.push_back(s);
    ret.push_back(x);
    ret.push_back(y);
    ret.push_back(r);
    return ret;
}

void GraphHandler::visualizeGraph(const Graph &G, std::string path) {
    Mat dst(GRAPH_DISPLAY_PATCH_SIZE * 20, GRAPH_DISPLAY_PATCH_SIZE * 20, CV_8UC3, Scalar(0));
    IndexMap index = get(vertex_index, g);

    typename graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei){               // for all edges in the graph
        Vertex src = source(*ei, G), targ = target(*ei, G);
        vector<int> sxyr1 = this->vIndex2sxyrIndex(index[src]);
        Point pt1;
//        pt1.x = int((float(sxyr1[0]) / 2 + sxyr1[1] + float(sxyr1[3]) / 3) * GRAPH_DISPLAY_PATCH_SIZE);
//        pt1.y = int((float(sxyr1[0]) / 2 + sxyr1[2] + float(sxyr1[3]) / 3) * GRAPH_DISPLAY_PATCH_SIZE);
        vector<int> sxyr2 = this->vIndex2sxyrIndex(index[targ]);
        Point pt2;
//        pt2.x = int((float(sxyr2[0]) / 2 + sxyr2[1] + float(sxyr2[3]) / 3) * GRAPH_DISPLAY_PATCH_SIZE);
//        pt2.y = int((float(sxyr2[0]) / 2 + sxyr2[2] + float(sxyr2[3]) / 3) * GRAPH_DISPLAY_PATCH_SIZE);
//        pt1.x = G[src].x_pixel * 4;
//        pt1.y = G[src].y_pixel * 4;
//        pt2.x = G[targ].x_pixel * 4;
//        pt2.y = G[targ].y_pixel * 4;
        pt1.x = G[src].x_pixel;
        pt1.y = G[src].y_pixel;
        pt2.x = G[targ].x_pixel;
        pt2.y = G[targ].y_pixel;

//        std::stringstream ss;
        circle(dst, pt1, 4, Scalar(255, 255, 0));
//        ss << "(" << sxyr1[0] << " " << sxyr1[1] << " " << sxyr1[2] << " " << sxyr1[3] << ")";
//        putText(dst, ss.str(), pt1, FONT_HERSHEY_PLAIN, 2.0, Scalar(255, 102, 0));
//
//        ss.str(std::string()); // clear the buffer
        circle(dst, pt2, 4, Scalar(255, 255, 0));
//        ss << "(" << sxyr2[0] << " " << sxyr2[1] << " " << sxyr2[2] << " " << sxyr2[3] << ")";
//        putText(dst, ss.str(), pt2, FONT_HERSHEY_PLAIN, 2.0, Scalar(255, 102, 0));

        line(dst, pt1, pt2, Scalar(0, 255, 0));


    }


    imwrite(path, dst);
    //imshow(winNameGraph, dst);
}
