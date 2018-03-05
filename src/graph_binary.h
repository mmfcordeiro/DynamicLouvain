// File: graph_binary.h
// -- graph handling header file
//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#ifndef GRAPH_H
#define GRAPH_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;

class Graph {
public:
    unsigned int nb_nodes;
    unsigned long nb_links;
    double total_weight;  

    unordered_map<unsigned int, unsigned long> degrees;
    unordered_multimap<unsigned int, unsigned int> links;
    unordered_multimap<unsigned int, float> weights;

    Graph();

    virtual ~Graph()
    {
        degrees.max_load_factor (9999999999);
        links.max_load_factor (9999999999);
        weights.max_load_factor (9999999999);

        degrees.clear();
        links.clear();
        weights.clear();
    };

    // binary file format is
    // 4 bytes for the number of nodes in the graph
    // 8*(nb_nodes) bytes for the cumulative degree for each node:
    //    deg(0)=degrees[0]
    //    deg(k)=degrees[k]-degrees[k-1]
    // 4*(sum_degrees) bytes for the links
    // IF WEIGHTED 4*(sum_degrees) bytes for the weights in a separate file
    Graph(char *filename, char *filename_w, int type);

    Graph(int nb_nodes, int nb_links, double total_weight, unordered_map<unsigned int, unsigned long>& degrees, unordered_multimap<unsigned int, unsigned int>& links, unordered_multimap<unsigned int, float>& weights);

    void display(void);
    void display_reverse(void);
    void display_binary(char *outfile);
    bool check_symmetry();


    // return the number of neighbors (degree) of the node
    inline unsigned int nb_neighbors(unsigned int node);

    // return the number of self loops of the node
    inline double nb_selfloops(unsigned int node);

    // return the weighted degree of the node
    inline double weighted_degree(unsigned int node);

    // return pointers to the first neighbor and first weight of the node
    inline pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, float>::iterator > neighbors(unsigned int node);

    void add(unsigned int src, unsigned int dest, float weight);
    void remove(unsigned int src, unsigned int dest, float weight);
    bool exist(unsigned int node);
};


inline unsigned int
    Graph::nb_neighbors(unsigned int node) {
        //assert(node>=0 && node<nb_nodes);

        unordered_map<unsigned int, unsigned long>::iterator degrees_it = degrees.find(node);

        if (degrees_it != degrees.end())
        {
            return degrees_it->second;
        }
		else
		{
			return 0;
		}
        //else
        //{
        //    long current_degree = degrees_it->second;
        //    degrees_it--;
        //    long previous_degree = degrees_it->second;
        //    return current_degree - previous_degree;
        //}
}

inline double
    Graph::nb_selfloops(unsigned int node) {
        //assert(node>=0 && node<nb_nodes);

        pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, float>::iterator > p = neighbors(node);
        for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
            if (p.first->second==node) {
                if (weights.size() != 0 && p.second != weights.end())
                    return (double)(p.second->second);
                else 
                    return 1.;
            }
            p.first++;
            if (weights.size() != 0 && p.second != weights.end())
                p.second++;
        }
        return 0.;
}

inline double
    Graph::weighted_degree(unsigned int node) {
        //assert(node>=0 && node<nb_nodes);

        if (weights.size()==0)
            return (double)nb_neighbors(node);
        else {
            pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, float>::iterator > p = neighbors(node);
            double res = 0;
            int a = nb_neighbors(node);
            for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
                if(p.second == weights.end()) {
                    return res;
                }

                res += (double)(p.second->second);
                p.second++;
            }
            return res;
        }
}

inline pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, float>::iterator >
    Graph::neighbors(unsigned int node) {
        //assert(node>=0 && node<nb_nodes);

        if (node==0)
            return make_pair(links.begin(), weights.begin());
        else if (weights.size()!=0)
            return make_pair(links.find(node), weights.find(node));
        else
            return make_pair(links.find(node), weights.begin());
}


#endif // GRAPH_H
