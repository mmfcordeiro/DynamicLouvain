// File: graph_binary.cpp
// -- graph handling source
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

# if !defined(WIN32)
#include <sys/mman.h>
# endif

#include <fstream>
#include "graph_binary.h"
#include "math.h"

Graph::Graph() {
    nb_nodes     = 0;
    nb_links     = 0;
    total_weight = 0;
}

Graph::Graph(char *filename, char *filename_w, int type) {
    ifstream finput;
    finput.open(filename,fstream::in | fstream::binary);

    nb_nodes     = 0;
    nb_links     = 0;
    total_weight = 0;

    // Read number of nodes on 4 bytes
    finput.read((char *)&nb_nodes, sizeof(unsigned int));
    assert(finput.rdstate() == ios::goodbit);

    // Read cumulative degree sequence: 8 bytes for each node
    // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
    degrees.rehash(nb_nodes);
    for (int i = 0; i < nb_nodes; i++)
    {
        unsigned int src;
        finput.read((char *)&src, sizeof(unsigned int));

        unsigned long degree;
        finput.read((char *)&degree, sizeof(unsigned int));

        nb_links += degree;

        degrees.insert(pair<unsigned int, unsigned long>(src, degree));
    }

    links.rehash(nb_links);
    for (int i = 0; i < nb_links; i++)
    {
        unsigned int src;
        finput.read((char *)&src, sizeof(unsigned int));

        unsigned long dest;
        finput.read((char *)&dest, sizeof(long));

        links.insert(pair<unsigned int, unsigned long>(src, dest));
    }

    // IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
    //weights.resize(0);
    total_weight=0;
    if (type==WEIGHTED) {
        ifstream finput_w;
        finput_w.open(filename_w,fstream::in | fstream::binary);
        weights.rehash(nb_links);
        for (int i = 0; i < nb_links; i++)
        {
            unsigned int src;
            finput.read((char *)&src, sizeof(unsigned int));

            long weight;
            finput.read((char *)&weight, sizeof(long));

            weights.insert(pair<unsigned int, float>(src, weight));
        } 
    }    

    // Compute total weight
    for(unordered_multimap<unsigned int, unsigned int>::iterator it = links.begin(), end = links.end(); it != end; it = links.upper_bound(it->first)) {
        total_weight += (double)weighted_degree(it->first);
        //total_weight += (double)nb_selfloops(it->first);
    }
}

Graph::Graph(int n, int m, double t, unordered_map<unsigned int, unsigned long>& d, unordered_multimap<unsigned int, unsigned int>& l, unordered_multimap<unsigned int, float>& w) {
    nb_nodes     = n;
    nb_links     = m;
    total_weight = t;

    degrees      = d;
    links        = l;
    weights      = w;
}

void Graph::add(unsigned int src, unsigned int dest, float weight)
{
    if (!exist(src))
    {
        nb_nodes++;
    }

    if (!exist(dest))
    {
        nb_nodes++;
    }

    nb_links+=2;
    degrees[src]++;
    degrees[dest]++;

    links.insert(pair<unsigned int, unsigned long>(src, dest));
    links.insert(pair<unsigned int, unsigned long>(dest, src));
    
    if (weight != 0) {
        weights.insert(pair<unsigned int, float>(src, weight));
    }    

    // Update total weight
    //total_weight += (double)weighted_degree(src);
    // Compute total weight
    total_weight = 0;
    for(unordered_multimap<unsigned int, unsigned int>::iterator it = links.begin(), end = links.end(); it != end; it = links.upper_bound(it->first)) {
        total_weight += (double)weighted_degree(it->first);
        //total_weight += (double)nb_selfloops(it->first);
    }
}

void Graph::remove(unsigned int src, unsigned int dest, float weight)
{

	bool found = false;
    std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair;
    
    iterpair = links.equal_range(src);
    for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
    {
        unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
        if (erase_iter->second == dest)
        {
            links.erase(erase_iter);
			found = true;
			nb_links-=1;
			degrees[src]--;

			if (degrees[src] == 0)
			{
				degrees.erase(degrees.find(src));
				nb_nodes--;
			}
        }
    }

    iterpair = links.equal_range(dest);
    for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
    {
        unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
        if (erase_iter->second == src)
        {
            links.erase(erase_iter);
			found = true;
			nb_links-=1;
			degrees[dest]--;

			if (degrees[dest] == 0)
			{
				degrees.erase(degrees.find(dest));
				nb_nodes--;
			}
        }
    }

	if (found)
	{
		if (weight != 0) {        
			std::pair<unordered_multimap<unsigned int, float>::iterator, unordered_multimap<unsigned int, float>::iterator> weight_iterpair;

			weight_iterpair = weights.equal_range(src);
			for (unordered_multimap<unsigned int, float>::iterator iter = weight_iterpair.first; iter != weight_iterpair.second;)
			{
				unordered_multimap<unsigned int, float>::iterator erase_iter = iter++;
				if (erase_iter->second == weight)
				{
					weights.erase(erase_iter);
				}
			}
		}    

		// Update total weight
		//total_weight += (double)weighted_degree(src);
		// Compute total weight
		total_weight = 0;
		for(unordered_multimap<unsigned int, unsigned int>::iterator it = links.begin(), end = links.end(); it != end; it = links.upper_bound(it->first)) {
			total_weight += (double)weighted_degree(it->first);
			//total_weight += (double)nb_selfloops(it->first);
		}
	}
}

bool Graph::exist(unsigned int node)
{
    return degrees.find(node) != degrees.end();
}

void
    Graph::display() {
        /*  for (unsigned int node=0 ; node<nb_nodes ; node++) {
        pair<vector<unsigned int>::iterator, vector<float>::iterator > p = neighbors(node);
        for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
        if (node<=*(p.first+i)) {
        if (weights.size()!=0)
        cout << node << " " << *(p.first+i) << " " << *(p.second+i) << endl;
        else
        cout << node << " " << *(p.first+i) << endl;
        }
        }   
        }*/
        for (unsigned int node=0 ; node<nb_nodes ; node++) {
            pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, float>::iterator > p = neighbors(node);
            cout << node << ":" ;
            for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
                if (true) {
                    if (weights.size()!=0)
                        cout << " (" << p.first->second << " " << p.second->second << ")";
                    else
                        cout << " " << p.first->second;
                }
                p.first++;

                if (weights.size() != 0)
                    p.second++;
            }
            cout << endl;
        }
}

void
    Graph::display_reverse() {
        for (unsigned int node=0 ; node<nb_nodes ; node++) {
            pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, float>::iterator > p = neighbors(node);
            for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
                if (node> p.first->second) {
                    if (weights.size()!=0)
                        cout << p.first->second << " " << node << " " << p.second->second << endl;
                    else
                        cout << p.first->second << " " << node << endl;
                }
                p.first++;

                if (weights.size() != 0)
                    p.second++;
            }   
        }
}


bool
    Graph::check_symmetry() {
        int error=0;
        for (unsigned int node=0 ; node<nb_nodes ; node++) {
            pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, float>::iterator > p = neighbors(node);
            for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
                unsigned int neigh = p.first->second;
                float weight = p.second->second;

                pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, float>::iterator > p_neigh = neighbors(neigh);
                for (unsigned int j=0 ; j<nb_neighbors(neigh) ; j++) {
                    unsigned int neigh_neigh = p_neigh.first->second;
                    float neigh_weight = p_neigh.second->second;

                    if (node==neigh_neigh && weight!=neigh_weight) {
                        cout << node << " " << neigh << " " << weight << " " << neigh_weight << endl;
                        if (error++==10)
                            exit(0);
                    }

                    p_neigh.first++;
                    p_neigh.second++;
                }

                p.first++;
                
                if (weights.size() != 0)
                    p.second++;
            }
        }
        return (error==0);
}


void
    Graph::display_binary(char *outfile) {
        //ofstream foutput;
        //foutput.open(outfile ,fstream::out | fstream::binary);

        //foutput.write((char *)(&nb_nodes),sizeof(unsigned int));
        //foutput.write((char *)(&degrees[0]),sizeof(pair<unsigned int, unsigned long>)*nb_nodes);
        //foutput.write((char *)(&links[0]),sizeof(pair<const unsigned int, unsigned int>)*nb_links);
}
