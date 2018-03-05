// File: community.h
// -- community detection header file
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

#ifndef COMMUNITY_H
#define COMMUNITY_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_map>
#include <set>

#include "graph_binary.h"

using namespace std;

class Community {
 public:
  unordered_map<unsigned int, double> neigh_weight;
  unordered_map<unsigned int, unsigned int> neigh_pos;
  unsigned int neigh_last;

  Graph* g; // network to compute communities for
  int size; // nummber of nodes in the network and size of all vectors
  unordered_map<unsigned int, int> n2c; // community to which each node belongs
  unordered_multimap<unsigned int, unsigned int> c2n; // nodes that belong to the community
  unordered_map<unsigned int, double> in,tot; // used to compute the modularity participation of each community

  // number of pass for one level computation
  // if -1, compute as many pass as needed to increase modularity
  int nb_pass;

  // a new pass is computed if the last one has generated an increase 
  // greater than min_modularity
  // if 0. even a minor increase is enough to go for one more pass
  double min_modularity;

  // constructors:
  // reads graph from file using graph constructor
  // type defined the weighted/unweighted status of the graph file
  Community (char *filename, char *filename_w, int type, int nb_pass, double min_modularity);
  // copy graph
  Community (Graph* g, int nb_pass, double min_modularity);

  Community (Graph* g, Community& c, int nb_pass, double min_modularity);

  virtual ~Community()
  {
      neigh_weight.clear();
      neigh_pos.clear();
      n2c.clear();
      c2n.clear();
      in.clear();
      tot.clear();
      delete g;
  }

  // initiliazes the partition with something else than all nodes alone
  void init_partition(char *filename_part);

  // display the community of each node
  void display();

  // remove the node from its current community with which it has dnodecomm links
  inline void remove(int node, int comm, double dnodecomm);

  // insert the node in comm with which it shares dnodecomm links
  inline void insert(int node, int comm, double dnodecomm, std::unordered_map<unsigned int, int>& new_n2c, bool community_change);

  // compute the gain of modularity if node where inserted in comm
  // given that node has dnodecomm links to comm.  The formula is:
  // [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
  // [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
  // where In(comm)    = number of half-links strictly inside comm
  //       Tot(comm)   = number of half-links inside or outside comm (sum(degrees))
  //       d(node,com) = number of links from node to comm
  //       deg(node)   = node degree
  //       m           = number of links
  inline double modularity_gain(int node, int comm, double dnodecomm, double w_degree);

  // compute the set of neighboring communities of node
  // for each community, gives the number of links from node to comm
  void neigh_comm(unsigned int node);

  // compute the modularity of the current partition
  double modularity();

  // displays the graph of communities as computed by one_level
  void partition2graph();
  // displays the current partition (with communities renumbered from 0 to k-1)
  void display_partition();

  // generates the binary graph of communities as computed by one_level
  Graph* partition2graph_binary();

  // compute communities of the graph for one level
  // return true if some nodes have been moved
  bool one_level(std::unordered_map<unsigned int, int>& new_n2c);

  void add(unsigned int src, unsigned int dest, float weight, std::set<unsigned int>& affected_nodes);
  void rem(unsigned int src, unsigned int dest, float weight, std::set<unsigned int>& affected_nodes);

  void update_add(Community* cp, std::set<unsigned int>& affected_nodes);
  void update_rem(Community* cp, std::set<unsigned int>& affected_nodes);
  void update_community (Community* cp, std::unordered_map<unsigned int, int> new_n2c);

  pair<double, float> print_data(const std::string& msg)
  {
  //    double m2 = (double)g -> total_weight;

  //    {

  //        for(unordered_multimap<unsigned int, unsigned int>::iterator it = c2n.begin(), end = c2n.end(); it != end; it = c2n.upper_bound(it->first))
  //        {
  //            bool found = false;
  //            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = c2n.equal_range(it->first);

  //            for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second; iter++)
  //            {
  //               if (iter->first == iter->second)
  //               {
  //                   found = true;
  //               }
  //            }

  //            if(found == false)
  //            {
  //               for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second; iter++)
  //               {
  //                   cerr << "\t\t check: \tc2n (" << iter->first << ") \t = " << iter->second << endl;
  //               }

  //               int a = 10;
  //            }
  //        }


  //        for(unordered_multimap<unsigned int, unsigned int>::iterator iter = c2n.begin(); iter != c2n.end(); iter++)
  //        {
  //            cerr << "\t\t c2n (" << iter->first << ") \t = " << iter->second << endl;
  //        }

  //        for(unordered_multimap<unsigned int, unsigned int>::iterator iter = c2n.begin(); iter != c2n.end(); iter++)
  //        {
  //            cerr << "\t\t c2n (" << iter->first << ") \t = " << iter->second << endl;
  //        }


  //        cerr << endl;

  //        for(unordered_map<unsigned int, int>::iterator iter = n2c.begin(); iter != n2c.end(); iter++)
  //        {
  //            cerr << "\t\t n2c (" << iter->first << ") \t = " << iter->second << endl;
  //        }
  //    }

  //    double a = 0.;
  //    for (unordered_map<unsigned int, double>::iterator it = tot.begin(); it != tot.end(); it++)
  //    {
  //        a += it->second;

  //        //if (it->second != 0)
  //        {
  //            double aa = 0.;

  //            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = c2n.equal_range(it->first);
  //            
  //            for(unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second; iter++)
  //            {
  //                std::pair<unordered_multimap<unsigned int, float>::iterator, unordered_multimap<unsigned int, float>::iterator> w_iterpair = this -> g -> weights.equal_range(iter->second);

  //                for (unordered_multimap<unsigned int, float>::iterator w_iter = w_iterpair.first; w_iter != w_iterpair.second; w_iter++)
  //                {
  //                    cerr << "\t\t weight (" << w_iter->first << ") = " << w_iter->second << endl;
  //                    aa += w_iter->second;
  //                }
  //            }         


  //            if (aa != it->second && aa != 0)
  //            {
  //                cerr << "\t\t tot [" << it->first << "] = " << it->second << " with problems: SUM of this -> g -> weights (" << aa << ")" << endl;
  //            }
  //            else
  //            {
  //                cerr << "\t\t tot [" << it->first << "] = " << it->second << " ok: SUM of this -> g -> weights (" << aa << ")" << endl;
  //            }
  //        }         

  //    }

  //    cerr << "\t\t---" << msg.c_str() << "-------" << endl;
  //    cerr << "\t\ttot: " << a << endl;

  //    double b = 0.;
  //    for (unordered_map<unsigned int, double>::iterator it = in.begin(); it != in.end(); it++)
  //    {
  //        b += it->second;
  //    }

  //    cerr << "\t\tin: " << b << endl;

  //    float c = 0.;
  //    for (unordered_multimap<unsigned int, float>::iterator it = g->weights.begin(); it != g->weights.end(); it++)
  //    {
  //        c += it->second;
  //    }

  //    cerr << "\t\tg->weights: " << c << endl;



  //    cerr << "\t\ttotal_weight: " << m2 << endl;

  //    return pair<double, float>(a,c);
    return pair<double, float>(1,1);
  }



};

inline void
Community::remove(int node, int comm, double dnodecomm) {
  //assert(node>=0 && node<size);

  //cerr << "tot[" << comm <<"] " << tot[comm] << endl; 
  tot[comm] -= g -> weighted_degree(node);
  //cerr << "weighted_degree(" << node <<") " << g -> weighted_degree(node) << endl; 
  //cerr << "tot[" << comm <<"] " << tot[comm] << endl; 
  in[comm]  -= 2*dnodecomm + g -> nb_selfloops(node);
  n2c[node]  = -1;

  std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = c2n.equal_range(comm);

  unordered_multimap<unsigned int, unsigned int>::iterator it = iterpair.first;
  for (; it != iterpair.second; ++it) {
      if (it->second == node) { 
          c2n.erase(it);
          break;
      }
  }
}

inline void
    Community::insert(int node, int comm, double dnodecomm, std::unordered_map<unsigned int, int>& new_n2c, bool community_change) {
  //assert(node>=0 && node<size);
    
    //cerr << "tot[" << comm <<"] " << tot[comm] << endl; 
    tot[comm] += g -> weighted_degree(node);
    //cerr << "weighted_degree(" << node <<") " << g -> weighted_degree(node) << endl; 
    //cerr << "tot[" << comm <<"] " << tot[comm] << endl; 
    in[comm]  += 2*dnodecomm + g -> nb_selfloops(node);
    n2c[node]=comm;

    if (community_change == true)
    {
        new_n2c[node]=comm;
		new_n2c[comm]=comm;
    }
    c2n.insert(make_pair(comm, node)); 

    if (node > comm)
    {
        //cerr << "tot[" << node << "]: " << tot[node] << " = tot[" << comm << "]: " << tot[comm] << endl; 
        tot[node] += tot[comm];
        tot[comm] = 0;

        in[node] += in[comm];
        in[comm] = 0;
        
        std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = c2n.equal_range(comm);

        vector<unsigned int> pairs;
        unordered_multimap<unsigned int, unsigned int>::iterator it = iterpair.first;
        for (; it != iterpair.second; ++it) {        
           n2c[it->second]=node;
           if (community_change == true)
           {
              new_n2c[it->second]=node;
			  new_n2c[node]=node;
           }
           pairs.push_back(it->second);
        }

        for (vector<unsigned int>::iterator it = pairs.begin(); it != pairs.end(); it++)
        {
            c2n.insert(make_pair(node, *it)); 
        }        

        iterpair = c2n.equal_range(comm);

        c2n.erase(iterpair.first, iterpair.second);
    }
}

inline double
Community::modularity_gain(int node, int comm, double dnodecomm, double w_degree) {
  //assert(node>=0 && node<size);

  double totc = (double)tot[comm];
  double degc = (double)w_degree;
  double m2   = (double)g -> total_weight;
  double dnc  = (double)dnodecomm;
  
  return (dnc - totc*degc/m2);
}


#endif // COMMUNITY_H
