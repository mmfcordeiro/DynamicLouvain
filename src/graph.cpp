// File: graph.cpp
// -- simple graph handling source file
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

#include "graph.h"

using namespace std;

//std::unordered_map<int,int> numbers_;

Graph::Graph(char *filename, int type) {
  ifstream finput;
  finput.open(filename,fstream::in);
  
  while (!finput.eof()) {
    unsigned int src, dest;
    double weight=1.;

    if (type==WEIGHTED) {
      finput >> src >> dest >> weight;
    } else {
      finput >> src >> dest;
    }
    
    if (finput) {

    /*if(numbers_.find(src) == numbers_.end())
        numbers_[src] = numbers_.size();
    src = numbers_[src];

    if(numbers_.find(dest) == numbers_.end())
        numbers_[dest] = numbers_.size();
    dest = numbers_[dest];*/

      links[src][dest] = (float)weight;

      if (src!=dest)
        links[dest][src] = (float)weight;
    }
  }

  finput.close();
}

void
Graph::renumber(int type) {
 // vector<int> linked(links.size(),-1);
 // vector<int> renum(links.size(),-1);
 // int nb=0;
 // 
 // for (unsigned int i=0 ; i<links.size() ; i++) {
 //   for (unsigned int j=0 ; j<links[i].size() ; j++) {
 //     linked[i]=1;
 //     linked[links[i][j].first]=1;
 //   }
 // }
 // 
 // for (unsigned int i=0 ; i<links.size() ; i++) {
 //   if (linked[i]==1)
 //     renum[i]=nb++;
 // }

 // for (unsigned int i=0 ; i<links.size() ; i++) {
 //   if (linked[i]==1) {
 //     for (unsigned int j=0 ; j<links[i].size() ; j++) {
	//links[i][j].first = renum[links[i][j].first];
 //     }
 //     links[renum[i]]=links[i];
 //   }
 // }
 // links.resize(nb);
}

void
Graph::clean(int type) {


  //for (unordered_map<int, unordered_map<int, float> >::iterator link_it_i = links.begin(); link_it_i != links.end(); link_it_i++) {

  //  unordered_map<int, float> m;
  //  unordered_map<int, float>::iterator it;

  //  for (unordered_map<int, float>::iterator link_it_j = link_it_i->second.begin(); link_it_j != link_it_i->second.end(); link_it_j++) {
  //    it = m.find(link_it_j->first);
  //    if (it==m.end())
	 //   m.insert(make_pair(link_it_j->first, link_it_j->second));
  //    else if (type==WEIGHTED)
  //    	it->second+=link_it_j->second;
  //  }
  //  
  //  link_it_i->second.clear();

  //  for (it = m.begin() ; it!=m.end() ; it++)
  //    link_it_i->second.insert(make_pair(it->first, it->second));
  //}
}

void
Graph::display(int type) {
  for (unordered_map<int, unordered_map<int, float> >::iterator link_it_i = links.begin(); link_it_i != links.end(); link_it_i++) {
    
    for (unordered_map<int, float>::iterator link_it_j = link_it_i->second.begin(); link_it_j != link_it_i->second.end(); link_it_j++) {

      int dest   = link_it_j->first;
      float weight = link_it_j->second;
      if (type==WEIGHTED)
	cout << link_it_i->first << " " << dest << " " << weight << endl;
      else
	cout << link_it_i->first << " " << dest << endl;
    }
  }
}

void
Graph::display_binary(char *filename, char *filename_w, int type) {
  ofstream foutput;
  foutput.open(filename, fstream::out | fstream::binary);

  unsigned int s = links.size();

  // outputs number of nodes
  foutput.write((char *)(&s),sizeof(unsigned int));

  // outputs cumulative degree sequence
  long tot=0;
  for (unordered_map<int, unordered_map<int, float> >::iterator link_it_i = links.begin(); link_it_i != links.end(); link_it_i++) {
    tot=(long)link_it_i->second.size();
    int src = link_it_i->first;
    foutput.write((char *)(&src),sizeof(int));
    foutput.write((char *)(&tot),sizeof(long));
  }

  // outputs links
  for (unordered_map<int, unordered_map<int, float> >::iterator link_it_i = links.begin(); link_it_i != links.end(); link_it_i++) {
    for (unordered_map<int, float>::iterator link_it_j = link_it_i->second.begin(); link_it_j != link_it_i->second.end(); link_it_j++) {
      int src = link_it_i->first;
      foutput.write((char *)(&src),sizeof(int));
      int dest = link_it_j->first;
      foutput.write((char *)(&dest),sizeof(int));
    }
  }
  foutput.close();

  // outputs weights in a separate file
  if (type==WEIGHTED) {
    ofstream foutput_w;
    foutput_w.open(filename_w,fstream::out | fstream::binary);
    for (unordered_map<int, unordered_map<int, float> >::iterator link_it_i = links.begin(); link_it_i != links.end(); link_it_i++) {
        for (unordered_map<int, float>::iterator link_it_j = link_it_i->second.begin(); link_it_j != link_it_i->second.end(); link_it_j++) {
            int src = link_it_i->first;
            foutput.write((char *)(&src),sizeof(int));
            float weight = link_it_j->second;
	        foutput_w.write((char *)(&weight),sizeof(float));
        }
    }
    foutput_w.close();
  }
}

