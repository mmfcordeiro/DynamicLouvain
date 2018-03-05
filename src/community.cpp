// File: community.h
// -- community detection source file
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

#include "community.h"
#include <set>

using namespace std;

Community::Community(char * filename, char * filename_w, int type, int nbp, double minm)
{
    g = new Graph(filename, filename_w, type);
    neigh_last=0;
    size = 0;

    n2c.rehash(g -> nb_nodes);
    c2n.rehash(g -> nb_nodes);
    tot.rehash(g -> nb_nodes);
    in.rehash(g -> nb_nodes);
    for(unordered_multimap<unsigned int, unsigned int>::iterator it = g -> links.begin(), end = g -> links.end(); it != end; it = g -> links.upper_bound(it->first)) {
        neigh_weight[it->first] = -1;

        n2c[it->first] = it->first;
        c2n.insert(make_pair(it->first, it->first));
        tot[it->first] = g -> weighted_degree(it->first);
        in[it->first] = g -> nb_selfloops(it->first);
        size++;
    }

    nb_pass = nbp;
    min_modularity = minm;
}

Community::Community(Graph* gc, int nbp, double minm)
{
    g = gc;
    neigh_last=0;
    size = 0;

    long degrees =  g -> degrees.size();
    n2c.rehash(degrees);
    c2n.rehash(degrees);
    tot.rehash(degrees);
    in.rehash(degrees);
    for(unordered_multimap<unsigned int, unsigned int>::iterator it = g -> links.begin(), end = g -> links.end(); it != end; it = g -> links.upper_bound(it->first)) {
        neigh_weight[it->first] = -1;

        n2c[it->first] = it->first;
        c2n.insert(make_pair(it->first, it->first));
        tot[it->first] = g -> weighted_degree(it->first);
        in[it->first] = g -> nb_selfloops(it->first);
        size++;
    }

    nb_pass = nbp;
    min_modularity = minm;
}


Community::Community(Graph* gc, Community& c, int nbp, double minm) 
{
    g = gc;
    neigh_last=c.neigh_last;
    size = c.size;

    n2c = c.n2c;
    c2n = c.c2n;
    tot = c.tot;
    in = c.in;

    nb_pass = nbp;
    min_modularity = minm;
}

void Community::add(unsigned int src, unsigned int dest, float weight, std::set<unsigned int>& affected_nodes)
{   
    unordered_map<unsigned int, int>::iterator n2c_src_it = n2c.find(src);
    unordered_map<unsigned int, int>::iterator n2c_dest_it = n2c.find(dest);

    if (n2c_src_it != n2c.end() && n2c_dest_it != n2c.end() &&
        n2c_src_it->second == n2c_dest_it->second)
    {
        //inner community edges
        tot[n2c_src_it->second] ++;
        in[n2c_src_it->second] ++;
        tot[n2c_dest_it->second] ++;
        in[n2c_dest_it->second] ++;

        g -> add(src, dest, weight);

        std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair;
        unordered_multimap<unsigned int, unsigned int>::iterator it;

        iterpair = c2n.equal_range(n2c_src_it->second);
        // src node is a community
        if (iterpair.first != iterpair.second)
        {
            it = iterpair.first;
            for (; it != iterpair.second; ++it) {

                affected_nodes.insert(it->first);
                affected_nodes.insert(it->second);
            }
        }
        else
        {
            affected_nodes.insert(n2c_src_it->second);
            affected_nodes.insert(n2c_dest_it->second);
        }

    }
    else
    {
        bool added = false;
        if (n2c_src_it == n2c.end() || n2c_dest_it == n2c.end())
        {
            g -> add(src, dest, weight);
            added = true;
        }

        vector<unsigned int> pairs;

        std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair;
        unordered_multimap<unsigned int, unsigned int>::iterator it;

        iterpair = c2n.equal_range(src);
        // src node is a community
        if (iterpair.first != iterpair.second)
        {
            it = iterpair.first;
            for (; it != iterpair.second; ++it) { 

                n2c[it->second] = it->second;
                pairs.push_back(it->second);

                tot[it->second] = g -> weighted_degree(it->second);
                double selfloops = g -> nb_selfloops(it->second);
                in[it->second] = selfloops;
                in[it->first] -= selfloops;

                affected_nodes.insert(it->first);
                affected_nodes.insert(it->second);
            }
        }
        else
        {
            //src node is not a community
            n2c_src_it = n2c.find(src);
            if (n2c_src_it == n2c.end())
            {
                // src is a new node
                n2c[src] = src;
                c2n.insert(make_pair(src, src));
                tot[src] = g -> weighted_degree(src);
                in[src] = g -> nb_selfloops(src);
                size++;

                affected_nodes.insert(src);
            }
            else
            {
                iterpair = c2n.equal_range(n2c[src]);

                for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
                {
                    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;

                    double degree = g -> weighted_degree(erase_iter->second);
                    tot[erase_iter->first] -= degree;
                    tot[erase_iter->second] = degree;
                    double selfloops = g -> nb_selfloops(erase_iter->second);
                    in[erase_iter->first] -= selfloops;
                    in[erase_iter->second] = selfloops;

                    affected_nodes.insert(erase_iter->first);
                    affected_nodes.insert(erase_iter->second);

                    n2c[erase_iter->second] = erase_iter->second;
                    c2n.insert(make_pair(erase_iter->second, erase_iter->second));

                    c2n.erase(erase_iter);
                }

                //double degree = g -> weighted_degree(src);
                //tot[src] = degree;
                //tot[n2c_src_it->second] -= degree;
                //double selfloops = g -> nb_selfloops(src);
                //in[src] = selfloops;
                //in[n2c_src_it->second] -= selfloops + 2;

                //affected_nodes.insert(src);
                //affected_nodes.insert(n2c_src_it->second);

                //iterpair = c2n.equal_range(n2c_src_it->second);

                //for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
                //{
                //    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                //    if (erase_iter->second == src)
                //    {
                //        c2n.erase(erase_iter);
                //    }
                //}

                //c2n.insert(make_pair(src, src)); 

                //n2c_src_it->second = src;
            }        
        }       

        iterpair = c2n.equal_range(dest);
        // dest node is a community
        if (iterpair.first != iterpair.second)
        {
            it = iterpair.first;
            for (; it != iterpair.second; ++it) { 

                n2c[it->second] = it->second;
                pairs.push_back(it->second);

                tot[it->second] = g -> weighted_degree(it->second);
                double selfloops = g -> nb_selfloops(it->second);
                in[it->second] = selfloops;
                in[it->first] -= selfloops;

                affected_nodes.insert(it->first);
                affected_nodes.insert(it->second);
            }
        }
        else
        {
            //dest node is not a community
            n2c_dest_it = n2c.find(dest);
            if (n2c_dest_it == n2c.end())
            {
                //dest in a new node
                n2c[dest] = dest;
                c2n.insert(make_pair(dest, dest));
                tot[dest] = g -> weighted_degree(dest);
                in[dest] = g -> nb_selfloops(dest);
                size++;

                affected_nodes.insert(dest);
            }
            else
            {
                iterpair = c2n.equal_range(n2c[dest]);

                for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
                {
                    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;

                    double degree = g -> weighted_degree(erase_iter->second);
                    tot[erase_iter->first] -= degree;
                    tot[erase_iter->second] = degree;
                    double selfloops = g -> nb_selfloops(erase_iter->second);
                    in[erase_iter->first] -= selfloops;
                    in[erase_iter->second] = selfloops;

                    affected_nodes.insert(erase_iter->first);
                    affected_nodes.insert(erase_iter->second);

                    n2c[erase_iter->second] = erase_iter->second;
                    c2n.insert(make_pair(erase_iter->second, erase_iter->second));

                    c2n.erase(erase_iter);
                }

                //double degree = g -> weighted_degree(dest);
                //tot[dest] = degree;
                //tot[n2c_dest_it->second] -= degree;
                //double selfloops = g -> nb_selfloops(dest);
                //in[dest] = selfloops;
                //in[n2c_dest_it->second] -= selfloops + 2;

                //affected_nodes.insert(dest);
                //affected_nodes.insert(n2c_dest_it->second);

                //iterpair = c2n.equal_range(n2c_dest_it->second);

                //for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
                //{
                //    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                //    if (erase_iter->second == dest)
                //    {
                //        c2n.erase(erase_iter);
                //    }
                //}

                //c2n.insert(make_pair(dest, dest)); 

                //n2c_dest_it->second = dest;
            }
        }

        for (vector<unsigned int>::iterator it = pairs.begin(); it != pairs.end(); it++)
        {
            iterpair = c2n.equal_range(*it);

            for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
            {
                unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                c2n.erase(erase_iter);
            }

            c2n.insert(make_pair(*it, *it)); 
        } 

        if (!added)
        {
            g -> add(src, dest, weight);

            tot[src] ++;

            tot[dest] ++;

            affected_nodes.insert(src);
            affected_nodes.insert(dest);
        }
    }

    pair<double, float> xxx = print_data("xxx");
}

void Community::rem(unsigned int src, unsigned int dest, float weight, std::set<unsigned int>& affected_nodes)
{
    unordered_map<unsigned int, int>::iterator n2c_src_it = n2c.find(src);
    unordered_map<unsigned int, int>::iterator n2c_dest_it = n2c.find(dest);

    if (n2c_src_it == n2c.end() || n2c_dest_it == n2c.end())
    {
        return;
    }

	if (src == dest)
    {
        return;
    }

    if (n2c_src_it != n2c.end() && n2c_dest_it != n2c.end() &&
        n2c_src_it->second == n2c_dest_it->second)
    {
        vector<unsigned int> pairs;

        std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = c2n.equal_range(n2c_src_it->second);
        unordered_multimap<unsigned int, unsigned int>::iterator it = iterpair.first;

        for (; it != iterpair.second; ++it) { 

            n2c[it->second] = it->second;
            pairs.push_back(it->second);

            tot[it->second] = g -> weighted_degree(it->second);
            in[it->second] = g -> nb_selfloops(it->second);
            in[it->first] -= g -> nb_selfloops(it->second);

            affected_nodes.insert(it->first);
            affected_nodes.insert(it->second);
        }

        for (vector<unsigned int>::iterator it = pairs.begin(); it != pairs.end(); it++)
        {
            iterpair = c2n.equal_range(*it);

            for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
            {
                unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                c2n.erase(erase_iter);
            }

            c2n.insert(make_pair(*it, *it)); 
        }

        tot[src] --;
        tot[dest] --;

        affected_nodes.insert(src);
        affected_nodes.insert(dest);

        g -> remove(src, dest, weight);
    }
	else if (n2c_src_it != n2c.end() && n2c_dest_it != n2c.end() &&
			 n2c_src_it->second != n2c_dest_it->second)
	{
		tot[n2c_src_it->second] --;
        tot[n2c_dest_it->second] --;
		
        affected_nodes.insert(src);
        affected_nodes.insert(dest);

        g -> remove(src, dest, weight);

	}
    else
    {
        tot[n2c_src_it->first] --;
        tot[n2c_dest_it->first] --;

        g -> remove(src, dest, weight);
    }

    if (!g -> exist(src))
    {
        unordered_map<unsigned int, int>::iterator n2c_it = n2c.find(src);
		unordered_map<unsigned int, double>::iterator tot_it = tot.find(src);
		unordered_map<unsigned int, double>::iterator in_it = in.find(src);

		if (n2c_it != n2c.end()) {
			n2c.erase(n2c_it);
		}

		if (tot_it != tot.end()) {
			tot.erase(tot_it);
		}

		if (in_it != in.end()) {
			in.erase(in_it);
		}

        //affected_nodes.erase(affected_nodes.find(src));

        std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = c2n.equal_range(src);

        for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
        {
            unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
            c2n.erase(erase_iter);
        }

        size--;
    }

    if (!g -> exist(dest))
    {
        unordered_map<unsigned int, int>::iterator n2c_it = n2c.find(dest);
		unordered_map<unsigned int, double>::iterator tot_it = tot.find(dest);
		unordered_map<unsigned int, double>::iterator in_it = in.find(dest);

		if (n2c_it != n2c.end()) {
			n2c.erase(n2c_it);
		}

		if (tot_it != tot.end()) {
			tot.erase(tot_it);
		}

		if (in_it != in.end()) {
			in.erase(in_it);
		}

        //affected_nodes.erase(affected_nodes.find(dest));

        std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = c2n.equal_range(dest);

        for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
        {
            unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
            c2n.erase(erase_iter);
        }

        size--;
    }
    
}

void Community::update_community (Community* cp, std::unordered_map<unsigned int, int> new_n2c)
{
    //c2n.clear();
    //c2n.rehash(g -> nb_nodes);
    for (unordered_map<unsigned int, int>::iterator new_n2c_it = new_n2c.begin(); new_n2c_it != new_n2c.end(); new_n2c_it++)
    {
		if(new_n2c_it->first == new_n2c_it->second)
		{
			this->in[new_n2c_it->second] = cp->in[new_n2c_it->second];
			this->tot[new_n2c_it->second] = cp->tot[new_n2c_it->second];
		}
		else
		{
			this->in[new_n2c_it->first] = 0;
			this->tot[new_n2c_it->first] = 0;
		}

        unordered_map<unsigned int, int>::iterator f = cp->n2c.find(new_n2c_it->first);
        if (f == cp->n2c.end())
        {
            unsigned int old_n2c_src = new_n2c_it->first;
            int old_n2c_dest = n2c [new_n2c_it->first];

            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = this -> c2n.equal_range(old_n2c_dest);
            
            unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first;

            std::set<unsigned int> other;
            for (; iter != iterpair.second;)
            {
                if (iter->first == old_n2c_dest && iter->second == old_n2c_src)
                {
                    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                    this -> c2n.erase(erase_iter);  
                }
                else
                {
                    other.insert(iter->second);
                    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                    this -> c2n.erase(erase_iter);  
                    //iter++;
                }
            }

            n2c [new_n2c_it->first] = cp->n2c [new_n2c_it->second];
            c2n.insert(make_pair(cp->n2c [new_n2c_it->second], new_n2c_it->first));

            for (std::set<unsigned int>::iterator it = other.begin(); it != other.end(); it++)
            {
                n2c [*it] = cp->n2c [new_n2c_it->second];
                c2n.insert(make_pair(cp->n2c [new_n2c_it->second], *it));
            }
        }
        else
        {
            unsigned int old_n2c_src = new_n2c_it->first;
            int old_n2c_dest = n2c [new_n2c_it->first];

            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = this -> c2n.equal_range(old_n2c_dest);

            unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first;

            std::set<unsigned int> other;
            for (; iter != iterpair.second;)
            {
                if (iter->first == old_n2c_dest && iter->second == old_n2c_src)
                {
                    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                    this -> c2n.erase(erase_iter);  
                }
                else
                {
                    other.insert(iter->second);
                    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                    this -> c2n.erase(erase_iter);  
                    //iter++;
                }
            }

            n2c[new_n2c_it->first] = new_n2c_it->second;
            c2n.insert(make_pair(new_n2c_it->second, new_n2c_it->first));

            for (std::set<unsigned int>::iterator it = other.begin(); it != other.end(); it++)
            {
                n2c[*it] = new_n2c_it->second;
                c2n.insert(make_pair(new_n2c_it->second, *it));
            }
        }        
    }
}

void Community::update_add(Community* cp1, std::set<unsigned int>& affected_nodes)
{
    pair<double, float> xxx = print_data("xxx");

    if (affected_nodes.size() == 0)
    {
        return;
    }

    std::set<unsigned int> affected_nodes_communities;
    for (std::set<unsigned int>::iterator it = affected_nodes.begin(); it != affected_nodes.end(); it++)
    {
        affected_nodes_communities.insert(cp1 -> n2c[*it]);
        if (affected_nodes_communities.size() > 1)
        {
            break;
        }
    }

    if (affected_nodes_communities.size() > 1)
    {
        //change this method to 
        std::set<unsigned int> new_affected_nodes;

        //float degrees_z = this -> g -> degrees.max_load_factor();
        //float links_z = this -> g -> links.max_load_factor();
        //float weights_z = this -> g -> weights.max_load_factor();

        //this -> g -> degrees.max_load_factor (9999999999);
        //this -> g -> links.max_load_factor (9999999999);
        //this -> g -> weights.max_load_factor (9999999999);

        for (std::set<unsigned int>::iterator it = affected_nodes.begin(); it != affected_nodes.end(); it++)
        {
            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = this -> g -> links.equal_range(*it);

            std::pair<unordered_multimap<unsigned int, float>::iterator, unordered_multimap<unsigned int, float>::iterator> w_iterpair = this -> g -> weights.equal_range(*it);

            unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first;
            unordered_multimap<unsigned int, float>::iterator w_iter = w_iterpair.first;

            for (; iter != iterpair.second;)
            {
                if (affected_nodes.find(iter->first) != affected_nodes.end())
                {
                    if (affected_nodes.find(iter->second) == affected_nodes.end() /*&& cp1->n2c[iter->first] == cp1->n2c[iter->second]*/)
                    {
                        new_affected_nodes.insert(iter->second);
                    }

                    if (this -> g -> degrees.find(iter->first) != this -> g -> degrees.end())
                    {
                        this -> g -> degrees.erase(iter->first);
                    }

                    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                    this -> g -> links.erase(erase_iter);

                    unordered_multimap<unsigned int, float>::iterator w_erase_iter = w_iter++;
                    this -> g -> weights.erase(w_erase_iter);
                }
                else
                {
                    iter++;
                    w_iter++;  
                }
            }
        }

        for (std::set<unsigned int>::iterator it = new_affected_nodes.begin(); it != new_affected_nodes.end(); it++)
        {
            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = this -> g -> links.equal_range(*it);

            std::pair<unordered_multimap<unsigned int, float>::iterator, unordered_multimap<unsigned int, float>::iterator> w_iterpair = this -> g -> weights.equal_range(*it);

            unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first;
            unordered_multimap<unsigned int, float>::iterator w_iter = w_iterpair.first;

            for (; iter != iterpair.second;)
            {
                if (affected_nodes.find(iter->second) != affected_nodes.end() || affected_nodes.find(iter->first) != affected_nodes.end())
                {
                    if (this -> g -> degrees.find(iter->first) != this -> g -> degrees.end())
                    {
                        unsigned long degree = this -> g -> degrees [iter->first];

                        if (degree != 0)
                        {
                            this -> g -> degrees [iter->first] = degree - 1;
                        }
                    }

                    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                    this -> g -> links.erase(erase_iter);  

                    unordered_multimap<unsigned int, float>::iterator w_erase_iter = w_iter++;
                    this -> g -> weights.erase(w_erase_iter);              
                }
                else
                {
                    iter++;
                    w_iter++; 
                }
            }
        }

        for (std::set<unsigned int>::iterator it = affected_nodes.begin(); it != affected_nodes.end(); it++)
        {
            // Update communities of the affected nodes
            if (this -> c2n.find(*it) == this -> c2n.end() && cp1 -> n2c.find(*it) != cp1 -> n2c.end())
            {
                this -> c2n.insert(make_pair(cp1 -> n2c [*it], *it));
            }

            if (cp1 -> n2c.find(*it) != cp1 -> n2c.end())
            {
                this -> n2c [*it] = cp1 -> n2c [*it];

                //std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = this -> c2n.equal_range(this -> n2c [*it]);

                //double in_tot = 0.;
                //double tot_tot = 0.;
                //double neigh_weight_tot = 0.;
                //for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second; iter++)
                //{
                //    in_tot += this -> in [iter->second];
                //    this -> in [iter->second] = 0;

                //    tot_tot += this -> tot [iter->second];
                //    this -> tot [iter->second] = 0;

                //    neigh_weight_tot += this -> neigh_weight [iter->second];
                //    this -> neigh_weight[iter->second] = 0;
                //}            


                //in_tot += this -> in [*it];
                this -> in [*it] = cp1 -> in [*it];
                //cerr  << "\t\tin -> " << this -> in [*it] << "(" << in_tot << ")" << endl;

                //tot_tot += this -> tot [*it];
                this -> tot [*it] = cp1 -> tot [*it];
                //cerr << "\t\ttot -> " << this -> tot [*it] << "(" << tot_tot << ")" <<endl;

                //neigh_weight_tot += this -> neigh_weight [*it];
                this -> neigh_weight[*it] = cp1 -> neigh_weight[*it];
                //cerr << "\t\tneigh_weight -> " << this -> neigh_weight [*it]  << "(" << neigh_weight_tot << ")" << endl;

            }
            else
            {
                if (this -> n2c.find(*it) != this -> n2c.end())
                {
                    this -> in [*it] = 0;
                    this -> tot [*it] = 0;
                    this -> neigh_weight[*it] = 0;
                }
            }

            //Update corresponding graph

            //unordered_map<unsigned int, unsigned long>::iterator degree_it = cp1 -> g -> degrees.find(*it);

            //if(degree_it != cp1 -> g -> degrees.end())
            //{
            //    this -> g -> degrees [degree_it->first] = degree_it->second;
            //}

            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = cp1 -> g -> links.equal_range(*it);

            for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;iter++)
            {
                int src_link = cp1 -> n2c [iter->first];
                int dest_link = cp1 -> n2c [iter->second];

                std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> link_iterpair = this -> g -> links.equal_range(src_link);

                std::pair<unordered_multimap<unsigned int, float>::iterator, unordered_multimap<unsigned int, float>::iterator> weights_iterpair = this -> g -> weights.equal_range(src_link);

                unordered_multimap<unsigned int, unsigned int>::iterator link_iter = link_iterpair.first;
                unordered_multimap<unsigned int, float>::iterator weights_iter = weights_iterpair.first;

                bool found = false;
                for (; link_iter != link_iterpair.second; link_iter++, weights_iter++)
                {
                    if(link_iter->second == dest_link)
                    {
                        weights_iter->second += 1;
                        found = true;
                    }
                }

                if (found == false)
                {
                    this -> g -> links.insert(make_pair(src_link, dest_link));

                    //if (new_affected_nodes.find(src_link) != new_affected_nodes.end())
                    {
                        unordered_map<unsigned int, unsigned long>::iterator degree_it = this -> g -> degrees.find(src_link);

                        if(degree_it != this -> g -> degrees.end())
                        {
                            this -> g -> degrees [src_link] = degree_it->second + 1;
                        }
                        else
                        {
                            this -> g -> degrees [src_link] = 1;
                        }
                    }

                    this -> g -> weights.insert(make_pair(src_link, 1));
                }

                if (affected_nodes.find(iter->second) == affected_nodes.end())
                {
                    link_iterpair = this -> g -> links.equal_range(dest_link);

                    weights_iterpair = this -> g -> weights.equal_range(dest_link);

                    link_iter = link_iterpair.first;
                    weights_iter = weights_iterpair.first;

                    bool found = false;
                    for (; link_iter != link_iterpair.second; link_iter++, weights_iter++)
                    {
                        if(link_iter->second == src_link)
                        {
                            weights_iter->second += 1;
                            found = true;
                        }
                    }

                    if (found == false)
                    {
                        this -> g -> links.insert(make_pair(dest_link, src_link));

                        //if (new_affected_nodes.find(dest_link) != new_affected_nodes.end())
                        {
                            unordered_map<unsigned int, unsigned long>::iterator degree_it = this -> g -> degrees.find(dest_link);

                            if(degree_it != this -> g -> degrees.end())
                            {
                                this -> g -> degrees [dest_link] = degree_it->second + 1;
                            }
                            else
                            {
                                this -> g -> degrees [dest_link] = 1;
                            }
                        }

                        this -> g -> weights.insert(make_pair(dest_link, 1));
                    }
                }

            }

            //if (this -> g -> degrees.find(*it) == this -> g -> degrees.end())
            //{
            //    unordered_map<unsigned int, unsigned long>::iterator degree_it = cp1 -> g -> degrees.find(*it);

            //    if(degree_it != cp1 -> g -> degrees.end())
            //    {
            //        this -> g -> degrees [degree_it->first] = degree_it->second;
            //    }
            //}

        }

        this -> g -> nb_nodes = this -> g -> degrees.size();
        this -> g -> nb_links = this -> g -> links.size();
        this -> g -> total_weight = cp1 -> g ->total_weight;

        //this -> g -> degrees.max_load_factor (degrees_z);
        //this -> g -> links.max_load_factor (links_z);
        //this -> g -> weights.max_load_factor (weights_z);

    }
    else
    {
        this -> g -> total_weight += 2;

        if (this -> g -> weights.find(*(affected_nodes_communities.begin())) == this -> g -> weights.end())
        {
            this -> g -> weights.insert(make_pair(*(affected_nodes_communities.begin()), 2));
        }
        else
        {
            this -> g -> weights.find(*(affected_nodes_communities.begin()))->second += 2;
        }

        this -> in [*(affected_nodes_communities.begin())] += 2;
        this -> tot [*(affected_nodes_communities.begin())] += 2;
    }  

    for (std::set<unsigned int>::iterator it = affected_nodes.begin(); it != affected_nodes.end();)
    {
        if (!cp1 -> g -> exist(*it))
        {
			unordered_map<unsigned int, int>::iterator n2c_it = n2c.find(*it);
			unordered_map<unsigned int, double>::iterator tot_it = tot.find(*it);
			unordered_map<unsigned int, double>::iterator in_it = in.find(*it);

			if (n2c_it != n2c.end()) {
				n2c.erase(n2c_it);
			}

			if (tot_it != tot.end()) {
				tot.erase(tot_it);
			}

			if (in_it != in.end()) {
				in.erase(in_it);
			}

            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = c2n.equal_range(*it);

            for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
            {
                unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                c2n.erase(erase_iter);
            }

            size--;

            std::set<unsigned int>::iterator erase_it = it++;
            affected_nodes.erase(erase_it);
        }
        else
        {
            it++;
        }
    }

    xxx = print_data("xxx");
}

void Community::update_rem(Community* cp1, std::set<unsigned int>& affected_nodes)
{
    pair<double, float> xxx = print_data("xxx");

    if (affected_nodes.size() == 0)
    {
        return;
    }

    std::set<unsigned int> affected_nodes_communities;
    for (std::set<unsigned int>::iterator it = affected_nodes.begin(); it != affected_nodes.end(); it++)
    {
		if (cp1 -> n2c.find(*it) != cp1 -> n2c.end())
		{
			affected_nodes_communities.insert(cp1 -> n2c[*it]);
			if (affected_nodes_communities.size() > 1)
			{
				break;
			}
		}
		else
		{
			affected_nodes_communities.insert(*it);
		}
    }

	//if (affected_nodes_communities.size() == 2 && affected_nodes.size() == 2)
 //   {	
	//	for (std::set<unsigned int>::iterator it = affected_nodes.begin(); it != affected_nodes.end(); it++)
 //       {
	//		if (cp1 -> tot.find (*it) != cp1 -> tot.end())
	//		{
	//			this -> tot [*it] = cp1-> tot [*it];
	//		}	
	//	}

	//	//change this method to 
 //       std::set<unsigned int> new_affected_nodes;

	//	for (std::set<unsigned int>::iterator it = affected_nodes.begin(); it != affected_nodes.end(); it++)
 //       {
 //           std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = this -> g -> links.equal_range(*it);

 //           std::pair<unordered_multimap<unsigned int, float>::iterator, unordered_multimap<unsigned int, float>::iterator> w_iterpair = this -> g -> weights.equal_range(*it);

 //           unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first;
 //           unordered_multimap<unsigned int, float>::iterator w_iter = w_iterpair.first;

 //           for (; iter != iterpair.second;)
 //           {
 //               if (affected_nodes.find(iter->first) != affected_nodes.end())
 //               {
 //                   if (affected_nodes.find(iter->second) == affected_nodes.end() /*&& cp1->n2c[iter->first] == cp1->n2c[iter->second]*/)
 //                   {
 //                       new_affected_nodes.insert(iter->second);
 //                   }

 //                   if (this -> g -> degrees.find(iter->first) != this -> g -> degrees.end())
 //                   {
 //                       this -> g -> degrees.erase(iter->first);
 //                   }

 //                   unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
 //                   this -> g -> links.erase(erase_iter);

 //                   unordered_multimap<unsigned int, float>::iterator w_erase_iter = w_iter++;
 //                   this -> g -> weights.erase(w_erase_iter);
 //               }
 //               else
 //               {
 //                   iter++;
 //                   w_iter++;  
 //               }
 //           }
 //       }

	//	for (std::set<unsigned int>::iterator it = new_affected_nodes.begin(); it != new_affected_nodes.end(); it++)
 //       {
 //           std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = this -> g -> links.equal_range(*it);

 //           std::pair<unordered_multimap<unsigned int, float>::iterator, unordered_multimap<unsigned int, float>::iterator> w_iterpair = this -> g -> weights.equal_range(*it);

 //           unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first;
 //           unordered_multimap<unsigned int, float>::iterator w_iter = w_iterpair.first;

 //           for (; iter != iterpair.second;)
 //           {
 //               if (affected_nodes.find(iter->second) != affected_nodes.end() || affected_nodes.find(iter->first) != affected_nodes.end())
 //               {
 //                   if (this -> g -> degrees.find(iter->first) != this -> g -> degrees.end())
 //                   {
 //                       unsigned long degree = this -> g -> degrees [iter->first];

 //                       if (degree != 0)
 //                       {
 //                           this -> g -> degrees [iter->first] = degree - 1;
 //                       }
 //                   }

 //                   unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
 //                   this -> g -> links.erase(erase_iter);  

 //                   unordered_multimap<unsigned int, float>::iterator w_erase_iter = w_iter++;
 //                   this -> g -> weights.erase(w_erase_iter);              
 //               }
 //               else
 //               {
 //                   iter++;
 //                   w_iter++; 
 //               }
 //           }
 //       }

	//	for (std::set<unsigned int>::iterator it = affected_nodes.begin(); it != affected_nodes.end(); it++)
 //       {
 //           //Update corresponding graph

 //           //unordered_map<unsigned int, unsigned long>::iterator degree_it = cp1 -> g -> degrees.find(*it);

 //           //if(degree_it != cp1 -> g -> degrees.end())
 //           //{
 //           //    this -> g -> degrees [degree_it->first] = degree_it->second;
 //           //}

 //           std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = cp1 -> g -> links.equal_range(*it);

 //           for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;iter++)
 //           {
 //               int src_link = cp1 -> n2c [iter->first];
 //               int dest_link = cp1 -> n2c [iter->second];

 //               std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> link_iterpair = this -> g -> links.equal_range(src_link);

 //               std::pair<unordered_multimap<unsigned int, float>::iterator, unordered_multimap<unsigned int, float>::iterator> weights_iterpair = this -> g -> weights.equal_range(src_link);

 //               unordered_multimap<unsigned int, unsigned int>::iterator link_iter = link_iterpair.first;
 //               unordered_multimap<unsigned int, float>::iterator weights_iter = weights_iterpair.first;

 //               bool found = false;
 //               for (; link_iter != link_iterpair.second; link_iter++, weights_iter++)
 //               {
 //                   if(link_iter->second == dest_link)
 //                   {
 //                       weights_iter->second += 1;
 //                       found = true;
 //                   }
 //               }

 //               if (found == false)
 //               {
 //                   this -> g -> links.insert(make_pair(src_link, dest_link));

 //                   //if (new_affected_nodes.find(src_link) != new_affected_nodes.end())
 //                   {
 //                       unordered_map<unsigned int, unsigned long>::iterator degree_it = this -> g -> degrees.find(src_link);

 //                       if(degree_it != this -> g -> degrees.end())
 //                       {
 //                           this -> g -> degrees [src_link] = degree_it->second + 1;
 //                       }
 //                       else
 //                       {
 //                           this -> g -> degrees [src_link] = 1;
 //                       }
 //                   }

 //                   this -> g -> weights.insert(make_pair(src_link, 1));
 //               }

 //               if (affected_nodes.find(iter->second) == affected_nodes.end())
 //               {
 //                   link_iterpair = this -> g -> links.equal_range(dest_link);

 //                   weights_iterpair = this -> g -> weights.equal_range(dest_link);

 //                   link_iter = link_iterpair.first;
 //                   weights_iter = weights_iterpair.first;

 //                   bool found = false;
 //                   for (; link_iter != link_iterpair.second; link_iter++, weights_iter++)
 //                   {
 //                       if(link_iter->second == src_link)
 //                       {
 //                           weights_iter->second += 1;
 //                           found = true;
 //                       }
 //                   }

 //                   if (found == false)
 //                   {
 //                       this -> g -> links.insert(make_pair(dest_link, src_link));

 //                       //if (new_affected_nodes.find(dest_link) != new_affected_nodes.end())
 //                       {
 //                           unordered_map<unsigned int, unsigned long>::iterator degree_it = this -> g -> degrees.find(dest_link);

 //                           if(degree_it != this -> g -> degrees.end())
 //                           {
 //                               this -> g -> degrees [dest_link] = degree_it->second + 1;
 //                           }
 //                           else
 //                           {
 //                               this -> g -> degrees [dest_link] = 1;
 //                           }
 //                       }

 //                       this -> g -> weights.insert(make_pair(dest_link, 1));
 //                   }
 //               }

 //           }

 //           //if (this -> g -> degrees.find(*it) == this -> g -> degrees.end())
 //           //{
 //           //    unordered_map<unsigned int, unsigned long>::iterator degree_it = cp1 -> g -> degrees.find(*it);

 //           //    if(degree_it != cp1 -> g -> degrees.end())
 //           //    {
 //           //        this -> g -> degrees [degree_it->first] = degree_it->second;
 //           //    }
 //           //}

 //       }

 //       this -> g -> nb_nodes = this -> g -> degrees.size();
 //       this -> g -> nb_links = this -> g -> links.size();
 //       this -> g -> total_weight = cp1 -> g ->total_weight;

 //       //this -> g -> degrees.max_load_factor (degrees_z);
 //       //this -> g -> links.max_load_factor (links_z);
 //       //this -> g -> weights.max_load_factor (weights_z);

	//}
 //   else 
	if (affected_nodes_communities.size() > 1)
    {
        //change this method to 
        std::set<unsigned int> new_affected_nodes;

        //float degrees_z = this -> g -> degrees.max_load_factor();
        //float links_z = this -> g -> links.max_load_factor();
        //float weights_z = this -> g -> weights.max_load_factor();

        //this -> g -> degrees.max_load_factor (9999999999);
        //this -> g -> links.max_load_factor (9999999999);
        //this -> g -> weights.max_load_factor (9999999999);

        for (std::set<unsigned int>::iterator it = affected_nodes.begin(); it != affected_nodes.end(); it++)
        {
            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = this -> g -> links.equal_range(*it);

            std::pair<unordered_multimap<unsigned int, float>::iterator, unordered_multimap<unsigned int, float>::iterator> w_iterpair = this -> g -> weights.equal_range(*it);

            unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first;
            unordered_multimap<unsigned int, float>::iterator w_iter = w_iterpair.first;

            for (; iter != iterpair.second;)
            {
                if (affected_nodes.find(iter->first) != affected_nodes.end())
                {
                    if (affected_nodes.find(iter->second) == affected_nodes.end() /*&& cp1->n2c[iter->first] == cp1->n2c[iter->second]*/)
                    {
                        new_affected_nodes.insert(iter->second);
                    }

                    if (this -> g -> degrees.find(iter->first) != this -> g -> degrees.end())
                    {
                        this -> g -> degrees.erase(iter->first);
                    }

                    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                    this -> g -> links.erase(erase_iter);

                    unordered_multimap<unsigned int, float>::iterator w_erase_iter = w_iter++;
                    this -> g -> weights.erase(w_erase_iter);
                }
                else
                {
                    iter++;
                    w_iter++;  
                }
            }
        }

        for (std::set<unsigned int>::iterator it = new_affected_nodes.begin(); it != new_affected_nodes.end(); it++)
        {
            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = this -> g -> links.equal_range(*it);

            std::pair<unordered_multimap<unsigned int, float>::iterator, unordered_multimap<unsigned int, float>::iterator> w_iterpair = this -> g -> weights.equal_range(*it);

            unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first;
            unordered_multimap<unsigned int, float>::iterator w_iter = w_iterpair.first;

            for (; iter != iterpair.second;)
            {
                if (affected_nodes.find(iter->second) != affected_nodes.end() || affected_nodes.find(iter->first) != affected_nodes.end())
                {
                    if (this -> g -> degrees.find(iter->first) != this -> g -> degrees.end())
                    {
                        unsigned long degree = this -> g -> degrees [iter->first];

                        if (degree != 0)
                        {
                            this -> g -> degrees [iter->first] = degree - 1;
                        }
                    }

                    unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                    this -> g -> links.erase(erase_iter);  

                    unordered_multimap<unsigned int, float>::iterator w_erase_iter = w_iter++;
                    this -> g -> weights.erase(w_erase_iter);              
                }
                else
                {
                    iter++;
                    w_iter++; 
                }
            }
        }

        for (std::set<unsigned int>::iterator it = affected_nodes.begin(); it != affected_nodes.end(); it++)
        {
            // Update communities of the affected nodes
            if (this -> c2n.find(*it) == this -> c2n.end() && cp1 -> n2c.find(*it) != cp1 -> n2c.end())
            {
				unordered_map<unsigned int, double>::iterator in_it = cp1 -> in.find(*it);
				unordered_map<unsigned int, double>::iterator tot_it = cp1 -> tot.find(*it);
				if((in_it != cp1 -> in.end() && in_it->second != 0) || (tot_it != cp1 -> tot.end() && tot_it->second != 0))
				{
					this -> c2n.insert(make_pair(cp1 -> n2c [*it], *it));
				}
            }

            if (cp1 -> n2c.find(*it) != cp1 -> n2c.end())
            {
				unordered_map<unsigned int, double>::iterator in_it = cp1 -> in.find(*it);
				unordered_map<unsigned int, double>::iterator tot_it = cp1 -> tot.find(*it);
				if((in_it != cp1 -> in.end() && in_it->second != 0) || (tot_it != cp1 -> tot.end() && tot_it->second != 0))
				{
					this -> n2c [*it] = cp1 -> n2c [*it];
				

					//std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = this -> c2n.equal_range(this -> n2c [*it]);

					//double in_tot = 0.;
					//double tot_tot = 0.;
					//double neigh_weight_tot = 0.;
					//for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second; iter++)
					//{
					//    in_tot += this -> in [iter->second];
					//    this -> in [iter->second] = 0;

					//    tot_tot += this -> tot [iter->second];
					//    this -> tot [iter->second] = 0;

					//    neigh_weight_tot += this -> neigh_weight [iter->second];
					//    this -> neigh_weight[iter->second] = 0;
					//}            


					//in_tot += this -> in [*it];
					this -> in [*it] = cp1 -> in [*it];
					//cerr  << "\t\tin -> " << this -> in [*it] << "(" << in_tot << ")" << endl;

					//tot_tot += this -> tot [*it];
					this -> tot [*it] = cp1 -> tot [*it];
					//cerr << "\t\ttot -> " << this -> tot [*it] << "(" << tot_tot << ")" <<endl;

					//neigh_weight_tot += this -> neigh_weight [*it];
					this -> neigh_weight[*it] = cp1 -> neigh_weight[*it];
					//cerr << "\t\tneigh_weight -> " << this -> neigh_weight [*it]  << "(" << neigh_weight_tot << ")" << endl;
				}
				else if (affected_nodes_communities.size() == 2 && affected_nodes.size() == 2)
				{
					std::set<unsigned int> affected_nodes_communities;
					for (std::set<unsigned int>::iterator affected_nodes_communities_it = affected_nodes_communities.begin(); affected_nodes_communities_it != affected_nodes_communities.end(); affected_nodes_communities_it++)
					{
						this -> in [*affected_nodes_communities_it] = cp1 -> in [*affected_nodes_communities_it];
						this -> tot [*affected_nodes_communities_it] = cp1 -> tot [*affected_nodes_communities_it];
						this -> neigh_weight[*affected_nodes_communities_it] = cp1 -> neigh_weight[*affected_nodes_communities_it];
					}
				}
            }
            else
            {
                if (this -> n2c.find(*it) != this -> n2c.end())
                {
                    this -> in [*it] = 0;
                    this -> tot [*it] = 0;
                    this -> neigh_weight[*it] = 0;
                }
            }

            //Update corresponding graph

            //unordered_map<unsigned int, unsigned long>::iterator degree_it = cp1 -> g -> degrees.find(*it);

            //if(degree_it != cp1 -> g -> degrees.end())
            //{
            //    this -> g -> degrees [degree_it->first] = degree_it->second;
            //}

            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = cp1 -> g -> links.equal_range(*it);

            for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;iter++)
            {
                int src_link = cp1 -> n2c [iter->first];
                int dest_link = cp1 -> n2c [iter->second];

                std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> link_iterpair = this -> g -> links.equal_range(src_link);

                std::pair<unordered_multimap<unsigned int, float>::iterator, unordered_multimap<unsigned int, float>::iterator> weights_iterpair = this -> g -> weights.equal_range(src_link);

                unordered_multimap<unsigned int, unsigned int>::iterator link_iter = link_iterpair.first;
                unordered_multimap<unsigned int, float>::iterator weights_iter = weights_iterpair.first;

                bool found = false;
                for (; link_iter != link_iterpair.second; link_iter++, weights_iter++)
                {
                    if(link_iter->second == dest_link)
                    {
                        weights_iter->second += 1;
                        found = true;
                    }
                }

                if (found == false)
                {
                    this -> g -> links.insert(make_pair(src_link, dest_link));

                    //if (new_affected_nodes.find(src_link) != new_affected_nodes.end())
                    {
                        unordered_map<unsigned int, unsigned long>::iterator degree_it = this -> g -> degrees.find(src_link);

                        if(degree_it != this -> g -> degrees.end())
                        {
                            this -> g -> degrees [src_link] = degree_it->second + 1;
                        }
                        else
                        {
                            this -> g -> degrees [src_link] = 1;
                        }
                    }

                    this -> g -> weights.insert(make_pair(src_link, 1));
                }

                if (affected_nodes.find(iter->second) == affected_nodes.end())
                {
                    link_iterpair = this -> g -> links.equal_range(dest_link);

                    weights_iterpair = this -> g -> weights.equal_range(dest_link);

                    link_iter = link_iterpair.first;
                    weights_iter = weights_iterpair.first;

                    bool found = false;
                    for (; link_iter != link_iterpair.second; link_iter++, weights_iter++)
                    {
                        if(link_iter->second == src_link)
                        {
                            weights_iter->second += 1;
                            found = true;
                        }
                    }

                    if (found == false)
                    {
                        this -> g -> links.insert(make_pair(dest_link, src_link));

                        //if (new_affected_nodes.find(dest_link) != new_affected_nodes.end())
                        {
                            unordered_map<unsigned int, unsigned long>::iterator degree_it = this -> g -> degrees.find(dest_link);

                            if(degree_it != this -> g -> degrees.end())
                            {
                                this -> g -> degrees [dest_link] = degree_it->second + 1;
                            }
                            else
                            {
                                this -> g -> degrees [dest_link] = 1;
                            }
                        }

                        this -> g -> weights.insert(make_pair(dest_link, 1));
                    }
                }

            }

            //if (this -> g -> degrees.find(*it) == this -> g -> degrees.end())
            //{
            //    unordered_map<unsigned int, unsigned long>::iterator degree_it = cp1 -> g -> degrees.find(*it);

            //    if(degree_it != cp1 -> g -> degrees.end())
            //    {
            //        this -> g -> degrees [degree_it->first] = degree_it->second;
            //    }
            //}

        }

        this -> g -> nb_nodes = this -> g -> degrees.size();
        this -> g -> nb_links = this -> g -> links.size();
        this -> g -> total_weight = cp1 -> g ->total_weight;

        //this -> g -> degrees.max_load_factor (degrees_z);
        //this -> g -> links.max_load_factor (links_z);
        //this -> g -> weights.max_load_factor (weights_z);

    }
    else if (affected_nodes_communities.size() == 0)
    {
        this -> g -> total_weight -= 2;

        if (this -> g -> weights.find(*(affected_nodes_communities.begin())) == this -> g -> weights.end())
        {
            this -> g -> weights.insert(make_pair(*(affected_nodes_communities.begin()), 2));
        }
        else if (this -> g -> weights.find(*(affected_nodes_communities.begin())) != this -> g -> weights.end())
		{
            this -> g -> weights.find(*(affected_nodes_communities.begin()))->second -= 2;
			this -> in [*(affected_nodes_communities.begin())] -= 2;
			this -> tot [*(affected_nodes_communities.begin())] -= 2;
        }

    }  

    for (std::set<unsigned int>::iterator it = affected_nodes.begin(); it != affected_nodes.end();)
    {
        if (!cp1 -> g -> exist(*it))
        {
			unordered_map<unsigned int, int>::iterator n2c_it = n2c.find(*it);
			unordered_map<unsigned int, double>::iterator tot_it = tot.find(*it);
			unordered_map<unsigned int, double>::iterator in_it = in.find(*it);

			if (n2c_it != n2c.end()) {
				n2c.erase(n2c_it);
			}

			if (tot_it != tot.end()) {
				tot.erase(tot_it);
			}

			if (in_it != in.end()) {
				in.erase(in_it);
			}

            std::pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, unsigned int>::iterator> iterpair = c2n.equal_range(*it);

            for (unordered_multimap<unsigned int, unsigned int>::iterator iter = iterpair.first; iter != iterpair.second;)
            {
                unordered_multimap<unsigned int, unsigned int>::iterator erase_iter = iter++;
                c2n.erase(erase_iter);
            }

            size--;

            std::set<unsigned int>::iterator erase_it = it++;
            affected_nodes.erase(erase_it);
        }
        else
        {
            it++;
        }
    }

    xxx = print_data("xxx");
}

void
    Community::init_partition(char * filename) {
        ifstream finput;
        finput.open(filename,fstream::in);
        
        std::unordered_map<unsigned int, int> new_n2c;

        // read partition
        while (!finput.eof()) {
            unsigned int node, comm;
            finput >> node >> comm;

            if (finput) {
                int old_comm = n2c[node];
                neigh_comm(node);

                remove(node, old_comm, neigh_weight[old_comm]);

                unsigned int i=0;
                for ( i=0 ; i<neigh_last ; i++) {
                    unsigned int best_comm     = neigh_pos[i];
                    float best_nblinks  = neigh_weight[neigh_pos[i]];
                    if (best_comm==comm) {
                        insert(node, best_comm, best_nblinks, new_n2c, false);
                        break;
                    }
                }
                if (i==neigh_last)
                    insert(node, comm, 0, new_n2c, false);
            }
        }
        finput.close();
}

// inline void
// Community::remove(int node, int comm, double dnodecomm) {
//   assert(node>=0 && node<size);

//   tot[comm] -= g -> weighted_degree(node);
//   in[comm]  -= 2*dnodecomm + g -> nb_selfloops(node);
//   n2c[node]  = -1;
// }

// inline void
// Community::insert(int node, int comm, double dnodecomm) {
//   assert(node>=0 && node<size);

//   tot[comm] += g -> weighted_degree(node);
//   in[comm]  += 2*dnodecomm + g -> nb_selfloops(node);
//   n2c[node]=comm;
// }

void
    Community::display() {
         for(unordered_multimap<unsigned int, unsigned int>::iterator it = g -> links.begin(), end = g -> links.end(); it != end; it = g -> links.upper_bound(it->first)) {
            cerr << " " << it->first << "/" << n2c[it->first] << "/" << in[it->first] << "/" << tot[it->first];
         }
        cerr << endl;
}


double
    Community::modularity() {
        double q  = 0.;
        double m2 = (double)g -> total_weight;

        //for(unordered_multimap<unsigned int, unsigned int>::iterator it = g -> links.begin(), end = g -> links.end(); it != end; it = g -> links.upper_bound(it->first)) {
        for(unordered_multimap<unsigned int, unsigned int>::iterator it = c2n.begin(), end = c2n.end(); it != end; it = c2n.upper_bound(it->first)){
            //if (tot[it->first]>0)
                q += (double)in[it->first]/m2 - ((double)tot[it->first]/m2)*((double)tot[it->first]/m2);           
        }

        //if(q < 0)
        

        return q;
}

void
    Community::neigh_comm(unsigned int node) {
        for (unsigned int i=0 ; i<neigh_last ; i++)
            neigh_weight[neigh_pos[i]]=-1;
        neigh_last=0;

        pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, float>::iterator > p = g -> neighbors(node);

        unsigned int deg = g -> nb_neighbors(node);

        neigh_pos[0]=n2c[node];
        neigh_weight[neigh_pos[0]]=0;
        neigh_last=1;

        for (unsigned int i=0 ; i<deg ; i++) {
            unsigned int neigh        = p.first->second;

            if (n2c.find(neigh) != n2c.end())
            {
                unsigned int neigh_comm = n2c[neigh];

                double neigh_w = 1.;

                if (g -> weights.size() != 0 && p.second != g -> weights.end()) {
                    neigh_w = p.second->second;
                }

                if (neigh!=node) {
                    if (neigh_weight[neigh_comm]==-1) {
                        neigh_weight[neigh_comm]=0.;
                        neigh_pos[neigh_last++]=neigh_comm;
                    }
                    neigh_weight[neigh_comm]+=neigh_w;
                }
                p.first++;

                if (g -> weights.size() != 0 && p.second != g -> weights.end())
                    p.second++;
            }            
        }
}

void
    Community::partition2graph() {
        vector<int> renumber(size, -1);
        for (int node=0 ; node<size ; node++) {
            renumber[n2c[node]]++;
        }

        int final=0;
        for (int i=0 ; i<size ; i++)
            if (renumber[i]!=-1)
                renumber[i]=final++;


        for (int i=0 ; i<size ; i++) {
            pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, float>::iterator > p = g -> neighbors(i);

            int deg = g -> nb_neighbors(i);
            for (int j=0 ; j<deg ; j++) {
                int neigh = p.first->second;
                cout << renumber[n2c[i]] << " " << renumber[n2c[neigh]] << endl;
                p.first++;
            }
        }
}

void
    Community::display_partition() {

        for(unordered_multimap<unsigned int, unsigned int>::iterator it = g -> links.begin(), end = g -> links.end(); it != end; it = g -> links.upper_bound(it->first)) {            
            cout << it->first << " " << n2c[it->first] << endl;
        }
}


Graph*
    Community::partition2graph_binary() {
        // Renumber communities
        //set<int> renumber;
        //for(unordered_multimap<unsigned int, unsigned int>::iterator it = g -> links.begin(), end = g -> links.end(); it != end; it = g -> links.upper_bound(it->first)) {
        //    renumber.insert(n2c[it->first]);
        //}

        //int final=0;
        //for (int i=0 ; i<size ; i++)
        //    if (renumber[i]!=-1)
        //        renumber[i]=final++;

        // Compute communities
        unordered_map<unsigned int, set<int> > comm_nodes;

        for(unordered_multimap<unsigned int, unsigned int>::iterator it = g -> links.begin(), end = g -> links.end(); it != end; it = g -> links.upper_bound(it->first)) {            
            comm_nodes [n2c[it->first]].insert(it->first);
        }

        // Compute weighted graph
        Graph* g2 = new Graph();
        g2 -> nb_nodes = comm_nodes.size();
        //g2 -> degrees.resize(comm_nodes.size());
        g2 -> degrees.rehash(g2 -> nb_nodes);
        g2 -> links.rehash(g2 -> nb_nodes);
        g2 -> weights.rehash(g2 -> nb_nodes);

        size_t previous_size = 0;
        for (unordered_map<unsigned int, set<int> >::iterator comm_nodes_it = comm_nodes.begin(); comm_nodes_it != comm_nodes.end(); comm_nodes_it++) {

            unordered_map<int,float> m;
            unordered_map<int,float>::iterator it;

            for (set<int>::iterator comm_nodes_it2 = comm_nodes_it->second.begin(); comm_nodes_it2 != comm_nodes_it->second.end(); comm_nodes_it2++) {
                pair<unordered_multimap<unsigned int, unsigned int>::iterator, unordered_multimap<unsigned int, float>::iterator > p = g -> neighbors(*comm_nodes_it2);


                int deg = g -> nb_neighbors(*comm_nodes_it2);
                for (int i=0 ; i<deg ; i++) {
                    int neigh        = p.first->second;
                    int neigh_comm   = n2c[neigh];
                    
                    double neigh_weight = 1.;

                    if (g -> weights.size() != 0 && p.second != g -> weights.end()) {
                        neigh_weight = p.second->second;
                    }

                    it = m.find(neigh_comm);
                    if (it==m.end())
                        m.insert(make_pair(neigh_comm, (float)neigh_weight));
                    else
                        it->second+=neigh_weight;

                    p.first++;

                    if (g -> weights.size() != 0 && p.second != g -> weights.end())
                        p.second++;
                }


            }

            //g2 -> degrees[comm_nodes_it->first]=(comm_nodes_it==comm_nodes.begin())?m.size():previous_size+m.size();
            g2 -> degrees[comm_nodes_it->first]=m.size();

            previous_size += m.size();
            g2 -> nb_links+=m.size();


            for (it = m.begin() ; it!=m.end() ; it++) {
                g2 -> total_weight  += it->second;
                g2 -> links.insert(make_pair(comm_nodes_it->first, it->first));
                g2 -> weights.insert(make_pair(comm_nodes_it->first, it->second));
            }
        }

        return g2;
}


bool
    Community::one_level(std::unordered_map<unsigned int, int>& new_n2c) {
        bool improvement=false ;
        unsigned int nb_moves;
        int nb_pass_done = 0;
        double new_mod   = modularity();
        double cur_mod   = new_mod;

        vector<int> random_order;
		if (size > 0)
			random_order.reserve(size);
        
        //for(unordered_multimap<unsigned int, unsigned int>::iterator it = g -> links.begin(), end = g -> links.end(); it != end; it = g -> links.upper_bound(it->first)) {
        //    random_order.push_back(it->first);
        //}

        for(unordered_multimap<unsigned int, unsigned int>::iterator it = c2n.begin(), end = c2n.end(); it != end; it = c2n.upper_bound(it->first))
        {
            random_order.push_back(it->first);
        }

        for (int i=0 ; i<random_order.size()-1 ; i++) {
            int rand_pos = rand()%(random_order.size()-i)+i;
            int tmp      = random_order[i];
            random_order[i] = random_order[rand_pos];
            random_order[rand_pos] = tmp;
        }

        // repeat while 
        //   there is an improvement of modularity
        //   or there is an improvement of modularity greater than a given epsilon 
        //   or a predefined number of pass have been done
        unsigned int previous_nb_moves = 9999999999999999999;
        nb_moves = 9999999999999999999;
        do {
            cur_mod = new_mod;
            previous_nb_moves = nb_moves;
            nb_moves = 0;
            nb_pass_done++;

            // for each node: remove the node from its community and insert it in the best community
            for (int node_tmp=0 ; node_tmp<random_order.size() ; node_tmp++) {
                //      int node = node_tmp;
                int node = random_order[node_tmp];
                //int node_comm     = n2c[max (n2c[node], node)];
                int node_comm     = n2c[node];
                double w_degree = g -> weighted_degree(node);

                // computation of all neighboring communities of current node
                neigh_comm(node);
                // remove node from its current community
                //this->print_data("before remove");
                remove(node, node_comm, neigh_weight[node_comm]);
                //this->print_data("after remove");

                // compute the nearest community for node
                // default choice for future insertion is the former community
                int best_comm        = node_comm;
                double best_nblinks  = 0.;
                double best_increase = 0.;
                for (unsigned int i=0 ; i<neigh_last ; i++) {
                    double increase = modularity_gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]], w_degree);
                    if (increase>best_increase) {
                        best_comm     = neigh_pos[i];
                        best_nblinks  = neigh_weight[neigh_pos[i]];
                        best_increase = increase;
                    }
                }

                // insert node in the nearest community
                //this->print_data("before insert");
                insert(node, best_comm, best_nblinks, new_n2c, best_comm != node_comm);
                //this->print_data("after insert");

                if (best_comm!=node_comm)
                    nb_moves++;
            }

            double total_tot=0;
            double total_in=0;
            for (unordered_map<unsigned int, double>::iterator it = tot.begin(); it != tot.end(); it++) {
                total_tot+=it->second;
            }

            for (unordered_map<unsigned int, double>::iterator it = in.begin(); it != in.end(); it++) {
                total_in+=it->second;
            }


            new_mod = modularity();
            if (nb_moves>0)
                improvement=true;

        } while (nb_moves < previous_nb_moves && nb_moves>0 && new_mod-cur_mod>min_modularity);

        return improvement;
}

