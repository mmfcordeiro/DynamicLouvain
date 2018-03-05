// File: main_community.cpp
// -- community detection, sample main file
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
// Time	    : February 2008Community
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

# if defined(WIN32)
#include <stdio.h>
#include <process.h>
# endif

#include "graph_binary.h"
#include "community.h"

using namespace std;

char *filename = NULL;
char *filename_w = NULL;
char *filename_part = NULL;
int type       = UNWEIGHTED;
int nb_pass    = 0;
double precision = 0.000001;
int display_level = -2;
int k1 = 16;

std::string directory = ".";

bool verbose = false;

void
usage(char *prog_name, const char *more) {
  cerr << more;
  cerr << "usage: " << prog_name << " input_file [-w weight_file] [-p part_file] [-q epsilon] [-l display_level] [-v] [-h]" << endl << endl;
  cerr << "input_file: file containing the graph to decompose in communities." << endl;
  cerr << "-w file\tread the graph as a weighted one (weights are set to 1 otherwise)." << endl;
  cerr << "-p file\tstart the computation with a given partition instead of the trivial partition." << endl;
  cerr << "\tfile must contain lines \"node community\"." << endl;
  cerr << "-q eps\ta given pass stops when the modularity is increased by less than epsilon." << endl;
  cerr << "-l k\tdisplays the graph of level k rather than the hierachical structure." << endl;
  cerr << "\tif k=-1 then displays the hierarchical structure rather than the graph at a given level." << endl;
  cerr << "-v\tverbose mode: gives computation time, information about the hierarchy and modularity." << endl;
  cerr << "-s\tsequence directory: indicates the directory where the sequence files are placed." << endl;
  cerr << "-h\tshow this usage message." << endl;
  exit(0);
}

void
parse_args(int argc, char **argv) {
  if (argc<2)
    usage(argv[0], "Bad arguments number\n");

  for (int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
      case 'w':
    type = WEIGHTED;
        filename_w = argv[i+1];
    i++;
    break;
      case 'p':
        filename_part = argv[i+1];
    i++;
    break;
      case 'q':
    precision = atof(argv[i+1]);
    i++;
    break;
      case 'l':
    display_level = atoi(argv[i+1]);
    i++;
    break;
      case 'k':
    k1 = atoi(argv[i+1]);
    i++;
    break;
    case 's':
    directory = argv[i+1];
    i++;
    break;
    case 'v':
        verbose=true;
    break;
      default:
    usage(argv[0], "Unknown option\n");
      }
    } else {
      if (filename==NULL)
        filename = argv[i];
      else
        usage(argv[0], "More than one filename\n");
    }
  }
}

void
display_time(const char *str) {
    clock_t rawtime;
    rawtime = clock();
    cerr << str << ": " << ((float)rawtime)/CLOCKS_PER_SEC*1000 << " milliseconds" << std::endl;
}

struct Edge
{
    Edge (unsigned int src, unsigned int dest, double weight)
        : src (src), dest (dest), weight (weight)
    {
    }

    Edge (unsigned int src, unsigned int dest)
        : src (src), dest (dest), weight (0.)
    {
    }

    const bool operator< (const Edge& rhs) const
    {
        if (this->src == rhs.src)
        {
            return this->dest < rhs.dest;
        }
        else
        {
            return this->src < rhs.src;
        }
    }

    unsigned int src;
    unsigned int dest;
    double weight;
};

std::set<Edge> read_from_file (const char *filename, bool& file_exists, int type = UNWEIGHTED)
{
    std::set<Edge> edges;

    ifstream finput;
    finput.open(filename, fstream::in);

    if (finput.good())
    {
        file_exists = true;

        while (!finput.eof()) {
            unsigned int src, dest;
            double weight=1.;

            if (type==WEIGHTED) {
                finput >> src >> dest >> weight;
                edges.insert(Edge(src, dest, weight));
            } else {
                finput >> src >> dest;
                edges.insert(Edge(src, dest));
            }
        }

        finput.close();
    }
    else
    {
        file_exists = false;
    }
    
    return edges;
}

bool check_g (Community* c)
{
#ifdef _DEBUG
	/*if (c -> n2c.size() != c -> c2n.size())
	{
		return false;
	}

	int sum_degress = 0;
	for (unordered_map<unsigned int, unsigned long>::iterator it = c -> g -> degrees.begin();
		it != c -> g -> degrees.end();
		it++)
	{
		sum_degress += it->second;
	}

	if (c -> g -> links.size() != sum_degress)
	{
		return false;
	}*/
#endif

	return true;
}


int
main(int argc, char **argv) {
# if defined(WIN32)
    srand(time(NULL)+_getpid());
# else
    srand(time(NULL)+getpid());
# endif

  parse_args(argc, argv);
  time_t time_begin, time_end;
  time(&time_begin);
  if (verbose)
    display_time("Begin");

  std::vector <Community*> iterations;

  Community* cp1 = new Community(filename, filename_w, type, -1, precision);
  Community* cp2;
  Community* cp0 = NULL;
  if (filename_part!=NULL)
      cp1 -> init_partition(filename_part);
  Graph* g;
  bool improvement=true;
  double mod=cp1 -> modularity(), new_mod;
  int level=0;
  
  bool has_add_remove_files = true;
  int index = 0;
  while (has_add_remove_files)
  {
      if (index > 0)
      {
          std::cout << "begin: " << std::setw(10) << std::setfill('0') << index - 1 << std::endl;
      }      


      /*do*/ {
          if (verbose) {
              cerr << "level " << level << ":\n";
              display_time("  start computation");
              cerr << "  network size: " 
                  << cp1 -> g ->nb_nodes << " nodes, " 
                  << cp1 -> g ->nb_links << " links, "
                  << cp1 -> g ->total_weight << " weight." << endl;
          }

          std::unordered_map<unsigned int, int> new_n2c;
          cp1 -> print_data("start one_level");
          improvement = cp1 -> one_level(new_n2c);
		  check_g (cp1); //mmc
          cp1 -> print_data("end one_level");

          //std::set<int> new_n2c_communities;
          //for (std::unordered_map<unsigned int, int>::iterator new_n2c_it = new_n2c.begin(); new_n2c_it != new_n2c.end(); new_n2c_it++)
          //{
          //    new_n2c_communities.insert(new_n2c_it->second);
          //    if (new_n2c_communities.size() > 1)
          //    {
          //        break;
          //    }
          //}
          
          if (cp0 != NULL /*&& new_n2c_communities.size() > 1*/)
          {
              //update cp0 final communities
              cp0 -> print_data("start update_community");
              cp0->update_community(cp1, new_n2c);
			  check_g (cp0); //mmc
              cp0 -> print_data("end update_community");

          }

          new_n2c.clear();

          new_mod = cp1 -> modularity();
          if (++level==display_level)
              g ->display();
          if (display_level==-1) {
			  if (cp0 != NULL) {
				cp0 -> display_partition();
			  } 
			  else {
				cp1 -> display_partition();
			  }

		  }

          cp1 -> print_data("start partition2graph_binary");
          g = cp1 -> partition2graph_binary();
          cp1 -> print_data("end partition2graph_binary");

          cp2 = new Community(g, -1, precision);
		  check_g (cp2); //mmc

          if (verbose)
              cerr << "  modularity increased from " << mod << " to " << new_mod << endl;

          mod=new_mod;
          if (verbose)
              display_time("  end computation");

          if (filename_part!=NULL && level==1) // do at least one more computation if partition is provided
              improvement=true;
      } //while(improvement);

    bool has_add_file = false;

    std::stringstream add;
    add << directory.c_str() << "/a" << std::setw(10) << std::setfill('0') << index << ".txt";
    //add << "a" << std::setw(10) << std::setfill('0') << index << ".txt";

    std::set<Edge> add_edges = read_from_file (add.str().c_str(), has_add_file, type);

    std::set<unsigned int> final_affected_nodes;

    int added_nodes = 0;

    for (std::set<Edge>::iterator it = add_edges.begin(); it != add_edges.end(); it++)
    {
        added_nodes++;
        std::set<unsigned int> affected_nodes;
        if (type == WEIGHTED)
        {
            cp1 -> add(it->src, it->dest, it->weight, affected_nodes);
        }
        else
        {
            if (cp0 == NULL)
            {
                std::stringstream added_nodes_s;
                added_nodes_s << added_nodes;

                //display_time(added_nodes_s.str().c_str());

                //display_time("  cp1 -> add () begin");
                std::stringstream mm;
                mm << "start add(" << it->src << "<->" << it->dest << ")" << endl;
                cp1 -> print_data(mm.str().c_str());
                cp1 -> add(it->src, it->dest, 0, affected_nodes);
                cp1 ->print_data("end add");
                //display_time("  cp1 -> add () end");

                //display_time("  cp2 -> update () begin");
                cp2 -> print_data("start update");
                cp2 -> update_add(cp1, affected_nodes);
                cp2 -> print_data("end update");
                //display_time("  cp2 -> update () end");
            }
            else
            {
                //display_time("  cp0 -> update () begin");
                std::stringstream mm;
                mm << "start add(" << it->src << "<->" << it->dest << ")" << endl;
                cp0 -> print_data(mm.str().c_str());
                cp0 -> add(it->src, it->dest, 0, affected_nodes);
                cp0 -> print_data("end add");
                //display_time("  cp0 -> update () end");

                //display_time("  cp2 -> update () begin");
                cp2 -> print_data("start update");
                cp2 -> update_add(cp0, affected_nodes);
                cp2 -> print_data("end update");
                //display_time("  cp2 -> update () end");
            }
           
        }

        final_affected_nodes.insert(affected_nodes.begin(), affected_nodes.end());
    }

    bool has_rem_file = false;

    std::stringstream rem;
    rem << directory.c_str() << "/r" << std::setw(10) << std::setfill('0') << index << ".txt";
    //rem << "r" << std::setw(10) << std::setfill('0') << index << ".txt";

    std::set<Edge> rem_edges = read_from_file (rem.str().c_str(), has_rem_file, type);
    
    int removed_nodes = 0;

    for (std::set<Edge>::iterator it = rem_edges.begin(); it != rem_edges.end(); it++)
    {
        removed_nodes++;
        std::set<unsigned int> affected_nodes;
        if (type == WEIGHTED)
        {
            cp1 -> rem(it->src, it->dest, it->weight, affected_nodes);
        }
        else
        {
            if (cp0 == NULL)
            {
                std::stringstream removed_nodes_s;
                removed_nodes_s << removed_nodes;

                //display_time(removed_nodes_s.str().c_str());

                //display_time("  cp1 -> add () begin");
                cp1 -> rem(it->src, it->dest, 0, affected_nodes);
				check_g (cp1); //mmc
                //display_time("  cp1 -> add () end");

                //display_time("  cp2 -> update () begin");
                cp2 -> update_rem(cp1, affected_nodes);
				check_g (cp2); //mmc
                //display_time("  cp2 -> update () end");
            }
            else
            {
                //display_time("  cp0 -> update () begin");
                cp0 -> rem(it->src, it->dest, 0, affected_nodes);
				check_g (cp0); //mmc
                //display_time("  cp0 -> update () end");

                //display_time("  cp2 -> update () begin");
                cp2 -> update_rem(cp0, affected_nodes);
				check_g (cp2); //mmc
                //display_time("  cp2 -> update () end");
            }

        }

        final_affected_nodes.insert(affected_nodes.begin(), affected_nodes.end());
    }

    //iterations[0] = cp1;

    index++; 

    has_add_remove_files = has_add_file || has_rem_file;
    
    if (cp0 == NULL)
    {
        cp0 = new Community(new Graph(cp1->g->nb_nodes, cp1->g->nb_links, cp1->g->total_weight, cp1->g->degrees, cp1->g->links, cp1->g->weights), *cp1, cp1->nb_pass, cp1->min_modularity);
		check_g (cp0); //mmc
    }
    else
    {
        /*for (std::set<unsigned int>::iterator it = final_affected_nodes.begin(); it != final_affected_nodes.end(); it++)
        {
            cp0 -> n2c [*it] = cp1 -> n2c [*it];
            cp0 -> in [*it] = cp1 -> in [*it];
            cp0 -> tot [*it] = cp1 -> tot [*it];
        }*/

        delete cp1;
    }
    
    cp1 = cp2;


  }
    
  //do {
  //    if (verbose) {
  //        cerr << "level " << level << ":\n";
  //        display_time("  start computation");
  //        cerr << "  network size: " 
  //            << cp1 -> g ->nb_nodes << " nodes, " 
  //            << cp1 -> g ->nb_links << " links, "
  //            << cp1 -> g ->total_weight << " weight." << endl;
  //    }

  //    std::unordered_map<unsigned int, int> new_n2c;
  //    improvement = cp1 -> one_level(new_n2c);

  //    std::set<int> new_n2c_communities;
  //    for (std::unordered_map<unsigned int, int>::iterator new_n2c_it = new_n2c.begin(); new_n2c_it != new_n2c.end(); new_n2c_it++)
  //    {
  //        new_n2c_communities.insert(new_n2c_it->second);
  //        if (new_n2c_communities.size() > 1)
  //        {
  //            break;
  //        }
  //    }

  //    if (cp0 != NULL && new_n2c_communities.size() > 1)
  //    {
  //        //update cp0 final communities
  //        cp0->update_community(cp1, new_n2c);

  //    }

  //    new_mod = cp1 -> modularity();
  //    if (++level==display_level)
  //        g ->display();
  //    if (display_level==-1)
  //        cp1 -> display_partition();

  //    g = cp1 -> partition2graph_binary();

  //    cp2 = new Community(g, -1, precision);

  //    if (verbose)
  //        cerr << "  modularity increased from " << mod << " to " << new_mod << endl;

  //    mod=new_mod;
  //    if (verbose)
  //        display_time("  end computation");

  //    if (filename_part!=NULL && level==1) // do at least one more computation if partition is provided
  //        improvement=true;
  //} while(improvement);

  time(&time_end);
  if (verbose) {
    display_time("End");
    cerr << "Total duration: " << (time_end-time_begin) << " sec." << endl;
  }
  cerr << new_mod << endl;
}

