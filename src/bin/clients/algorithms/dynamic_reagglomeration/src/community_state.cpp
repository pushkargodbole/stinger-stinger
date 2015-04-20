#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>
#include <utility>      // std::pair, std:: make_pair
#include <queue>
# include <set>
#include <cmath>
#include "stinger_core/stinger.h"
#include "stinger_core/stinger_atomics.h"
#include "stinger_core/xmalloc.h"
#include "stinger_core/stinger_error.h"
#include "stinger_net/stinger_alg.h"

#define ASSERT(left,operator,right) { if(!((left) operator (right))){ std::cerr << "ASSERT FAILED: " << #left << #operator << #right << " @ " << __FILE__ << " (" << __LINE__ << "). " << #left << "=" << (left) << "; " << #right << "=" << (right) << std::endl; } }

using namespace std;

class community_state
{
    public:
        vector<vector<int64_t> > parents; // Chronological list of parental history for each node
        set<int64_t> active; // Set of presently active communities
        vector<int64_t> size; // Present number of nodes inside each community
        vector<int64_t> inedges; // Present number of edges inside each community
        vector<int64_t> exedges; // Present number of edges incident on each community
        vector<vector<int64_t> > children; // Chronological list of children merged into each community (recursive)
        //map<pair<int64_t, int64_t>, int64_t> adjacency; // Upper triangular sparse matrix of the present number of edges between each community pair (including self edges)
        vector<map<int64_t, int64_t> > nbrs; // Present set of neighbors for each community, based on the adjacency
        int64_t nv; // Number of vertices
        double ne; // Number of edges
        double modularity = 0;
        int64_t time = 0;
        
        int64_t get_adjacency(int64_t u, int64_t v)
        {
            if(nbrs[u].find(v) == nbrs[u].end()) return 0;
            else return nbrs[u][v];
        }
        
        void add_adjacency(int64_t u, int64_t v, int64_t w)
        {
            exedges[u] += w;
            exedges[v] += w;
            if(nbrs[u].find(v) == nbrs[u].end())
            {
                nbrs[u][v] = w;
                if(u != v) nbrs[v][u] = w;
            }
            else
            {
                if(nbrs[u][v] + w == 0)
                {
                    nbrs[u].erase(v);
                    if(u != v) nbrs[v].erase(u);
                }
                else
                {
                    nbrs[u][v] += w;
                    if(u != v) nbrs[v][u] += w;
                }
            }
        }
        
        void change_parent(int64_t v, int64_t p) // Recursively change parent of vertex v and all its children to p
        {
            children[p].push_back(v);
            //children[p].erase(-1); // Mark children[p] set as not empty (by removing -1 from its set)
            parents[v].push_back(p);
            queue<int64_t> child_queue;
            if(children[v].size() > 0) child_queue.push(v);
            while(child_queue.size() > 0)
            {
                int64_t child = child_queue.front();
                child_queue.pop();
                for(int64_t i=0;i<children[child].size();i++)
                {
                    int64_t newchild = children[child][i];
                    parents[newchild].push_back(p);
                    if(children[newchild].size() > 0) child_queue.push(newchild);
                }
            }
        }
        
        void revert_parent(int64_t v) // Recursively revert parent of vertex v and all its children to previous
        {
            //children[p].erase(-1); // Mark children[p] set as not empty (by removing -1 from its set)
            parents[v].pop_back();
            queue<int64_t> child_queue;
            if(children[v].size() > 0) child_queue.push(v);
            while(child_queue.size() > 0)
            {
                int64_t child = child_queue.front();
                child_queue.pop();
                for(int64_t i=0;i<children[child].size();i++)
                {
                    int64_t newchild = children[child][i];
                    parents[newchild].pop_back();
                    if(children[newchild].size() > 0) child_queue.push(newchild);
                }
            }
        }
        
        int64_t volume(int64_t u) // Volume: Sumation over degrees of all vertices inside community u
        {
            int64_t parent_u = find(u);
            ASSERT(u, ==, parent_u);
            int64_t vol = 0;
            // External edges
            for(map<int64_t, int64_t>::iterator it=nbrs[u].begin(); it!=nbrs[u].end(); ++it)
            {
                if(active.find(it->first)!=active.end()) vol += it->second;
            }
            // Internal edges
            vol += 2*get_adjacency(u, u);
            return vol;
        }
        
    //public:
        void init(const struct stinger * S) // Initializer
        {
            nv = stinger_max_active_vertex(S)+1; // +1, for indices begin at 0 and end at (nv-1)
            time = 0;
            vector<int64_t> empty_vector;
            map<int64_t, int64_t> empty_map;
            //cout << "woohoo " << nv << endl;
            for(int64_t v=0;v<nv;++v)
            {
                cout << v << " ";
                parents.push_back(empty_vector);
                children.push_back(empty_vector);
                nbrs.push_back(empty_map);
                inedges.push_back(0);
                exedges.push_back(0);
            }
            cout << endl << "-----------------------------" << endl;
            for(int64_t v=0;v<nv;++v)
            {
                parents[v].push_back(v);
                active.insert(v);
                size.push_back(1);
                //children[v].insert(-1); // -1 represents no children (empty set)
                //cout << "woohoo1 " << v << endl;
                cout << v << " : ";
                STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, v)
                {
                    int64_t dest = STINGER_EDGE_DEST;
                    cout << dest << " ";
                    int64_t w = STINGER_EDGE_WEIGHT;
                    if(nbrs[v].find(dest) == nbrs[v].end())
                    {
                        add_adjacency(v, dest, 1);
                        ne++;
                    }
                    //cout << "woohoo5" << endl;
                } STINGER_FORALL_EDGES_OF_VTX_END();
                cout << endl;
                //cout << "woohoo2" << endl;
            }
        }
        
        int64_t find(int64_t v) // Find parent of vertex v
        {
            return parents[v].back();
        }
        
        void merge(int64_t u, int64_t v) // Merge communities u and v
        {
            int64_t parent_u = find(u);
            int64_t parent_v = find(v);
            ASSERT(u, ==, parent_u);
            ASSERT(v, ==, parent_v);
            if(u != v)
            {
                int64_t parent, child;
                if(size[u] >= size[v])
                {
                    parent = u;
                    child = v;
                }
                else
                {
                    parent = v;
                    child = u;
                }
                change_parent(child, parent);
                size[parent] += size[child];
                active.erase(child);
                inedges[parent] += inedges[child] + nbrs[parent][child];
                for(map<int64_t, int64_t>::iterator nbr_it=nbrs[child].begin(); nbr_it!=nbrs[child].end(); ++nbr_it)
                {
                    int64_t nbr = nbr_it->first;
                    int64_t wt = nbr_it->second;
                    if(active.find(nbr)!=active.end()) add_adjacency(parent, nbr, wt);
                }
                add_adjacency(parent, parent, get_adjacency(child, child)); // Add community internal edges
            }
            time++;
        }
        
        void split(int64_t u)
        {
            int64_t parent_u = find(u);
            ASSERT(u, ==, parent_u);
            if(children[u].size()>0)
            {
                cout << "woohoo1" << endl;
                int64_t child = children[u].back();
                children[u].pop_back();
                revert_parent(child);
                size[u] -= size[child];
                inedges[u] -= inedges[child] + nbrs[u][child];
                cout << "woohoo2" << endl;
                for(map<int64_t, int64_t>::iterator nbr_it=nbrs[child].begin(); nbr_it!=nbrs[child].end(); ++nbr_it)
                {
                    int64_t nbr = nbr_it->first;
                    int64_t wt = nbr_it->second;
                    if(active.find(nbr)!=active.end()) add_adjacency(u, nbr, -wt);
                }
                add_adjacency(u, u, -get_adjacency(child, child)); // Remove community internal edges
                active.insert(child);
                cout << "woohoo3" << endl;
            }
            time++;
        }
        
        double mod(int64_t u) // Modularity contribution by community node u
        {
            double inedges = (double) get_adjacency(u, u);
            double vol = (double) volume(u);
            double mod = inedges/ne - pow(vol, 2)/(4*pow(ne, 2));
            return mod;
        }
        
        double dmod(int64_t u, int64_t v) // Differential change in modularity on merging community nodes u and v
        {
            double inedges_u = (double) get_adjacency(u, u);
            double vol_u = (double) volume(u);
            double mod_u = inedges_u/ne - pow(vol_u, 2)/(4*pow(ne, 2));
            //cout << u << " " << mod_u << endl;
            double inedges_v = (double) get_adjacency(v, v);
            double vol_v = (double) volume(v);
            double mod_v = inedges_v/ne - pow(vol_v, 2)/(4*pow(ne, 2));
            //cout << v << " " << mod_v << endl;
            double inedges_uv = inedges_u + inedges_v + get_adjacency(u, v);
            double vol_uv = vol_u + vol_v;
            double mod_uv = inedges_uv/ne - pow(vol_uv, 2)/(4*pow(ne, 2));
            //cout << "mod_uv " << mod_uv << endl;
            return mod_uv - mod_u - mod_v;
        }
        
        double update_mod() // Update and return modularity of the cstate
        {
            modularity = 0;
            for(set<int64_t>::iterator it=active.begin(); it!=active.end(); ++it)
            {
                modularity += mod(*it);
            }
            return modularity;
        }
        
        void agglomerate()
        {
            double delta_mod;
            do
            {
                delta_mod = 0;
                set<int64_t>::iterator it=active.begin();
                while(it!=active.end())
                {
                    int n = *it;
                    for(map<int64_t, int64_t>::iterator nbr_it=nbrs[n].begin(); nbr_it!=nbrs[n].end(); ++nbr_it)
                    {
                        int64_t nbr = nbr_it->first;
                        if(active.find(nbr)!=active.end() && n!=nbr)
                        {
                            double d_mod = dmod(n, nbr);
                            if(d_mod > 0)
                            {
                                cout << n << " " << nbr << " " << d_mod << endl;
                                delta_mod += d_mod;
                                merge(n, nbr);
                            }
                            if(active.find(n)==active.end()) break;
                        }
                    }
                    it++;
                }
            } while(delta_mod > 0);
        }
        
        void add_batch(const stinger_registered_alg * alg)
        {
            const struct stinger * S = alg->stinger;
            const int64_t nins = alg->num_insertions;
            const stinger_edge_update * restrict ins = alg->insertions;
            const int64_t nrem = alg->num_deletions;
            const stinger_edge_update * restrict rem = alg->deletions;
            
            for(int64_t i = 0; i < nrem; ++i)
            {
                const int64_t i = rem[k].source;
                const int64_t j = rem[k].destination;
                const int64_t w = rem[k].weight;
                ne--;
                if(find(i)==find(j))
                {
                    int64_t parent = find(i);
                    while(find(i)==find(j))
                    {
                        split(parent);
                        parent = find(i);
                    }
                }
            }
            
            
};
                    

double
cstate_preproc_alg (struct community_state * restrict cstate,
                   const stinger_registered_alg * alg)
{
  //printf("cstatepreproc1\n");
  const struct stinger * S = alg->stinger;
  const int64_t nincr = alg->num_insertions;
  const stinger_edge_update * restrict incr = alg->insertions;
  const int64_t nrem = alg->num_deletions;
  const stinger_edge_update * restrict rem = alg->deletions;
  const int64_t * restrict cmap = cstate->cmap;
  int64_t nvlist = cstate->nvlist;
  int64_t * restrict vlist = cstate->vlist;
  int64_t * restrict mark = cstate->mark;
  intvtx_t * restrict d = cstate->cg.d;
  int64_t n_new_edges = 0;
  //printf("cstatepreproc2\n");
  tic ();
  OMP("omp parallel") {
    struct insqueue q;
    q.n = 0;

    OMP("omp for")
      for (int64_t k = 0; k < nvlist; ++k) mark[vlist[k]] = -1;
    OMP("omp single") nvlist = 0;
    //printf("cstatepreproc3\n");
    OMP("omp for reduction(+: n_new_edges)")
      for (int64_t k = 0; k < nincr; ++k) {
        const int64_t i = incr[k].source;
        const int64_t j = incr[k].destination;
        const int64_t ci = cmap[i];
        const int64_t cj = cmap[j];
        if (ci != cj)
          ++n_new_edges;
      }

    OMP("omp for reduction(+: n_new_edges)")
      for (int64_t k = 0; k < nrem; ++k) {
        const int64_t i = rem[k].source;
        const int64_t j = rem[k].destination;
        const int64_t ci = cmap[i];
        const int64_t cj = cmap[j];
        //if (ci != cj)
          ++n_new_edges;
      }
    //printf("cstatepreproc4\n");
    OMP("omp single") {
      if (realloc_graph (&cstate->cg, cstate->cg.nv, cstate->cg.ne + n_new_edges))
       abort ();
    }

    OMP("omp for")
      for (int64_t k = 0; k < nincr; ++k) {
        const int64_t i = incr[k].source;
        const int64_t j = incr[k].destination;
        const int64_t w = incr[k].weight;
        const int64_t ci = cmap[i];
        const int64_t cj = cmap[j];
        if (ci != cj) {
          append_to_vlist (&nvlist, vlist, mark, i);
          //append_to_vlist (&nvlist, vlist, mark, j);
          STINGER_FORALL_EDGES_OF_VTX_BEGIN (S, i) {
            int64_t i_j = STINGER_EDGE_DEST;
            append_to_vlist (&nvlist, vlist, mark, i_j);
          } STINGER_FORALL_EDGES_OF_VTX_END ();
          STINGER_FORALL_EDGES_OF_VTX_BEGIN (S, j) {
            int64_t j_j = STINGER_EDGE_DEST;
            if(j_j != i) {
              append_to_vlist (&nvlist, vlist, mark, j_j);
            }
          } STINGER_FORALL_EDGES_OF_VTX_END ();
          enqueue (&q, ci, cj, w, &cstate->cg);
        } else
          OMP("omp atomic") d[ci] += w;
      }

    OMP("omp for")
      for (int64_t k = 0; k < nrem; ++k) {
        const int64_t i = rem[k].source;
        const int64_t j = rem[k].destination;
        const int64_t ci = cmap[i];
        const int64_t cj = cmap[j];
        //if (ci != cj) {
        append_to_vlist (&nvlist, vlist, mark, i);
        //append_to_vlist (&nvlist, vlist, mark, j);
        //}
        STINGER_FORALL_EDGES_OF_VTX_BEGIN (S, i) {
            int64_t i_j = STINGER_EDGE_DEST;
            append_to_vlist (&nvlist, vlist, mark, i_j);
          } STINGER_FORALL_EDGES_OF_VTX_END ();
          STINGER_FORALL_EDGES_OF_VTX_BEGIN (S, j) {
            int64_t j_j = STINGER_EDGE_DEST;
            if(j_j != i) {
              append_to_vlist (&nvlist, vlist, mark, j_j);
            }
          } STINGER_FORALL_EDGES_OF_VTX_END ();
        /* Find the weight to remove.  ugh. */
        STINGER_FORALL_EDGES_OF_VTX_BEGIN (S, i) {
          if (STINGER_EDGE_DEST == j) {
            if (ci != cj)
              enqueue (&q, ci, cj, -STINGER_EDGE_WEIGHT, &cstate->cg);
            else
              OMP("omp atomic") d[ci] -= STINGER_EDGE_WEIGHT;
            break;
            /* XXX: Technically, could have many of different types. */
          }
        } STINGER_FORALL_EDGES_OF_VTX_END ();
      }

    qflush (&q, &cstate->cg);
  }
  cstate->nvlist = nvlist;
  if (n_new_edges) {
    realloc_ws (&cstate->ws, &cstate->wslen, cstate->cg.nv_orig, cstate->cg.ne_orig);
    contract_self (&cstate->cg, cstate->ws);
  }
  return toc();
}
