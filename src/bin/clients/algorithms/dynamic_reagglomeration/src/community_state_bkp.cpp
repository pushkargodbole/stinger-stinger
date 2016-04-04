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
//#include "utils.h"

#define ASSERT(left,operator,right) { if(!((left) operator (right))){ std::cerr << "ASSERT FAILED: " << #left << #operator << #right << " @ " << __FILE__ << " (" << __LINE__ << "). " << #left << "=" << (left) << "; " << #right << "=" << (right) << std::endl; throw;} }

using namespace std;

class community_state
{
    public:
        vector<vector<pair<int64_t, int64_t> > > parents; // Chronological list of parental history for each node with timestamps (parent, timestamp)
        set<int64_t> active; // Set of presently active communities
        vector<int64_t> size; // Present number of nodes inside each community
        //vector<int64_t> inedges; // Present number of edges inside each community
        //vector<int64_t> exedges; // Present number of edges incident on each community
        vector<vector<int64_t> > children; // Chronological list of children merged into each community (recursive)
        //map<pair<int64_t, int64_t>, int64_t> adjacency; // Upper triangular sparse matrix of the present number of edges between each community pair (including self edges)
        vector<map<int64_t, int64_t> > nbrs; // Present set of neighbors for each community, based on the adjacency
        int64_t nv; // Number of vertices
        double ne; // Number of edges
        double modularity;
        int64_t time;
        
        int64_t get_adjacency(int64_t u, int64_t v)
        {
            if(nbrs[u].find(v) == nbrs[u].end()) return 0;
            else return nbrs[u][v];
        }
        
        void add_adjacency(int64_t u, int64_t v, int64_t w, int64_t up=-1, int64_t vp=-1)
        {
            //exedges[u] += w;
            //exedges[v] += w;
            if(nbrs[u].find(v) == nbrs[u].end())
            {
                if(w > 0)
                {
                    nbrs[u][v] = w;
                    if(u != v) nbrs[v][u] = w;
                }
                else
                {
                    cout << "WARNING: Attempt to delete non-existent edge (" << u << ", " << v << ", " << up << ", " << vp << ")" << endl;
                    throw;
                }
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
            parents[v].push_back(make_pair(p, time));
            queue<int64_t> child_queue;
            if(children[v].size() > 0) child_queue.push(v);
            while(child_queue.size() > 0)
            {
                int64_t child = child_queue.front();
                child_queue.pop();
                for(int64_t i=0;i<children[child].size();i++)
                {
                    int64_t newchild = children[child][i];
                    parents[newchild].push_back(make_pair(p, time));
                    if(children[newchild].size() > 0) child_queue.push(newchild);
                }
            }
        }
        
        void revert_parent(int64_t v) // Recursively revert parent of vertex v and all its children to previous
        {
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
            int64_t parent_u = find(u).first;
            if(u!=parent_u)
            {
                cout << u << endl;
                for(int64_t i=0;i<parents[u].size();i++) cout << parents[u][i].first << " ";
                cout << endl;
                if(active.find(u)!=active.end()) cout << "ACTIVE" << endl;
                throw;
            }
            
            //ASSERT(u, ==, parent_u);
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
        
        void edit_edge(int64_t u, int64_t v, int64_t w) // Recursively add an edge of weight w, between vertices u and v from the community graph (w can be negative for deletions)
        {
            //ASSERT(find(u).first, !=, find(v).first);
            int64_t i = parents[u].size()-1;
            int64_t j = parents[v].size()-1;
            while(i >= 0 && j >= 0)
            {
                int64_t pu = parents[u][i].first;
                int64_t pv = parents[v][j].first;
                int64_t tu = parents[u][i].second;
                int64_t tv = parents[v][j].second;
                //cout << "E: " << pu << " " << pv << endl;
                add_adjacency(pu, pv, w);
                if(tu >= tv) i--;
                else j--;
                
            }
        }
        
    //public:
        void init(const struct stinger * S) // Initializer
        {
            nv = stinger_max_active_vertex(S)+1; // +1, for indices begin at 0 and end at (nv-1)
            time = 0;
            ne = 0;
            modularity = 0;
            parents.clear();
            size.clear();
            active.clear();
            children.clear();
            nbrs.clear();
            vector<pair<int64_t, int64_t> > empty_pair_vector;
            vector<int64_t> empty_vector;
            map<int64_t, int64_t> empty_map;
            //cout << "woohoo " << nv << endl;
            for(int64_t v=0;v<nv;++v)
            {
                //cout << v << " ";
                parents.push_back(empty_pair_vector);
                children.push_back(empty_vector);
                nbrs.push_back(empty_map);
                //inedges.push_back(0);
                //exedges.push_back(0);
            }
            //cout << endl << "-----------------------------" << endl;
            for(int64_t v=0;v<nv;++v)
            {
                parents[v].push_back(make_pair(v, time));
                active.insert(v);
                size.push_back(1);
                //children[v].insert(-1); // -1 represents no children (empty set)
                //cout << "woohoo1 " << v << endl;
                //cout << v << " : ";
                STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, v)
                {
                    int64_t dest = STINGER_EDGE_DEST;
                    //cout << dest << " ";
                    int64_t w = STINGER_EDGE_WEIGHT;
                    if(nbrs[v].find(dest) == nbrs[v].end())
                    {
                        add_adjacency(v, dest, 1);
                        ne++;
                    }
                    //cout << "woohoo5" << endl;
                } STINGER_FORALL_EDGES_OF_VTX_END();
                //cout << endl;
                //cout << "woohoo2" << endl;
            }
        }
        
        pair<int64_t, int64_t> find(int64_t v) // Find parent of vertex v
        {
            return parents[v].back();
        }
        
        void merge(int64_t u, int64_t v) // Merge communities u and v
        {
            int64_t parent_u = find(u).first;
            int64_t parent_v = find(v).first;
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
                //inedges[parent] += inedges[child] + nbrs[parent][child];
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
            int64_t parent_u = find(u).first;
            ASSERT(u, ==, parent_u);
            if(children[u].size()>0)
            {
                //cout << "woohoo1" << endl;
                int64_t child = children[u].back();
                children[u].pop_back();
                revert_parent(child);
                size[u] -= size[child];
                //inedges[u] -= inedges[child] + nbrs[u][child];
                //cout << "woohoo2" << endl;
                for(map<int64_t, int64_t>::iterator nbr_it=nbrs[child].begin(); nbr_it!=nbrs[child].end(); ++nbr_it)
                {
                    int64_t nbr = nbr_it->first;
                    int64_t wt = nbr_it->second;
                    if(active.find(nbr)!=active.end()) add_adjacency(u, nbr, -wt);
                }
                add_adjacency(u, u, -get_adjacency(child, child)); // Remove community internal edges
                active.insert(child);
                //cout << "woohoo3" << endl;
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
            //do
            //{
                delta_mod = 0;
                set<int64_t>::iterator it=active.begin();
                set<int64_t>::iterator current;
                while(it!=active.end())
                {
                    cout << "woohoo" << endl;
                    while(active.find(*it)==active.end() && it!=active.end()) it++;
                    if(it==active.end()) break;
                    current = it++;
                    cout << "woohoo1" << endl;
                    int n = *current;
                    cout << "woohoo2" << endl;
                    for(map<int64_t, int64_t>::iterator nbr_it=nbrs[n].begin(); nbr_it!=nbrs[n].end(); ++nbr_it)
                    {
                        int64_t nbr = nbr_it->first;
                        if(active.find(nbr)!=active.end() && n!=nbr)
                        {
                            //cout << "s" << endl;
                            if(active.find(n)!=active.end()) cout << "ACTIVE n " << n << endl;
                            else cout << "INACTIVE n " << n << endl;
                            if(active.find(nbr)!=active.end()) cout << "ACTIVE nbr " << nbr << endl;
                            double d_mod = dmod(n, nbr);
                            //cout << "e" << endl;
                            if(d_mod > 0)
                            {
                                cout << n << " " << nbr << " " << d_mod << endl;
                                delta_mod += d_mod;
                                merge(n, nbr);
                                cout << "woogoo" << endl;
                            }
                            if(active.find(n)==active.end())
                            {
                                cout << "break" << endl;
                                break;
                            }
                        }
                    }
                    /*
                    do
                    {
                        it++;
                        cout << *it << endl;
                    } while(active.find(*it)==active.end() && it!=active.end());
                    //it++;
                    */
                }
            //} while(delta_mod > 0);
        }
        
        void add_batch(const stinger_registered_alg * alg)
        {
            //const struct stinger * S = alg->stinger;
            const int64_t nins = alg->num_insertions;
            const stinger_edge_update * ins = alg->insertions;
            const int64_t nrem = alg->num_deletions;
            const stinger_edge_update * rem = alg->deletions;
            
            for(int64_t i = 0; i < nrem; ++i)
            {
                //cout << "werdf" << endl;
                const int64_t u = rem[i].source;
                const int64_t v = rem[i].destination;
                const int64_t w = rem[i].weight;
                //cout << u << " " << v << " " << w << " " << endl;
                if(nbrs[u].find(v) != nbrs[u].end())
                {
                    //cout << i << endl;
                    ne--;
                    if(find(u).first==find(v).first)
                    {
                        int64_t parent = find(u).first;
                        while(find(u).first==find(v).first)
                        {
                            split(parent);
                            parent = find(u).first;
                        }
                    }
                    edit_edge(u, v, -1);
                }
            }
            
            for(int64_t i = 0; i < nins; ++i)
            {
                const int64_t u = ins[i].source;
                const int64_t v = ins[i].destination;
                const int64_t w = ins[i].weight;
                if(nbrs[u].find(v) == nbrs[u].end())
                {
                    ne++;
                    edit_edge(u, v, 1);
                }
            }
        }
};
