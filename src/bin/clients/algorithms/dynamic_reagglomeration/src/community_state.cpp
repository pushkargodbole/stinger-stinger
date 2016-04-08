#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <utility>      // std::pair, std:: make_pair
#include <queue>
# include <set>
#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>
#include "stinger_core/stinger.h"
#include "stinger_core/stinger_atomics.h"
#include "stinger_core/xmalloc.h"
#include "stinger_core/stinger_error.h"
#include "stinger_utils/timer.h"
#include "stinger_net/stinger_alg.h"

//#include "utils.h"

#define ASSERT(left,operator,right) { if(!((left) operator (right))){ std::cerr << "ASSERT FAILED: " << #left << #operator << #right << " @ " << __FILE__ << " (" << __LINE__ << "). " << #left << "=" << (left) << "; " << #right << "=" << (right) << std::endl; throw;} }

using namespace std;

inline bool edge_compare(const pair<double, pair<int64_t, int64_t> >& a, const pair<double, pair<int64_t, int64_t> >& b)
{
    if(a > b) return true;
    else return false;
}

class community_state
{
    public:
        vector<vector<pair<int64_t, int64_t> > > parents; // Chronological list of parental history for each node with timestamps (parent, timestamp)
        vector<int64_t> LastParents; // Parenthood before the last batch addition
        vector<int64_t> active; // Presently active communities marked with 1 and inactive marked with 0
        vector<int64_t> size; // Present number of nodes inside each community
        //vector<int64_t> inedges; // Present number of edges inside each community
        //vector<int64_t> exedges; // Present number of edges incident on each community
        vector<vector<int64_t> > children; // Chronological list of children merged into each community (recursive)
        //map<pair<int64_t, int64_t>, int64_t> adjacency; // Upper triangular sparse matrix of the present number of edges between each community pair (including self edges)
        vector<unordered_map<int64_t, int64_t> > nbrs; // Present set of neighbors for each community, based on the adjacency (neighborid, weight)
        int64_t nv; // Number of vertices
        double ne; // Number of edges
        double modularity;
        int64_t time;
        string inter_ins_bt = "ee"; // Inter community edge insertion backtracking policy [Default: ee (Edit Edge)]
        string inter_rem_bt = "ee"; // Inter community edge removal backtracking policy [Default: ee (Edit Edge)]
        string intra_ins_bt = "ee"; // Intra community edge insertion backtracking policy [Default: ee (Edit Edge)]
        string intra_rem_bt = "sep"; // Intra community edge removal backtracking policy [Default: sep (Split to Separation)]
        
        string merge_policy = "match"; //Merge policy for agglomeration [Default: match (Greedy Edge Matching)]
        
        unordered_map<string, double> Diagnostic;
        void ResetDiagnostic(string flag = "init")
        {
            if(Diagnostic.size() == 0)
            {
                Diagnostic["nedit_edge"] = 0; //Number of calls to edit_edge
                Diagnostic["nedit_edgeedits"] = 0; //Number of calls to add_adjacency due to edit_edge
                Diagnostic["nmerge"] = 0; //Number calls to merge
                Diagnostic["nmergeedits"] = 0; //Number of calls to add_adjacency due to merge
                Diagnostic["nsplit"] = 0; //Number of calls to split
                Diagnostic["nsplitedits"] = 0; //Number of calls to add_adjacency due to split
                Diagnostic["nedges"] = 0; //Total number of edges in dendogram
                Diagnostic["split_time"] = 0; //Total time required for split
                Diagnostic["merge_time"] = 0; //Total time required for merge
                Diagnostic["editedge_time"] = 0; //Total time required for edit_edge
                Diagnostic["agglomerate_time"] = 0; //Total time required for agglomerate
                Diagnostic["match_time"] = 0; //Total time required for matching in agglomerate_match
                Diagnostic["addbatch_time"] = 0; //Total time required for add_batch
                Diagnostic["init_time"] = 0; //Total time required for init
            }
            else
            {
                for (unordered_map<string, double>::iterator it=Diagnostic.begin(); it!=Diagnostic.end(); ++it)
                {
                    if(flag == "init" || (flag == "add_batch" && it->first != "nedges")) it->second = 0;
                }
            }
        }
        
        int64_t get_adjacency(int64_t u, int64_t v)
        {
            if(nbrs[u].find(v) == nbrs[u].end()) return 0;
            else return nbrs[u][v];
        }
        
        void add_adjacency(int64_t u, int64_t v, int64_t w, string flag="base")
        {
            if(nbrs[u].find(v) == nbrs[u].end())
            {
                if(w > 0)
                {
                    nbrs[u][v] = w;
                    if(u != v) nbrs[v][u] = w;
                    Diagnostic["nedges"]++;
                }
                else if(w < 0)
                {
                    cout << "WARNING: Attempt to delete non-existent edge (" << u << ", " << v << ", " << w << ", " << flag << ")" << endl;
                    for(int i=0; i<children[u].size(); i++) cout << children[u][i] << " ";
                    cout << endl;
                    
                    throw;
                }
            }
            else
            {
                if(nbrs[u][v] + w == 0)
                {
                    nbrs[u].erase(v);
                    if(u != v) nbrs[v].erase(u);
                    Diagnostic["nedges"]--;
                }
                else if(nbrs[u][v] + w < 0)
                {
                    cout << "WARNING: Negative edge weight addition (" << nbrs[u][v] << ", " << u << ", " << v << ", " << w << ", " << flag << ")" << endl;
                    throw;
                }
                else
                {
                    nbrs[u][v] += w;
                    if(u != v) nbrs[v][u] += w;
                }
            }
            
            if(w != 0)
            {
                if(flag == "edit_edge") Diagnostic["nedit_edgeedits"]++;
                if(flag == "split") Diagnostic["nsplitedits"]++;
                if(flag == "merge") Diagnostic["nmergeedits"]++;
            }
        }
        
        void change_parent(int64_t v, int64_t p) // Recursively change parent of vertex v and all its children to p
        {
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
                if(active[u]==1) cout << "ACTIVE" << endl;
                throw;
            }
            int64_t vol = 0;
            // External edges
            for(unordered_map<int64_t, int64_t>::iterator it=nbrs[u].begin(); it!=nbrs[u].end(); ++it)
            {
                if(it->first != u && active[it->first]==1) vol += it->second;
            }
            // Internal edges
            vol += 2*get_adjacency(u, u);
            return vol;
        }
        
        void edit_edge(int64_t u, int64_t v, int64_t w) // Recursively add an edge of weight w, between vertices u and v from the community graph (w can be negative for deletions)
        {
            double tstart = timer();
            int64_t upos = parents[u].size()-1;
            int64_t vpos = parents[v].size()-1;
            while(upos >= 0 && vpos >= 0)
            {
                int64_t pu = parents[u][upos].first;
                int64_t pv = parents[v][vpos].first;
                int64_t tu = parents[u][upos].second;
                int64_t tv = parents[v][vpos].second;
                if(tu == tv)
                {
                    if(upos != 0 && vpos != 0 && pu != pv)
                    {
                        cout << "Simultaneous merge with different parents detected!" << endl;
                        throw;
                    }
                    add_adjacency(pu, pv, w, "edit_edge");
                    upos--;
                    vpos--;
                }
                if(tu > tv)
                {
                    int64_t i = vpos+1;
                    do
                    {
                        i--;
                        pv = parents[v][i].first;
                        add_adjacency(pu, pv, w, "edit_edge");
                    }while(i>0 && pv!=pu);
                    upos--;
                }
                if(tu < tv)
                {
                    int64_t i = upos+1;
                    do
                    {
                        i--;
                        pu = parents[u][i].first;
                        add_adjacency(pu, pv, w, "edit_edge");
                    }while(i>0 && pv!=pu);
                    vpos--;
                }
            }
            Diagnostic["nedit_edge"]++;
            Diagnostic["editedge_time"] += timer() - tstart;
        }
        
        void init(const struct stinger * S) // Initializer
        {
            double tstart = timer();
            nv = stinger_max_active_vertex(S)+1; // +1, because indices begin at 0 and end at (nv-1)
            time = 0;
            ne = 0;
            modularity = 0;
            LastParents.clear();
            ResetDiagnostic("init");
            if(parents.size() > 0)
            {
                for(int64_t v=0; v<nv; v++) LastParents.push_back(find(v).first);
            }
            else
            {
                for(int64_t v=0; v<nv; v++) LastParents.push_back(v);
            }
            parents.clear();
            size.clear();
            active.clear();
            children.clear();
            nbrs.clear();
            vector<pair<int64_t, int64_t> > empty_pair_vector;
            vector<int64_t> empty_vector;
            unordered_map<int64_t, int64_t> empty_map;
            for(int64_t v=0;v<nv;++v)
            {
                parents.push_back(empty_pair_vector);
                children.push_back(empty_vector);
                nbrs.push_back(empty_map);
            }
            for(int64_t v=0;v<nv;++v)
            {
                parents[v].push_back(make_pair(v, time));
                active.push_back(1);
                size.push_back(1);
                STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, v)
                {
                    int64_t dest = STINGER_EDGE_DEST;
                    int64_t w = STINGER_EDGE_WEIGHT;
                    if(dest > v)
                    {
                        if(w<0)
                        {
                            cout << "w = " << w << endl;
                            throw;
                        }
                        add_adjacency(v, dest, w);
                        ne = ne + w;
                    }
                } STINGER_FORALL_EDGES_OF_VTX_END();
            }
            time++;
            Diagnostic["init_time"] += timer() - tstart;
        }
        
        pair<int64_t, int64_t> find(int64_t v) // Find parent of vertex v
        {
            return parents[v].back();
        }
        
        int64_t merge(int64_t u, int64_t v, string flag = "dynamic") // Merge communities u and v. Returns id of parent.
        {
            double tstart = timer();
            int64_t parent_u = find(u).first;
            int64_t parent_v = find(v).first;
            int64_t parent, child;
            ASSERT(u, ==, parent_u);
            ASSERT(v, ==, parent_v);
            if(u != v)
            {
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
                children[parent].push_back(child);
                change_parent(child, parent);
                size[parent] += size[child];
                active[child] = 0;
                for(unordered_map<int64_t, int64_t>::iterator nbr_it=nbrs[child].begin(); nbr_it!=nbrs[child].end(); ++nbr_it)
                {
                    int64_t nbr = nbr_it->first;
                    int64_t wt = nbr_it->second;
                    if(flag == "dynamic" || (flag == "static" && active[nbr] == 1))
                    {
                        if(find(nbr).first!=child && find(nbr).first!=parent)
                        {
                            add_adjacency(parent, nbr, wt, "merge");
                        }
                    }
                }
                add_adjacency(parent, parent, get_adjacency(child, child), "merge"); // Add community internal edges
                add_adjacency(parent, parent, get_adjacency(parent, child), "merge"); // Add joining community edge
            }
            time++;
            Diagnostic["nmerge"]++;
            Diagnostic["merge_time"] += timer() - tstart;
            return parent;
        }
        
        void split(int64_t u)
        {
            double tstart = timer();
            int64_t parent_u = find(u).first;
            ASSERT(u, ==, parent_u);
            if(children[u].size()>0)
            {
                int64_t child = children[u].back();
                children[u].pop_back();
                revert_parent(child);
                size[u] -= size[child];
                for(unordered_map<int64_t, int64_t>::iterator nbr_it=nbrs[child].begin(); nbr_it!=nbrs[child].end(); ++nbr_it)
                {
                    int64_t nbr = nbr_it->first;
                    int64_t wt = nbr_it->second;
                    if(find(nbr).first!=child && find(nbr).first!=u) 
                    {
                        add_adjacency(u, nbr, -wt, "split");
                    }
                }
                add_adjacency(u, u, -get_adjacency(child, child), "split"); // Remove community internal edges
                add_adjacency(u, u, -get_adjacency(u, child), "split"); // Remove joining community edge
                active[child] = 1;
            }
            Diagnostic["nsplit"]++;
            Diagnostic["split_time"] += timer() - tstart;
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
            double inedges_v = (double) get_adjacency(v, v);
            double vol_v = (double) volume(v);
            double mod_v = inedges_v/ne - pow(vol_v, 2)/(4*pow(ne, 2));
            double inedges_uv = inedges_u + inedges_v + get_adjacency(u, v);
            double vol_uv = vol_u + vol_v;
            double mod_uv = inedges_uv/ne - pow(vol_uv, 2)/(4*pow(ne, 2));
            return mod_uv - mod_u - mod_v;
        }
        
        double update_mod() // Update and return modularity of the cstate
        {
            modularity = 0;
            for(int64_t i=0;i<nv;i++)
            {
                if(active[i]==1)
                {
                    modularity += mod(i);
                }
            }
            return modularity;
        }
        
        template<typename T>
        double Mean(vector<T> x)
        {
            double sum = accumulate(x.begin(), x.end(), 0);
            return (double) sum/x.size();
        }
        
        template<typename T>
        double Stddev(vector<T> x)
        {
            double mean = Mean(x);
            vector<double> diff(x.size());
            transform(x.begin(), x.end(), diff.begin(), bind2nd(minus<double>(), mean));
            double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
            double stddev = sqrt(sq_sum / x.size());
            
            return stddev;
        }

        int SizeofChange(const struct stinger * S)
        {
            int64_t soc = 0;
            vector<double> Cj;
            vector<double> Cl;
            for(int64_t v=0; v<nv; v++)
            {
                int64_t Join = 0;
                int64_t Leave = 0;
                int64_t Stay = 0;
                STINGER_FORALL_EDGES_OF_VTX_BEGIN(S, v)
                {
                    int64_t dest = STINGER_EDGE_DEST;
                    if(find(v).first == find(dest).first)
                    {
                        if(LastParents[v] == LastParents[dest]) Stay++;
                        if(LastParents[v] != LastParents[dest]) Join++;
                    }
                    if(find(v).first != find(dest).first && LastParents[v] == LastParents[dest]) Leave++;
                    
                } STINGER_FORALL_EDGES_OF_VTX_END();
                
                double cj, cl;
                if(Join+Stay > 0) cj = (double) Join/(Join+Stay);
                else cj = 0;
                if(Leave+Stay > 0) cl = (double) Leave/(Leave+Stay);
                else cl = 0;
                Cj.push_back(cj);
                Cl.push_back(cl);
            }
            
            double meanCj = Mean(Cj);
            double meanCl = Mean(Cl);
            double stddevCj = Stddev(Cj);
            double stddevCl = Stddev(Cl);
            
            for(int64_t v=0; v<nv; v++)
            {
                if(Cj[v] > meanCj + 2*stddevCj || Cl[v] > meanCl + 2*stddevCl) soc++;
            }
            
            return soc;
        }
        
        int64_t ncoms()
        {
            int64_t n = 0;
            for(int64_t v=0; v<nv; v++)
            {
                if(active[v] == 1) n++;
            }
            return n;
        }
        
        int64_t max_comm_size()
        {
            vector<int64_t> active_sizes;
            for(int64_t v=0; v<nv; v++)
            {
                if(active[v] == 1) active_sizes.push_back(size[v]);
            }
            int64_t max_csize = *max_element(active_sizes.begin(), active_sizes.end());
            return max_csize;
        }
        
        int64_t min_comm_size()
        {
            vector<int64_t> active_sizes;
            for(int64_t v=0; v<nv; v++)
            {
                if(active[v] == 1) active_sizes.push_back(size[v]);
            }
            int64_t min_csize = *min_element(active_sizes.begin(), active_sizes.end());
            return min_csize;
        }
        
        double mean_comm_size()
        {
            vector<int64_t> active_sizes;
            for(int64_t v=0; v<nv; v++)
            {
                if(active[v] == 1) active_sizes.push_back(size[v]);
            }
            double mean_csize = Mean(active_sizes);
            return mean_csize;
        }
        
        double stddev_comm_size()
        {
            vector<int64_t> active_sizes;
            for(int64_t v=0; v<nv; v++)
            {
                if(active[v] == 1) active_sizes.push_back(size[v]);
            }
            double stddev_csize = Stddev(active_sizes);
            return stddev_csize;
        }
        
        void set_merge_policy(string mergepolicy)
        {
            if(mergepolicy!="nodespan" && mergepolicy!="match" && mergepolicy!="bestfirst")
            {
                cout << "Possible Merge Policies:" << endl
                     << "1. nodespan" << endl
                     << "2. match" << endl
                     << "3. bestfirst" << endl
                     << mergepolicy << " given." << endl;
                throw;
            }
            
            merge_policy = mergepolicy;
            cout << "Merge Policy: " << merge_policy << endl;
        }
        
        void agglomerate_nodespan(string flag = "dynamic")
        {
            double delta_mod = 0;
            for(int64_t i=0;i<nv;i++)
            {
                while(i<nv && active[i]==0) i++;
                if(i==nv) break;
                for(unordered_map<int64_t, int64_t>::iterator nbr_it=nbrs[i].begin(); nbr_it!=nbrs[i].end(); ++nbr_it)
                {
                    int64_t nbr = nbr_it->first;
                    if(active[nbr]==1 && i!=nbr)
                    {
                        double d_mod = dmod(i, nbr);
                        if(d_mod > 0)
                        {
                            delta_mod += d_mod;
                            int64_t parent = merge(i, nbr, flag);
                            if(parent==nbr) break;
                        }
                    }
                }
            }
        }
        
        void agglomerate_match(string flag = "dynamic")
        {
            bool(*edge_compare_ptr)(const pair<double, pair<int64_t, int64_t> >&, const pair<double, pair<int64_t, int64_t> >&) = edge_compare;
            while(1)
            {
                double tstart = timer();
                set<pair<double, pair<int64_t, int64_t> >, bool(*)(const pair<double, pair<int64_t, int64_t> >&, const pair<double, pair<int64_t, int64_t> >&) > edges (edge_compare_ptr);
                for(int64_t v=0; v<nv; v++)
                {
                    if(active[v] == 1)
                    {
                        for(unordered_map<int64_t, int64_t>::iterator nbr_it=nbrs[v].begin(); nbr_it!=nbrs[v].end(); ++nbr_it)
                        {
                            int64_t nbr = nbr_it->first;
                            if(active[nbr]==1 && v<nbr)
                            {
                                double d_mod = dmod(v, nbr);
                                if(d_mod > 0) edges.insert(make_pair(d_mod, make_pair(v, nbr)));
                            }
                        }
                    }
                }
                Diagnostic["match_time"] += timer() - tstart;
                if(edges.size() == 0) return;
                
                set<int64_t> visited_nodes;
                for (set<pair<double, pair<int64_t, int64_t> > >::iterator it=edges.begin(); it!=edges.end(); ++it)
                {
                    if(visited_nodes.find(it->second.first) == visited_nodes.end() &&
                       visited_nodes.find(it->second.second) == visited_nodes.end())
                    {
                        merge(it->second.first, it->second.second, flag);
                        visited_nodes.insert(it->second.first);
                        visited_nodes.insert(it->second.second);
                    }
                }
            }
        }
        
        void agglomerate_bestfirst(string flag = "dynamic")
        {
            while(1)
            {
                double best_dmod = 0;
                pair<int64_t, int64_t> best_edge;
                for(int64_t v=0; v<nv; v++)
                {
                    if(active[v] == 1)
                    {
                        for(unordered_map<int64_t, int64_t>::iterator nbr_it=nbrs[v].begin(); nbr_it!=nbrs[v].end(); ++nbr_it)
                        {
                            int64_t nbr = nbr_it->first;
                            if(active[nbr]==1 && v<nbr)
                            {
                                double d_mod = dmod(v, nbr);
                                if(d_mod > best_dmod)
                                {
                                    best_dmod = d_mod;
                                    best_edge = make_pair(v, nbr);
                                }
                            }
                        }
                    }
                }
                
                if(best_dmod == 0) return;
                
                merge(best_edge.first, best_edge.second, flag);
            }
        }
        
        void agglomerate(string flag = "dynamic")
        {
            double tstart = timer();
            if(merge_policy == "nodespan") agglomerate_nodespan(flag);
            if(merge_policy == "match") agglomerate_match(flag);
            if(merge_policy == "bestfirst") agglomerate_bestfirst(flag);
            Diagnostic["agglomerate_time"] += timer() - tstart;
        }
        
        void set_bt_policy(string interinsbt="ee", string interrembt="ee", string intrainsbt="ee", string intrarembt="sep")
        {
            //Inter community edge insertion backtracking policy: Edit Edge or Split to Singletons
            if(interinsbt !="ee" && interinsbt!="sing")
            {
                cout << "Possible strategies for inter community edge insertions:" << endl
                     << "1. ee:   Edit Edge" << endl
                     << "2. sing: Split to Singletons" << endl
                     << interinsbt << " given." << endl;
                throw;
            }
            
            inter_ins_bt = interinsbt;
            
            //Inter community edge removal backtracking policy: Edit Edge
            if(interrembt !="ee")
            {
                cout << "Possible strategies for inter community edge removals:" << endl
                     << "1. ee:   Edit Edge" << endl
                     << interrembt << " given." << endl;
                throw;
            }
            
            inter_rem_bt = interrembt;
            
            //Intra community edge insertion backtracking policy: Edit Edge or Split to Separation or Split to Singletons
            if(intrainsbt !="ee" && intrainsbt!="sep" && intrainsbt!="sing")
            {
                cout << "Possible strategies for intra community edge insertions:" << endl
                     << "1. ee:   Edit Edge" << endl
                     << "2. sep:  Split to Separation" << endl
                     << "3. sing: Split to Singletons" << endl
                     << intrainsbt << " given." << endl;
                throw;
            }
            
            intra_ins_bt = intrainsbt;
            
            //Intra community edge removal backtracking policy: Edit Edge or Split to Separation or Split to Singletons
            if(intrarembt !="ee" && intrarembt!="sep" && intrarembt!="sing")
            {
                cout << "Possible strategies for intra community edge removals:" << endl
                     << "1. ee:   Edit Edge" << endl
                     << "2. sep:  Split to Separation" << endl
                     << "3. sing: Split to Singletons" << endl
                     << intrarembt << " given." << endl;
                throw;
            }
            
            intra_rem_bt = intrarembt;
            
            cout << "Backtracking Policies:" << endl
                 << "Inter Edge Insertion: " << inter_ins_bt << endl
                 << "Inter Edge Removal: " << inter_rem_bt << endl
                 << "Intra Edge Insertion: " << intra_ins_bt << endl
                 << "Intra Edge Removal: " << intra_rem_bt << endl;
        }
        
        void add_batch(const stinger_registered_alg * alg)
        {
            double tstart = timer();
            ResetDiagnostic("add_batch");
            for(int64_t v=0; v<nv; v++)
            {
                LastParents[v] = find(v).first;
            }
            
            const int64_t nins = alg->num_insertions;
            const stinger_edge_update * ins = alg->insertions;
            const int64_t nrem = alg->num_deletions;
            const stinger_edge_update * rem = alg->deletions;
            
            for(int64_t i = 0; i < nrem; ++i)
            {
                const int64_t u = rem[i].source;
                const int64_t v = rem[i].destination;
                const int64_t w = rem[i].weight;
                if(nbrs[u].find(v) != nbrs[u].end())
                {
                    ne--;
                    if(find(u).first==find(v).first) //Intra removal
                    {
                        if(intra_rem_bt=="sep" || intra_rem_bt=="sing")
                        {
                            int64_t parent = find(u).first;
                            while(find(u).first==find(v).first)
                            {
                                split(parent);
                                parent = find(u).first;
                            }
                        }
                        
                        if(intra_rem_bt=="sing")
                        {
                            int64_t parent = find(u).first;
                            while(parent!=u)
                            {
                                split(parent);
                                parent = find(u).first;
                            }
                            
                            parent = find(v).first;
                            while(parent!=v)
                            {
                                split(parent);
                                parent = find(v).first;
                            }
                        }
                    }
                    edit_edge(u, v, -1);
                }
                else
                {
                    cout << "Remove stream edge not found: " << u << " " << v << endl;
                    throw;
                }
            }
            
            for(int64_t i = 0; i < nins; ++i)
            {
                const int64_t u = ins[i].source;
                const int64_t v = ins[i].destination;
                const int64_t w = ins[i].weight;
                ne++;
                
                bool intra_edge = false;
                if(find(u).first==find(v).first) //Intra insertion
                {
                    intra_edge = true;
                    if(intra_ins_bt=="sep" || intra_ins_bt=="sing")
                    {
                        int64_t parent = find(u).first;
                        while(find(u).first==find(v).first)
                        {
                            split(parent);
                            parent = find(u).first;
                        }
                    }
                }
                
                if(intra_ins_bt=="sing" && intra_edge==true || inter_ins_bt=="sing" && intra_edge==false)
                {
                    int64_t parent = find(u).first;
                    while(parent!=u)
                    {
                        split(parent);
                        parent = find(u).first;
                    }
                    
                    parent = find(v).first;
                    while(parent!=v)
                    {
                        split(parent);
                        parent = find(v).first;
                    }
                }
                
                edit_edge(u, v, 1);
            }
            Diagnostic["addbatch_time"] += timer() - tstart;
        }
};
