// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// This is the main() file
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
//#include <utility>      // std::pair, std:: make_pair
#include <queue>
# include <set>
#include "stinger_core/stinger.h"
#include "stinger_core/stinger_atomics.h"
#include "stinger_core/xmalloc.h"
#include "stinger_core/stinger_error.h"
#include "stinger_net/stinger_alg.h"
#include "stinger_utils/timer.h"
#include "community_state.cpp"


bool do_static = false;
bool do_dynamic = false;
int nbatch;

int main(int argc, char *argv[])
{
    community_state cstate_s, cstate_d;
    double modularity_s, modularity_d;
    vector<double> modularities_s;
    vector<double> modularities_d;
    int64_t soc_s, soc_d;
    int64_t ncoms_s, ncoms_d;
    int64_t max_csize_s, max_csize_d;
    int64_t min_csize_s, min_csize_d;
    double mean_csize_s, mean_csize_d;
    double stddev_csize_s, stddev_csize_d;
    double time_s, time_d;
    int64_t nedges_s, nedges_d;
    int64_t nedit_edge_s, nedit_edge_d;
    int64_t nmerge_s, nmerge_d;
    int64_t nsplit_s, nsplit_d;
    int64_t nedit_edgeedits_s, nedit_edgeedits_d;
    int64_t nmergeedits_s, nmergeedits_d;
    int64_t nsplitedits_s, nsplitedits_d;
    double split_time_d;
    double merge_time_s, merge_time_d;
    double editedge_time_d;
    double agglomerate_time_s, agglomerate_time_d;
    double match_time_s, match_time_d;
    double init_time_s;
    double addbatch_time_d;
    vector<int64_t> SOC_s;
    vector<int64_t> SOC_d;
    vector<int64_t> Ncoms_s;
    vector<int64_t> Ncoms_d;
    vector<int64_t> Max_csize_s;
    vector<int64_t> Max_csize_d;
    vector<int64_t> Min_csize_s;
    vector<int64_t> Min_csize_d;
    vector<double> Mean_csize_s;
    vector<double> Mean_csize_d;
    vector<double> Stddev_csize_s;
    vector<double> Stddev_csize_d;
    vector<double> Time_s;
    vector<double> Time_d;
    vector<int64_t> Nedges_s;
    vector<int64_t> Nedges_d;
    vector<int64_t> Nedit_edge_s;
    vector<int64_t> Nedit_edge_d;
    vector<int64_t> Nmerge_s;
    vector<int64_t> Nmerge_d;
    vector<int64_t> Nsplit_s;
    vector<int64_t> Nsplit_d;
    vector<int64_t> Nedit_edgeedits_s;
    vector<int64_t> Nedit_edgeedits_d;
    vector<int64_t> Nmergeedits_s;
    vector<int64_t> Nmergeedits_d;
    vector<int64_t> Nsplitedits_s;
    vector<int64_t> Nsplitedits_d;
    vector<double> Split_time_d;
    vector<double> Merge_time_s;
    vector<double> Merge_time_d;
    vector<double> Editedge_time_d;
    vector<double> Agglomerate_time_s;
    vector<double> Agglomerate_time_d;
    vector<double> Match_time_s;
    vector<double> Match_time_d;
    vector<double> Init_time_s;
    vector<double> Addbatch_time_d;
    const struct stinger * S;
    int port = 10103;
    string static_merge_policy = "match";
    string dynamic_merge_policy = "match";
    string inter_ins_bt = "ee";
    string inter_rem_bt = "ee";
    string intra_ins_bt = "ee";
    string intra_rem_bt = "sep";
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "--static") == 0) do_static = true;
        else if (strcmp(argv[i], "--dynamic") == 0) do_dynamic = true;
        else if (0 == strcmp(argv[i], "-nb")) nbatch = atoi(argv[i+1]);
        else if (0 == strcmp(argv[i], "--port") || 0 == strcmp(argv[i], "-p")) port = atoi(argv[i+1]);
        else if (0 == strcmp(argv[i], "-smp")) static_merge_policy = argv[i+1];
        else if (0 == strcmp(argv[i], "-dmp")) dynamic_merge_policy = argv[i+1];
        else if (0 == strcmp(argv[i], "--interins")) inter_ins_bt = argv[i+1];
        else if (0 == strcmp(argv[i], "--interrem")) inter_rem_bt = argv[i+1];
        else if (0 == strcmp(argv[i], "--intrains")) intra_ins_bt = argv[i+1];
        else if (0 == strcmp(argv[i], "--intrarem")) intra_rem_bt = argv[i+1];
        else if (0 == strcmp(argv[i], "--help") || 0 == strcmp(argv[i], "-h"))
        {
            cout << endl; // Write help!
        }
    }
    
    if(do_static == false && do_dynamic == false) do_dynamic = true;
    cstate_s.set_merge_policy(static_merge_policy);
    cstate_d.set_merge_policy(dynamic_merge_policy);
    cstate_d.set_bt_policy(inter_ins_bt, inter_rem_bt, intra_ins_bt, intra_rem_bt);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Setup and register algorithm with the server
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    string alg_name_str = "reagglomeration";
    char * alg_name = new char[alg_name_str.length() + 1];
    strcpy(alg_name, alg_name_str.c_str());
    
    stinger_registered_alg * alg = stinger_register_alg(
        alg_name, // Algorithm name
        "localhost", // Host
        port, // Port
        0, // is.remote
        0, // map_private
        sizeof(int64_t), // data_per_vertex
        "l community_label", // data_description
        NULL, // dependencies
        0 // num_dependencies
    );

    if(!alg)
    {
        LOG_E("Registering algorithm failed.  Exiting");
        return -1;
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initial static computation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    S = alg->stinger;
    cout << "--------------- Static init ------------------" << endl;
    tic();
    cstate_s.init(S);
    cstate_s.agglomerate("static");
    time_s = toc();
    cout << "--------------- Dynamic init ------------------" << endl;
    tic();
    cstate_d.init(S);
    cstate_d.agglomerate();
    time_d = toc();
    modularity_s = cstate_s.update_mod();
    modularity_d = cstate_d.update_mod();
    soc_s = cstate_s.SizeofChange(S);
    soc_d = cstate_d.SizeofChange(S);
    ncoms_s = cstate_s.ncoms();
    ncoms_d = cstate_d.ncoms();
    max_csize_s = cstate_s.max_comm_size();
    max_csize_d = cstate_d.max_comm_size();
    min_csize_s = cstate_s.min_comm_size();
    min_csize_d = cstate_d.min_comm_size();
    mean_csize_s = cstate_s.mean_comm_size();
    mean_csize_d = cstate_d.mean_comm_size();
    stddev_csize_s = cstate_s.stddev_comm_size();
    stddev_csize_d = cstate_d.stddev_comm_size();
    nedges_s = cstate_s.Diagnostic["nedges"];
    nedges_d = cstate_d.Diagnostic["nedges"];
    nedit_edge_s = cstate_s.Diagnostic["nedit_edge"];
    nedit_edge_d = cstate_d.Diagnostic["nedit_edge"];
    nmerge_s = cstate_s.Diagnostic["nmerge"];
    nmerge_d = cstate_d.Diagnostic["nmerge"];
    nsplit_s = cstate_s.Diagnostic["nsplit"];
    nsplit_d = cstate_d.Diagnostic["nsplit"];
    nedit_edgeedits_s = cstate_s.Diagnostic["nedit_edgeedits"];
    nedit_edgeedits_d = cstate_d.Diagnostic["nedit_edgeedits"];
    nmergeedits_s = cstate_s.Diagnostic["nmergeedits"];
    nmergeedits_d = cstate_d.Diagnostic["nmergeedits"];
    nsplitedits_s = cstate_s.Diagnostic["nsplitedits"];
    nsplitedits_d = cstate_d.Diagnostic["nsplitedits"];
    split_time_d = cstate_d.Diagnostic["split_time"];
    merge_time_s = cstate_s.Diagnostic["merge_time"];
    merge_time_d = cstate_d.Diagnostic["merge_time"];
    editedge_time_d = cstate_d.Diagnostic["editedge_time"];
    agglomerate_time_s = cstate_s.Diagnostic["agglomerate_time"];
    agglomerate_time_d = cstate_d.Diagnostic["agglomerate_time"];
    match_time_s = cstate_s.Diagnostic["match_time"];
    match_time_d = cstate_d.Diagnostic["match_time"];
    addbatch_time_d = cstate_d.Diagnostic["addbatch_time"];
    init_time_s = cstate_s.Diagnostic["init_time"];
    cout << "Iter: 0" << endl;
    cout << "Modularity (Static): " << modularity_s << endl;
    cout << "Modularity (Dynamic): " << modularity_d << endl;
    cout << "SOC (Static): " << soc_s << endl;
    cout << "SOC (Dynamic): " << soc_d << endl;
    cout << "Ncoms (Static): " << ncoms_s << endl;
    cout << "Ncoms (Dynamic): " << ncoms_d << endl;
    cout << "Max csize (Static): " << max_csize_s << endl;
    cout << "Max csize (Dynamic): " << max_csize_d << endl;
    cout << "Min csize (Static): " << min_csize_s << endl;
    cout << "Min csize (Dynamic): " << min_csize_d << endl; 
    cout << "Mean csize (Static): " << mean_csize_s << endl;
    cout << "Mean csize (Dynamic): " << mean_csize_d << endl;
    cout << "Stddev csize (Static): " << stddev_csize_s << endl;
    cout << "Stddev csize (Dynamic): " << stddev_csize_d << endl;
    cout << "Nedges (Static): " << nedges_s << endl;
    cout << "Nedges (Dynamic): " << nedges_d << endl;
    cout << "Nedit_edge (Static): " << nedit_edge_s << endl;
    cout << "Nedit_edge (Dynamic): " << nedit_edge_d << endl;
    cout << "Nmerge (Static): " << nmerge_s << endl;
    cout << "Nmerge (Dynamic): " << nmerge_d << endl;
    cout << "Nsplit (Static): " << nsplit_s << endl;
    cout << "Nsplit (Dynamic): " << nsplit_d << endl;
    cout << "Nedit_edge edits (Static): " << nedit_edgeedits_s << endl;
    cout << "Nedit_edge edits (Dynamic): " << nedit_edgeedits_d << endl;
    cout << "Nmerge edits (Static): " << nmergeedits_s << endl;
    cout << "Nmerge edits (Dynamic): " << nmergeedits_d << endl;
    cout << "Nsplit edits (Static): " << nsplitedits_s << endl;
    cout << "Nsplit edits (Dynamic): " << nsplitedits_d << endl;
    cout << "Time (Static): " << time_s << endl;
    cout << "Time (Dynamic): " << time_d << endl;
    cout << "Init time (Static): " << init_time_s << endl;
    cout << "Add Batch time (Dynamic): " << addbatch_time_d << endl;
    cout << "Split time (Dynamic): " << split_time_d << endl;
    cout << "Edit Edge time (Dynamic): " << editedge_time_d << endl;
    cout << "Agglomerate time (Static): " << agglomerate_time_s << endl;
    cout << "Agglomerate time (Dynamic): " << agglomerate_time_d << endl;
    cout << "Match time (Static): " << match_time_s << endl;
    cout << "Match time (Dynamic): " << match_time_d << endl;
    cout << "Merge time (Static): " << merge_time_s << endl;
    cout << "Merge time (Dynamic): " << merge_time_d << endl;
    modularities_s.push_back(modularity_s);
    modularities_d.push_back(modularity_d);
    SOC_s.push_back(soc_s);
    SOC_d.push_back(soc_d);
    Ncoms_s.push_back(ncoms_s);
    Ncoms_d.push_back(ncoms_d);
    Max_csize_s.push_back(max_csize_s);
    Max_csize_d.push_back(max_csize_d);
    Min_csize_s.push_back(min_csize_s);
    Min_csize_d.push_back(min_csize_d);
    Mean_csize_s.push_back(mean_csize_s);
    Mean_csize_d.push_back(mean_csize_d);
    Stddev_csize_s.push_back(stddev_csize_s);
    Stddev_csize_d.push_back(stddev_csize_d);
    Nedges_s.push_back(nedges_s);
    Nedges_d.push_back(nedges_d);
    Nedit_edge_s.push_back(nedit_edge_s);
    Nedit_edge_d.push_back(nedit_edge_d);
    Nmerge_s.push_back(nmerge_s);
    Nmerge_d.push_back(nmerge_d);
    Nsplit_s.push_back(nsplit_s);
    Nsplit_d.push_back(nsplit_d);
    Nedit_edgeedits_s.push_back(nedit_edgeedits_s);
    Nedit_edgeedits_d.push_back(nedit_edgeedits_d);
    Nmergeedits_s.push_back(nmergeedits_s);
    Nmergeedits_d.push_back(nmergeedits_d);
    Nsplitedits_s.push_back(nsplitedits_s);
    Nsplitedits_d.push_back(nsplitedits_d);
    Time_s.push_back(time_s);
    Time_d.push_back(time_d);
    Init_time_s.push_back(init_time_s);
    Addbatch_time_d.push_back(addbatch_time_d);
    Split_time_d.push_back(split_time_d);
    Editedge_time_d.push_back(editedge_time_d);
    Agglomerate_time_s.push_back(agglomerate_time_s);
    Agglomerate_time_d.push_back(agglomerate_time_d);
    Match_time_s.push_back(match_time_s);
    Match_time_d.push_back(match_time_d);
    Merge_time_s.push_back(merge_time_s);
    Merge_time_d.push_back(merge_time_d);
    stinger_alg_begin_init(alg);
    {
    } stinger_alg_end_init(alg);

    //LOG_V_A("simple_communities: init comm time = %g\n", comm_time);

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
     * Streaming Phase
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    int64_t n = 0;
    while(n<nbatch)
    {
        if(alg->enabled)
        {
            // Pre processing
            if(stinger_alg_begin_pre(alg))
            {
                stinger_alg_end_pre(alg);
            }
            S = alg->stinger;
            cout << "--------------- Static ------------------" << endl;
            tic();
            cstate_s.init(S);
            cstate_s.agglomerate("static");
            time_s = toc();
            cout << "--------------- Dynamic ------------------" << endl;
            tic();
            cstate_d.add_batch(alg);
            cstate_d.agglomerate();
            time_d = toc();
            modularity_s = cstate_s.update_mod();
            modularity_d = cstate_d.update_mod();
            soc_s = cstate_s.SizeofChange(S);
            soc_d = cstate_d.SizeofChange(S);
            ncoms_s = cstate_s.ncoms();
            ncoms_d = cstate_d.ncoms();
            max_csize_s = cstate_s.max_comm_size();
            max_csize_d = cstate_d.max_comm_size();
            min_csize_s = cstate_s.min_comm_size();
            min_csize_d = cstate_d.min_comm_size();
            mean_csize_s = cstate_s.mean_comm_size();
            mean_csize_d = cstate_d.mean_comm_size();
            stddev_csize_s = cstate_s.stddev_comm_size();
            stddev_csize_d = cstate_d.stddev_comm_size();
            nedges_s = cstate_s.Diagnostic["nedges"];
            nedges_d = cstate_d.Diagnostic["nedges"];
            nedit_edge_s = cstate_s.Diagnostic["nedit_edge"];
            nedit_edge_d = cstate_d.Diagnostic["nedit_edge"];
            nmerge_s = cstate_s.Diagnostic["nmerge"];
            nmerge_d = cstate_d.Diagnostic["nmerge"];
            nsplit_s = cstate_s.Diagnostic["nsplit"];
            nsplit_d = cstate_d.Diagnostic["nsplit"];
            nedit_edgeedits_s = cstate_s.Diagnostic["nedit_edgeedits"];
            nedit_edgeedits_d = cstate_d.Diagnostic["nedit_edgeedits"];
            nmergeedits_s = cstate_s.Diagnostic["nmergeedits"];
            nmergeedits_d = cstate_d.Diagnostic["nmergeedits"];
            nsplitedits_s = cstate_s.Diagnostic["nsplitedits"];
            nsplitedits_d = cstate_d.Diagnostic["nsplitedits"];
            split_time_d = cstate_d.Diagnostic["split_time"];
            merge_time_s = cstate_s.Diagnostic["merge_time"];
            merge_time_d = cstate_d.Diagnostic["merge_time"];
            editedge_time_d = cstate_d.Diagnostic["editedge_time"];
            agglomerate_time_s = cstate_s.Diagnostic["agglomerate_time"];
            agglomerate_time_d = cstate_d.Diagnostic["agglomerate_time"];
            match_time_s = cstate_s.Diagnostic["match_time"];
            match_time_d = cstate_d.Diagnostic["match_time"];
            addbatch_time_d = cstate_d.Diagnostic["addbatch_time"];
            init_time_s = cstate_s.Diagnostic["init_time"];
            cout << "Iter: " << n+1 << endl;
            cout << "Modularity (Static): " << modularity_s << endl;
            cout << "Modularity (Dynamic): " << modularity_d << endl;
            cout << "SOC (Static): " << soc_s << endl;
            cout << "SOC (Dynamic): " << soc_d << endl;
            cout << "Ncoms (Static): " << ncoms_s << endl;
            cout << "Ncoms (Dynamic): " << ncoms_d << endl;
            cout << "Max csize (Static): " << max_csize_s << endl;
            cout << "Max csize (Dynamic): " << max_csize_d << endl;
            cout << "Min csize (Static): " << min_csize_s << endl;
            cout << "Min csize (Dynamic): " << min_csize_d << endl; 
            cout << "Mean csize (Static): " << mean_csize_s << endl;
            cout << "Mean csize (Dynamic): " << mean_csize_d << endl;
            cout << "Stddev csize (Static): " << stddev_csize_s << endl;
            cout << "Stddev csize (Dynamic): " << stddev_csize_d << endl;
            cout << "Nedges (Static): " << nedges_s << endl;
            cout << "Nedges (Dynamic): " << nedges_d << endl;
            cout << "Nedit_edge (Static): " << nedit_edge_s << endl;
            cout << "Nedit_edge (Dynamic): " << nedit_edge_d << endl;
            cout << "Nmerge (Static): " << nmerge_s << endl;
            cout << "Nmerge (Dynamic): " << nmerge_d << endl;
            cout << "Nsplit (Static): " << nsplit_s << endl;
            cout << "Nsplit (Dynamic): " << nsplit_d << endl;
            cout << "Nedit_edge edits (Static): " << nedit_edgeedits_s << endl;
            cout << "Nedit_edge edits (Dynamic): " << nedit_edgeedits_d << endl;
            cout << "Nmerge edits (Static): " << nmergeedits_s << endl;
            cout << "Nmerge edits (Dynamic): " << nmergeedits_d << endl;
            cout << "Nsplit edits (Static): " << nsplitedits_s << endl;
            cout << "Nsplit edits (Dynamic): " << nsplitedits_d << endl;
            cout << "Time (Static): " << time_s << endl;
            cout << "Time (Dynamic): " << time_d << endl;
            cout << "Init time (Static): " << init_time_s << endl;
            cout << "Add Batch time (Dynamic): " << addbatch_time_d << endl;
            cout << "Split time (Dynamic): " << split_time_d << endl;
            cout << "Edit Edge time (Dynamic): " << editedge_time_d << endl;
            cout << "Agglomerate time (Static): " << agglomerate_time_s << endl;
            cout << "Agglomerate time (Dynamic): " << agglomerate_time_d << endl;
            cout << "Match time (Static): " << match_time_s << endl;
            cout << "Match time (Dynamic): " << match_time_d << endl;
            cout << "Merge time (Static): " << merge_time_s << endl;
            cout << "Merge time (Dynamic): " << merge_time_d << endl;
            modularities_s.push_back(modularity_s);
            modularities_d.push_back(modularity_d);
            SOC_s.push_back(soc_s);
            SOC_d.push_back(soc_d);
            Ncoms_s.push_back(ncoms_s);
            Ncoms_d.push_back(ncoms_d);
            Max_csize_s.push_back(max_csize_s);
            Max_csize_d.push_back(max_csize_d);
            Min_csize_s.push_back(min_csize_s);
            Min_csize_d.push_back(min_csize_d);
            Mean_csize_s.push_back(mean_csize_s);
            Mean_csize_d.push_back(mean_csize_d);
            Stddev_csize_s.push_back(stddev_csize_s);
            Stddev_csize_d.push_back(stddev_csize_d);
            Nedges_s.push_back(nedges_s);
            Nedges_d.push_back(nedges_d);
            Nedit_edge_s.push_back(nedit_edge_s);
            Nedit_edge_d.push_back(nedit_edge_d);
            Nmerge_s.push_back(nmerge_s);
            Nmerge_d.push_back(nmerge_d);
            Nsplit_s.push_back(nsplit_s);
            Nsplit_d.push_back(nsplit_d);
            Nedit_edgeedits_s.push_back(nedit_edgeedits_s);
            Nedit_edgeedits_d.push_back(nedit_edgeedits_d);
            Nmergeedits_s.push_back(nmergeedits_s);
            Nmergeedits_d.push_back(nmergeedits_d);
            Nsplitedits_s.push_back(nsplitedits_s);
            Nsplitedits_d.push_back(nsplitedits_d);
            Time_s.push_back(time_s);
            Time_d.push_back(time_d);
            Init_time_s.push_back(init_time_s);
            Addbatch_time_d.push_back(addbatch_time_d);
            Split_time_d.push_back(split_time_d);
            Editedge_time_d.push_back(editedge_time_d);
            Agglomerate_time_s.push_back(agglomerate_time_s);
            Agglomerate_time_d.push_back(agglomerate_time_d);
            Match_time_s.push_back(match_time_s);
            Match_time_d.push_back(match_time_d);
            Merge_time_s.push_back(merge_time_s);
            Merge_time_d.push_back(merge_time_d);
            
            if(stinger_alg_begin_post(alg))
            {
                stinger_alg_end_post(alg);
            }
            n++;
        }
        else
        {
            cout << "Alg not enabled!" << endl;
            break;
        }
    }
    cout << "Bno. Mod(s) Mod(d) SOC(s) SOC(d) Ncoms(s) Ncoms(d) MaxCs(s) MaxCs(d) MinCs(s) MinCs(d) MeanCs(s) MeanCs(d) StddevCs(s) StddevCs(d) Time(s) Time(d) InitTime(s) AddBatchTime(d) SplitTime(d) EditEdgeTime(d) AgglomerateTime(s) AgglomerateTime(d) MatchTime(s) MatchTime(d) MergeTime(s) MergeTime(d) Nedges(s) Nedges(d) Nedit_edge(s) Nedit_edge(d) Nmerge(s) Nmerge(d) Nsplit(s) Nsplit(d) Nedit_edgeedits(s) Nedit_edgeedits(d) Nmergeedits(s) Nmergeedits(d) Nsplitedits(s) Nsplitedits(d)" << endl;
    for(int i=0; i<modularities_s.size(); i++)
    {
        cout << i << " " << modularities_s[i] << " " << modularities_d[i] << " " << SOC_s[i] << " " << SOC_d[i] << " " << Ncoms_s[i] << " " << Ncoms_d[i] << " " << Max_csize_s[i] << " " << Max_csize_d[i] << " " << Min_csize_s[i] << " " << Min_csize_d[i] << " " << Mean_csize_s[i] << " " << Mean_csize_d[i] << " " << Stddev_csize_s[i] << " " << Stddev_csize_d[i] << " " << Time_s[i] << " " << Time_d[i] << " " << Init_time_s[i] << " " << Addbatch_time_d[i] << " " << Split_time_d[i] << " " << Editedge_time_d[i] << " " << Agglomerate_time_s[i] << " " << Agglomerate_time_d[i] << " " << Match_time_s[i] << " " << Match_time_d[i] << " " << Merge_time_s[i] << " " << Merge_time_d[i] << " " << Nedges_s[i] << " " << Nedges_d[i] << " " << Nedit_edge_s[i] << " " << Nedit_edge_d[i] << " " << Nmerge_s[i] << " " << Nmerge_d[i] << " " << Nsplit_s[i] << " " << Nsplit_d[i] << " " << Nedit_edgeedits_s[i] << " " << Nedit_edgeedits_d[i] << " " << Nmergeedits_s[i] << " " << Nmergeedits_d[i] << " " << Nsplitedits_s[i] << " " << Nsplitedits_d[i] << endl;
    }
}
