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

#include "community_state.cpp"
#include "utils.h"

bool do_static = false;
bool do_dynamic = false;
int nbatch;

int
main(int argc, char *argv[])
{
    community_state cstate;
    const struct stinger * S;
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "--static") == 0) do_static = true;
        else if (strcmp(argv[i], "--dynamic") == 0) do_dynamic = true;
        else if (0 == strcmp(argv[i], "-nb")) nbatch = atoi(argv[i+1]);
        else if (0 == strcmp(argv[i], "--help") || 0 == strcmp(argv[i], "-h"))
        {
            cout << endl; // Write help!
        }
    }
    
    if(do_static == false && do_dynamic == false) do_dynamic = true;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Setup and register algorithm with the server
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    string alg_name_str = "reagglomeration";
    char * alg_name = new char[alg_name_str.length() + 1];
    strcpy(alg_name, alg_name_str.c_str());
    
    stinger_registered_alg * alg = stinger_register_alg(
        alg_name, // Algorithm name
        "localhost", // Host
        10103, // Port
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
    stinger_alg_begin_init(alg);
    {
        S = alg->stinger;
        cstate.init(S);
    } stinger_alg_end_init(alg);

    //LOG_V_A("simple_communities: init comm time = %g\n", comm_time);

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
     * Streaming Phase
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    while(alg->enabled)
    {
        // Pre processing
        if(stinger_alg_begin_pre(alg))
        {
            stinger_alg_end_pre(alg);
        }
        cout << cstate.dmod(0, 1) << endl;
        cstate.agglomerate();
        
        for(int64_t i=0;i<cstate.nv;i++)
        {
            cout << cstate.find(i) << " ";
        }
        cout << endl;
        double modularity = cstate.update_mod();
        cout << "Modularity: " << modularity << endl;
        if(stinger_alg_begin_post(alg))
        {
            stinger_alg_end_post(alg);
        }
    }
    
}
