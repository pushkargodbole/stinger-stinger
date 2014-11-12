#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <assert.h>

#include "stinger_core/stinger.h"
#include "stinger_core/xmalloc.h"
#include "stinger_core/stinger_error.h"
#include "stinger_net/stinger_alg.h"
#include "stinger_utils/timer.h"

#include "compat.h"
#include "pagerank.h"
#include "spmspv.h"
#include "spmspv_ompcas_batch.h"

static inline int append_to_vlist (int64_t * restrict nvlist,
                                   int64_t * restrict vlist,
                                   int64_t * restrict mark /* |V|, init to -1 */,
                                   const int64_t i);
void clear_vlist (int64_t * restrict nvlist,
                  int64_t * restrict vlist,
                  int64_t * restrict mark);

static inline void normalize_pr (const int64_t nv, double * restrict pr_val);
static int64_t find_max_pr (const int64_t nv, const double * restrict pr_val);

static int nonunit_weights = 1;

struct spvect {
  int64_t nv;
  int64_t * idx;
  double * val;
};
static inline struct spvect alloc_spvect (int64_t nvmax);

static void gather_pre (const stinger_registered_alg * alg, struct spvect * x,
                        int64_t * restrict mark);
static void dpr_pre (const int64_t nv, struct stinger * S, struct spvect * restrict b,
                     struct spvect x, const double * restrict pr_val,
                     int64_t * restrict mark, double * restrict dzero_workspace,
                     int64_t * restrict pr_vol);
static void dpr_held_pre (const int64_t nv, struct stinger * S, struct spvect * restrict b,
                          struct spvect * restrict x_held,
                          struct spvect x, const double * restrict pr_val,
                          int64_t * restrict mark, double * restrict dzero_workspace,
                          int64_t * restrict pr_vol);
static void dpr_update (const int64_t nv, struct stinger * S, struct spvect * x, struct spvect * b,
                        struct spvect * dpr, double * restrict pr_val, const double alpha,
                        const int maxiter, int64_t * restrict mark, int64_t * restrict iworkspace,
                        double * workspace, double * dzero_workspace, int * niter, int64_t * pr_vol);
static void dpr_held_update (const int64_t nv, struct stinger * S, struct spvect * x, struct spvect * b,
                             struct spvect * dpr, double * restrict pr_val, const double alpha,
                             const int maxiter, int64_t * restrict mark, int64_t * restrict iworkspace,
                             double * workspace, double * dzero_workspace, int * niter, int64_t * pr_vol);
static void pr_update (const int64_t nv, const int64_t NE, struct stinger * S, double * restrict pr_val, double * restrict v, const double alpha, const int maxiter, double * workspace, int * niter, int64_t * pr_vol);
static void rpr_update (const int64_t nv, const int64_t NE, struct stinger * S, double * restrict pr_val, double * restrict v, const double alpha, const int maxiter, double * workspace, int * niter, int64_t * pr_vol);


#define BASELINE 0
#define RESTART 1
#define DPR 2
#define DPR_HELD 3
#define NPR_ALG 4

static int n_alg_run = 0;
static int run_alg[NPR_ALG] = {0};
static const char *desc_alg[NPR_ALG] = {
  "pr delta-pr",
  "rpr delta-rpr",
  "dpr delta-dpr",
  "dprheld delta-dprheld",
};
static const char *pr_name[NPR_ALG] = {"pr", "pr_restart", "dpr", "dpr_held"};

int
main(int argc, char *argv[])
{
  size_t data_desc_len;
  char * data_desc;

  double init_time, gather_time;
  double cwise_err;

  int iter = 0;
  int64_t nv;
  int64_t max_nv;

  double * pr_val[NPR_ALG] = { NULL };
  int pr_val_init = NPR_ALG - 1; /* initialize into this one... */
  double * old_pr_val[NPR_ALG] = { NULL }; /* Only used for BASELINE and RESTART */
  double * pr_val_delta[NPR_ALG] = { NULL };
  int niter[NPR_ALG] = { 0 };
  double pr_pre_time[NPR_ALG] = { 0.0 };
  double pr_time[NPR_ALG] = { 0.0 };
  int64_t pr_vol[NPR_ALG] = { 0 };
  int64_t pr_nupd[NPR_ALG] = { 0 };

  double alpha = 0.15;
  int maxiter = 100;

  struct spvect x;
  struct spvect x_held;
  struct spvect b;
  struct spvect b_held;
  struct spvect dpr;
  struct spvect dpr_held;

  int64_t * mark;
  double * v;
  double * dzero_workspace;
  double * workspace;
  int64_t * iworkspace;

  for (int k = 1; k < argc; ++k) {
    if (0 == strcmp(argv[k], "--help") || 0 == strcmp(argv[k], "-h")) {
      fprintf (stderr,
               "pagerank_updating [--run-pr] [--run-rpr] [--run-dpr] [--run-dprheld]\n");
      return EXIT_SUCCESS;
    } else if (0 == strcmp(argv[k], "--run-pr")) {
      run_alg[BASELINE] = 1;
    } else if (0 == strcmp(argv[k], "--run-rpr")) {
      run_alg[RESTART] = 1;
    } else if (0 == strcmp(argv[k], "--run-dpr")) {
      run_alg[DPR] = 1;
    } else if (0 == strcmp(argv[k], "--run-dprheld")) {
      run_alg[DPR_HELD] = 1;
    }
  }

  /* Build the data description */
  n_alg_run = 0;
  for (int k = 0; k < NPR_ALG; ++k) {
    if (run_alg[k] && k < pr_val_init) pr_val_init = k; /* use the first one... */
    n_alg_run += run_alg[k];
  }
  if (!n_alg_run) {
    n_alg_run = NPR_ALG;
    /* If none specified, run them all... */
    for (int k = 0; k < NPR_ALG; ++k) run_alg[k] = 1;
  }
  data_desc_len = 1; /* first space */
  for (int k = 0; k < NPR_ALG; ++k)
    if (run_alg[k]) /* Null byte counts as joining space */
      data_desc_len += 3 + strlen (desc_alg[k]);
  data_desc = calloc (data_desc_len+1, sizeof (*data_desc));
  {
    size_t off = 0;
    for (int k = 0; k < NPR_ALG; ++k)
      if (run_alg[k]) { data_desc[off++] = 'd'; data_desc[off++] = 'd'; }

    for (int k = 0; k < NPR_ALG; ++k)
      if (run_alg[k]) {
        data_desc[off++] = ' ';
        strcpy (&data_desc[off], desc_alg[k]);
        off += strlen (desc_alg[k]);
      }
  }
  /* fprintf (stderr, "n_alg_run = %d\n", n_alg_run); */
  /* fprintf (stderr, "data_desc = \"%s\"\n", data_desc); */
  /* fprintf (stderr, "initial into %d \"%s\"\n", pr_val_init, desc_alg[pr_val_init]); */

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
   * Setup and register algorithm with the server
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  stinger_registered_alg * alg =
    stinger_register_alg(
      .name="pagerank_updating",
      .data_per_vertex=2*n_alg_run*sizeof(double),
      .data_description=data_desc,
      .host="localhost",
    );

  if(!alg) {
    LOG_E("Registering algorithm failed.  Exiting");
    return -1;
  }

  /* Assume the number of vertices does not change. */
  //nv = stinger_max_active_vertex (alg->stinger) + 1;
  max_nv = stinger_max_nv (alg->stinger);
  nv = max_nv;
  if (getenv ("NV")) {
    long nv_env = strtol (getenv ("NV"), NULL, 10);
    if (nv_env > 0 && nv_env < nv)
      nv = nv_env;
  }
  fprintf (stderr, "NV: %ld   of %ld\n", (long)nv, (long)max_nv);

  {
    int off = 0;
    for (int k = 0; k < NPR_ALG; ++k)
      if (run_alg[k]) {
        double * alg_data = (double*)alg->alg_data;
        pr_val[k] = &alg_data[2 * off * max_nv];
        pr_val_delta[k] = &alg_data[(1 + 2 * off) * max_nv];
        if (k <= RESTART) old_pr_val[k] = calloc (nv, sizeof (*old_pr_val[k]));
        ++off;
      }
  }

  mark = xmalloc (nv * sizeof (*mark));
  v = xmalloc (nv * sizeof (*v));
  dzero_workspace = xcalloc (nv, sizeof (*dzero_workspace));
  workspace = xmalloc (2 * nv * sizeof (*workspace));
  iworkspace = xmalloc (nv * sizeof (*iworkspace));
  x = alloc_spvect (nv);
  x_held = alloc_spvect (nv);
  b = alloc_spvect (nv);
  b_held = alloc_spvect (nv);
  dpr = alloc_spvect (nv);
  dpr_held = alloc_spvect (nv);

  OMP(omp parallel) {
    for (int k = 0; k < NPR_ALG; ++k)
      if (run_alg[k]) {
        OMP(omp for OMP_SIMD nowait)
          for (int64_t i = 0; i < nv; ++i) {
            pr_val[k][i] = 0.0;
            pr_val_delta[k][i] = 0.0;
          }
      }
    OMP(omp for OMP_SIMD nowait)
      for (int64_t i = 0; i < nv; ++i)
        mark[i] = -1;
    OMP(omp for OMP_SIMD nowait)
      for (int64_t i = 0; i < nv; ++i)
        v[i] = 1.0;
  }

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
   * Initial static computation
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  stinger_alg_begin_init(alg); {
    if (stinger_total_edges (alg->stinger)) {
      tic ();
      niter[pr_val_init] = pagerank (nv, alg->stinger, pr_val[pr_val_init], v, alpha, maxiter, workspace);
      pr_time[pr_val_init] = toc ();
#if !defined(NDEBUG)
      for (int64_t k = 0; k < nv; ++k)
        assert (!isnan(pr_val[pr_val_init][k]));
#endif
      normalize_pr (nv, pr_val[pr_val_init]);
#if !defined(NDEBUG)
      for (int64_t k = 0; k < nv; ++k)
        assert (!isnan(pr_val[pr_val_init][k]));
#endif
      /* Copy the starting vector */
      for (int k = 0; k < NPR_ALG; ++k)
        if (run_alg[k] && k != pr_val_init) {
          OMP(omp parallel for OMP_SIMD)
            for (int64_t i = 0; i < nv; ++i)
              pr_val[k][i] = pr_val[pr_val_init][i];
        }
    }
  } stinger_alg_end_init(alg);


  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
   * Streaming Phase
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  while(alg->enabled) {
    ++iter;
    /* Batch pre-processing, only needed for DPR and DPR_HELD but computed for all */
    if(stinger_alg_begin_pre(alg)) {
      clear_vlist (&x.nv, x.idx, mark);
      OMP(omp parallel for OMP_SIMD)
        for (int64_t k = 0; k < nv; ++k) mark[k] = -1;
      /* Gather vertices that are affected by the batch. */
      /* (Really should be done in the framework, but ends up being handy here...) */
      for (int k = 1; k < NPR_ALG; ++k) {
        pr_time[k] = -1;
        pr_vol[k] = 0;
        pr_nupd[k] = 0;
      }

      tic ();
      gather_pre (alg, &x, mark);
      gather_time = toc ();

      if (run_alg[DPR]) {
        tic ();
        dpr_pre (nv, alg->stinger, &b, x, pr_val[DPR], mark, dzero_workspace, &pr_vol[DPR]);
        pr_pre_time[DPR] = toc ();
      }

      /* Compute x, b0 for dpr_held */
      if (run_alg[DPR_HELD]) {
        tic ();
        dpr_held_pre (nv, alg->stinger, &b_held, &x_held, x, pr_val[DPR_HELD], mark, dzero_workspace, &pr_vol[DPR_HELD]);
        pr_pre_time[DPR_HELD] = toc ();
      }

      stinger_alg_end_pre(alg);
    }

    fprintf (stderr, "%ld: NV CHANGED %ld    %ld %ld\n", (long)iter, (long)x.nv,
             (long)alg->num_insertions, (long)alg->num_deletions);

    /* Post processing */
    if(stinger_alg_begin_post(alg)) {
      struct stinger * S = alg->stinger;
      const int64_t NE = stinger_total_edges (S);

      /* Compute the change in PageRank */
      if (run_alg[DPR]) {
        tic ();
        dpr_update (nv, S, &x, &b, &dpr, pr_val[DPR], alpha, maxiter, mark, iworkspace, workspace, dzero_workspace, &niter[DPR], &pr_vol[DPR]);
        pr_time[DPR] = toc () + pr_pre_time[DPR];
        pr_nupd[DPR] = dpr.nv;

        OMP(omp parallel for OMP_SIMD)
          for (int64_t k = 0; k < dpr.nv; ++k) {
            assert (mark[dpr.idx[k]] != -1);
            mark[dpr.idx[k]] = -1;
          }
      }

      if (run_alg[DPR_HELD]) {
        /* Compute the change in PageRank, holding back small changes */
        tic ();
        dpr_update (nv, S, &x_held, &b_held, &dpr_held, pr_val[DPR_HELD], alpha, maxiter, mark, iworkspace, workspace, dzero_workspace, &niter[DPR_HELD], &pr_vol[DPR_HELD]);
        pr_time[DPR_HELD] = toc () + pr_pre_time[DPR_HELD];
        pr_nupd[DPR_HELD] = dpr_held.nv;
      }

      if (run_alg[BASELINE]) {
        /* Compute PageRank from scratch */
        tic ();
        pr_update (nv, NE, S, pr_val[BASELINE], v, alpha, maxiter, workspace, &niter[BASELINE], &pr_vol[BASELINE]);
        pr_time[BASELINE] = toc ();
        pr_nupd[BASELINE] = nv;
      }

      if (run_alg[RESTART]) {
        /* Restart using the previous vector, cumulative... */
        tic ();
        rpr_update (nv, NE, S, pr_val[RESTART], v, alpha, maxiter, workspace, &niter[RESTART], &pr_vol[RESTART]);
        pr_time[RESTART] = toc ();
        pr_nupd[RESTART] = nv;
      }

      stinger_alg_end_post(alg);
    }

    int64_t loc_baseline;
    for (int alg = 0; alg < NPR_ALG; ++alg)
      if (run_alg[alg]) {
        double err = (run_alg[BASELINE]? 0.0 : NAN), mxerr = 0.0;
        int64_t where = -1, nerr = 0;
        const int64_t loc = find_max_pr (nv, pr_val[alg]);
        if (alg == BASELINE) loc_baseline = loc;
        if (alg > 0 && run_alg[BASELINE])
          OMP(omp parallel for OMP_SIMD reduction(+: err))
            for (int64_t i = 0; i < nv; ++i) {
              const double ei = fabs (pr_val[alg][i] - pr_val[BASELINE][i]);
              if (ei > mxerr) { where = i; mxerr = ei; }
              err += ei;
              if (ei > 0.0) ++nerr;
            }
        fprintf (stderr, "%ld: %12s %d\t%d %10g %8ld %10e %8ld %8ld\n",
                 (long)iter, pr_name[alg], alg,
                 niter[alg], pr_time[alg], pr_vol[alg], err, (long)nerr, (long)pr_nupd[alg]);
      }
  }

  LOG_I("Algorithm complete... shutting down");
}

/* Utility functions */

MTA("mta inline") MTA("mta expect parallel context")
int
append_to_vlist (int64_t * restrict nvlist,
                 int64_t * restrict vlist,
                 int64_t * restrict mark /* |V|, init to -1 */,
                 const int64_t i)
{
  int out = 0;

  /* for (int64_t k = 0; k < *nvlist; ++k) { */
  /*   assert (vlist[k] >= 0); */
  /*   if (k != mark[vlist[k]]) */
  /*     fprintf (stderr, "wtf mark[vlist[%ld] == %ld] == %ld\n", */
  /*              (long)k, (long)vlist[k], (long)mark[vlist[k]]); */
  /*   assert (mark[vlist[k]] == k); */
  /* } */

  if (mark[i] < 0) {
    int64_t where;
    if (bool_int64_compare_and_swap (&mark[i], -1, -2)) {
      assert (mark[i] == -2);
      /* Own it. */
      where = int64_fetch_add (nvlist, 1);
      mark[i] = where;
      vlist[where] = i;
      out = 1;
    }
  }

  return out;
}

void
clear_vlist (int64_t * restrict nvlist,
             int64_t * restrict vlist,
             int64_t * restrict mark)
{
  const int64_t nvl = *nvlist;

  OMP(omp parallel for)
    for (int64_t k = 0; k < nvl; ++k)
      mark[vlist[k]] = -1;

  *nvlist = 0;
}

struct spvect
alloc_spvect (int64_t nvmax)
{
  struct spvect out;
  out.nv = 0;
  out.idx = xmalloc (nvmax * sizeof (*out.idx));
  out.val = xmalloc (nvmax * sizeof (*out.val));
  return out;
}

void
normalize_pr (const int64_t nv, double * restrict pr_val)
{
  double n1 = 0.0;
  OMP(omp parallel) {
    OMP(omp for OMP_SIMD reduction(+: n1))
      for (int64_t k = 0; k < nv; ++k)
        n1 += fabs (pr_val[k]);
    OMP(omp for OMP_SIMD)
      for (int64_t k = 0; k < nv; ++k)
        pr_val[k] /= n1;
  }
}

int64_t
find_max_pr (const int64_t nv, const double * restrict pr_val)
{
  double max_val = -1.0;
  int64_t loc = -1;
  OMP(omp parallel) {
    double t_max_val = -1.0;
    int64_t t_loc = -1;
    OMP(omp for)
      for (int64_t k = 0; k < nv; ++k) {
        const double v = pr_val[k];
        if (v > t_max_val) {
          t_max_val = v;
          t_loc = k;
        }
      }
    OMP(omp critical) {
      if (t_max_val > max_val) {
        max_val = t_max_val;
        loc = t_loc;
      }
    }
  }
  return loc;
}

void
gather_pre (const stinger_registered_alg * alg, struct spvect * x, int64_t * restrict mark)
{
  const int64_t nins = alg->num_insertions;
  const int64_t nrem = alg->num_deletions;
  const stinger_edge_update * restrict ins = alg->insertions;
  const stinger_edge_update * restrict rem = alg->deletions;

  OMP(omp parallel) {
    OMP(omp for nowait)
      for (int64_t k = 0; k < nins; ++k) {
        append_to_vlist (&x->nv, x->idx, mark, ins[k].source);
        append_to_vlist (&x->nv, x->idx, mark, ins[k].destination);
      }
    OMP(omp for)
      for (int64_t k = 0; k < nrem; ++k) {
        append_to_vlist (&x->nv, x->idx, mark, rem[k].source);
        append_to_vlist (&x->nv, x->idx, mark, rem[k].destination);
      }
    OMP(omp for OMP_SIMD)
      for (int64_t k = 0; k < x->nv; ++k)
        mark[x->idx[k]] = -1;
  }
}

void
dpr_pre (const int64_t nv, struct stinger * S, struct spvect * restrict b, struct spvect x, const double * restrict pr_val, int64_t * restrict mark, double * restrict dzero_workspace, int64_t * restrict pr_vol)
{
  int64_t b_nv = 0;
  /* Compute b0 in b */
  OMP(omp parallel) {
    OMP(omp for OMP_SIMD)
      for (int64_t k = 0; k < x.nv; ++k) {
        assert(!isnan(pr_val[x.idx[k]]));
        x.val[k] = pr_val[x.idx[k]];
      }
    stinger_unit_dspmTspv_degscaled_ompcas_batch (nv,
                                                  1.0, S,
                                                  x.nv, x.idx, x.val,
                                                  0.0,
                                                  &b_nv, b->idx, b->val,
                                                  mark, dzero_workspace,
                                                  pr_vol);
    OMP(omp for OMP_SIMD)
      for (int64_t k = 0; k < b_nv; ++k) mark[b->idx[k]] = -1;
  }
  b->nv = b_nv;
}

void
dpr_held_pre (const int64_t nv, struct stinger * S, struct spvect * restrict b, struct spvect * restrict x_held, const struct spvect x, const double * restrict pr_val, int64_t * restrict mark, double * restrict dzero_workspace, int64_t * restrict pr_vol)
{
  int64_t b_nv = 0;
  int64_t x_held_nv = x.nv;
  /* Compute b0 in b */
  OMP(omp parallel) {
    OMP(omp for OMP_SIMD)
      for (int64_t k = 0; k < x.nv; ++k) {
        assert(!isnan(pr_val[x.idx[k]]));
        x_held->idx[k] = x.idx[k];
        x_held->val[k] = pr_val[x.idx[k]];
      }
    stinger_unit_dspmTspv_degscaled_ompcas_batch (nv,
                                                  1.0, S,
                                                  x_held_nv, x_held->idx, x_held->val,
                                                  0.0,
                                                  &b_nv, b->idx, b->val,
                                                  mark, dzero_workspace,
                                                  pr_vol);
    OMP(omp for OMP_SIMD)
      for (int64_t k = 0; k < b_nv; ++k) mark[b->idx[k]] = -1;
  }
  b->nv = b_nv;
  x_held->nv = x_held_nv;
}

void
dpr_update (const int64_t nv, struct stinger * S, struct spvect * x, struct spvect * b, struct spvect * dpr, double * restrict pr_val, const double alpha, const int maxiter, int64_t * restrict mark, int64_t * restrict iworkspace, double * workspace, double * dzero_workspace, int * niter, int64_t * pr_vol)
{
  *niter = pagerank_dpr (nv, S,
                         &x->nv, x->idx, x->val,
                         alpha, maxiter,
                         &b->nv, b->idx, b->val,
                         &dpr->nv, dpr->idx, dpr->val,
                         mark, iworkspace, workspace, dzero_workspace,
                         pr_vol);
  const int64_t dpr_nv = dpr->nv;
  OMP(omp parallel for OMP_SIMD)
    for (int64_t k = 0; k < dpr_nv; ++k)
      pr_val[dpr->idx[k]] += dpr->val[k];
  normalize_pr (nv, pr_val);
}

/* Nearly identical, but useful for profiling. */
void
dpr_held_update (const int64_t nv, struct stinger * S, struct spvect * x, struct spvect * b, struct spvect * dpr, double * restrict pr_val, const double alpha, const int maxiter, int64_t * restrict mark, int64_t * restrict iworkspace, double * workspace, double * dzero_workspace, int * niter, int64_t * pr_vol)
{
  *niter = pagerank_dpr_held (nv, S,
                              &x->nv, x->idx, x->val,
                              alpha, maxiter,
                              &b->nv, b->idx, b->val,
                              &dpr->nv, dpr->idx, dpr->val,
                              mark, iworkspace, workspace, dzero_workspace,
                              pr_vol);
  const int64_t dpr_nv = dpr->nv;
  OMP(omp parallel for OMP_SIMD)
    for (int64_t k = 0; k < dpr_nv; ++k)
      pr_val[dpr->idx[k]] += dpr->val[k];
  normalize_pr (nv, pr_val);
}

void
pr_update (const int64_t nv, const int64_t NE, struct stinger * S, double * restrict pr_val, double * restrict v, const double alpha, const int maxiter, double * workspace, int * niter, int64_t * pr_vol)
{
  *niter = pagerank (nv, S, pr_val, v, alpha, maxiter, workspace);
  normalize_pr (nv, pr_val);
  *pr_vol = *niter * NE;
}

void
rpr_update (const int64_t nv, const int64_t NE, struct stinger * S, double * restrict pr_val, double * restrict v, const double alpha, const int maxiter, double * workspace, int * niter, int64_t * pr_vol)
{
  *niter = pagerank_restart (nv, S, pr_val, v, alpha, maxiter, workspace);
  normalize_pr (nv, pr_val);
  *pr_vol = *niter * NE;
}
