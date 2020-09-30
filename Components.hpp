/*
* Copyright (c) 2020 Marcus Ritt, Artur Brum
*
* Permission is hereby granted, free of charge, to any person (the "Person")
* obtaining a copy of this software and associated documentation files (the
* "Software"), to deal in the Software, including the rights to use, copy, modify,
* merge, publish, distribute the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* 1. The above copyright notice and this permission notice shall be included in
*    all copies or substantial portions of the Software.
* 2. Under no circumstances shall the Person be permitted, allowed or authorized
*    to commercially exploit the Software.
* 3. Changes made to the original Software shall be labeled, demarcated or
*    otherwise identified and attributed to the Person.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
* FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
* IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
* CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#pragma once

#include <vector>
#include <random>
#include <limits>
#include "hfs/instance.hpp"
#include "hfs/solution.hpp"
#include "hfs/parameters.hpp"
#include "hfs/constructive.hpp"
#include "hfs/ls.hpp"
#include "hfs/rules.hpp"
#include "hfs/perturbations.hpp"
#include "hfs/beam_search.hpp"

using namespace std::placeholders;

std::function<unsigned(unsigned,unsigned)> getTiebreak(const Instance& I, const unsigned option) {
   if (option == 0) {
      std::vector<char> tiebreak(I.n+1,0);
      getKK1(I, tiebreak);
      return bind(usebool,_1,_2,tiebreak);
   }
   else if (option == 1) {
      std::vector<char> tiebreak(I.n+1,0);
      getKK2(I, tiebreak);
      return bind(usebool,_1,_2,tiebreak);
   }
   else if (option == 2) {
      return firstposition;
   }
   else if (option == 3) {
      return lastposition;
   }
   else {
      return randomtiebreak;
   }
}

std::vector<char> getTiebreakAsVector(const Instance& I, const unsigned option) {
   if (option == 0) {
      std::vector<char> tiebreak(I.n+1,0);
      getKK1(I, tiebreak);
      return tiebreak;
   }
   else if (option == 1) {
      std::vector<char> tiebreak(I.n+1,0);
      getKK2(I, tiebreak);
      return tiebreak;
   }
   else if (option == 2) {
      std::vector<char> tiebreak(I.n+1,1);
      return tiebreak;
   }
   else {
      std::vector<char> tiebreak(I.n+1,0);
      return tiebreak;
   }
}

struct SGEDecider {
   unsigned dv_initial_sol,
            dv_order,
            dv_tiebreak_initial,
            dv_tiebreak_perturb,
            dv_tiebreak_ls1,
            dv_tiebreak_ls2,
            dv_perturb,
            dv_local_search,
            dv_ls_procedure1,
            dv_ls_procedure2,
            dv_acceptance,
            dv_tiebreak_perturb_np,
            dv_tiebreak_ls1_np,
            dv_tiebreak_ls2_np,
            dv_local_search_np,
            dv_ls_procedure1_np,
            dv_ls_procedure2_np,
            dv_acceptance_np,
            dv_perturb_np;
   std::function<unsigned(unsigned,unsigned)> tiebreak_initial,
                                              tiebreak_perturb,
                                              tiebreak_ls1,
                                              tiebreak_ls2,
                                              tiebreak_perturb_np, 
                                              tiebreak_ls1_np, 
                                              tiebreak_ls2_np;

   SGEDecider(const Instance& I, const std::vector<unsigned>& v, const char obj, const bool np) {
      dv_initial_sol = v[0];
      dv_order = v[1];
      dv_tiebreak_initial = v[2];
      dv_tiebreak_perturb = v[3];
      dv_tiebreak_ls1 = v[4];
      dv_tiebreak_ls2 = v[5];
      dv_perturb = v[6];
      dv_local_search = v[7];
      dv_ls_procedure1 = v[8];
      dv_ls_procedure2 = v[9];
      dv_acceptance = v[10];
      dv_tiebreak_perturb_np = v[11];
      dv_tiebreak_ls1_np = v[12];
      dv_tiebreak_ls2_np = v[13];
      dv_local_search_np = v[14];
      dv_ls_procedure1_np = v[15];
      dv_ls_procedure2_np = v[16];
      dv_acceptance_np = v[17];
      dv_perturb_np = v[18];
      tiebreak_initial = getTiebreak(I, dv_tiebreak_initial);
      tiebreak_perturb = getTiebreak(I, dv_tiebreak_perturb);
      tiebreak_ls1 = getTiebreak(I, dv_tiebreak_ls1);
      tiebreak_ls2 = getTiebreak(I, dv_tiebreak_ls2);
      tiebreak_perturb_np = getTiebreak(I, dv_tiebreak_perturb_np);
      tiebreak_ls1_np = getTiebreak(I, dv_tiebreak_ls1_np);
      tiebreak_ls2_np = getTiebreak(I, dv_tiebreak_ls2_np);

   }
};

void job_order(const Instance& I, std::vector<Job>& pi, SGEDecider& d, const char obj){
   if (obj == 'm' || obj =='f') {
      if (d.dv_order == 0) {
         pi = getTotalTimeOrder(I);
      }
      else if (d.dv_order == 1) {
         pi = getTotalTimeOrder(I);
         std::reverse(pi.begin()+1,pi.end());
      }
      else if (d.dv_order == 2) {
         std::vector<char> first(I.n+1);
         pi = getKK1(I,first);
      }
      else if (d.dv_order == 3) {
         std::vector<char> first(I.n+1);
         pi = getKK2(I,first);
      }
   }
   else if (obj == 't') {
      if (d.dv_order == 0) {
         pi = getEDD(I);
      }
   }
   else {
      std::cout << "Unrecognized objective!" << std::endl;
   }
}

PFSSolution initial_sol(Instance& I, SGEDecider& d, Params& params, std::mt19937& rng, builder& construct, const std::chrono::system_clock::time_point& start) {
   PFSSolution cs(I);

   std::vector<Job> pi;
   job_order(I, pi, d, params.obj);
   
   if (params.obj == 'm') {
      switch (d.dv_initial_sol) {
         case 0:
            cs = construct.construct(I,pi,d.tiebreak_initial);
            break;
         case 1:
            cs = construct.FRB5(I,pi,params.obj,d.tiebreak_initial,params.n_ls,start,params.timelimit);
            break;
         case 2:
            cs = construct.repeatedConstr(I,rng,params.obj,params.rcon,start,params.timelimit);
            break;
         case 3:
            cs = construct.nehFF(I);
            break;
         case 4:
            cs = construct.FRB5_FF(I,pi,params.obj,params.n_ls,start,params.timelimit);
            break;
         case 5:
            cs = construct.G8(I,pi,params.obj,d.tiebreak_initial,start,params.timelimit);
            break;
         case 6:
            cs = construct.G8_FF(I,pi,params.obj,start,params.timelimit);
            break;
      }
      cs.computeMakespan(I);
      if (params.both) {
         I.reverse();
         PFSSolution csr(I);
         switch (d.dv_initial_sol) {
            case 0:
               csr = construct.construct(I,pi,d.tiebreak_initial);
               break;
            case 1:
               csr = construct.FRB5(I,pi,params.obj,d.tiebreak_initial,params.n_ls,start,params.timelimit);
               break;
            case 2:
               csr = construct.repeatedConstr(I,rng,params.obj,params.rcon,start,params.timelimit);
               break;
            case 3:
               csr = construct.nehFF(I);
               break;
            case 4:
               cs = construct.FRB5_FF(I,pi,params.obj,params.n_ls,start,params.timelimit);
               break;
            case 5:
               cs = construct.G8(I,pi,params.obj,d.tiebreak_initial,start,params.timelimit);
               break;
            case 6:
               cs = construct.G8_FF(I,pi,params.obj,start,params.timelimit);
               break;
         }
         csr.computeMakespan(I);
         if (csr.ms < cs.ms) {
            csr.reverse();
            cs = csr;
         }
         I.reverse();
      }
   }
   
   else if (params.obj == 'f') {
      switch(d.dv_initial_sol) {
         case 0:
            cs = LR(I,params.x);
            break;
         case 1:
            cs = construct.constructCsum(I,pi,d.tiebreak_initial); 
            break;
         case 2:
            cs = construct.FRB5(I,pi,params.obj,d.tiebreak_initial,params.n_ls,start,params.timelimit);
            break;
         case 3:
            beam_search(I, cs, params.w);
            break;
         case 4:
            cs = construct.repeatedConstr(I,rng,params.obj,params.rcon,start,params.timelimit);
            break;

      }   
      cs.computeFlowtime(I);   
   }
   else {
      switch (d.dv_initial_sol) {
         case 0:
            cs = construct.constructTT(I,pi,d.tiebreak_initial);
            break;
         case 1:
            cs = construct.constructTT_IT1(I,pi);
            break;
         case 2:
            cs = construct.constructTT_IT2(I,pi);
            break;
         case 3:
            cs = construct.FRB5(I,pi,params.obj,d.tiebreak_initial,params.n_ls,start,params.timelimit);
            break;
         case 4:
            cs = construct.FRB5_IT1(I,pi,params.obj,params.n_ls,start,params.timelimit);
            break;
         case 5:
            cs = construct.FRB5_IT2(I,pi,params.obj,params.n_ls,start,params.timelimit);
            break;
         case 6:
            cs = construct.repeatedConstrTT(I,rng,params.obj,params.rcon,start,params.timelimit);
            break;
         case 7:
            beam_search_tt(I, cs, params.w);
            break;
      }
   }
   
   return cs;
}

inline void iga_step(Instance& I, SGEDecider& d, Params& params, PFSSolution& n, std::mt19937& rng, builder& construct, const std::chrono::system_clock::time_point& start) {
   switch (d.dv_perturb) {
      case 0:
         ri(I, n, rng, construct, d.tiebreak_perturb, params);
         break;
      case 1:
         gi(I, n, rng, construct, d.tiebreak_perturb, params);
         break;
      case 2:
         rs(I, n, rng, construct, d.tiebreak_perturb, params);
         break;
      case 3:
         gs(I, n, rng, construct, d.tiebreak_perturb, params);
         break;
      case 4:
         ras(I, n, rng, construct, d.tiebreak_perturb, params);
         break;
      case 5:
         gi_asls(I, n, rng, construct, d.tiebreak_perturb, params);
         break;
      case 6:
         ils_gi(I, n, rng, construct, d.tiebreak_perturb, params);
         break;
      case 7:
         gi_ils(I, n, rng, construct, d.tiebreak_perturb, params);
         break;
      case 8:
         giFF(I, n, rng, construct, params);
         break;
      case 9:
         gi_aslsFF(I, n, rng, construct, params);
         break;
      case 10:
         ils_giFF(I, n, rng, construct, params);
         break;
      case 11:
         gi_ilsFF(I, n, rng, construct, params);
         break;
      case 12:
         giIT1(I, n, rng, construct, params);
         break;
      case 13:
         gi_aslsIT1(I, n, rng, construct, params);
         break;
      case 14:
         ils_giIT1(I, n, rng, construct, params);
         break;
      case 15:
         gi_ilsIT1(I, n, rng, construct, params);
         break;
      case 16:
         giIT2(I, n, rng, construct, params);
         break;
      case 17:
         gi_aslsIT2(I, n, rng, construct, params);
         break;
      case 18:
         ils_giIT2(I, n, rng, construct, params);
         break;
      case 19:
         gi_ilsIT2(I, n, rng, construct, params);
         break;
   }
}

inline void local_search_procedure_cmax(Instance& I, SGEDecider& d, const unsigned d_index, Params& params, PFSSolution& n, std::mt19937& rng, builder& construct, const std::chrono::system_clock::time_point& start) {
   unsigned v;
   std::function<unsigned(unsigned,unsigned)> tiebreak_func;
   if (d_index == 0) {
      v = d.dv_ls_procedure1;
      tiebreak_func = d.tiebreak_ls1;
   }
   else {
      v = d.dv_ls_procedure2;
      tiebreak_func = d.tiebreak_ls2;
   }

   switch (v) {
      case 0: 
         insertion(I, n, construct, tiebreak_func, params);
         break;
      case 1:
         N5(I, n, params);
         break;
      case 2:
         Pc(I, n, construct, tiebreak_func, params);
         break;
      case 3:
         insertionFF(I, n, construct, params);
         break;
   }
}

inline void local_search_procedure_csum(Instance& I, SGEDecider& d, const unsigned d_index, Params& params, PFSSolution& n, std::mt19937& rng, builder& construct, const std::chrono::system_clock::time_point& start) {
   unsigned v;
   std::function<unsigned(unsigned,unsigned)> tiebreak_func;
   if (d_index == 0) {
      v = d.dv_ls_procedure1;
      tiebreak_func = d.tiebreak_ls1;
   }
   else {
      v = d.dv_ls_procedure2;
      tiebreak_func = d.tiebreak_ls2;
   }

   switch (v) {
      case 0:
         insertion(I, n, construct, tiebreak_func, params);
         break;   
      case 1:
         fpe(I, n, start, params.timelimit, params.x);
         break;
      case 2:
         swapTasgetiren(I, n, start, params.timelimit);
         break;
      case 3:
         swapExp(I, n, start, params.timelimit);
         break;
      case 4:
         insertTasgetiren(I, n, start, params.timelimit);
         break;
      case 5:
         lsTasgetiren(I, n, start, params.timelimit);
         break;
      case 6:
         bpe(I, n, start, params.timelimit, params.x);
         break;
      case 7:
         iRZ(I, n, start, params.timelimit, params.n_ls);
         break;
      case 8:
         riRZ(I, n, start, params.timelimit, params.n_ls);
         break;
      case 9:
         raiRZ(I, n, start, params.timelimit, rng, params.n_ls);
         break;
      case 10:
         iSwapFirst(I, n, start, params.timelimit, params.n_ls);
         break;
      case 11:
         iSwapBest(I, n, start, params.timelimit, params.n_ls);
         break;
      case 12:
         viRZ(I, n, start, params.timelimit);
         break;
      case 13:
         powerSwap(I, n, start, params.timelimit);
         break;
      case 14:
         insertJarboui(I, n, start, params.timelimit);
         break;
   }     
}

inline void local_search_procedure_tt(Instance& I, SGEDecider& d, const unsigned d_index, Params& params, PFSSolution& n, std::mt19937& rng, builder& construct, const std::chrono::system_clock::time_point& start) {
   unsigned v;
   std::function<unsigned(unsigned,unsigned)> tiebreak_func;
   if (d_index == 0) {
      v = d.dv_ls_procedure1;
      tiebreak_func = d.tiebreak_ls1;
   }
   else {
      v = d.dv_ls_procedure2;
      tiebreak_func = d.tiebreak_ls2;
   }

   switch (v) {
      case 0: 
         insertion(I, n, construct, tiebreak_func, params);
         break;
      case 1:
         insertionIT1(I, n, construct, params);
         break;
      case 2:
         insertionIT2(I, n, construct, params);
         break;
      case 3:
         ASLS_TT(I, n, params);
         break;
      case 4:
         BISLS_TT(I, n, params);
         break;
      case 5:
         CH3(I, n, construct, tiebreak_func, params);
         break;
      case 6:
         CH3_IT1(I, n, construct, params);
         break;
      case 7:
         CH3_IT2(I, n, construct, params);
         break;
      case 8:
         CH4(I, n, construct, tiebreak_func, params);
         break;
      case 9:
         CH4_IT1(I, n, construct, params);
         break;
      case 10:
         CH4_IT2(I, n, construct, params);
         break;
      case 11:
         CH5(I, n, construct, tiebreak_func, params);
         break;
      case 12:
         CH5_IT1(I, n, construct, params);
         break;
      case 13:
         CH5_IT2(I, n, construct, params);
         break;
      case 14:
         CH6(I, n, construct, tiebreak_func, params);
         break;
      case 15:
         CH6_IT1(I, n, construct, params);
         break;
      case 16:
         CH6_IT2(I, n, construct, params);
         break;
   }
}

inline void local_search(Instance& I, SGEDecider& d, Params& params, PFSSolution& n, std::mt19937& rng, builder& construct, unsigned steps, const std::chrono::system_clock::time_point& start) {
   if (d.dv_local_search == 0)
      return;
   else if (d.dv_local_search == 1 || steps % params.alternate != 0)
      if (params.obj == 'm')
         local_search_procedure_cmax(I, d, 0, params, n, rng, construct, start);
      else if (params.obj == 'f')
         local_search_procedure_csum(I, d, 0, params, n, rng, construct, start);
      else
         local_search_procedure_tt(I, d, 0, params, n, rng, construct, start);
   else
      if (params.obj == 'm')
         local_search_procedure_cmax(I, d, 1, params, n, rng, construct, start);
      else if (params.obj == 'f')
         local_search_procedure_csum(I, d, 1, params, n, rng, construct, start);
      else
         local_search_procedure_tt(I, d, 1, params, n, rng, construct, start);
}

inline bool acceptance(PFSSolution& ns, PFSSolution& cs, const unsigned bs, SGEDecider& d, Params& params, std::vector<unsigned>& past, const unsigned pos) {
   if (ns.ms < cs.ms) {
      cs = ns;
      return true;
   }
   else{
      switch(d.dv_acceptance){
         case 0:
            if (drand48() < exp(-double(ns.ms-cs.ms)/params.T)) {
               cs = ns;
               return true;
            }
            break;
         case 1:
            if (ns.ms < past[pos] || cs.ms == ns.ms) {
               cs = ns;
               return true;
            }
            break;
         case 2:
            if (ns.ms < bs + params.rrtd) {
               cs = ns;
               return true;
            }
            break;
         case 3:
            if (ns.ms < bs + (bs*params.thres)) {
               cs = ns;
               return true;
            }
            break;
      }

      return false;
   }
}

inline void iga_step_np(Instance& I, SGEDecider& d, Params& params, NPFSSolution& n, std::mt19937& rng, nBuilder& construct, const std::chrono::system_clock::time_point& start) {
   switch (d.dv_perturb_np) {
      case 0:
         gi_np(I, n, rng, construct, d.tiebreak_perturb_np, params);
         break;
      case 1:
         gi_IT1_np(I, n, rng, construct, params);
         break;
      case 2:
         gi_IT2_np(I, n, rng, construct, params);
         break;
   }
}

inline void local_search_procedure_cmax_np(Instance& I, SGEDecider& d, const unsigned d_index, Params& params, NPFSSolution& n, std::mt19937& rng, nBuilder& construct, const std::chrono::system_clock::time_point& start) {
   unsigned v;
   std::function<unsigned(unsigned,unsigned)> tiebreak_func;
   if (d_index == 0) {
      v = d.dv_ls_procedure1_np;
      tiebreak_func = d.tiebreak_ls1_np;
   }
   else {
      v = d.dv_ls_procedure2_np;
      tiebreak_func = d.tiebreak_ls2_np;
   }

   switch (v) {
      case 0:
         insertion_np(I, n, params, tiebreak_func);
         break;
      case 1:
         CPASLS_Cmax_np(I, n, params);
         break;
      case 2:
         ASLS_Cmax_np(I, n, params);
         break;
      case 3:
         Pc_np(I, n, construct, tiebreak_func, rng, params);
         break;
   }   
}

inline void local_search_procedure_csum_np(Instance& I, SGEDecider& d, const unsigned d_index, Params& params, NPFSSolution& n, std::mt19937& rng, nBuilder& construct, const std::chrono::system_clock::time_point& start) {
   unsigned v;
   std::function<unsigned(unsigned,unsigned)> tiebreak_func;
   if (d_index == 0) {
      v = d.dv_ls_procedure1_np;
      tiebreak_func = d.tiebreak_ls1_np;
   }
   else {
      v = d.dv_ls_procedure2_np;
      tiebreak_func = d.tiebreak_ls2_np;
   }

   switch (v) {
      case 0:
         AASLS(I, n, params);
         break;   
      case 1:
         insertion_np(I, n, params, tiebreak_func);
         break;
      case 2:
         ARNASLS(I, n, params);
         break;
      case 3:
         AASLS_G8(I, n, params);
         break;
      case 4:
         AASLS_r(I, n, params);
         break;
   }   
}

inline void local_search_procedure_tt_np(Instance& I, SGEDecider& d, const unsigned d_index, Params& params, NPFSSolution& n, std::mt19937& rng, nBuilder& construct, const std::chrono::system_clock::time_point& start) {
   unsigned v;
   std::function<unsigned(unsigned,unsigned)> tiebreak_func;
   if (d_index == 0) {
      v = d.dv_ls_procedure1_np;
      tiebreak_func = d.tiebreak_ls1_np;
   }
   else {
      v = d.dv_ls_procedure2_np;
      tiebreak_func = d.tiebreak_ls2_np;
   }

   switch (v) {
      case 0:
         AASLS_TT(I, n, params);
         break;
      case 1:
         AASLS_r_TT(I, n, params);
         break;
      case 2:
         AASLS_G8_TT(I, n, params);
         break;
      case 3:
         insertion_np(I, n, params, tiebreak_func);
         break;
      case 4:
         insertionTT_IT1_np(I, n, params);
         break;
      case 5:
         insertionTT_IT2_np(I, n, params);
         break;
   }
}

inline void local_search_np(Instance& I, SGEDecider& d, Params& params, NPFSSolution& n, std::mt19937& rng, nBuilder& construct, unsigned steps, const std::chrono::system_clock::time_point& start) {
   if (d.dv_local_search_np == 0)
      return;
   else if (d.dv_local_search_np == 1 || steps % params.alternate_np != 0)
      if (params.obj == 'm')
         local_search_procedure_cmax_np(I, d, 0, params, n, rng, construct, start);
      else if (params.obj == 'f')
         local_search_procedure_csum_np(I, d, 0, params, n, rng, construct, start);
      else
         local_search_procedure_tt_np(I, d, 0, params, n, rng, construct, start);
   else
      if (params.obj == 'm')
         local_search_procedure_cmax_np(I, d, 1, params, n, rng, construct, start);
      else if (params.obj == 'f')
         local_search_procedure_csum_np(I, d, 1, params, n, rng, construct, start);
      else
         local_search_procedure_tt_np(I, d, 1, params, n, rng, construct, start);
}

inline bool acceptance_np(NPFSSolution& ns, NPFSSolution& cs, const unsigned bs, SGEDecider& d, Params& params, std::vector<unsigned>& past, const unsigned pos) {
   if (ns.ms < cs.ms) {
      cs = ns;
      return true;
   }
   else{
      switch(d.dv_acceptance_np){
         case 0:
            if (drand48() < exp(-double(ns.ms-cs.ms)/params.T_np)) {
               cs = ns;
               return true;
            }
            break;
         case 1:
            if (ns.ms < past[pos] || cs.ms == ns.ms) {
               cs = ns;
               return true;
            }
            break;
         case 2:
            if (ns.ms < bs + params.rrtd_np) {
               cs = ns;
               return true;
            }
            break;
         case 3:
            if (ns.ms < bs + (bs*params.thres_np)) {
               cs = ns;
               return true;
            }
            break;
      }

      return false;
   }
}