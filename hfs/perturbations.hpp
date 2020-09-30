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
#include <algorithm>
#include "instance.hpp"
#include "solution.hpp"
#include "constructive.hpp"
#include "ls.hpp"
#include "rules.hpp"
#include "parameters.hpp"
template<typename Tiebreak>
inline void ri(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Tiebreak tiebreak, Params& p) {
   std::uniform_int_distribution<unsigned> dis(1, I.n);
   std::vector<unsigned> available(n.pi.begin()+1, n.pi.end());
   std::shuffle(available.begin(), available.end(), rng);
   unsigned index = 0;
   for(unsigned i = 0; i < p.dc; i++){
      unsigned job = available[index];
      index++;
      auto pos = std::find(n.pi.begin()+1, n.pi.end(), job);
      assert(pos != n.pi.end());
      unsigned rp1 = pos - n.pi.begin();
      unsigned rp2 = dis(rng);
      while (rp2 == rp1) {
         rp2 = dis(rng);
      }
      if (rp1 < rp2) {
         for (unsigned j = rp1; j < rp2; j++) {
            n.pi[j] = n.pi[j+1];
         }
      }
      else {
         for (unsigned j = rp1; j > rp2; j--) {
            n.pi[j] = n.pi[j-1];
         }
      }
      n.pi[rp2] = job;
   }
   if (p.obj == 'm') {
      n.computeMakespan(I);
   }
   else if (p.obj == 'f') {
      n.computeFlowtime(I);
   }
   else {
      n.computeTotalTardiness(I);
   }
}
template <typename Tiebreak>
inline void gi(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Tiebreak tiebreak, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i=0; i<p.dc; i++, dst++, left--)
      std::swap(*dst,*(dst+(dis(rng)%left)));
   remove.erase(dst,remove.end());
   auto nend = n.pi.end();
   for(unsigned i=1; i<=p.dc; i++)
      nend=std::remove(n.pi.begin()+kl,nend,remove[i]);
   std::fill(nend,n.pi.end(),0);
   if (p.obj == 'm') {
      construct.constructp(I,n,I.n-p.dc+1,kl,ku-p.dc,remove,tiebreak);
   }
   else if (p.obj == 'f') {
      construct.constructpCsum(I,n,I.n-p.dc+1,kl,ku-p.dc,remove,tiebreak);
   }
   else {
      construct.constructpTT(I,n,I.n-p.dc+1,kl,ku-p.dc,remove,tiebreak);
   }
}
template<typename Tiebreak>
inline void rs(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Tiebreak tiebreak, Params& p) {
   std::uniform_int_distribution<unsigned> dis(1, I.n);
   std::vector<unsigned> available(n.pi.begin()+1, n.pi.end());
   std::shuffle(available.begin(), available.end(), rng);
   unsigned index = 0;
   for(unsigned i = 0; i < p.dc; i++){
      unsigned job = available[index];
      index++;
      auto pos = std::find(n.pi.begin()+1, n.pi.end(), job);
      assert(pos != n.pi.end());
      unsigned rp1 = pos - n.pi.begin(); 
      unsigned rp2 = dis(rng); 
      while (rp2 == rp1) {
         rp2 = dis(rng);
      }
      std::swap(n.pi[rp1], n.pi[rp2]);
   }
   if (p.obj == 'm') {
      n.computeMakespan(I);
   }
   else if (p.obj == 'f') {
      n.computeFlowtime(I);
   }
   else {
      n.computeTotalTardiness(I);
   }
}
template<typename Tiebreak>
inline void gs(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Tiebreak tiebreak, Params& p) {
   std::vector<unsigned> available(n.pi.begin()+1, n.pi.end());
   std::shuffle(available.begin(), available.end(), rng);
   unsigned index = 0;
   for(unsigned i = 0; i < p.dc; i++) {
      unsigned rp1 = available[index]; 
      index++;
      unsigned best = std::numeric_limits<unsigned>::max();
      unsigned bestj = 0;
      for (unsigned j = 1; j <= I.n; j++) {
         if (j != rp1) {
            std::swap(n.pi[rp1], n.pi[j]);
            if (p.obj == 'm') {
               n.computeMakespan(I);
            }
            else if (p.obj == 'f') {
               n.computeFlowtime(I);
            }
            else {
               n.computeTotalTardiness(I);
            }
            if (n.ms < best) {
               best = n.ms;
               bestj = j;
            }
            std::swap(n.pi[rp1], n.pi[j]);
         }
      }
      assert(bestj != 0);
      std::swap(n.pi[rp1], n.pi[bestj]);
   }
   if (p.obj == 'm') {
      n.computeMakespan(I);
   }
   else if (p.obj == 'f') {
      n.computeFlowtime(I);
   }
   else {
      n.computeTotalTardiness(I);
   }
}
template<typename Tiebreak>
inline void ras(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Tiebreak tiebreak, Params& p) {
   std::uniform_int_distribution<unsigned> dis(1, I.n - 1);
   for(unsigned i = 0; i < p.dc; i++){
      unsigned rp1 = dis(rng); 
      std::swap(n.pi[rp1], n.pi[rp1+1]);
   }
   if (p.obj == 'm') {
      n.computeMakespan(I);
   }
   else if (p.obj == 'f') {
      n.computeFlowtime(I);
   }
   else {
      n.computeTotalTardiness(I);
   }
}
template <typename Tiebreak>
inline void gi_asls(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Tiebreak tiebreak, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i = 0; i < p.dc; i++, dst++, left--)
      std::swap(*dst, *(dst+(dis(rng)%left)));
   remove.erase(dst, remove.end());
   auto nend = n.pi.end();
   for (unsigned i = 1; i <= p.dc; i++)
      nend = std::remove(n.pi.begin()+kl, nend, remove[i]);
   std::fill(nend, n.pi.end(), 0);
   for (unsigned i = 1; i <= p.dc; i++) {
      std::vector<Job> insert = { 0, remove[i] };
      if (p.obj == 'm') {
         construct.constructp(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert,tiebreak);
      }
      else if (p.obj == 'f') {
         construct.constructpCsum(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert,tiebreak);
      }
      else {
         construct.constructpTT(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert,tiebreak);
      }
      auto b = std::find(n.pi.begin()+1, n.pi.end(), remove[i]); 
      assert(b != n.pi.end());
      unsigned best = n.ms;
      unsigned bestj = 0;
      unsigned index = (b - n.pi.begin()) + 1; 
      for (unsigned j = I.n-p.dc+i; j > index; j--) {
         std::swap(n.pi[j], n.pi[j-1]);
         if (p.obj == 'm') {
            n.computeMakespan(I,n.pi.size()-p.dc+i);
         }
         else if (p.obj == 'f') {
            n.computeFlowtime(I,n.pi.size()-p.dc+i);
         }
         else {
            n.computeTotalTardiness(I,n.pi.size()-p.dc+i);
         }
         if (n.ms < best) {
            best = n.ms;
            bestj = j;
         }
         std::swap(n.pi[j], n.pi[j-1]);
      }
      if (bestj > 0) { 
         std::swap(n.pi[bestj], n.pi[bestj-1]);
      }
   }
   if (p.obj == 'm') {
      n.computeMakespan(I);
   }
   else if (p.obj == 'f') {
      n.computeFlowtime(I);
   }
   else {
      n.computeTotalTardiness(I);
   }
}
template <typename Tiebreak>
inline void ils_gi(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Tiebreak tiebreak, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i=0; i<p.dc; i++, dst++, left--)
      std::swap(*dst,*(dst+(dis(rng)%left)));
   remove.erase(dst,remove.end());
   auto nend = n.pi.end();
   for(unsigned i=1; i<=p.dc; i++)
      nend=std::remove(n.pi.begin()+kl,nend,remove[i]);
   std::fill(nend,n.pi.end(),0);
   if (p.obj == 'm') {
      n.computeMakespan(I,n.pi.size()-p.dc);
   }
   else if (p.obj == 'f') {
      n.computeFlowtime(I,n.pi.size()-p.dc);   
   }
   else {
      n.computeTotalTardiness(I,n.pi.size()-p.dc);
   }
   lsps_i(I,n,construct,p.obj,kl,ku-p.dc,p.dc,p.n_ls,p.start,p.timelimit,tiebreak);
   if (p.obj == 'm') {
      construct.constructp(I,n,I.n-p.dc+1,kl,ku-p.dc,remove,tiebreak);
   }
   else if (p.obj == 'f') {
      construct.constructpCsum(I,n,I.n-p.dc+1,kl,ku-p.dc,remove,tiebreak);
   }
   else {
      construct.constructpTT(I,n,I.n-p.dc+1,kl,ku-p.dc,remove,tiebreak);
   }
}
template <typename Tiebreak>
inline void gi_ils(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Tiebreak tiebreak, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i = 0; i < p.dc; i++, dst++, left--)
      std::swap(*dst, *(dst+(dis(rng)%left)));
   remove.erase(dst, remove.end());
   auto nend = n.pi.end();
   for (unsigned i = 1; i <= p.dc; i++)
      nend = std::remove(n.pi.begin()+kl, nend, remove[i]);
   std::fill(nend, n.pi.end(), 0);
   for (unsigned i = 1; i <= p.dc; i++) {
      std::vector<Job> insert = { 0, remove[i] };
      if (p.obj == 'm') {
         construct.constructp(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert,tiebreak);
      }
      else if (p.obj == 'f') {
         construct.constructpCsum(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert,tiebreak);
      }
      else {
         construct.constructpTT(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert,tiebreak);
      }
      auto b = std::find(n.pi.begin()+1, n.pi.end(), remove[i]);
      assert(b != n.pi.end());
      insert.erase(insert.begin()+1);
      nend = n.pi.end()-p.dc+i;
      unsigned index = b - n.pi.begin();
      if (index > 1) {
         insert.push_back(n.pi[index-1]);
      }
      if (index < I.n-p.dc+i) {
         insert.push_back(n.pi[index+1]);
      }
      for (unsigned j = 1; j < insert.size(); j++) {
         nend = std::remove(n.pi.begin()+kl, nend, insert[j]);
      }
      std::fill(nend, n.pi.end(), 0);
      if (p.obj == 'm') {
         construct.constructp(I,n,I.n-p.dc+i-(insert.size()-1)+1,kl,ku-p.dc+i-(insert.size()-1),insert,tiebreak);
      }
      else if (p.obj == 'f') {
         construct.constructpCsum(I,n,I.n-p.dc+i-(insert.size()-1)+1,kl,ku-p.dc+i-(insert.size()-1),insert,tiebreak);
      }
      else {
         construct.constructpTT(I,n,I.n-p.dc+i-(insert.size()-1)+1,kl,ku-p.dc+i-(insert.size()-1),insert,tiebreak);
      }
   }
}
inline void giFF(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i=0; i<p.dc; i++, dst++, left--)
      std::swap(*dst,*(dst+(dis(rng)%left)));
   remove.erase(dst,remove.end());
   auto nend = n.pi.end();
   for(unsigned i=1; i<=p.dc; i++)
      nend=std::remove(n.pi.begin()+kl,nend,remove[i]);
   std::fill(nend,n.pi.end(),0);
   construct.constructpFF(I,n,I.n-p.dc+1,kl,ku-p.dc,remove);
}
inline void gi_aslsFF(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i = 0; i < p.dc; i++, dst++, left--)
      std::swap(*dst, *(dst+(dis(rng)%left)));
   remove.erase(dst, remove.end());
   auto nend = n.pi.end();
   for (unsigned i = 1; i <= p.dc; i++)
      nend = std::remove(n.pi.begin()+kl, nend, remove[i]);
   std::fill(nend, n.pi.end(), 0);
   for (unsigned i = 1; i <= p.dc; i++) {
      std::vector<Job> insert = { 0, remove[i] };
      construct.constructpFF(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert);
      auto b = std::find(n.pi.begin()+1, n.pi.end(), remove[i]); 
      assert(b != n.pi.end());
      unsigned best = n.ms;
      unsigned bestj = 0;
      unsigned index = (b - n.pi.begin()) + 1; 
      for (unsigned j = I.n-p.dc+i; j > index; j--) {
         std::swap(n.pi[j], n.pi[j-1]);
         n.computeMakespan(I,n.pi.size()-p.dc+i);
         if (n.ms < best) {
            best = n.ms;
            bestj = j;
         }
         std::swap(n.pi[j], n.pi[j-1]);
      }
      if (bestj > 0) { 
         std::swap(n.pi[bestj], n.pi[bestj-1]);
      }
   }
   n.computeMakespan(I);
}
inline void ils_giFF(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i=0; i<p.dc; i++, dst++, left--)
      std::swap(*dst,*(dst+(dis(rng)%left)));
   remove.erase(dst,remove.end());
   auto nend = n.pi.end();
   for(unsigned i=1; i<=p.dc; i++)
      nend=std::remove(n.pi.begin()+kl,nend,remove[i]);
   std::fill(nend,n.pi.end(),0);
   n.computeMakespan(I,n.pi.size()-p.dc);
   lsps_iFF(I,n,construct,p.obj,kl,ku-p.dc,p.dc,p.n_ls,p.start,p.timelimit);
   construct.constructpFF(I,n,I.n-p.dc+1,kl,ku-p.dc,remove);
}
inline void gi_ilsFF(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i = 0; i < p.dc; i++, dst++, left--)
      std::swap(*dst, *(dst+(dis(rng)%left)));
   remove.erase(dst, remove.end());
   auto nend = n.pi.end();
   for (unsigned i = 1; i <= p.dc; i++)
      nend = std::remove(n.pi.begin()+kl, nend, remove[i]);
   std::fill(nend, n.pi.end(), 0);
   for (unsigned i = 1; i <= p.dc; i++) {
      std::vector<Job> insert = { 0, remove[i] };
      construct.constructpFF(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert);
      auto b = std::find(n.pi.begin()+1, n.pi.end(), remove[i]);
      assert(b != n.pi.end());
      insert.erase(insert.begin()+1);
      nend = n.pi.end()-p.dc+i;
      unsigned index = b - n.pi.begin();
      if (index > 1) {
         insert.push_back(n.pi[index-1]);
      }
      if (index < I.n-p.dc+i) {
         insert.push_back(n.pi[index+1]);
      }
      for (unsigned j = 1; j < insert.size(); j++) {
         nend = std::remove(n.pi.begin()+kl, nend, insert[j]);
      }
      std::fill(nend, n.pi.end(), 0);
      construct.constructpFF(I,n,I.n-p.dc+i-(insert.size()-1)+1,kl,ku-p.dc+i-(insert.size()-1),insert);
   }
}
inline void giIT1(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i=0; i<p.dc; i++, dst++, left--)
      std::swap(*dst,*(dst+(dis(rng)%left)));
   remove.erase(dst,remove.end());
   auto nend = n.pi.end();
   for(unsigned i=1; i<=p.dc; i++)
      nend=std::remove(n.pi.begin()+kl,nend,remove[i]);
   std::fill(nend,n.pi.end(),0);
   construct.constructpTT_IT1(I,n,I.n-p.dc+1,kl,ku-p.dc,remove);
}
inline void giIT2(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i=0; i<p.dc; i++, dst++, left--)
      std::swap(*dst,*(dst+(dis(rng)%left)));
   remove.erase(dst,remove.end());
   auto nend = n.pi.end();
   for(unsigned i=1; i<=p.dc; i++)
      nend=std::remove(n.pi.begin()+kl,nend,remove[i]);
   std::fill(nend,n.pi.end(),0);
   construct.constructpTT_IT2(I,n,I.n-p.dc+1,kl,ku-p.dc,remove);
}
inline void gi_aslsIT1(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i = 0; i < p.dc; i++, dst++, left--)
      std::swap(*dst, *(dst+(dis(rng)%left)));
   remove.erase(dst, remove.end());
   auto nend = n.pi.end();
   for (unsigned i = 1; i <= p.dc; i++)
      nend = std::remove(n.pi.begin()+kl, nend, remove[i]);
   std::fill(nend, n.pi.end(), 0);
   for (unsigned i = 1; i <= p.dc; i++) {
      std::vector<Job> insert = { 0, remove[i] };
      construct.constructpTT_IT1(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert);
      auto b = std::find(n.pi.begin()+1, n.pi.end(), remove[i]); 
      assert(b != n.pi.end());
      unsigned best = n.ms;
      unsigned bestj = 0;
      unsigned index = (b - n.pi.begin()) + 1; 
      for (unsigned j = I.n-p.dc+i; j > index; j--) {
         std::swap(n.pi[j], n.pi[j-1]);
         n.computeTotalTardiness(I,n.pi.size()-p.dc+i);
         if (n.ms < best) {
            best = n.ms;
            bestj = j;
         }
         std::swap(n.pi[j], n.pi[j-1]);
      }
      if (bestj > 0) { 
         std::swap(n.pi[bestj], n.pi[bestj-1]);
      }
   }
   n.computeTotalTardiness(I);
}
inline void gi_aslsIT2(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i = 0; i < p.dc; i++, dst++, left--)
      std::swap(*dst, *(dst+(dis(rng)%left)));
   remove.erase(dst, remove.end());
   auto nend = n.pi.end();
   for (unsigned i = 1; i <= p.dc; i++)
      nend = std::remove(n.pi.begin()+kl, nend, remove[i]);
   std::fill(nend, n.pi.end(), 0);
   for (unsigned i = 1; i <= p.dc; i++) {
      std::vector<Job> insert = { 0, remove[i] };
      construct.constructpTT_IT2(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert);
      auto b = std::find(n.pi.begin()+1, n.pi.end(), remove[i]); 
      assert(b != n.pi.end());
      unsigned best = n.ms;
      unsigned bestj = 0;
      unsigned index = (b - n.pi.begin()) + 1; 
      for (unsigned j = I.n-p.dc+i; j > index; j--) {
         std::swap(n.pi[j], n.pi[j-1]);
         n.computeTotalTardiness(I,n.pi.size()-p.dc+i);
         if (n.ms < best) {
            best = n.ms;
            bestj = j;
         }
         std::swap(n.pi[j], n.pi[j-1]);
      }
      if (bestj > 0) { 
         std::swap(n.pi[bestj], n.pi[bestj-1]);
      }
   }
   n.computeTotalTardiness(I);
}
inline void ils_giIT1(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i=0; i<p.dc; i++, dst++, left--)
      std::swap(*dst,*(dst+(dis(rng)%left)));
   remove.erase(dst,remove.end());
   auto nend = n.pi.end();
   for(unsigned i=1; i<=p.dc; i++)
      nend=std::remove(n.pi.begin()+kl,nend,remove[i]);
   std::fill(nend,n.pi.end(),0);
   n.computeTotalTardiness(I,n.pi.size()-p.dc);
   lsps_iIT1(I,n,construct,p.obj,kl,ku-p.dc,p.dc,p.n_ls,p.start,p.timelimit);
   construct.constructpTT_IT1(I,n,I.n-p.dc+1,kl,ku-p.dc,remove);
}
inline void ils_giIT2(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i=0; i<p.dc; i++, dst++, left--)
      std::swap(*dst,*(dst+(dis(rng)%left)));
   remove.erase(dst,remove.end());
   auto nend = n.pi.end();
   for(unsigned i=1; i<=p.dc; i++)
      nend=std::remove(n.pi.begin()+kl,nend,remove[i]);
   std::fill(nend,n.pi.end(),0);
   n.computeTotalTardiness(I,n.pi.size()-p.dc);
   lsps_iIT2(I,n,construct,p.obj,kl,ku-p.dc,p.dc,p.n_ls,p.start,p.timelimit);
   construct.constructpTT_IT2(I,n,I.n-p.dc+1,kl,ku-p.dc,remove);
}
inline void gi_ilsIT1(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i = 0; i < p.dc; i++, dst++, left--)
      std::swap(*dst, *(dst+(dis(rng)%left)));
   remove.erase(dst, remove.end());
   auto nend = n.pi.end();
   for (unsigned i = 1; i <= p.dc; i++)
      nend = std::remove(n.pi.begin()+kl, nend, remove[i]);
   std::fill(nend, n.pi.end(), 0);
   for (unsigned i = 1; i <= p.dc; i++) {
      std::vector<Job> insert = { 0, remove[i] };
      construct.constructpTT_IT1(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert);
      auto b = std::find(n.pi.begin()+1, n.pi.end(), remove[i]);
      assert(b != n.pi.end());
      insert.erase(insert.begin()+1);
      nend = n.pi.end()-p.dc+i;
      unsigned index = b - n.pi.begin();
      if (index > 1) {
         insert.push_back(n.pi[index-1]);
      }
      if (index < I.n-p.dc+i) {
         insert.push_back(n.pi[index+1]);
      }
      for (unsigned j = 1; j < insert.size(); j++) {
         nend = std::remove(n.pi.begin()+kl, nend, insert[j]);
      }
      std::fill(nend, n.pi.end(), 0);
      construct.constructpTT_IT1(I,n,I.n-p.dc+i-(insert.size()-1)+1,kl,ku-p.dc+i-(insert.size()-1),insert);
   }
}
inline void gi_ilsIT2(Instance& I, PFSSolution& n, std::mt19937& rng, builder& construct, Params& p) {
   unsigned kl = 1;
   unsigned ku = I.n+1;
   std::uniform_int_distribution<> dis(0);
   std::vector<Job> remove(n.pi.begin()+kl-1,n.pi.begin()+ku);
   unsigned left = ku-kl;
   auto dst = remove.begin()+1;
   for(unsigned i = 0; i < p.dc; i++, dst++, left--)
      std::swap(*dst, *(dst+(dis(rng)%left)));
   remove.erase(dst, remove.end());
   auto nend = n.pi.end();
   for (unsigned i = 1; i <= p.dc; i++)
      nend = std::remove(n.pi.begin()+kl, nend, remove[i]);
   std::fill(nend, n.pi.end(), 0);
   for (unsigned i = 1; i <= p.dc; i++) {
      std::vector<Job> insert = { 0, remove[i] };
      construct.constructpTT_IT2(I,n,I.n-p.dc+i,kl,ku-p.dc+i-1,insert);
      auto b = std::find(n.pi.begin()+1, n.pi.end(), remove[i]);
      assert(b != n.pi.end());
      insert.erase(insert.begin()+1);
      nend = n.pi.end()-p.dc+i;
      unsigned index = b - n.pi.begin();
      if (index > 1) {
         insert.push_back(n.pi[index-1]);
      }
      if (index < I.n-p.dc+i) {
         insert.push_back(n.pi[index+1]);
      }
      for (unsigned j = 1; j < insert.size(); j++) {
         nend = std::remove(n.pi.begin()+kl, nend, insert[j]);
      }
      std::fill(nend, n.pi.end(), 0);
      construct.constructpTT_IT2(I,n,I.n-p.dc+i-(insert.size()-1)+1,kl,ku-p.dc+i-(insert.size()-1),insert);
   }
}
template <typename Tiebreak>
inline void gi_np(Instance& I, NPFSSolution& n, std::mt19937& rng, nBuilder& construct, Tiebreak tiebreak, Params& params) {
   unsigned dc = params.dc_np;
   std::vector<Job> remove(I.n+1, 0);
   std::iota(remove.begin(), remove.end(), 0);
   std::shuffle(remove.begin()+1, remove.end(), rng);
   remove.erase(remove.begin()+1+dc, remove.end());
   for(unsigned r = 1; r < remove.size(); r++) {
      for(unsigned i = 1; i <= I.m; i++) {
         unsigned dj = 1;
         for(unsigned j = 1; j <= I.n; j++) {
            if (n.pi[i][j] != remove[r]) {
               n.pi[i][dj] = n.pi[i][j];
               dj++;
            }
         }
         while (dj <= I.n) {
            n.pi[i][dj] = 0;
            dj++;
         }
      }
   }
   if (params.obj == 'm')
      construct.constructnpNEW(n, I.n-dc+1, remove, tiebreak);
   else if (params.obj == 'f')
      construct.constructnpNEWCsum(n, I.n-dc+1, remove, tiebreak);
   else
      construct.constructnpNEWTT(n, I.n-dc+1, remove, tiebreak);
}
inline void gi_IT1_np(Instance& I, NPFSSolution& n, std::mt19937& rng, nBuilder& construct, Params& params) {
   unsigned dc = params.dc_np;
   std::vector<Job> remove(I.n+1, 0);
   std::iota(remove.begin(), remove.end(), 0);
   std::shuffle(remove.begin()+1, remove.end(), rng);
   remove.erase(remove.begin()+1+dc, remove.end());
   for(unsigned r = 1; r < remove.size(); r++) {
      for(unsigned i = 1; i <= I.m; i++) {
         unsigned dj = 1;
         for(unsigned j = 1; j <= I.n; j++) {
            if (n.pi[i][j] != remove[r]) {
               n.pi[i][dj] = n.pi[i][j];
               dj++;
            }
         }
         while (dj <= I.n) {
            n.pi[i][dj] = 0;
            dj++;
         }
      }
   }
   construct.constructnpNEWTT_IT1(n, I.n-dc+1, remove);
}
inline void gi_IT2_np(Instance& I, NPFSSolution& n, std::mt19937& rng, nBuilder& construct, Params& params) {
   unsigned dc = params.dc_np;
   std::vector<Job> remove(I.n+1, 0);
   std::iota(remove.begin(), remove.end(), 0);
   std::shuffle(remove.begin()+1, remove.end(), rng);
   remove.erase(remove.begin()+1+dc, remove.end());
   for(unsigned r = 1; r < remove.size(); r++) {
      for(unsigned i = 1; i <= I.m; i++) {
         unsigned dj = 1;
         for(unsigned j = 1; j <= I.n; j++) {
            if (n.pi[i][j] != remove[r]) {
               n.pi[i][dj] = n.pi[i][j];
               dj++;
            }
         }
         while (dj <= I.n) {
            n.pi[i][dj] = 0;
            dj++;
         }
      }
   }
   construct.constructnpNEWTT_IT2(n, I.n-dc+1, remove);
}
