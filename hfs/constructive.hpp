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
#include <utility>
#include <set>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "instance.hpp"
#include "solution.hpp"
#include "rules.hpp"
#include "units.hpp"
typedef boost::multi_array_types::index_range irange;
const std::tuple<unsigned,unsigned,unsigned> nomove { 0,0,0 };
struct builder {
  unsigned n,m;
  const Instance& I;
  boost::multi_array<Time,2> h, t, f;  
  boost::multi_array<Time,2> p;        
  boost::multi_array<Time,2> p_;       
  boost::multi_array<unsigned,2> mup;  
  std::vector<Time> Cm, Cpsum;
  Time Csum;
  builder(const Instance& _I) : n(_I.n), m(_I.m), I(_I),
        h(boost::extents[m+1][n+1]), t(boost::extents[m+2][n+2]), f(boost::extents[m+1][n+1]),
        p(boost::extents[m+1][n+1]), p_(boost::extents[n+1][m+1]),
        mup(boost::extents[m+2][n+2]), Cm(m+1), Cpsum(n+1) {
    for(unsigned i=0; i<=m; i++)
      for(unsigned j=0; j<=n; j++) {
        p[i][j] = _I.p[j][i];
        p_[j][i] = _I.p[j][i];
      }
  }
  template <typename Tiebreak>
  inline void constructp(const Instance& I, PFSSolution& r, unsigned k, unsigned kl, unsigned ku, const std::vector<Job>& pi, Tiebreak tiebreak, Time ub = infinite_time, unsigned minPos = 1, unsigned maxPos = 0) {
    const unsigned n = I.n;
    const unsigned m = I.m;
    Time FC = infinite_time;
    if (maxPos==0)
      maxPos=k-1;
    for(unsigned i=m; i>=1; i--)
      for(unsigned j=maxPos; j>=1; j--)
        t[i][j+n+1-k] = std::max(t[i][j+1+n+1-k],t[i+1][j+n+1-k])+I.p[r.pi[j]][i];
    for(unsigned i=1; i<=m; i++)
      for(unsigned j=minPos; j<k; j++)
        h[i][j] = std::max(h[i-1][j],h[i][j-1])+I.p[r.pi[j]][i];
    for(unsigned kk=1; kk<pi.size(); kk++, k++, ku++) { 
      const unsigned job = pi[kk];
      unsigned bp = 0, bt = std::numeric_limits<unsigned>::max(); 
      FC = infinite_time;
      for(unsigned j=kl; j<=ku; j++) { 
        Time Fk = 0, ct = 0;
        for(unsigned i=1; i<=m; i++) {
          ct = std::max(ct,h[i][j-1])+I.p[job][i];
          Fk = std::max(Fk,ct+t[i][j+n+1-k]);
          if (Fk>FC)
            break;
        }
        unsigned tb = tiebreak(job,j);
        if (Fk < FC || (Fk == FC && tb < bt)) {
          FC = Fk;
          bp = j;
          bt = tb;
        }
      }
      assert(bp > 0);
      std::copy_backward(r.pi.begin()+bp,r.pi.begin()+k,r.pi.begin()+k+1);
      r.pi[bp]=job;
      if (FC>=ub) {
        r.n = k;
        return; 
      }
      if (kk+1<pi.size()) {
        for(unsigned i=1; i<=m; i++)
          for(unsigned j=bp; j<k+1; j++)
            h[i][j] = std::max(h[i-1][j],h[i][j-1])+I.p[r.pi[j]][i];
        for(unsigned i=m; i>=1; i--)
          for(unsigned j=bp; j>=1; j--)
            t[i][j+n-k] = std::max(t[i][j+1+n-k],t[i+1][j+n-k])+I.p[r.pi[j]][i];
      }
    } 
    r.ms = FC;
  }
  inline void constructpFF(const Instance& I, PFSSolution& r, unsigned k, unsigned kl, unsigned ku, const std::vector<Job>& pi, Time ub = infinite_time, unsigned minPos = 1, unsigned maxPos = 0) {
    const unsigned n = I.n;
    const unsigned m = I.m;
    std::vector<unsigned> ptb(I.n); 
    Time FC = infinite_time;
    if (maxPos==0)
      maxPos=k-1;
    for(unsigned i=m; i>=1; i--)
      for(unsigned j=maxPos; j>=1; j--)
        t[i][j+n+1-k] = std::max(t[i][j+1+n+1-k],t[i+1][j+n+1-k])+I.p[r.pi[j]][i];
    for(unsigned i=1; i<=m; i++)
      for(unsigned j=minPos; j<k; j++)
        h[i][j] = std::max(h[i-1][j],h[i][j-1])+I.p[r.pi[j]][i];
    for(unsigned kk=1; kk<pi.size(); kk++, k++, ku++) { 
      const unsigned job = pi[kk];
      unsigned bp = 0, tb = 0; 
      int it = std::numeric_limits<int>::max(); 
      std::fill(ptb.begin(), ptb.end(), 0);
      FC = infinite_time;
      for(unsigned j=1; j<=k; j++) { 
        Time Fk = 0, ct = 0;
        for(unsigned i=1; i<=m; i++) {
          ct = std::max(ct, h[i][j-1]) + I.p[job][i];
          f[i][j] = ct; 
          Fk = std::max(Fk, ct + t[i][j+n+1-k]);
          if (Fk > FC)
            break;
        }
        if (Fk < FC) {
          FC = Fk;
          bp = j;
          tb = 1;
          ptb[tb-1] = bp;
        } else if (Fk == FC) {
          tb++;
          ptb[tb-1] = j;
        }
      }
      assert(bp > 0 && bp <= k && tb > 0);
      if (tb > 1) {
        for (unsigned l=0; l<tb; l++) {
          int F = 0;
          if (ptb[l] == k) {
            for (unsigned i=2; i<=m; i++)
              F += f[i][k] - h[i][k-1];
          } else {
            unsigned ff = f[1][ptb[l]] + I.p[r.pi[ptb[l]]][1];
            for (unsigned i=2; i<=m; i++) {
              int diff = ff - f[i][ptb[l]];
              F += f[i][ptb[l]] - h[i][ptb[l]] + I.p[r.pi[ptb[l]]][i] + std::max(0, diff);
              ff = std::max(ff, f[i][ptb[l]]) + I.p[r.pi[ptb[l]]][i];
            }
          }
          if (F < it) {
            bp = ptb[l];
            it = F;
          }
        }
      }
      std::copy_backward(r.pi.begin()+bp,r.pi.begin()+k,r.pi.begin()+k+1);
      r.pi[bp]=job;
      if (FC>=ub) {
        r.n = k;
        return; 
      }
      if (kk+1<pi.size()) {
        for(unsigned i=1; i<=m; i++)
          for(unsigned j=bp; j<k+1; j++)
            h[i][j] = std::max(h[i-1][j],h[i][j-1])+I.p[r.pi[j]][i];
        for(unsigned i=m; i>=1; i--)
          for(unsigned j=bp; j>=1; j--)
            t[i][j+n-k] = std::max(t[i][j+1+n-k],t[i+1][j+n-k])+I.p[r.pi[j]][i];
      }
    } 
    r.ms = FC;
  }
  template <typename Tiebreak>
  inline void constructpMaxMakes(const Instance& I, PFSSolution& r, unsigned k, unsigned kl, unsigned ku, const std::vector<Job>& pi, Tiebreak tiebreak, Time ub = infinite_time, unsigned minPos = 1, unsigned maxPos = 0) {
    const unsigned n = I.n;
    const unsigned m = I.m;
    Time FC = infinite_time;
    if (maxPos==0)
      maxPos=k-1;
    for(unsigned i=m; i>=1; i--)
      for(unsigned j=maxPos; j>=1; j--)
        t[i][j+n+1-k] = std::max(t[i][j+1+n+1-k],t[i+1][j+n+1-k])+I.p[r.pi[j]][i];
    for(unsigned i=1; i<=m; i++)
      for(unsigned j=minPos; j<k; j++)
        h[i][j] = std::max(h[i-1][j],h[i][j-1])+I.p[r.pi[j]][i];
    for(unsigned kk=1; kk<pi.size(); kk++, k++, ku++) { 
      const unsigned job = pi[kk];
      unsigned bp = 0, bt = std::numeric_limits<unsigned>::max(); 
      FC = infinite_time;
      for(unsigned j=kl; j<=ku; j++) { 
        Time Fk = 0, ct = 0;
        for(unsigned i=1; i<=m; i++) {
          ct = std::max(ct,h[i][j-1])+I.p[job][i];
          Fk = std::max(Fk,ct+t[i][j+n+1-k]);
          if (Fk>FC)
            break;
        }
        unsigned tb = tiebreak(job,j);
        if (Fk > FC || (Fk == FC && tb < bt)) {
          FC = Fk;
          bp = j;
          bt = tb;
        }
      }
      assert(bp > 0);
      std::copy_backward(r.pi.begin()+bp,r.pi.begin()+k,r.pi.begin()+k+1);
      r.pi[bp]=job;
      if (FC>=ub) {
        r.n = k;
        return; 
      }
      if (kk+1<pi.size()) {
        for(unsigned i=1; i<=m; i++)
          for(unsigned j=bp; j<k+1; j++)
            h[i][j] = std::max(h[i-1][j],h[i][j-1])+I.p[r.pi[j]][i];
        for(unsigned i=m; i>=1; i--)
          for(unsigned j=bp; j>=1; j--)
            t[i][j+n-k] = std::max(t[i][j+1+n-k],t[i+1][j+n-k])+I.p[r.pi[j]][i];
      }
    } 
    r.ms = FC;
  }
  template <typename Tiebreak>
  inline void constructpCsum(const Instance& I, PFSSolution& r, unsigned k, unsigned kl, unsigned ku, const std::vector<Job>& pi, Tiebreak tiebreak, Time ub = infinite_time) {
    Time FC = infinite_time;
    for(unsigned i=m; i>=1; i--)
      for(unsigned j=k-1; j>=1; j--) {
        if (t[i][j+1+n+1-k]<=t[i+1][j+n+1-k]) {
          mup[i][j+n+1-k] = mup[i+1][j+n+1-k];
        } else {
          mup[i][j+n+1-k] = i;
        }
        t[i][j+n+1-k] = std::max(t[i][j+1+n+1-k], t[i+1][j+n+1-k])+p[i][r.pi[j]];
      }
    for(unsigned i=1; i<=m; i++)
      for(unsigned j=1; j<k; j++)
        h[i][j] = std::max(h[i-1][j],h[i][j-1])+p[i][r.pi[j]];
    for(unsigned j=1; j<k; j++)
      Cpsum[j]=Cpsum[j-1]+h[m][j];
    for(unsigned kk=1; kk<pi.size(); kk++, k++, ku++) { 
      const unsigned job = pi[kk];
      unsigned bp = 0, bt = std::numeric_limits<unsigned>::max(); 
      FC = infinite_time;
      for(unsigned j=kl; j<=ku; j++) { 
        Time Fk = 0;
        unsigned im = 0;
        Time pt = 0;
        auto cp = p_[job].begin()+1;
        for(unsigned i=1; i<=m; i++, cp++) {
          pt = Cm[i] = std::max(pt,h[i][j-1])+*cp;
          Time Fi = pt+t[i][j+n+1-k];
          if (Fi >= Fk) {
            Fk = Fi;
            im = i;
          }
        }
        assert(im>0);
        Csum = Cpsum[j-1] + Cm[m];
        auto pcj = r.pi.begin()+j;
        if (j<k) {
          cp = p_[*pcj].begin()+im;
          Cm[im] += *cp++;
          pt = Cm[im];
          for(unsigned i=m-im, ii=im+1; i>0; i--,ii++,cp++)
            pt = Cm[ii] = std::max(pt,Cm[ii])+*cp;
          Csum  += pt;
        }
        pcj++;
        for(unsigned jj=j+1+n+1-k-1, je=k+n+1-k-1; jj<je; jj++, pcj++) {
          im = mup[im][jj];
          assert(im <= m);
          cp = p_[*pcj].begin()+im;
          Cm[im] += *cp++;
          pt = Cm[im];
          for(unsigned i=m-im, ii=im+1; i>0; i--,ii++,cp++)
            pt = Cm[ii] = std::max(pt,Cm[ii])+*cp;
          Csum  += pt;
        }
        {
          unsigned tb = tiebreak(job,j);
          if (Csum < FC || (Csum == FC && tb < bt)) {
            FC = Csum;
            bp = j;
            bt = tb;
          }
        }
      }
      assert(bp > 0);
      std::copy_backward(r.pi.begin()+bp,r.pi.begin()+k,r.pi.begin()+k+1);
      r.pi[bp]=job;
      if (FC>=ub) {
        r.n = k;
        return; 
      }
      if (kk+1<pi.size()) {
        for(unsigned i=1; i<=m; i++)
          for(unsigned j=bp; j<k+1; j++)
            h[i][j] = std::max(h[i-1][j],h[i][j-1])+p[i][r.pi[j]];
        for(unsigned j=bp; j<k+1; j++)
          Cpsum[j] = Cpsum[j-1]+h[m][j];
        for(unsigned i=m; i>=1; i--) {
          for(unsigned j=bp+n-k, je=1+n-k; j>=je; j--) {
            if (t[i][j+1]<=t[i+1][j]) {
              mup[i][j] = mup[i+1][j];
            } else {
              mup[i][j] = i;
            }
            t[i][j] = std::max(t[i][j+1], t[i+1][j])+p[i][r.pi[j+k-n]];
          }
        }
      }
    } 
    r.ms = FC;
  }
  template <typename Tiebreak>
  inline void constructpTT(const Instance& I, PFSSolution& r, unsigned k, unsigned kl, unsigned ku, const std::vector<Job>& pi, Tiebreak tiebreak, Time ub = infinite_time) {
    Time bestTT = infinite_time;
    unsigned tt = 0;
    std::vector<unsigned> Cptt(I.n+1);
    for(unsigned i=m; i>=1; i--)
      for(unsigned j=k-1; j>=1; j--) {
        if (t[i][j+1+n+1-k]<=t[i+1][j+n+1-k]) {
          mup[i][j+n+1-k] = mup[i+1][j+n+1-k];
        } else {
          mup[i][j+n+1-k] = i;
        }
        t[i][j+n+1-k] = std::max(t[i][j+1+n+1-k], t[i+1][j+n+1-k])+p[i][r.pi[j]];
      }
    for(unsigned i=1; i<=m; i++)
      for(unsigned j=1; j<k; j++)
        h[i][j] = std::max(h[i-1][j],h[i][j-1])+p[i][r.pi[j]];
    for(unsigned j=1; j<k; j++) {
      Cpsum[j]=Cpsum[j-1]+h[m][j];
      int tard = h[m][j]-I.dd[r.pi[j]];
      Cptt[j]=Cptt[j-1]+std::max(0,tard);
    }
    for(unsigned kk=1; kk<pi.size(); kk++, k++, ku++) { 
      const unsigned job = pi[kk];
      unsigned bp = 0, bt = std::numeric_limits<unsigned>::max(); 
      bestTT = infinite_time;
      for(unsigned j=kl; j<=ku; j++) { 
        Time Fk = 0;
        unsigned im = 0;
        Time pt = 0;
        auto cp = p_[job].begin()+1;
        for(unsigned i=1; i<=m; i++, cp++) {
          pt = Cm[i] = std::max(pt,h[i][j-1])+*cp;
          Time Fi = pt+t[i][j+n+1-k];
          if (Fi >= Fk) {
            Fk = Fi;
            im = i;
          }
        }
        assert(im>0);
        Csum = Cpsum[j-1] + Cm[m];
        int tard = Cm[m]-I.dd[job];
        tt = Cptt[j-1] + std::max(0, tard);
        auto pcj = r.pi.begin()+j;
        if (j<k) {
          cp = p_[*pcj].begin()+im;
          Cm[im] += *cp++;
          pt = Cm[im];
          for(unsigned i=m-im, ii=im+1; i>0; i--,ii++,cp++)
            pt = Cm[ii] = std::max(pt,Cm[ii])+*cp;
          Csum  += pt;
          int tard = pt-I.dd[*pcj];
          tt += std::max(0, tard);
        }
        pcj++;
        for(unsigned jj=j+1+n+1-k-1, je=k+n+1-k-1; jj<je; jj++, pcj++) {
          im = mup[im][jj];
          assert(im <= m);
          cp = p_[*pcj].begin()+im;
          Cm[im] += *cp++;
          pt = Cm[im];
          for(unsigned i=m-im, ii=im+1; i>0; i--,ii++,cp++)
            pt = Cm[ii] = std::max(pt,Cm[ii])+*cp;
          Csum  += pt;
          int tard = pt-I.dd[*pcj];
          tt += std::max(0, tard);
        }
        {
          unsigned tb = tiebreak(job,j);
          if (tt < bestTT || (tt == bestTT && tb < bt)) {
            bestTT = tt;
            bp = j;
            bt = tb;
          } 
        }
      } 
      assert(bp > 0);
      std::copy_backward(r.pi.begin()+bp,r.pi.begin()+k,r.pi.begin()+k+1);
      r.pi[bp]=job;
      if (bestTT>=ub) {
        r.n = k;
        return; 
      }
      if (kk+1<pi.size()) {
        for(unsigned i=1; i<=m; i++)
          for(unsigned j=bp; j<k+1; j++)
            h[i][j] = std::max(h[i-1][j],h[i][j-1])+p[i][r.pi[j]];
        for(unsigned j=bp; j<k+1; j++) {
          Cpsum[j] = Cpsum[j-1]+h[m][j];
          int tard = h[m][j]-I.dd[r.pi[j]];
          Cptt[j] = Cptt[j-1]+std::max(0,tard);
        }
        for(unsigned i=m; i>=1; i--) {
          for(unsigned j=bp+n-k, je=1+n-k; j>=je; j--) {
            if (t[i][j+1]<=t[i+1][j]) {
              mup[i][j] = mup[i+1][j];
            } else {
              mup[i][j] = i;
            }
            t[i][j] = std::max(t[i][j+1], t[i+1][j])+p[i][r.pi[j+k-n]];
          }
        }
      }
    } 
    r.ms = bestTT;
  }
  inline void constructpTT_IT1(const Instance& I, PFSSolution& r, unsigned k, unsigned kl, unsigned ku, const std::vector<Job>& pi, Time ub = infinite_time) {
    Time bestTT = infinite_time;
    boost::multi_array<unsigned,2> C(boost::extents[m+1][n+1]); 
    std::vector<unsigned> pt(I.n+1);  
    std::vector<unsigned> tpt(I.m+1); 
    for(unsigned j=1; j<k; j++) {
      unsigned Cj = 0;
      for(unsigned i=1; i<=m; i++) {
        C[i][j] = Cj = std::max(Cj, C[i][j-1]) + p[i][r.pi[j]];
        tpt[i] += p[i][r.pi[j]];
      }
      pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[r.pi[j]]));
    }
    for(unsigned kk=1; kk<pi.size(); kk++, k++, ku++) { 
      const unsigned job = pi[kk];
      for(unsigned i=1; i <=m; i++)
        tpt[i] += p[i][job];
      unsigned bp = 0, bt = std::numeric_limits<unsigned>::max(); 
      bestTT = infinite_time;
      for(unsigned j=kl; j<=ku; j++) { 
        unsigned tt = pt[j-1];
        std::vector<unsigned> Cm(I.m+1);
        for(unsigned i=1; i<=m; i++)
          Cm[i] = std::max(Cm[i-1], C[i][j-1]) + p[i][job];
        tt += std::max(0,int(Cm[m]-I.dd[job]));
        for(unsigned jj=j; jj<k; jj++) {
          for(unsigned i=1; i<=m; i++)
            Cm[i] = std::max(Cm[i-1], Cm[i]) + p[i][r.pi[jj]];
          tt += std::max(0,int(Cm[m]-I.dd[r.pi[jj]]));
        }
        unsigned tb = 0;
        for(unsigned i=1; i<=m; i++) {
          assert(Cm[i] >= tpt[i]);
          tb += Cm[i]-tpt[i];
        }
        if (tt < bestTT || (tt == bestTT && tb < bt)) {
          bestTT = tt;
          bp = j;
          bt = tb;
        } 
      } 
      assert(bp > 0);
      std::copy_backward(r.pi.begin()+bp,r.pi.begin()+k,r.pi.begin()+k+1);
      r.pi[bp]=job;
      if (bestTT>=ub) {
        r.n = k;
        return; 
      }
      if (kk+1<pi.size()) {
        for(unsigned j=bp; j<k+1; j++) {
          for(unsigned i=1; i<=m; i++)
            C[i][j] = std::max(C[i-1][j],C[i][j-1]) + p[i][r.pi[j]];
          pt[j] = pt[j-1] + std::max(0,int(C[m][j]-I.dd[r.pi[j]]));
        }
      }
    } 
    r.ms = bestTT;
  }
  inline void constructpTT_IT2(const Instance& I, PFSSolution& r, unsigned k, unsigned kl, unsigned ku, const std::vector<Job>& pi, Time ub = infinite_time) {
    Time bestTT = infinite_time;
    boost::multi_array<unsigned,2> C(boost::extents[m+1][n+1]); 
    std::vector<unsigned> pt(I.n+1);  
    std::vector<unsigned> tit(I.n+1); 
    for(unsigned j=1; j<k; j++) {
      unsigned Cj = 0;
      tit[j] = tit[j-1];
      for(unsigned i=1; i<=m; i++) {
        if(i>1 && j>1)
          tit[j] += std::max(0,int(Cj-C[i][j-1]));
        C[i][j] = Cj = std::max(Cj, C[i][j-1]) + p[i][r.pi[j]];
      }
      pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[r.pi[j]]));
    }
    for(unsigned kk=1; kk<pi.size(); kk++, k++, ku++) { 
      const unsigned job = pi[kk];
      unsigned bp = 0, bt = std::numeric_limits<unsigned>::max(); 
      bestTT = infinite_time;
      for(unsigned j=kl; j<=ku; j++) { 
        unsigned tt = pt[j-1];
        unsigned it = tit[j-1];
        std::vector<unsigned> Cm(I.m+1);
        for(unsigned i=1; i<=m; i++) {
          if(i>1 && j>1)
            it += std::max(0,int(Cm[i-1]-C[i][j-1]));
          Cm[i] = std::max(Cm[i-1], C[i][j-1]) + p[i][job];
        }
        tt += std::max(0,int(Cm[m]-I.dd[job]));
        for(unsigned jj=j; jj<k; jj++) {
          for(unsigned i=1; i<=m; i++) {
            if(i>1 && jj>1)
              it += std::max(0,int(Cm[i-1]-Cm[i]));
            Cm[i] = std::max(Cm[i-1], Cm[i]) + p[i][r.pi[jj]];
          }
          tt += std::max(0,int(Cm[m]-I.dd[r.pi[jj]]));
        }
        unsigned tb = it;
        if (tt < bestTT || (tt == bestTT && tb < bt)) {
          bestTT = tt;
          bp = j;
          bt = tb;
        }
      } 
      assert(bp > 0);
      std::copy_backward(r.pi.begin()+bp,r.pi.begin()+k,r.pi.begin()+k+1);
      r.pi[bp]=job;
      if (bestTT>=ub) {
        r.n = k;
        return; 
      }
      if (kk+1<pi.size()) {
        for(unsigned j=bp; j<k+1; j++) {
          tit[j] = tit[j-1];
          for(unsigned i=1; i<=m; i++) {
            if(i>1 && j>1)
              tit[j] += std::max(0,int(C[i-1][j]-C[i][j-1]));
            C[i][j] = std::max(C[i-1][j],C[i][j-1]) + p[i][r.pi[j]];
          }
          pt[j] = pt[j-1] + std::max(0,int(C[m][j]-I.dd[r.pi[j]]));
        }
      }
    } 
    r.ms = bestTT;
  }
  template <typename Tiebreak>
  inline PFSSolution construct(const Instance& I, const std::vector<Job>& pi, Tiebreak tiebreak, Time ub = infinite_time) {
    PFSSolution r(I);
    constructp(I,r,1,1,1,pi,tiebreak,ub);
    return r;
  }
  inline PFSSolution constructFF(const Instance& I, const std::vector<Job>& pi, Time ub = infinite_time) {
    PFSSolution r(I);
    constructpFF(I,r,1,1,1,pi,ub);
    return r;
  }
  template <typename Tiebreak>
  inline PFSSolution constructCsum(const Instance& I, const std::vector<Job>& pi, Tiebreak tiebreak, Time ub = infinite_time) {
    PFSSolution r(I);
    constructpCsum(I,r,1,1,1,pi,tiebreak,ub);
    return r;
  }
  template <typename Tiebreak>
  inline PFSSolution constructTT(const Instance& I, const std::vector<Job>& pi, Tiebreak tiebreak, Time ub = infinite_time) {
    PFSSolution r(I);
    constructpTT(I,r,1,1,1,pi,tiebreak,ub);
    return r;
  }
  inline PFSSolution constructTT_IT1(const Instance& I, const std::vector<Job>& pi, Time ub = infinite_time) {
    PFSSolution r(I);
    constructpTT_IT1(I,r,1,1,1,pi,ub);
    return r;
  }
  inline PFSSolution constructTT_IT2(const Instance& I, const std::vector<Job>& pi, Time ub = infinite_time) {
    PFSSolution r(I);
    constructpTT_IT2(I,r,1,1,1,pi,ub);
    return r;
  }
  inline PFSSolution neh(const Instance& I, Time ub = infinite_time) {
    std::vector<Job> pi = getTotalTimeOrder(I);
    return construct(I,pi,firstposition,ub);
  }
  inline PFSSolution nehFF(const Instance& I, Time ub = infinite_time) {
    std::vector<Job> pi = getTotalTimeOrderFF(I);
    return constructFF(I,pi,ub);
  }
  inline PFSSolution nehCsum(const Instance& I, Time ub = infinite_time) {
    std::vector<Job> pi = getTotalTimeOrder(I);
    std::reverse(pi.begin()+1,pi.end());
    return constructCsum(I,pi,firstposition,ub);
  }
  inline PFSSolution nehEDD(const Instance& I, Time ub = infinite_time) {
    std::vector<Job> pi = getEDD(I);
    return constructTT(I,pi,firstposition,ub);
  }
  inline PFSSolution nehEDD_IT1(const Instance& I, Time ub = infinite_time) {
    std::vector<Job> pi = getEDD(I);
    return constructTT_IT1(I,pi,ub);
  }
  inline PFSSolution nehEDD_IT2(const Instance& I, Time ub = infinite_time) {
    std::vector<Job> pi = getEDD(I);
    return constructTT_IT2(I,pi,ub);
  }
  inline PFSSolution nehkk1(const Instance& I, Time ub = infinite_time) {
    using namespace std::placeholders;
    std::vector<char> first(I.n+1);
    std::vector<Job> pi = getKK1(I,first);
    return construct(I,pi,std::bind(usebool,_1,_2,first),ub);
  }
  inline PFSSolution nehkk2(const Instance& I, Time ub = infinite_time) {
    using namespace std::placeholders;
    std::vector<char> first(I.n+1);
    std::vector<Job> pi = getKK2(I,first);
    return construct(I,pi,std::bind(usebool,_1,_2,first),ub);
  }
  template <typename Rng>
  inline PFSSolution nehkk2_2r(const Instance& I, Rng& rng, Time ub = infinite_time) {
    using namespace std::placeholders;
    std::vector<char> first(I.n+1);
    std::vector<Job> pi = getKK2(I,first);
    std::uniform_int_distribution<int> dis(1,I.n);
    std::swap(pi[1],pi[dis(rng)]);
    std::swap(pi[2],pi[dis(rng)]);
    return construct(I,pi,std::bind(usebool,_1,_2,first),ub);
  }
  template <typename Tiebreak>
  inline PFSSolution FRB5(const Instance& I, std::vector<Job>& jobOrder, const char obj, Tiebreak tiebreak, unsigned n_ls, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
    PFSSolution r(I);
    for (unsigned i = 1; i <= I.n; i++) {
      Job j = jobOrder[i];
      std::vector<Job> insert = { 0, j };
      if (obj == 'm') 
        constructp(I,r,i,1,i,insert,tiebreak);
      else if (obj == 'f')
        constructpCsum(I,r,i,1,i,insert,tiebreak);
      else
        constructpTT(I,r,i,1,i,insert,tiebreak);
      if (i > 2) {
        bool improved = true;
        while (improved && --n_ls > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
           improved = false;
           unsigned cms = r.ms;
           std::vector<Job> pi = r.pi;
           for(unsigned ji = i; ji >= 1; ji--) {
              j = pi[ji];
              unsigned pos_current_job = std::find(r.pi.begin(),r.pi.end(),j)-r.pi.begin();
              auto nend = std::remove(r.pi.begin(),r.pi.end(),j);
              std::fill(nend, r.pi.end(), 0);
              std::vector<Job> remove = { 0, j };
              if (obj == 'm') {
                unsigned minPos, maxPos;
                if (ji<i) {
                  unsigned new_pos_last_job = std::find(r.pi.begin(),r.pi.end(),pi[ji+1])-r.pi.begin();
                  minPos = std::min(pos_current_job,new_pos_last_job);
                  maxPos = std::max(pos_current_job,new_pos_last_job);
                  if (maxPos == i) 
                    maxPos--;
                  assert(maxPos < i);
                } else {
                  minPos = 1;
                  maxPos = 0;
                }
                constructp(I,r,i,1,i,remove,tiebreak,infinite_time,minPos,maxPos);
              }
              else if (obj == 'f') {
                constructpCsum(I,r,i,1,i,remove,tiebreak);
              }
              else {
                constructpTT(I,r,i,1,i,remove,tiebreak);
              }
           }
           if(r.ms < cms)
            improved = true;
        }
      }
    }
    return r;
  }
  template <typename Tiebreak>
  inline PFSSolution FRB5(const Instance& I, const char obj, Tiebreak tiebreak, unsigned n_ls, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
    PFSSolution r(I);
    std::vector<Job> jobs = getTotalTimeOrder(I);
    if (obj != 'm'){ 
      std::reverse(jobs.begin()+1,jobs.end());
    }
    return FRB5(I,jobs,obj,tiebreak,n_ls,start,timelimit);
  }
  inline PFSSolution FRB5_FF(const Instance& I, std::vector<Job>& jobOrder, const char obj, unsigned n_ls, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
    PFSSolution r(I);
    for (unsigned i = 1; i <= I.n; i++) {
      Job j = jobOrder[i];
      std::vector<Job> insert = { 0, j };
      constructpFF(I,r,i,1,i,insert);
      if (i > 2) {
        bool improved = true;
        while (improved && --n_ls > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
           improved = false;
           unsigned cms = r.ms;
           std::vector<Job> pi = r.pi;
           for(unsigned ji = i; ji >= 1; ji--) {
              j = pi[ji];
              unsigned pos_current_job = std::find(r.pi.begin(),r.pi.end(),j)-r.pi.begin();
              auto nend = std::remove(r.pi.begin(),r.pi.end(),j);
              std::fill(nend, r.pi.end(), 0);
              std::vector<Job> remove = { 0, j };
              unsigned minPos, maxPos;
              if (ji<i) {
                unsigned new_pos_last_job = std::find(r.pi.begin(),r.pi.end(),pi[ji+1])-r.pi.begin();
                minPos = std::min(pos_current_job,new_pos_last_job);
                maxPos = std::max(pos_current_job,new_pos_last_job);
                if (maxPos == i) 
                  maxPos--;
                assert(maxPos < i);
              } else {
                minPos = 1;
                maxPos = 0;
              }
              constructpFF(I,r,i,1,i,remove,infinite_time,minPos,maxPos);
           }
           if(r.ms < cms)
            improved = true;
        }
      }
    }
    return r;
  }
  inline PFSSolution FRB5_IT1(const Instance& I, std::vector<Job>& jobOrder, const char obj, unsigned n_ls, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
    PFSSolution r(I);
    for (unsigned i = 1; i <= I.n; i++) {
      Job j = jobOrder[i];
      std::vector<Job> insert = { 0, j };
      constructpTT_IT1(I,r,i,1,i,insert);
      if (i > 2) {
        bool improved = true;
        while (improved && --n_ls > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
           improved = false;
           unsigned cms = r.ms;
           std::vector<Job> pi = r.pi;
           for(unsigned ji = i; ji >= 1; ji--) {
              j = pi[ji];
              auto nend = std::remove(r.pi.begin(),r.pi.end(),j);
              std::fill(nend, r.pi.end(), 0);
              std::vector<Job> remove = { 0, j };
              constructpTT_IT1(I,r,i,1,i,remove);
           }
           if(r.ms < cms)
            improved = true;
        }
      }
    }
    return r;
  }
  inline PFSSolution FRB5_IT2(const Instance& I, std::vector<Job>& jobOrder, const char obj, unsigned n_ls, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
    PFSSolution r(I);
    for (unsigned i = 1; i <= I.n; i++) {
      Job j = jobOrder[i];
      std::vector<Job> insert = { 0, j };
      constructpTT_IT2(I,r,i,1,i,insert);
      if (i > 2) {
        bool improved = true;
        while (improved && --n_ls > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
           improved = false;
           unsigned cms = r.ms;
           std::vector<Job> pi = r.pi;
           for(unsigned ji = i; ji >= 1; ji--) {
              j = pi[ji];
              auto nend = std::remove(r.pi.begin(),r.pi.end(),j);
              std::fill(nend, r.pi.end(), 0);
              std::vector<Job> remove = { 0, j };
              constructpTT_IT2(I,r,i,1,i,remove);
           }
           if(r.ms < cms)
            improved = true;
        }
      }
    }
    return r;
  }
  template <typename Rng>
  inline PFSSolution repeatedConstr(const Instance& I, Rng& rng, const char obj, const unsigned rlimit, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
    PFSSolution s(I), ns(I);
    std::uniform_real_distribution<double> dis(0.0,1.0);
    std::vector<double> randj(I.n+1);
    for(unsigned j = 1; j <= I.n; j++)
      randj[j] = dis(rng);
    std::vector<Time> T(I.n+1,0);
    for(unsigned j = 1; j <= I.n; j++)
      for(unsigned i = 1; i <= I.m; i++)
        T[j] += I.p[j][i];
    std::vector<Job> pi(I.n+1);
    std::iota(pi.begin(), pi.end(),0);
    sort(pi.begin()+1, pi.end(), [&T,&randj](Job i, Job j) {
      return T[i]>T[j] || (T[i]==T[j] && randj[i]<randj[j]);
    });
    if (obj == 'f')
      std::reverse(pi.begin()+1,pi.end());
    unsigned r = 0;
    do {
      unsigned j = 1;
      while (j < I.n) {
        unsigned k = j+1;
        while (k <= I.n && T[pi[j]] == T[pi[k]])
          k++;
        if (k-j > 1)
          std::shuffle(pi.begin()+j, pi.begin()+k, rng);
        j = k;
      }
      if(obj == 'm') 
        ns = construct(I, pi, randomtiebreak);
      else
        ns = constructCsum(I, pi, randomtiebreak);
      if (ns.ms < s.ms)
        s = ns;
      r++;
    } while (r < rlimit && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
    return s;
  }
  template <typename Rng>
  inline PFSSolution repeatedConstrTT(const Instance& I, Rng& rng, const char obj, const unsigned rlimit, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
    PFSSolution s(I), ns(I);
    std::uniform_real_distribution<double> dis(0.0,1.0);
    std::vector<double> randj(I.n+1);
    for(unsigned j = 1; j <= I.n; j++)
      randj[j] = dis(rng);
    std::vector<unsigned> dd = I.dd;
    std::vector<Job> pi(I.n+1);
    std::iota(pi.begin(), pi.end(),0);
    sort(pi.begin()+1, pi.end(), [&dd,&randj](Job i, Job j) {
      return dd[i]<dd[j] || (dd[i]==dd[j] && randj[i]<randj[j]);
    });
    unsigned r = 0;
    do {
      unsigned j = 1;
      while (j < I.n) {
        unsigned k = j+1;
        while (k <= I.n && I.dd[pi[j]] == I.dd[pi[k]])
          k++;
        if (k-j > 1)
          std::shuffle(pi.begin()+j, pi.begin()+k, rng);
        j = k;
      }
      ns = constructTT(I, pi, randomtiebreak);
      if (ns.ms < s.ms)
        s = ns;
      r++;
    } while (r < rlimit && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
    return s;
  }
  template <typename Tiebreak>
  PFSSolution G8(const Instance& I, const std::vector<Job>& jobOrder, const char obj, Tiebreak tiebreak, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
    PFSSolution s(I);
    s.pi[1] = jobOrder[1];
    for(unsigned j = 2; j <= I.n; j++) {
      Job job = jobOrder[j];
      std::vector<Job> insert = { 0, job };
      constructp(I,s,j,1,j,insert,tiebreak);
      if(j > 4 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
        for(unsigned k = 1; k < j; k+=2) {
          Job job1 = s.pi[k];
          Job job2 = s.pi[k+1];
          s.removeJob(job1);
          s.removeJob(job2);
          std::vector<Job> insert = { 0, job1, job2 };
          constructp(I,s,j-1,1,j-1,insert,tiebreak);
        }
        for(unsigned k = 2; k < j-1; k+=2) {
          Job job1 = s.pi[k];
          Job job2 = s.pi[k+1];
          s.removeJob(job1);
          s.removeJob(job2);
          std::vector<Job> insert = { 0, job1, job2 };
          constructp(I,s,j-1,1,j-1,insert,tiebreak);
        }
      }
    }
    return s;
  }
  PFSSolution G8_FF(const Instance& I, const std::vector<Job>& jobOrder, const char obj, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
    PFSSolution s(I);
    s.pi[1] = jobOrder[1];
    for(unsigned j = 2; j <= I.n; j++) {
      Job job = jobOrder[j];
      std::vector<Job> insert = { 0, job };
      constructpFF(I,s,j,1,j,insert);
      if(j > 4 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
        for(unsigned k = 1; k < j; k+=2) {
          Job job1 = s.pi[k];
          Job job2 = s.pi[k+1];
          s.removeJob(job1);
          s.removeJob(job2);
          std::vector<Job> insert = { 0, job1, job2 };
          constructpFF(I,s,j-1,1,j-1,insert);
        }
        for(unsigned k = 2; k < j-1; k+=2) {
          Job job1 = s.pi[k];
          Job job2 = s.pi[k+1];
          s.removeJob(job1);
          s.removeJob(job2);
          std::vector<Job> insert = { 0, job1, job2 };
          constructpFF(I,s,j-1,1,j-1,insert);
        }
      }
    }
    return s;
  }
};
struct nBuilder {
  unsigned n,m;
  const Instance& I;
  boost::multi_array<Time,2> p; 
  boost::multi_array<Time,2> h, t;      
  boost::multi_array<bool,2> tup;       
  boost::multi_array<unsigned,2> mi, Mi;
  std::vector<unsigned> mu, nu;         
  std::vector<Time> C;                  
  nBuilder(const Instance& _I) : n(_I.n), m(_I.m), I(_I), p(boost::extents[n+1][m+1]),
         h(boost::extents[m+1][n+1]), t(boost::extents[m+2][n+2]), tup(boost::extents[m+2][n+2]), mi(boost::extents[m+1][n+2]), Mi(boost::extents[m+1][n+1]), mu(m+1,0), nu(m+1,0), C(n+1,0)
  {
    for(unsigned i=0; i<=m; i++)
      for(unsigned j=0; j<=n; j++)
        p[j][i] = _I.p[j][i];
  }
  unsigned gij(const boost::multi_array<unsigned,2>& pi, unsigned i, unsigned j, unsigned job, unsigned k) {
    if (j<k)
      return pi[i][j];
    if (j==k)
      return job;
    assert(j>0);
    return pi[i][j-1];
  }
  Time propagateNEW(std::vector<Time>& C, const boost::multi_array<unsigned,2>& npi, const boost::multi_array<unsigned,2>& invpi, const unsigned i, unsigned lj, unsigned uj) {
    Time ct = h[i][lj-1];
    for(unsigned ji=lj; ji<=uj; ji++) { 
      unsigned cj = npi[i][ji];
      if (C[cj]<infinite_time)
        ct = std::max(ct,C[cj])+p[cj][i];
      else
        ct = std::max(ct,h[i-1][invpi[i-1][cj]])+p[cj][i];
      C[cj] = ct;
    }
    ct += t[i][uj];
    return ct;
  }
  void propagateNEWFLOW(std::vector<Time>& C, const boost::multi_array<unsigned,2>& npi, const boost::multi_array<unsigned,2>& invpi, const unsigned i, unsigned lj, unsigned uj) {
    Time ct = h[i][lj-1];
    unsigned cj;
    if (uj<=nu[i]) {
      for(unsigned ji=lj; ji<=uj; ji++) { 
        cj = npi[i][ji];
        if (C[cj]<infinite_time)
          ct = std::max(ct,C[cj]);
        else
          ct = std::max(ct,h[i-1][invpi[i-1][cj]]);
        ct += p[cj][i];
        C[cj] = ct;
      }
    } else {
      for(unsigned ji=lj; ji<=nu[i]; ji++) {
        cj = npi[i][ji];
        if (C[cj]<infinite_time)
          ct = std::max(ct,C[cj]);
        else
          ct = std::max(ct,h[i-1][invpi[i-1][cj]]);
        ct += p[cj][i];
        C[cj] = ct;
      }
      for(unsigned ji=nu[i]+1; ji<=uj; ji++) {
        cj = npi[i][ji];
        if (C[cj]<infinite_time && C[cj]>ct)
          ct = C[cj];
        ct += p[cj][i];
        C[cj] = ct;
      }
    }
  }
  void constructnpNEW(NPFSSolution& r, unsigned k, const std::vector<Job>& pi, const std::vector<char>& first, Time ub = infinite_time, unsigned lbf = 1, unsigned lbb = 1) {
    assert(constructionPrecondition(I,k,pi,first));
    Time FC = infinite_time;
    for(unsigned j=0; j<=n; j++) {
      h[0][j] = 0;     
      t[m+1][j] = 0;   
    }
    for(unsigned i=0; i<=m; i++) {
      h[i][0] = 0;     
      for(unsigned j=k; j<=n+1; j++)
        t[i][j] = 0;   
    }
    boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
    unsigned bp; 
    unsigned bl; 
    bool     bf; 
    for(unsigned kk=1; kk<pi.size(); kk++, k++) { 
      const unsigned job = pi[kk];
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++) {
          npi[i][j]=r.pi[i][j];
          invpi[i][npi[i][j]]=j;
        }
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++)
          h[i][j] = std::max(h[i-1][invpi[i-1][npi[i][j]]],h[i][j-1])+p[r.pi[i][j]][i];
      for(unsigned i=m; i>=1; i--)
        for(unsigned j=k-1; j>=1; j--)
          t[i][j] = std::max(t[i+1][invpi[i+1][npi[i][j]]],t[i][j+1])+p[r.pi[i][j]][i];
      for(unsigned i=1; i<m; i++) {
        mi[i][k]=k; 
        for(unsigned j=k-1; j>0; j--)
          mi[i][j] = std::min(mi[i][j+1],invpi[i+1][npi[i][j]]);
      }
      for(unsigned i=m; i>1; i--)
        for(unsigned j=1; j<=k; j++)
          Mi[i][j] = std::max(Mi[i][j-1],invpi[i-1][npi[i][j-1]]);
      bp = bl = 0; 
      bf = true;   
      FC = infinite_time;    
      Time Fk;
      for(unsigned i=1; i<=m; i++)
        npi[i][k]=job;
      for(unsigned j=k; j>=1; j--) { 
        Fk = 0;
        fill(C.begin(), C.end(), infinite_time);
        mu[1] = j;
        for(unsigned i=1; i<m; i++)
          mu[i+1] = mi[i][mu[i]];
        nu[m] = j;
        for(unsigned i=m; i>1; i--)
          nu[i-1] = Mi[i][nu[i]]+1;
        for(unsigned i=1; i<=m; i++) {
          Time ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
          Fk = std::max(Fk,ct);
          if (Fk>FC) {
            break;
          }
        }
        if (Fk < FC || (Fk == FC && !first[job])) {
          FC = Fk;
          bp = j;
          bl = m;
          bf = false;
        }
        if (j>1)
          for(unsigned i=1; i<=m; i++)
            std::swap(npi[i][j],npi[i][j-1]);
      } 
#define BACKWARD_ENABLED
#if defined(BACKWARD_ENABLED)
      lbb = std::max(2u,lbb);
      if (k>1 && m>3) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        for(unsigned i=3; i<=m; i++)
          std::swap(npi[i][k-1],npi[i][k]);
        for(unsigned j=k; j>=lbb; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Fk = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          if (Fk>FC) {
            for(unsigned i=3; i<m; i++)
              std::swap(npi[i][j-1],npi[i][j]);
            goto nextj;
          }
          for(unsigned l=2; l<m-1; l++) {
            mu[l] = mi[l-1][mu[l-1]];
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            Fk = std::max(Fk,ct);
            if (Fk>FC) {
              for(unsigned i=l+1; i<m; i++)
                std::swap(npi[i][j-1],npi[i][j]);
              goto nextj;
            }
            C = Cl;
            Time Fi = Fk;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = std::min(mi[i-1][mu[i-1]],j-1);
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              Fi = std::max(Fi,ct);
              if (Fi>FC && i<m) {
                break;
              }
            }
            if (Fi < FC || (Fi == FC && ((!first[job] && j>=bp) || (first[job] && j<=bp)))) {  
              FC = Fi;
              bp = j;
              bl = l;
              bf = false;
            }
            std::swap(npi[l+1][j-1],npi[l+1][j]);
          }
          nextj:
          std::swap(npi[m][j-1],npi[m][j]);
          if (j>2) {
            for(unsigned i=1; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j]);
            for(unsigned i=3; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j-2]);
          }
        } 
      } 
#endif
#define FORWARD_ENABLED
#if defined(FORWARD_ENABLED)
      if (k>1 && m>3) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        std::swap(npi[1][k-1],npi[1][k]);
        std::swap(npi[2][k-1],npi[2][k]);
        for(unsigned j=k-1; j>=lbf; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j+1;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Fk = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          if (Fk>FC) {
            for(unsigned i=3; i<m; i++)
              std::swap(npi[i][j],npi[i][j+1]);
            goto nextj2;
          }
          for(unsigned l=2; l<m-1; l++) { 
            mu[l] = mi[l-1][mu[l-1]];
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            Fk = std::max(Fk,ct);
            if (Fk>FC) {
              for(unsigned i=l+1; i<m; i++)
                std::swap(npi[i][j],npi[i][j+1]);
              goto nextj2;
            }
            C = Cl;
            Time Fi = Fk;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = mi[i-1][mu[i-1]];
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              Fi = std::max(Fi,ct);
              if (Fi>FC && i<m)
                break;
            }
            if (Fi < FC || (Fi == FC && ((!first[job] && j>bp) || (first[job] && j<bp)))) {  
              FC = Fi;
              bp = j;
              bl = l;
              bf = true;
            }
            std::swap(npi[l+1][j],npi[l+1][j+1]);
          }
          nextj2:;
          std::swap(npi[m][j],npi[m][j+1]);
          if (j>1) {
            std::swap(npi[1][j-1],npi[1][j]);
            std::swap(npi[2][j-1],npi[2][j]);
          }
        } 
      } 
#endif
      assert(1<=bp && bp<=k);
      for(unsigned i=1; i<=m; i++)
        if (i<=bl) {
          for(unsigned p=k-1; p>=bp; p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][bp]=job;
        } 
        else {
          for(unsigned p=k-1; p>=(bf?bp+1:bp-1); p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][(bf?bp+1:bp-1)]=job;
        }
      if (FC>=ub) {
        r.n = k;
        return; 
      }
    } 
    r.ms = FC;
  }
  template <typename Tiebreak>
  void constructnpNEW(NPFSSolution& r, unsigned k, const std::vector<Job>& pi, Tiebreak tiebreak, Time ub = infinite_time, unsigned lbf = 1, unsigned lbb = 1) {
    Time FC = infinite_time;
    for(unsigned j=0; j<=n; j++) {
      h[0][j] = 0;     
      t[m+1][j] = 0;   
    }
    for(unsigned i=0; i<=m; i++) {
      h[i][0] = 0;     
      for(unsigned j=k; j<=n+1; j++)
        t[i][j] = 0;   
    }
    boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
    unsigned bp; 
    unsigned bt; 
    unsigned bl; 
    bool     bf; 
    for(unsigned kk=1; kk<pi.size(); kk++, k++) { 
      const unsigned job = pi[kk];
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++) {
          npi[i][j]=r.pi[i][j];
          invpi[i][npi[i][j]]=j;
        }
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++)
          h[i][j] = std::max(h[i-1][invpi[i-1][npi[i][j]]],h[i][j-1])+p[r.pi[i][j]][i];
      for(unsigned i=m; i>=1; i--)
        for(unsigned j=k-1; j>=1; j--)
          t[i][j] = std::max(t[i+1][invpi[i+1][npi[i][j]]],t[i][j+1])+p[r.pi[i][j]][i];
      for(unsigned i=1; i<m; i++) {
        mi[i][k]=k; 
        for(unsigned j=k-1; j>0; j--)
          mi[i][j] = std::min(mi[i][j+1],invpi[i+1][npi[i][j]]);
      }
      for(unsigned i=m; i>1; i--)
        for(unsigned j=1; j<=k; j++)
          Mi[i][j] = std::max(Mi[i][j-1],invpi[i-1][npi[i][j-1]]);
      bp = bl = 0; 
      bf = true;   
      bt = std::numeric_limits<unsigned>::max();
      FC = infinite_time;    
      Time Fk;
      for(unsigned i=1; i<=m; i++)
        npi[i][k]=job;
      for(unsigned j=k; j>=1; j--) { 
        Fk = 0;
        fill(C.begin(), C.end(), infinite_time);
        mu[1] = j;
        for(unsigned i=1; i<m; i++)
          mu[i+1] = mi[i][mu[i]];
        nu[m] = j;
        for(unsigned i=m; i>1; i--)
          nu[i-1] = Mi[i][nu[i]]+1;
        for(unsigned i=1; i<=m; i++) {
          Time ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
          Fk = std::max(Fk,ct);
          if (Fk>FC)
            break;
        }
        unsigned tb = tiebreak(job,j);
        if (Fk < FC || (Fk == FC && tb < bt)) {
          FC = Fk;
          bp = j;
          bl = m;
          bf = false;
          bt = tb;
        }
        if (j>1)
          for(unsigned i=1; i<=m; i++)
            std::swap(npi[i][j],npi[i][j-1]);
      } 
#define BACKWARD_ENABLED
#if defined(BACKWARD_ENABLED)
      lbb = std::max(2u,lbb);
      if (k>1 && m>3) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        for(unsigned i=3; i<=m; i++)
          std::swap(npi[i][k-1],npi[i][k]);
        for(unsigned j=k; j>=lbb; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Fk = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          if (Fk>FC) {
            for(unsigned i=3; i<m; i++)
              std::swap(npi[i][j-1],npi[i][j]);
            goto nextj;
          }
          for(unsigned l=2; l<m-1; l++) {
            mu[l] = mi[l-1][mu[l-1]];
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            Fk = std::max(Fk,ct);
            if (Fk>FC) {
              for(unsigned i=l+1; i<m; i++)
                std::swap(npi[i][j-1],npi[i][j]);
              goto nextj;
            }
            C = Cl;
            Time Fi = Fk;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = std::min(mi[i-1][mu[i-1]],j-1);
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              Fi = std::max(Fi,ct);
              if (Fi>FC && i<m)
                break;
            }
            unsigned tb = tiebreak(job,j);
            if (Fi < FC || (Fi == FC && tb<bt)) {
              FC = Fi;
              bp = j;
              bl = l;
              bf = false;
              bt = tb;
            }
            std::swap(npi[l+1][j-1],npi[l+1][j]);
          }
          nextj:
          std::swap(npi[m][j-1],npi[m][j]);
          if (j>2) {
            for(unsigned i=1; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j]);
            for(unsigned i=3; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j-2]);
          }
        } 
      } 
#endif
#define FORWARD_ENABLED
#if defined(FORWARD_ENABLED)
      if (k>1 && m>3) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        std::swap(npi[1][k-1],npi[1][k]);
        std::swap(npi[2][k-1],npi[2][k]);
        for(unsigned j=k-1; j>=lbf; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j+1;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Fk = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          if (Fk>FC) {
            for(unsigned i=3; i<m; i++)
              std::swap(npi[i][j],npi[i][j+1]);
            goto nextj2;
          }
          for(unsigned l=2; l<m-1; l++) { 
            mu[l] = mi[l-1][mu[l-1]];
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            Fk = std::max(Fk,ct);
            if (Fk>FC) {
              for(unsigned i=l+1; i<m; i++)
                std::swap(npi[i][j],npi[i][j+1]);
              goto nextj2;
            }
            C = Cl;
            Time Fi = Fk;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = mi[i-1][mu[i-1]];
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              Fi = std::max(Fi,ct);
              if (Fi>FC && i<m)
                break;
            }
            unsigned tb = tiebreak(job,j);
            if (Fi < FC || (Fi == FC && tb<bt)) {
              FC = Fi;
              bp = j;
              bl = l;
              bf = true;
              bt = tb;
            }
            std::swap(npi[l+1][j],npi[l+1][j+1]);
          }
          nextj2:;
          std::swap(npi[m][j],npi[m][j+1]);
          if (j>1) {
            std::swap(npi[1][j-1],npi[1][j]);
            std::swap(npi[2][j-1],npi[2][j]);
          }
        } 
      } 
#endif
      assert(1<=bp && bp<=k);
      for(unsigned i=1; i<=m; i++)
        if (i<=bl) {
          for(unsigned p=k-1; p>=bp; p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][bp]=job;
        } 
        else {
          for(unsigned p=k-1; p>=(bf?bp+1:bp-1); p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][(bf?bp+1:bp-1)]=job;
        }
      if (FC>=ub) {
        r.n = k;
        return; 
      }
    } 
    r.ms = FC;
  }
  Time computeFlowtime(unsigned k, unsigned imax, const std::vector<Time>& Catmax, std::vector<Time>& C, const boost::multi_array<unsigned,2>& npi, const boost::multi_array<unsigned,2>& invpi) {
    Time Fk = 0;
    unsigned cnu = nu[imax];
    if (imax<m || cnu<k) {
      C = Catmax;
      cnu++;
      for(unsigned i=imax; i<=m; i++) {
        while (!tup[i][cnu-1]) {
          assert(t[i][cnu-1] == t[i][cnu]+ p[npi[i][cnu]][i]);
          cnu++;
        }
        propagateNEWFLOW(C,npi,invpi,i,mu[i],cnu); 
      }
    }
    for(unsigned jj=1; jj<=k; jj++) {
      unsigned cj = npi[m][jj];
      if (C[cj]<infinite_time)
        Fk += C[cj];
      else
        Fk += h[m][jj];
    }
    return Fk;
  }
  Time computeTotalTardiness(unsigned k, unsigned imax, const std::vector<Time>& Catmax, std::vector<Time>& C, const boost::multi_array<unsigned,2>& npi, const boost::multi_array<unsigned,2>& invpi) {
    Time Fk = 0;
    unsigned cnu = nu[imax];
    if (imax<m || cnu<k) {
      C = Catmax;
      cnu++;
      for(unsigned i=imax; i<=m; i++) {
        while (!tup[i][cnu-1]) {
          assert(t[i][cnu-1] == t[i][cnu]+ p[npi[i][cnu]][i]);
          cnu++;
        }
        propagateNEWFLOW(C,npi,invpi,i,mu[i],cnu); 
      }
    }
    for(unsigned jj=1; jj<=k; jj++) {
      unsigned cj = npi[m][jj];
      if (C[cj]<infinite_time) {
        int tard = C[cj]-I.dd[cj];
        Fk += std::max(0, tard);
      } else {
        int tard = h[m][jj]-I.dd[cj];
        Fk += std::max(0,tard);
      }
    }
    return Fk;
  }
  unsigned computeIT2(unsigned k, const boost::multi_array<unsigned,2>& npi) {
    unsigned tit = 0;
    std::vector<unsigned> C(m+1);
    for(unsigned j=1; j<k; j++) {
      unsigned Cj = 0;
      for(unsigned i=1; i<=m; i++) {
        if(i>1 && j>1)
          tit += std::max(0,int(Cj-C[i]));
        C[i] = Cj = std::max(Cj, C[i-1]) + p[i][npi[i][j]];
      }
    }
    return tit;
  }
  void constructnpNEWCsum(NPFSSolution& r, unsigned k, const std::vector<Job>& pi, const std::vector<char>& first, Time ub = infinite_time, unsigned lbf = 1, unsigned lbb = 1) {
    assert(constructionPrecondition(I,k,pi,first));
    Time FC = infinite_time;
    for(unsigned j=0; j<=n; j++) {
      h[0][j] = 0;     
      t[m+1][j] = 0;   
    }
    for(unsigned i=0; i<=m; i++) {
      h[i][0] = 0;     
      for(unsigned j=k; j<=n+1; j++)
        t[i][j] = 0;   
    }
    boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
    unsigned bp; 
    unsigned bl; 
    bool     bf; 
    for(unsigned kk=1; kk<pi.size(); kk++, k++) { 
      const unsigned job = pi[kk];
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++) {
          npi[i][j]=r.pi[i][j];
          invpi[i][npi[i][j]]=j;
        }
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++)
          h[i][j] = std::max(h[i-1][invpi[i-1][npi[i][j]]],h[i][j-1])+p[r.pi[i][j]][i];
      for(unsigned i=m; i>=1; i--)
        for(unsigned j=k-1; j>=1; j--) {
          tup[i][j] = (t[i+1][invpi[i+1][npi[i][j]]]>=t[i][j+1]);
          t[i][j] = std::max(t[i+1][invpi[i+1][npi[i][j]]],t[i][j+1])+p[r.pi[i][j]][i];
        }
      for(unsigned i=1; i<m; i++) {
        mi[i][k]=k; 
        for(unsigned j=k-1; j>0; j--)
          mi[i][j] = std::min(mi[i][j+1],invpi[i+1][npi[i][j]]);
      }
      for(unsigned i=m; i>1; i--)
        for(unsigned j=1; j<=k; j++)
          Mi[i][j] = std::max(Mi[i][j-1],invpi[i-1][npi[i][j-1]]);
      bp = bl = 0; 
      bf = true;   
      FC = infinite_time;    
      Time Fk;
      Time Cmax;
      unsigned imax;
      std::vector<Time> Catmax(C), Cold(C), Catmaxi(C); 
      for(unsigned i=1; i<=m; i++)
        npi[i][k]=job;
      for(unsigned j=k; j>=1; j--) { 
        Cmax = 0;
        imax = 1;
        fill(C.begin(), C.end(), infinite_time);
        mu[1] = j;
        for(unsigned i=1; i<m; i++)
          mu[i+1] = mi[i][mu[i]];
        nu[m] = j;
        for(unsigned i=m; i>1; i--)
          nu[i-1] = Mi[i][nu[i]]+1;
        for(unsigned i=1; i<=m; i++) {
          Cold = C;
          Time ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
          if (ct>=Cmax) {
            Cmax = ct;
            Catmax = Cold;
            imax = i;
          }
        }
        Fk = computeFlowtime(k,imax,Catmax,C,npi,invpi);
        if (Fk < FC || (Fk == FC && !first[job])) {
          FC = Fk;
          bp = j;
          bl = m;
          bf = false;
        }
        if (j>1)
          for(unsigned i=1; i<=m; i++)
            std::swap(npi[i][j],npi[i][j-1]);
      } 
#define FLOW_BACKWARD_ENABLED
#if defined(FLOW_BACKWARD_ENABLED)
      lbb = std::max(2u,lbb);
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        for(unsigned i=3; i<=m; i++)
          std::swap(npi[i][k-1],npi[i][k]);
        for(unsigned j=k; j>=lbb; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          for(unsigned l=2; l<m; l++) {
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = std::min(mi[i-1][mu[i-1]],j-1);
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeFlowtime(k,imaxi,Catmaxi,C,npi,invpi);
            if (Fi < FC || (Fi == FC && ((!first[job] && j>=bp) || (first[job] && j<=bp)))) {  
              FC = Fi;
              bp = j;
              bl = l;
              bf = false;
            }
            std::swap(npi[l+1][j-1],npi[l+1][j]);
          }
          if (j>2) {
            for(unsigned i=1; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j]);
            for(unsigned i=3; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j-2]);
          }
        } 
      } 
#endif
#define FLOW_FORWARD_ENABLED
#if defined(FLOW_FORWARD_ENABLED)
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        std::swap(npi[1][k-1],npi[1][k]);
        std::swap(npi[2][k-1],npi[2][k]);
        for(unsigned j=k-1; j>=lbf; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j+1;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          for(unsigned l=2; l<m; l++) { 
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = mi[i-1][mu[i-1]];
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeFlowtime(k,imaxi,Catmaxi,C,npi,invpi);
            if (Fi < FC || (Fi == FC && ((!first[job] && j>bp) || (first[job] && j<bp)))) {
              FC = Fi;
              bp = j;
              bl = l;
              bf = true;
            }
            std::swap(npi[l+1][j],npi[l+1][j+1]);
          }
          if (j>1) {
            std::swap(npi[1][j-1],npi[1][j]);
            std::swap(npi[2][j-1],npi[2][j]);
          }
        } 
      } 
#endif
      assert(1<=bp && bp<=k);
      for(unsigned i=1; i<=m; i++) {
        if (i<=bl) {
          for(unsigned p=k-1; p>=bp; p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][bp]=job;
        } else {
          for(unsigned p=k-1; p>=(bf?bp+1:bp-1); p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][(bf?bp+1:bp-1)]=job;
        }
      }
      if (FC>=ub) {
        r.n = k;
        return; 
      }
      assert(flowtimeConsistent(I,r,FC,k+1));
    } 
    assert(flowtimeConsistent(I,r,FC,n+1));
    r.ms = FC;
  }
  template<typename Tiebreak>
  void constructnpNEWCsum(NPFSSolution& r, unsigned k, const std::vector<Job>& pi, Tiebreak tiebreak, Time ub = infinite_time, unsigned lbf = 1, unsigned lbb = 1) {
    Time FC = infinite_time;
    for(unsigned j=0; j<=n; j++) {
      h[0][j] = 0;     
      t[m+1][j] = 0;   
    }
    for(unsigned i=0; i<=m; i++) {
      h[i][0] = 0;     
      for(unsigned j=k; j<=n+1; j++)
        t[i][j] = 0;   
    }
    boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
    unsigned bp; 
    unsigned bt; 
    unsigned bl; 
    bool     bf; 
    for(unsigned kk=1; kk<pi.size(); kk++, k++) { 
      const unsigned job = pi[kk];
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++) {
          npi[i][j]=r.pi[i][j];
          invpi[i][npi[i][j]]=j;
        }
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++)
          h[i][j] = std::max(h[i-1][invpi[i-1][npi[i][j]]],h[i][j-1])+p[r.pi[i][j]][i];
      for(unsigned i=m; i>=1; i--)
        for(unsigned j=k-1; j>=1; j--) {
          tup[i][j] = (t[i+1][invpi[i+1][npi[i][j]]]>=t[i][j+1]);
          t[i][j] = std::max(t[i+1][invpi[i+1][npi[i][j]]],t[i][j+1])+p[r.pi[i][j]][i];
        }
      for(unsigned i=1; i<m; i++) {
        mi[i][k]=k; 
        for(unsigned j=k-1; j>0; j--)
          mi[i][j] = std::min(mi[i][j+1],invpi[i+1][npi[i][j]]);
      }
      for(unsigned i=m; i>1; i--)
        for(unsigned j=1; j<=k; j++)
          Mi[i][j] = std::max(Mi[i][j-1],invpi[i-1][npi[i][j-1]]);
      bp = bl = 0; 
      bf = true;   
      bt = std::numeric_limits<unsigned>::max();
      FC = infinite_time;    
      Time Fk;
      Time Cmax;
      unsigned imax;
      std::vector<Time> Catmax(C), Cold(C), Catmaxi(C); 
      for(unsigned i=1; i<=m; i++)
        npi[i][k]=job;
      for(unsigned j=k; j>=1; j--) { 
        Cmax = 0;
        imax = 1;
        fill(C.begin(), C.end(), infinite_time);
        mu[1] = j;
        for(unsigned i=1; i<m; i++)
          mu[i+1] = mi[i][mu[i]];
        nu[m] = j;
        for(unsigned i=m; i>1; i--)
          nu[i-1] = Mi[i][nu[i]]+1;
        for(unsigned i=1; i<=m; i++) {
          Cold = C;
          Time ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
          if (ct>=Cmax) {
            Cmax = ct;
            Catmax = Cold;
            imax = i;
          }
        }
        Fk = computeFlowtime(k,imax,Catmax,C,npi,invpi);
        unsigned tb = tiebreak(job,j);
        if (Fk < FC || (Fk == FC && tb<bt)) {
          FC = Fk;
          bp = j;
          bl = m;
          bf = false;
          bt = tb;
        }
        if (j>1)
          for(unsigned i=1; i<=m; i++)
            std::swap(npi[i][j],npi[i][j-1]);
      } 
#define FLOW_BACKWARD_ENABLED
#if defined(FLOW_BACKWARD_ENABLED)
      lbb = std::max(2u,lbb);
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        for(unsigned i=3; i<=m; i++)
          std::swap(npi[i][k-1],npi[i][k]);
        for(unsigned j=k; j>=lbb; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          for(unsigned l=2; l<m; l++) {
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = std::min(mi[i-1][mu[i-1]],j-1);
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeFlowtime(k,imaxi,Catmaxi,C,npi,invpi);
            unsigned tb = tiebreak(job,j);
            if (Fi < FC || (Fi == FC && tb<bt)) {  
              FC = Fi;
              bp = j;
              bl = l;
              bf = false;
              bt = tb;
            }
            std::swap(npi[l+1][j-1],npi[l+1][j]);
          }
          if (j>2) {
            for(unsigned i=1; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j]);
            for(unsigned i=3; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j-2]);
          }
        } 
      } 
#endif
#define FLOW_FORWARD_ENABLED
#if defined(FLOW_FORWARD_ENABLED)
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        std::swap(npi[1][k-1],npi[1][k]);
        std::swap(npi[2][k-1],npi[2][k]);
        for(unsigned j=k-1; j>=lbf; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j+1;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          for(unsigned l=2; l<m; l++) { 
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = mi[i-1][mu[i-1]];
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeFlowtime(k,imaxi,Catmaxi,C,npi,invpi);
            unsigned tb = tiebreak(job,j);
            if (Fi < FC || (Fi == FC && tb<bt)) {
              FC = Fi;
              bp = j;
              bl = l;
              bf = true;
              bt = tb;
            }
            std::swap(npi[l+1][j],npi[l+1][j+1]);
          }
          if (j>1) {
            std::swap(npi[1][j-1],npi[1][j]);
            std::swap(npi[2][j-1],npi[2][j]);
          }
        } 
      } 
#endif
      assert(1<=bp && bp<=k);
      for(unsigned i=1; i<=m; i++) {
        if (i<=bl) {
          for(unsigned p=k-1; p>=bp; p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][bp]=job;
        } else {
          for(unsigned p=k-1; p>=(bf?bp+1:bp-1); p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][(bf?bp+1:bp-1)]=job;
        }
      }
      if (FC>=ub) {
        r.n = k;
        return; 
      }
      assert(flowtimeConsistent(I,r,FC,k+1));
    } 
    assert(flowtimeConsistent(I,r,FC,n+1));
    r.ms = FC;
  }
  void constructnpNEWTT(NPFSSolution& r, unsigned k, const std::vector<Job>& pi, const std::vector<char>& first, Time ub = infinite_time, unsigned lbf = 1, unsigned lbb = 1) {
    assert(constructionPrecondition(I,k,pi,first));
    Time FC = infinite_time;
    for(unsigned j=0; j<=n; j++) {
      h[0][j] = 0;     
      t[m+1][j] = 0;   
    }
    for(unsigned i=0; i<=m; i++) {
      h[i][0] = 0;     
      for(unsigned j=k; j<=n+1; j++)
        t[i][j] = 0;   
    }
    boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
    unsigned bp; 
    unsigned bl; 
    bool     bf; 
    for(unsigned kk=1; kk<pi.size(); kk++, k++) { 
      const unsigned job = pi[kk];
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++) {
          npi[i][j]=r.pi[i][j];
          invpi[i][npi[i][j]]=j;
        }
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++)
          h[i][j] = std::max(h[i-1][invpi[i-1][npi[i][j]]],h[i][j-1])+p[r.pi[i][j]][i];
      for(unsigned i=m; i>=1; i--)
        for(unsigned j=k-1; j>=1; j--) {
          tup[i][j] = (t[i+1][invpi[i+1][npi[i][j]]]>=t[i][j+1]);
          t[i][j] = std::max(t[i+1][invpi[i+1][npi[i][j]]],t[i][j+1])+p[r.pi[i][j]][i];
        }
      for(unsigned i=1; i<m; i++) {
        mi[i][k]=k; 
        for(unsigned j=k-1; j>0; j--)
          mi[i][j] = std::min(mi[i][j+1],invpi[i+1][npi[i][j]]);
      }
      for(unsigned i=m; i>1; i--)
        for(unsigned j=1; j<=k; j++)
          Mi[i][j] = std::max(Mi[i][j-1],invpi[i-1][npi[i][j-1]]);
      bp = bl = 0; 
      bf = true;   
      FC = infinite_time;    
      Time Fk;
      Time Cmax;
      unsigned imax;
      std::vector<Time> Catmax(C), Cold(C), Catmaxi(C); 
      for(unsigned i=1; i<=m; i++)
        npi[i][k]=job;
      for(unsigned j=k; j>=1; j--) { 
        Cmax = 0;
        imax = 1;
        fill(C.begin(), C.end(), infinite_time);
        mu[1] = j;
        for(unsigned i=1; i<m; i++)
          mu[i+1] = mi[i][mu[i]];
        nu[m] = j;
        for(unsigned i=m; i>1; i--)
          nu[i-1] = Mi[i][nu[i]]+1;
        for(unsigned i=1; i<=m; i++) {
          Cold = C;
          Time ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
          if (ct>=Cmax) {
            Cmax = ct;
            Catmax = Cold;
            imax = i;
          }
        }
        Fk = computeTotalTardiness(k,imax,Catmax,C,npi,invpi);
        if (Fk < FC || (Fk == FC && !first[job])) {
          FC = Fk;
          bp = j;
          bl = m;
          bf = false;
        }
        if (j>1)
          for(unsigned i=1; i<=m; i++)
            std::swap(npi[i][j],npi[i][j-1]);
      } 
#define FLOW_BACKWARD_ENABLED
#if defined(FLOW_BACKWARD_ENABLED)
      lbb = std::max(2u,lbb);
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        for(unsigned i=3; i<=m; i++)
          std::swap(npi[i][k-1],npi[i][k]);
        for(unsigned j=k; j>=lbb; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          for(unsigned l=2; l<m; l++) {
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = std::min(mi[i-1][mu[i-1]],j-1);
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeTotalTardiness(k,imaxi,Catmaxi,C,npi,invpi);
            if (Fi < FC || (Fi == FC && ((!first[job] && j>=bp) || (first[job] && j<=bp)))) {  
              FC = Fi;
              bp = j;
              bl = l;
              bf = false;
            }
            std::swap(npi[l+1][j-1],npi[l+1][j]);
          }
          if (j>2) {
            for(unsigned i=1; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j]);
            for(unsigned i=3; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j-2]);
          }
        } 
      } 
#endif
#define FLOW_FORWARD_ENABLED
#if defined(FLOW_FORWARD_ENABLED)
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        std::swap(npi[1][k-1],npi[1][k]);
        std::swap(npi[2][k-1],npi[2][k]);
        for(unsigned j=k-1; j>=lbf; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j+1;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          for(unsigned l=2; l<m; l++) { 
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = mi[i-1][mu[i-1]];
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeTotalTardiness(k,imaxi,Catmaxi,C,npi,invpi);
            if (Fi < FC || (Fi == FC && ((!first[job] && j>bp) || (first[job] && j<bp)))) {
              FC = Fi;
              bp = j;
              bl = l;
              bf = true;
            }
            std::swap(npi[l+1][j],npi[l+1][j+1]);
          }
          if (j>1) {
            std::swap(npi[1][j-1],npi[1][j]);
            std::swap(npi[2][j-1],npi[2][j]);
          }
        } 
      } 
#endif
      assert(1<=bp && bp<=k);
      for(unsigned i=1; i<=m; i++) {
        if (i<=bl) {
          for(unsigned p=k-1; p>=bp; p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][bp]=job;
        } else {
          for(unsigned p=k-1; p>=(bf?bp+1:bp-1); p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][(bf?bp+1:bp-1)]=job;
        }
      }
      if (FC>=ub) {
        r.n = k;
        return; 
      }
    } 
    r.ms = FC;
  }
  template<typename Tiebreak>
  void constructnpNEWTT(NPFSSolution& r, unsigned k, const std::vector<Job>& pi, Tiebreak tiebreak, Time ub = infinite_time, unsigned lbf = 1, unsigned lbb = 1) {
    Time FC = infinite_time;
    for(unsigned j=0; j<=n; j++) {
      h[0][j] = 0;     
      t[m+1][j] = 0;   
    }
    for(unsigned i=0; i<=m; i++) {
      h[i][0] = 0;     
      for(unsigned j=k; j<=n+1; j++)
        t[i][j] = 0;   
    }
    boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
    unsigned bp; 
    unsigned bt; 
    unsigned bl; 
    bool     bf; 
    for(unsigned kk=1; kk<pi.size(); kk++, k++) { 
      const unsigned job = pi[kk];
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++) {
          npi[i][j]=r.pi[i][j];
          invpi[i][npi[i][j]]=j;
        }
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++)
          h[i][j] = std::max(h[i-1][invpi[i-1][npi[i][j]]],h[i][j-1])+p[r.pi[i][j]][i];
      for(unsigned i=m; i>=1; i--)
        for(unsigned j=k-1; j>=1; j--) {
          tup[i][j] = (t[i+1][invpi[i+1][npi[i][j]]]>=t[i][j+1]);
          t[i][j] = std::max(t[i+1][invpi[i+1][npi[i][j]]],t[i][j+1])+p[r.pi[i][j]][i];
        }
      for(unsigned i=1; i<m; i++) {
        mi[i][k]=k; 
        for(unsigned j=k-1; j>0; j--)
          mi[i][j] = std::min(mi[i][j+1],invpi[i+1][npi[i][j]]);
      }
      for(unsigned i=m; i>1; i--)
        for(unsigned j=1; j<=k; j++)
          Mi[i][j] = std::max(Mi[i][j-1],invpi[i-1][npi[i][j-1]]);
      bp = bl = 0; 
      bf = true;   
      bt = std::numeric_limits<unsigned>::max();
      FC = infinite_time;    
      Time Fk;
      Time Cmax;
      unsigned imax;
      std::vector<Time> Catmax(C), Cold(C), Catmaxi(C); 
      for(unsigned i=1; i<=m; i++)
        npi[i][k]=job;
      for(unsigned j=k; j>=1; j--) { 
        Cmax = 0;
        imax = 1;
        fill(C.begin(), C.end(), infinite_time);
        mu[1] = j;
        for(unsigned i=1; i<m; i++)
          mu[i+1] = mi[i][mu[i]];
        nu[m] = j;
        for(unsigned i=m; i>1; i--)
          nu[i-1] = Mi[i][nu[i]]+1;
        for(unsigned i=1; i<=m; i++) {
          Cold = C;
          Time ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
          if (ct>=Cmax) {
            Cmax = ct;
            Catmax = Cold;
            imax = i;
          }
        }
        Fk = computeTotalTardiness(k,imax,Catmax,C,npi,invpi);
        unsigned tb = tiebreak(job,j);
        if (Fk < FC || (Fk == FC && tb<bt)) {
          FC = Fk;
          bp = j;
          bl = m;
          bf = false;
          bt = tb;
        }
        if (j>1)
          for(unsigned i=1; i<=m; i++)
            std::swap(npi[i][j],npi[i][j-1]);
      } 
#define FLOW_BACKWARD_ENABLED
#if defined(FLOW_BACKWARD_ENABLED)
      lbb = std::max(2u,lbb);
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        for(unsigned i=3; i<=m; i++)
          std::swap(npi[i][k-1],npi[i][k]);
        for(unsigned j=k; j>=lbb; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          for(unsigned l=2; l<m; l++) {
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = std::min(mi[i-1][mu[i-1]],j-1);
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeTotalTardiness(k,imaxi,Catmaxi,C,npi,invpi);
            unsigned tb = tiebreak(job,j);
            if (Fi < FC || (Fi == FC && tb<bt)) {
              FC = Fi;
              bp = j;
              bl = l;
              bf = false;
              bt = tb;
            }
            std::swap(npi[l+1][j-1],npi[l+1][j]);
          }
          if (j>2) {
            for(unsigned i=1; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j]);
            for(unsigned i=3; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j-2]);
          }
        } 
      } 
#endif
#define FLOW_FORWARD_ENABLED
#if defined(FLOW_FORWARD_ENABLED)
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        std::swap(npi[1][k-1],npi[1][k]);
        std::swap(npi[2][k-1],npi[2][k]);
        for(unsigned j=k-1; j>=lbf; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j+1;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          for(unsigned l=2; l<m; l++) { 
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = mi[i-1][mu[i-1]];
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeTotalTardiness(k,imaxi,Catmaxi,C,npi,invpi);
            unsigned tb = tiebreak(job,j);
            if (Fi < FC || (Fi == FC && tb<bt)) {
              FC = Fi;
              bp = j;
              bl = l;
              bf = true;
              bt = tb;
            }
            std::swap(npi[l+1][j],npi[l+1][j+1]);
          }
          if (j>1) {
            std::swap(npi[1][j-1],npi[1][j]);
            std::swap(npi[2][j-1],npi[2][j]);
          }
        } 
      } 
#endif
      assert(1<=bp && bp<=k);
      for(unsigned i=1; i<=m; i++) {
        if (i<=bl) {
          for(unsigned p=k-1; p>=bp; p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][bp]=job;
        } else {
          for(unsigned p=k-1; p>=(bf?bp+1:bp-1); p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][(bf?bp+1:bp-1)]=job;
        }
      }
      if (FC>=ub) {
        r.n = k;
        return; 
      }
    } 
    r.ms = FC;
  }
  void constructnpNEWTT_IT1(NPFSSolution& r, unsigned k, const std::vector<Job>& pi, Time ub = infinite_time, unsigned lbf = 1, unsigned lbb = 1) {
    Time FC = infinite_time;
    for(unsigned j=0; j<=n; j++) {
      h[0][j] = 0;     
      t[m+1][j] = 0;   
    }
    for(unsigned i=0; i<=m; i++) {
      h[i][0] = 0;     
      for(unsigned j=k; j<=n+1; j++)
        t[i][j] = 0;   
    }
    boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
    unsigned bp; 
    unsigned bt; 
    unsigned bl; 
    bool     bf; 
    for(unsigned kk=1; kk<pi.size(); kk++, k++) { 
      const unsigned job = pi[kk];
      std::vector<unsigned> tpt(I.m+1); 
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++) {
          npi[i][j]=r.pi[i][j];
          invpi[i][npi[i][j]]=j;
          tpt[i]+=I.p[npi[i][j]][i];
        }
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++)
          h[i][j] = std::max(h[i-1][invpi[i-1][npi[i][j]]],h[i][j-1])+p[r.pi[i][j]][i];
      for(unsigned i=m; i>=1; i--)
        for(unsigned j=k-1; j>=1; j--) {
          tup[i][j] = (t[i+1][invpi[i+1][npi[i][j]]]>=t[i][j+1]);
          t[i][j] = std::max(t[i+1][invpi[i+1][npi[i][j]]],t[i][j+1])+p[r.pi[i][j]][i];
        }
      for(unsigned i=1; i<m; i++) {
        mi[i][k]=k; 
        for(unsigned j=k-1; j>0; j--)
          mi[i][j] = std::min(mi[i][j+1],invpi[i+1][npi[i][j]]);
      }
      for(unsigned i=m; i>1; i--)
        for(unsigned j=1; j<=k; j++)
          Mi[i][j] = std::max(Mi[i][j-1],invpi[i-1][npi[i][j-1]]);
      bp = bl = 0; 
      bf = true;   
      bt = std::numeric_limits<unsigned>::max();
      FC = infinite_time;    
      Time Fk;
      Time Cmax;
      unsigned imax;
      std::vector<Time> Catmax(C), Cold(C), Catmaxi(C); 
      for(unsigned i=1; i<=m; i++) {
        npi[i][k]=job;
        tpt[i]+=I.p[job][i];
      }
      for(unsigned j=k; j>=1; j--) { 
        Cmax = 0;
        imax = 1;
        fill(C.begin(), C.end(), infinite_time);
        mu[1] = j;
        for(unsigned i=1; i<m; i++)
          mu[i+1] = mi[i][mu[i]];
        nu[m] = j;
        for(unsigned i=m; i>1; i--)
          nu[i-1] = Mi[i][nu[i]]+1;
        unsigned tb = 0; 
        for(unsigned i=1; i<=m; i++) {
          Cold = C;
          Time ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
          assert(ct >= tpt[i]);
          tb += ct-tpt[i];
          if (ct>=Cmax) {
            Cmax = ct;
            Catmax = Cold;
            imax = i;
          }
        }
        Fk = computeTotalTardiness(k,imax,Catmax,C,npi,invpi);
        if (Fk < FC || (Fk == FC && tb<bt)) {
          FC = Fk;
          bp = j;
          bl = m;
          bf = false;
          bt = tb;
        }
        if (j>1)
          for(unsigned i=1; i<=m; i++)
            std::swap(npi[i][j],npi[i][j-1]);
      } 
#define FLOW_BACKWARD_ENABLED
#if defined(FLOW_BACKWARD_ENABLED)
      lbb = std::max(2u,lbb);
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j];
          npi[i][k]=job;
        }
        for(unsigned i=3; i<=m; i++)
          std::swap(npi[i][k-1],npi[i][k]);
        for(unsigned j=k; j>=lbb; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          unsigned tb = 0; 
          assert(Cmax >= tpt[1]);
          tb+=Cmax-tpt[1];
          for(unsigned l=2; l<m; l++) {
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            assert(ct >= tpt[l]);
            tb += ct-tpt[l];
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            unsigned tbb = 0;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = std::min(mi[i-1][mu[i-1]],j-1);
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              assert(ct >= tpt[i]);
              tbb += ct-tpt[i];
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeTotalTardiness(k,imaxi,Catmaxi,C,npi,invpi);
            if (Fi < FC || (Fi == FC && (tb+tbb)<bt)) {
              FC = Fi;
              bp = j;
              bl = l;
              bf = false;
              bt = tb;
            }
            std::swap(npi[l+1][j-1],npi[l+1][j]);
          }
          if (j>2) {
            for(unsigned i=1; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j]);
            for(unsigned i=3; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j-2]);
          }
        } 
      } 
#endif
#define FLOW_FORWARD_ENABLED
#if defined(FLOW_FORWARD_ENABLED)
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j];
          npi[i][k]=job;
        }
        std::swap(npi[1][k-1],npi[1][k]);
        std::swap(npi[2][k-1],npi[2][k]);
        for(unsigned j=k-1; j>=lbf; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j+1;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          unsigned tb = 0; 
          assert(Cmax >= tpt[1]);
          tb+=Cmax-tpt[1];
          for(unsigned l=2; l<m; l++) { 
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            assert(ct >= tpt[l]);
            tb += ct-tpt[l];
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            unsigned tbb = 0;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = mi[i-1][mu[i-1]];
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              assert(ct >= tpt[i]);
              tbb += ct-tpt[i];
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeTotalTardiness(k,imaxi,Catmaxi,C,npi,invpi);
            if (Fi < FC || (Fi == FC && (tb+tbb)<bt)) {
              FC = Fi;
              bp = j;
              bl = l;
              bf = true;
              bt = tb;
            }
            std::swap(npi[l+1][j],npi[l+1][j+1]);
          }
          if (j>1) {
            std::swap(npi[1][j-1],npi[1][j]);
            std::swap(npi[2][j-1],npi[2][j]);
          }
        } 
      } 
#endif
      assert(1<=bp && bp<=k);
      for(unsigned i=1; i<=m; i++) {
        if (i<=bl) {
          for(unsigned p=k-1; p>=bp; p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][bp]=job;
        } else {
          for(unsigned p=k-1; p>=(bf?bp+1:bp-1); p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][(bf?bp+1:bp-1)]=job;
        }
      }
      if (FC>=ub) {
        r.n = k;
        return; 
      }
    } 
    r.ms = FC;
  }
  void constructnpNEWTT_IT2(NPFSSolution& r, unsigned k, const std::vector<Job>& pi, Time ub = infinite_time, unsigned lbf = 1, unsigned lbb = 1) {
    Time FC = infinite_time;
    for(unsigned j=0; j<=n; j++) {
      h[0][j] = 0;     
      t[m+1][j] = 0;   
    }
    for(unsigned i=0; i<=m; i++) {
      h[i][0] = 0;     
      for(unsigned j=k; j<=n+1; j++)
        t[i][j] = 0;   
    }
    boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
    unsigned bp; 
    unsigned bt; 
    unsigned bl; 
    bool     bf; 
    for(unsigned kk=1; kk<pi.size(); kk++, k++) { 
      const unsigned job = pi[kk];
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++) {
          npi[i][j]=r.pi[i][j];
          invpi[i][npi[i][j]]=j;
        }
      for(unsigned i=1; i<=m; i++)
        for(unsigned j=1; j<k; j++)
          h[i][j] = std::max(h[i-1][invpi[i-1][npi[i][j]]],h[i][j-1])+p[r.pi[i][j]][i];
      for(unsigned i=m; i>=1; i--)
        for(unsigned j=k-1; j>=1; j--) {
          tup[i][j] = (t[i+1][invpi[i+1][npi[i][j]]]>=t[i][j+1]);
          t[i][j] = std::max(t[i+1][invpi[i+1][npi[i][j]]],t[i][j+1])+p[r.pi[i][j]][i];
        }
      for(unsigned i=1; i<m; i++) {
        mi[i][k]=k; 
        for(unsigned j=k-1; j>0; j--)
          mi[i][j] = std::min(mi[i][j+1],invpi[i+1][npi[i][j]]);
      }
      for(unsigned i=m; i>1; i--)
        for(unsigned j=1; j<=k; j++)
          Mi[i][j] = std::max(Mi[i][j-1],invpi[i-1][npi[i][j-1]]);
      bp = bl = 0; 
      bf = true;   
      bt = std::numeric_limits<unsigned>::max();
      FC = infinite_time;    
      Time Fk;
      Time Cmax;
      unsigned imax;
      std::vector<Time> Catmax(C), Cold(C), Catmaxi(C); 
      for(unsigned i=1; i<=m; i++) {
        npi[i][k]=job;
      }
      for(unsigned j=k; j>=1; j--) { 
        Cmax = 0;
        imax = 1;
        fill(C.begin(), C.end(), infinite_time);
        mu[1] = j;
        for(unsigned i=1; i<m; i++)
          mu[i+1] = mi[i][mu[i]];
        nu[m] = j;
        for(unsigned i=m; i>1; i--)
          nu[i-1] = Mi[i][nu[i]]+1;
        for(unsigned i=1; i<=m; i++) {
          Cold = C;
          Time ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
          if (ct>=Cmax) {
            Cmax = ct;
            Catmax = Cold;
            imax = i;
          }
        }
        Fk = computeTotalTardiness(k,imax,Catmax,C,npi,invpi);
        unsigned tb = computeIT2(k,npi);
        if (Fk < FC || (Fk == FC && tb<bt)) {
          FC = Fk;
          bp = j;
          bl = m;
          bf = false;
          bt = tb;
        }
        if (j>1)
          for(unsigned i=1; i<=m; i++)
            std::swap(npi[i][j],npi[i][j-1]);
      } 
#define FLOW_BACKWARD_ENABLED
#if defined(FLOW_BACKWARD_ENABLED)
      lbb = std::max(2u,lbb);
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j]; 
          npi[i][k]=job;
        }
        for(unsigned i=3; i<=m; i++)
          std::swap(npi[i][k-1],npi[i][k]);
        for(unsigned j=k; j>=lbb; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          for(unsigned l=2; l<m; l++) {
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = std::min(mi[i-1][mu[i-1]],j-1);
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeTotalTardiness(k,imaxi,Catmaxi,C,npi,invpi);
            unsigned tb = computeIT2(k,npi);
            if (Fi < FC || (Fi == FC && tb<bt)) {
              FC = Fi;
              bp = j;
              bl = l;
              bf = false;
              bt = tb;
            }
            std::swap(npi[l+1][j-1],npi[l+1][j]);
          }
          if (j>2) {
            for(unsigned i=1; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j]);
            for(unsigned i=3; i<=m; i++)
              std::swap(npi[i][j-1],npi[i][j-2]);
          }
        } 
      } 
#endif
#define FLOW_FORWARD_ENABLED
#if defined(FLOW_FORWARD_ENABLED)
      if (k>1 && m>2) {
        for(unsigned i=1; i<=m; i++) {
          for(unsigned j=1; j<k; j++)
            npi[i][j]=r.pi[i][j];
          npi[i][k]=job;
        }
        std::swap(npi[1][k-1],npi[1][k]);
        std::swap(npi[2][k-1],npi[2][k]);
        for(unsigned j=k-1; j>=lbf; j--) { 
          std::vector<Time> Cl(n+1,infinite_time); 
          mu[1] = j;
          nu[m] = j+1;
          for(unsigned i=m; i>1; i--)
            nu[i-1] = Mi[i][nu[i]]+1;
          Catmax = Cl;
          Cmax = propagateNEW(Cl,npi,invpi,1,mu[1],nu[1]);
          imax = 1;
          for(unsigned l=2; l<m; l++) { 
            mu[l] = mi[l-1][mu[l-1]];
            Cold = Cl;
            Time ct = propagateNEW(Cl,npi,invpi,l,mu[l],nu[l]);
            if (ct >= Cmax) {
              Cmax = ct;
              Catmax.swap(Cold);
              imax = l;
            }
            C = Cl;
            Time Cmaxi = Cmax;
            unsigned imaxi = imax;
            Catmaxi = Catmax;
            for(unsigned i=l+1; i<=m; i++) {
              mu[i] = mi[i-1][mu[i-1]];
              Cold = C;
              ct = propagateNEW(C,npi,invpi,i,mu[i],nu[i]);
              if (ct >= Cmaxi) {
                Cmaxi = ct;
                Catmaxi.swap(Cold);
                imaxi = i;
              }
            }
            Time Fi = computeTotalTardiness(k,imaxi,Catmaxi,C,npi,invpi);
            unsigned tb = computeIT2(k,npi);
            if (Fi < FC || (Fi == FC && tb<bt)) {
              FC = Fi;
              bp = j;
              bl = l;
              bf = true;
              bt = tb;
            }
            std::swap(npi[l+1][j],npi[l+1][j+1]);
          }
          if (j>1) {
            std::swap(npi[1][j-1],npi[1][j]);
            std::swap(npi[2][j-1],npi[2][j]);
          }
        } 
      } 
#endif
      assert(1<=bp && bp<=k);
      for(unsigned i=1; i<=m; i++) {
        if (i<=bl) {
          for(unsigned p=k-1; p>=bp; p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][bp]=job;
        } else {
          for(unsigned p=k-1; p>=(bf?bp+1:bp-1); p--)
            r.pi[i][p+1]=r.pi[i][p];
          r.pi[i][(bf?bp+1:bp-1)]=job;
        }
      }
      if (FC>=ub) {
        r.n = k;
        return; 
      }
    } 
    r.ms = FC;
  }
  NPFSSolution BR_R(const float p, Time ub = infinite_time, unsigned lbf = 1, unsigned lbb = 1) {
    assert(p >= 0.0 && p <= 1.0);
    std::vector<Job> queue = getTotalTimeOrder(I);
    PFSSolution s(I);
    s.pi[1] = queue[1];
    builder construct(I);
    unsigned jj = ceil(I.n-(p*I.n));
    for(unsigned j = 2; j < jj; j++) {
      Job job = queue[j];
      std::vector<Job> insert = { 0, job };
      construct.constructp(I,s,j,1,j,insert,firstposition);
      if(j > 4) {
        for(unsigned k = 1; k < j; k+=2) {
          Job job1 = s.pi[k];
          Job job2 = s.pi[k+1];
          s.removeJob(job1);
          s.removeJob(job2);
          std::vector<Job> insert = { 0, job1, job2 };
          construct.constructp(I,s,j-1,1,j-1,insert,firstposition);
        }
        for(unsigned k = 2; k < j-1; k+=2) {
          Job job1 = s.pi[k];
          Job job2 = s.pi[k+1];
          s.removeJob(job1);
          s.removeJob(job2);
          std::vector<Job> insert = { 0, job1, job2 };
          construct.constructp(I,s,j-1,1,j-1,insert,firstposition);
        }
      }
    }
    NPFSSolution ns(I,s);
    if(jj < 2) jj = 2;
    for(unsigned j = jj; j <= I.n; j++) {
      Job job = queue[j];
      std::vector<Job> insert = { 0, job };
      constructnpNEW(ns,j,insert,firstposition);
      if(j > 4) {
        std::vector<char> first(I.n+1,true);
        for(unsigned k = 1; k < j; k+=2) {
          Job job1 = ns.pi[1][k];
          Job job2 = ns.pi[1][k+1];
          ns.removeJob(job1);
          ns.removeJob(job2);
          std::vector<Job> insert = { 0, job1, job2 };
          constructnpNEW(ns,j-1,insert,firstposition);
        }
        for(unsigned k = 2; k < j-1; k+=2) {
          Job job1 = ns.pi[1][k];
          Job job2 = ns.pi[1][k+1];
          ns.removeJob(job1);
          ns.removeJob(job2);
          std::vector<Job> insert = { 0, job1, job2 };
          constructnpNEW(ns,j-1,insert,firstposition);
        }
      }
    }
    return ns;
  }
};
inline double IT(const Instance& I, const std::vector<Time>& C, std::vector<Time>& Cj, unsigned j, unsigned k) {
  double result = 0.0;
  Cj[1] = C[1]+I.p[j][1];
  for(unsigned i=2; i<=I.m; i++) {
    if (C[i]>Cj[i-1]) {
      Cj[i] = C[i] + I.p[j][i];
    } else {
      Cj[i] = Cj[i-1] + I.p[j][i];
      result += double(I.m*(I.n-2))/double(i*(I.n-k-2)+k*I.m)*(Cj[i-1]-C[i]);
    }
  }
  return result;
}
inline double AT(const Instance& I, const std::vector<Time>& Csum, const std::vector<Time>& Cj, unsigned j, unsigned k) {
  double ct = Cj[1] + double(Csum[1] - I.p[j][1])/(I.n-k-1);
  for(unsigned i=2; i<=I.m; i++)
    ct = std::max(ct,double(Cj[i])) + double(Csum[i] - I.p[j][i])/(I.n-k-1);
  return double(Cj[I.m])+ct;
}
inline PFSSolution LR(const Instance& I, unsigned x) {
  const unsigned n = I.n;
  const unsigned m = I.m;
  assert(1<=x && x<=n);
  std::vector<Time> C(I.m+1,0), Cj(I.m+1), Csum(I.m+1);
  for(unsigned i=1; i<=m; i++)
    for(unsigned j=1; j<=n; j++)
      Csum[i] += I.p[j][i];
  std::vector<std::pair<Job,double> > seed(I.n+1);
  for(unsigned j=1; j<=n; j++) {
    double f = (n-2)*IT(I,C,Cj,j,0);
    f += AT(I,Csum,Cj,j,0);
    seed[j]=std::make_pair(j,f);
  }
  std::nth_element(seed.begin()+1,seed.begin()+x,seed.end(),
       [](const std::pair<Job,double>& v1, const std::pair<Job,double>& v2) { return v1.second < v2.second; });
  PFSSolution bs(I);
  bs.ms = std::numeric_limits<unsigned>::max();
  for(unsigned ji=1; ji<=x; ji++) {
    Job fj = seed[ji].first;
    std::vector<Job> U(I.n+1,0);
    std::iota(U.begin(),U.end(),0);
    std::swap(U[fj],U[1]);
    fill(C.begin(),C.end(),0);
    C[1] = I.p[fj][1];
    for(unsigned i=2; i<=m; i++)
      C[i] = std::max(C[i],C[i-1])+I.p[fj][i];
    Time Fsum = C[m];
    fill(Csum.begin(),Csum.end(),0);
    for(unsigned i=1; i<=m; i++) {
      for(unsigned j=1; j<=n; j++)
        Csum[i] += I.p[j][i];
      Csum[i] -= I.p[fj][i];
    }
    for(unsigned k=1; k+1<n; k++) {
      double fm = std::numeric_limits<double>::max(), f0m = std::numeric_limits<double>::max();
      unsigned jm = 0;
      std::vector<Time> Cjm(I.m+1);
      for(unsigned nji=k+1; nji <= n; nji++) {
        double f0 = IT(I,C,Cj,U[nji],k);
        double f = (n-k-2)*f0+AT(I,Csum,Cj,U[nji],k);
        if (f < fm || ( f == fm && f0 < f0m) ) {
          fm = f;
          f0m = f0;
          jm = nji;
          Cjm.swap(Cj);
        }
      }
      assert(jm>0);
      std::swap(U[k+1],U[jm]);
      for(unsigned i=1; i<=m; i++)
        Csum[i] -= I.p[U[k+1]][i];
      C.swap(Cjm);
      Fsum += C[m];
    }
    Time ct = C[1]+I.p[U[n]][1];
    for(unsigned i=2; i<=I.m; i++)
      ct = std::max(ct,C[i]) + I.p[U[n]][i];
    Fsum += ct;
    if (Fsum<bs.ms) {
      bs.pi = U;
      bs.ms = Fsum;
    }
  }
  return bs;
}
template <typename Rng>
inline PFSSolution randP(const Instance& I, Rng& rng, unsigned rep = 1000, Time ub = infinite_time) {
  PFSSolution b(I);
  b.computeMakespan(I);
  while (rep-->0) {
    PFSSolution n(I);
    std::shuffle(n.pi.begin()+1,n.pi.end(),rng);
    n.computeMakespan(I);
    if (n.ms<b.ms)
      b = n;
  }
  return b;
}
