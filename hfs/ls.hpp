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
#include <utility>
#include <chrono>
#include <random>
#include <numeric>
#include "constructive.hpp"
#include "parameters.hpp"
template <typename Tiebreak>
inline void insertion(const Instance& I, PFSSolution &n, builder& construct, Tiebreak tiebreak, Params& params) {
  unsigned kl = 1;
  unsigned ku = I.n+1;
  unsigned n_ls = params.n_ls;
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned cms = n.ms;
    std::vector<Job> pi = n.pi;
    for (unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = pi[ji];
      unsigned pos_current_job = std::find(n.pi.begin(),n.pi.end(),j)-n.pi.begin();
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      if (params.obj == 'm') {
        unsigned minPos, maxPos;
        if (ji<ku-1) {
          unsigned new_pos_last_job = std::find(n.pi.begin(),n.pi.end(),pi[ji+1])-n.pi.begin();
          minPos = std::min(pos_current_job,new_pos_last_job);
          maxPos = std::max(pos_current_job,new_pos_last_job);
          if (maxPos == I.n) 
            maxPos--;
          assert(maxPos < I.n);
        } else {
          minPos = 1;
          maxPos = 0;
        }
        construct.constructp(I,n,I.n,kl,ku-1,remove,tiebreak,infinite_time,minPos,maxPos);
      }
      else if (params.obj == 'f')
        construct.constructpCsum(I,n,I.n,kl,ku-1,remove,tiebreak,infinite_time);
      else
        construct.constructpTT(I,n,I.n,kl,ku-1,remove,tiebreak,infinite_time);
    }
    if (n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
template<typename Tiebreak>
inline void lsps_i(const Instance& I, PFSSolution& n, builder& construct, const char obj, const unsigned kl, const unsigned ku, const unsigned dc, unsigned n_ls, const std::chrono::system_clock::time_point& start, const unsigned timelimit, Tiebreak tiebreak) {
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
    unsigned cms = n.ms;
    std::vector<Job> pi = n.pi;
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = pi[ji];
      unsigned pos_current_job = std::find(n.pi.begin(),n.pi.end(),j)-n.pi.begin();
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      if (obj == 'm') {
        unsigned minPos, maxPos;
        if (ji<ku-1) {
          unsigned new_pos_last_job = std::find(n.pi.begin(),n.pi.end(),pi[ji+1])-n.pi.begin();
          minPos = std::min(pos_current_job,new_pos_last_job);
          maxPos = std::max(pos_current_job,new_pos_last_job);
          if (maxPos == I.n-dc) 
            maxPos--;
          assert(maxPos < I.n-dc);
        } else {
          minPos = 1;
          maxPos = 0;
        }
        construct.constructp(I,n,I.n-dc,kl,ku-1,remove,tiebreak,infinite_time,minPos,maxPos);
      }
      else if (obj == 'f')
        construct.constructpCsum(I,n,I.n-dc,kl,ku-1,remove,tiebreak,infinite_time);
      else
        construct.constructpTT(I,n,I.n-dc,kl,ku-1,remove,tiebreak,infinite_time);
    }
    if (n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void insertionFF(const Instance& I, PFSSolution &n, builder& construct, Params& params) {
  unsigned kl = 1;
  unsigned ku = I.n+1;
  unsigned n_ls = params.n_ls;
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned cms = n.ms;
    std::vector<Job> pi = n.pi;
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = pi[ji];
      unsigned pos_current_job = std::find(n.pi.begin(),n.pi.end(),j)-n.pi.begin();
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      unsigned minPos, maxPos;
      if (ji<ku-1) {
        unsigned new_pos_last_job = std::find(n.pi.begin(),n.pi.end(),pi[ji+1])-n.pi.begin();
        minPos = std::min(pos_current_job,new_pos_last_job);
        maxPos = std::max(pos_current_job,new_pos_last_job);
        if (maxPos == I.n) 
          maxPos--;
        assert(maxPos < I.n);
      } else {
        minPos = 1;
        maxPos = 0;
      }
      construct.constructpFF(I,n,I.n,kl,ku-1,remove,infinite_time,minPos,maxPos);
    }
    if (n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void lsps_iFF(const Instance& I, PFSSolution& n, builder& construct, const char obj, const unsigned kl, const unsigned ku, const unsigned dc, unsigned n_ls, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
    unsigned cms = n.ms;
    std::vector<Job> pi = n.pi;
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = pi[ji];
      unsigned pos_current_job = std::find(n.pi.begin(),n.pi.end(),j)-n.pi.begin();
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      unsigned minPos, maxPos;
      if (ji<ku-1) {
        unsigned new_pos_last_job = std::find(n.pi.begin(),n.pi.end(),pi[ji+1])-n.pi.begin();
        minPos = std::min(pos_current_job,new_pos_last_job);
        maxPos = std::max(pos_current_job,new_pos_last_job);
        if (maxPos == I.n-dc) 
          maxPos--;
        assert(maxPos < I.n-dc);
      } else {
        minPos = 1;
        maxPos = 0;
      }
      construct.constructpFF(I,n,I.n-dc,kl,ku-1,remove,infinite_time,minPos,maxPos);
    }
    if (n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void N5(const Instance& I, PFSSolution& n, Params& params) {
  unsigned n_ls = params.n_ls;
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    boost::multi_array<Time,2> C;
    C.resize(boost::extents[I.m+1][I.n+1]);
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = 1; j <= I.n; j++) {
        C[i][n.pi[j]] = std::max(C[i-1][n.pi[j]],C[i][n.pi[j-1]]) + I.p[n.pi[j]][i];
      }
    }
    boost::multi_array<bool,2> c;
    c.resize(boost::extents[I.m+1][I.n+1]);
    c[I.m][n.pi[I.n]] = true;
    for(unsigned i = I.m; i >= 1; i--) {
      for(unsigned j = I.n; j >= 1; j--) {
        unsigned tj = n.pi[j];
        if (c[i][tj] && C[i][tj]-I.p[tj][i] == C[i][n.pi[j-1]])
          c[i][n.pi[j-1]]=true;
        if (c[i][tj] && C[i][tj]-I.p[tj][i] == C[i-1][tj])
          c[i-1][tj]=true;
      }
    }
    std::vector<std::pair<unsigned,unsigned>> r;
    c[1][0] = true; 
    unsigned col = 1;
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = col; j < I.n-1; j++) {
        if (c[i][n.pi[j]] && c[i][n.pi[j+1]] && !c[i][n.pi[j-1]]) { 
          r.push_back(std::make_pair(j, j+1));
        }
        else if (c[i][n.pi[j]] && c[i][n.pi[j+1]] && !c[i][n.pi[j+2]]) { 
          r.push_back(std::make_pair(j, j+1));
          col = j;
          break;
        }
      }
    }
    PFSSolution bs = n;  
    Time bm = n.ms;      
    Time pm = n.ms;      
    for (unsigned k = 0; k < r.size(); k++) {
      assert(r[k].first != r[k].second);
      std::swap(n.pi[r[k].first], n.pi[r[k].second]);
      n.computeMakespan(I);
      if (n.ms < bm) {
        bm = n.ms;
        bs = n;
      }
      std::swap(n.pi[r[k].first], n.pi[r[k].second]);
    }
    if (bm < pm) {
      n = bs;
      improved = true;
    }
    else {
      improved = false;
    } 
  }
}
template <typename Tiebreak>
inline void Pc(const Instance& I, PFSSolution& n, builder& construct, Tiebreak tiebreak, Params& params) {
  unsigned n_ls = params.n_ls;
  boost::multi_array<Time,2> h(boost::extents[I.m+1][I.n+1]), t(boost::extents[I.m+2][I.n+2]);
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    for(unsigned i = 1; i <= I.m; i++)
      for(unsigned j = 1; j <= I.n; j++)
        h[i][j] = std::max(h[i-1][j], h[i][j-1]) + I.p[n.pi[j]][i];
    for(unsigned i = I.m; i >= 1; i--)
      for(unsigned j = I.n; j >= 1; j--)
        t[i][j] = std::max(t[i+1][j], t[i][j+1]) + I.p[n.pi[j]][i];
    std::vector<std::pair<unsigned,unsigned>> r;
    std::vector<char> duplicate(I.n+1,0);
    unsigned cmax = n.ms;
    for(unsigned i = 1; i <= I.m; i++)
      for(unsigned j = 1; j < I.n; j++)
        if(!duplicate[j] && (h[i][j] + t[i+1][j] == cmax || h[i][j+1] + t[i+1][j+1] == cmax)) {
          r.push_back(std::make_pair(j,j+1));
          duplicate[j] = 1;
        }
    PFSSolution bs = n;  
    Time bm = n.ms;      
    Time pm = n.ms;      
    for (unsigned k = 0; k < r.size(); k++) {
      assert(r[k].first != r[k].second);
      PFSSolution ns = n;
      auto nend = ns.pi.end();
      nend = std::remove(ns.pi.begin()+1, nend, r[k].first);
      nend = std::remove(ns.pi.begin()+1, nend, r[k].second);
      std::fill(nend, ns.pi.end(), 0);
      std::vector<Job> insert = { 0, r[k].first, r[k].second };
      construct.constructp(I,ns,I.n+1-2,1,I.n+1-2,insert,tiebreak);
      if (ns.ms < bm) {
        bm = ns.ms;
        bs = ns;
      }
    }
    if (bm < pm) {
      n = bs;
      improved = true;
    }
    else {
      improved = false;
    }
  }
}
inline void fpe(const Instance& I, PFSSolution &s, const std::chrono::system_clock::time_point& start, unsigned timelimit, unsigned x, bool iterate=true) {
  PFSSolution cs(s),ns(s);
  do {
    cs = ns;
    for(unsigned ji=1; ji<=I.n; ji++) {
      unsigned j = cs.pi[ji];
      for(unsigned d=1; d<=x; d++) {
	      PFSSolution ts(ns);
	      unsigned jp = std::find(ts.pi.begin()+1,ts.pi.end(),j)-ts.pi.begin();
	      if (jp+d>I.n)
	        break;
	      std::swap(ts.pi[jp],ts.pi[jp+d]);
	      ts.computeFlowtime(I);
	      if (ts.ms<ns.ms)
	        ns = ts;
      }
    }
  } while (ns.ms < cs.ms && iterate && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
  s = ns;
}
inline void bpe(const Instance& I, PFSSolution &s, const std::chrono::system_clock::time_point& start, unsigned timelimit, unsigned x, bool iterate=true) {
  PFSSolution cs(s),ns(s);
  do {
    cs = ns;
    for(unsigned ji=I.n; ji>=1; ji--) {
      unsigned j = cs.pi[ji];
      for(unsigned d=1; d<=x; d++) {
        PFSSolution ts(ns);
        unsigned jp = std::find(ts.pi.begin()+1,ts.pi.end(),j)-ts.pi.begin();
        if (jp<d+1)
          break;
        std::swap(ts.pi[jp-d],ts.pi[jp]);
        ts.computeFlowtime(I);
        if (ts.ms<ns.ms) {
          ns = ts;
          break;
        }
      }
    }
  } while (ns.ms < cs.ms && iterate && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
  s = ns;
}
inline void swapTasgetiren(const Instance& I, PFSSolution &s, const std::chrono::system_clock::time_point& start, unsigned timelimit) {
  std::array<Time,64> Cp;
  Time ct, pms;
  const Time *ctp;
  PFSSolution cs(s);
  unsigned i = 1;
  std::fill(Cp.begin()+1,Cp.end(),0);
  pms = 0;
  do {
    unsigned j = i+1;
    do {
      PFSSolution ns(cs);
      std::swap(ns.pi[i],ns.pi[j]);
      ns.ms = pms;
      ns.computeFlowtime0(I,Cp,i,I.n+1,cs.ms);
      if (ns.ms<cs.ms) {
	      cs = ns;
	      i = 1;
	      j = i+1;
	      std::fill(Cp.begin()+1,Cp.end(),0);
	      pms = 0;
      }
      else j++;
    } while (j<=I.n);
    ct = 0;
    ctp = &I.p[cs.pi[i]][1];
    for(unsigned ii=I.m; ii>0; ii--, ctp++)
      ct = Cp[ii] = std::max(Cp[ii],ct)+*ctp;
    pms += ct;
    i++;
  } while (i<I.n && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
  s = cs;
}
inline std::pair<unsigned,unsigned> bestInsert(const Instance& I, PFSSolution ns, unsigned jp) {
  std::array<Time,64> Cp;
  Time ct, pms;
  const Time *ctp;
  std::pair<unsigned,unsigned> r = std::make_pair(0,ns.ms);
  PFSSolution ts(ns);
  std::rotate(ts.pi.begin()+1,ts.pi.begin()+jp,ts.pi.begin()+jp+1);
  std::fill(Cp.begin()+1,Cp.end(),0);
  ts.ms = 0;
  pms = 0;
  for(unsigned np=1; np<=I.n; np++) {
    if (np != jp) {
      ts.computeFlowtime0(I,Cp,np,I.n+1,r.second);
      if (ts.ms<r.second)
        r = std::make_pair(np,ts.ms);
    }
    if (np<I.n) {
      std::swap(ts.pi[np],ts.pi[np+1]);
      ct = 0;
      ctp = &I.p[ts.pi[np]][1];
      for(unsigned i=I.m; i>0; i--, ctp++)
        ct = Cp[i] = std::max(Cp[i],ct)+*ctp;
      pms += ct;
      ts.ms = pms;
    }
  }
  return r;
}
inline void swapExp(const Instance& I, PFSSolution &s, const std::chrono::system_clock::time_point& start, unsigned timelimit, unsigned min_d=1, unsigned max_d=0) {
  std::array<Time,64> Cp;
  Time ct, pms, cms;
  const Time *ctp;
  PFSSolution cs(s);
  unsigned d = min_d;
  unsigned touched = 3*I.n*I.n;
  if (max_d == 0)
    max_d = I.n;
  do {
    unsigned i = 1;
    std::fill(Cp.begin()+1,Cp.end(),0);
    pms = 0;
    cms = cs.ms;
    do {
      unsigned j = i+d;
      PFSSolution ns(cs);
      std::swap(ns.pi[i],ns.pi[j]);
      ns.ms = pms;
      ns.computeFlowtime0(I,Cp,i,I.n+1,cs.ms);
      touched--;
      if (ns.ms<cs.ms) {
	      cs = ns;
      }
      ct = 0;
      ctp = &I.p[cs.pi[i]][1];
      for(unsigned ii=I.m; ii>0; ii--, ctp++)
	      ct = Cp[ii] = std::max(Cp[ii],ct)+*ctp;
      pms += ct;
      i++;
    } while (i+d<=I.n && touched>0);
    if (cs.ms < cms)
      d = min_d;
    else
      d++;
  } while (d<max_d && touched>0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
  s = cs;
}
inline void insertTasgetiren(const Instance& I, PFSSolution &s, const std::chrono::system_clock::time_point& start, unsigned timelimit) {
  PFSSolution cs(s);
  unsigned i = 1;
  do {
    unsigned j = i+1;
    do {
      PFSSolution ns(cs);
      std::rotate(ns.pi.begin()+i,ns.pi.begin()+i+1,ns.pi.begin()+j+1);
      ns.computeFlowtime(I);
      if (ns.ms<cs.ms) {
        cs = ns;
        i = 1;
        j = i+1;
      }
      else j++;
    }
    while (j<=I.n);
    i++;
  }
  while (i<I.n && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
  s = cs;
}
inline void insertJarboui(const Instance& I, PFSSolution &s, const std::chrono::system_clock::time_point& start, unsigned timelimit) {
  PFSSolution cs(s);
  unsigned i = 1;
  do {
    unsigned j = 1;
    do {
      PFSSolution ns(cs);
      if (i==j) {
        j++;
        continue;
      }
      if (i<j)
        std::rotate(ns.pi.begin()+i,ns.pi.begin()+i+1,ns.pi.begin()+j+1);
      else
        std::rotate(ns.pi.begin()+j,ns.pi.begin()+i,ns.pi.begin()+i+1);
      ns.computeFlowtime(I);
      if (ns.ms<cs.ms) {
        cs = ns;
        i = i-1;  
        j = I.n + 1; 
      }
      else j++;
    }
    while (j<=I.n);
    i++;
  }
  while (i<I.n && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
  s = cs;
}
inline void lsTasgetiren(const Instance& I, PFSSolution &s, const std::chrono::system_clock::time_point& start, unsigned timelimit) {
  PFSSolution cs(s);
  do {
    s = cs;
    insertTasgetiren(I,cs,start,timelimit);
    swapTasgetiren(I,cs,start,timelimit);
  } while (cs.ms < s.ms && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
}
inline void shiftJob(PFSSolution& s, unsigned i, unsigned j) {
  if (i<j)
    std::rotate(s.pi.begin()+i,s.pi.begin()+i+1,s.pi.begin()+j+1);
  else if (i>j)
    std::rotate(s.pi.begin()+j,s.pi.begin()+i,s.pi.begin()+i+1);
}
inline std::pair<unsigned,unsigned> bestSwap(const Instance& I, PFSSolution ns) {
  std::array<Time,64> Cp;
  Time ct, pms;
  const Time *ctp;
  std::pair<unsigned,unsigned> r = std::make_pair(0,infinite_time);
  PFSSolution ts(ns);
  std::fill(Cp.begin(),Cp.end(),0);
  ts.ms = 0;
  pms = 0;
  for(unsigned np=1; np<I.n; np++) {
    ts.pi[np]=ns.pi[np+1];
    ts.pi[np+1]=ns.pi[np];
    ts.computeFlowtime0(I,Cp,np,I.n+1,r.second);
    if (ts.ms<r.second)
	    r = std::make_pair(np,ts.ms);
    if (np<I.n-1) {
      ts.pi[np]=ns.pi[np];
      ts.pi[np+1]=ns.pi[np+1];
      ct = 0;
      ctp = &I.p[ts.pi[np]][1];
      for(unsigned i=I.m; i>0; i--, ctp++)
	      ct = Cp[i] = std::max(Cp[i],ct)+*ctp;
      pms += ct;
      ts.ms = pms;
    }
  }
  return r;
}
inline void iRZ(const Instance& I, PFSSolution& s, const std::chrono::system_clock::time_point& start, unsigned timelimit, unsigned n_ls=3) {
  PFSSolution ns(s);
  unsigned touched = 0;
  do {
    s = ns;
    for(unsigned ji=1; ji<=I.n; ji++) {
      unsigned j = s.pi[ji];
      unsigned jp = std::find(ns.pi.begin()+1,ns.pi.end(),j)-ns.pi.begin();
      std::pair<unsigned,unsigned> r = bestInsert(I,ns,jp);
      if (r.second<ns.ms) {
        shiftJob(ns,jp,r.first);
        ns.ms = r.second;
        touched = 0;
      }
      if (touched>=I.n)
        break;
    }
  } while (ns.ms < s.ms && n_ls-->0 && touched<I.n && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
  s = ns; 
}
inline void riRZ(const Instance& I, PFSSolution& s, const std::chrono::system_clock::time_point& start, unsigned timelimit, unsigned n_ls=3) {
  PFSSolution ns(s);
  unsigned touched = 0;
  do {
    s = ns;
    for(unsigned ji=I.n; ji>=1; ji--) {
      unsigned j = s.pi[ji];
      unsigned jp = std::find(ns.pi.begin()+1,ns.pi.end(),j)-ns.pi.begin();
      std::pair<unsigned,unsigned> r = bestInsert(I,ns,jp);
      if (r.second<ns.ms) {
        shiftJob(ns,jp,r.first);
        ns.ms = r.second;
        touched = 0;
      }
      if (touched>=I.n)
        break;
    }
  } while (ns.ms < s.ms && n_ls-->0 && touched<I.n && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
  s = ns; 
}
inline void raiRZ(const Instance& I, PFSSolution& s, const std::chrono::system_clock::time_point& start, unsigned timelimit, std::mt19937& rng, unsigned n_ls=3) {
  PFSSolution ns(s);
  unsigned touched = 0;
  do {
    s = ns;
    std::vector<Job> order(I.n+1,0);
    std::iota(order.begin(),order.end(),0);
    std::shuffle(order.begin()+1,order.end(),rng);
    for(unsigned jii=I.n; jii>=1; jii--) { 
      unsigned ji=order[jii];
      unsigned j = s.pi[ji];
      unsigned jp = std::find(ns.pi.begin()+1,ns.pi.end(),j)-ns.pi.begin();
      std::pair<unsigned,unsigned> r = bestInsert(I,ns,jp);
      if (r.second<ns.ms) {
        shiftJob(ns,jp,r.first);
        ns.ms = r.second;
        touched = 0;
      }
      if (touched>=I.n)
        break;
    }
  } while (ns.ms < s.ms && n_ls-->0 && touched<I.n && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
  s = ns; 
}
inline void iSwapFirst(const Instance& I, PFSSolution& s, const std::chrono::system_clock::time_point& start, unsigned timelimit, unsigned n_ls=3) {
  std::array<Time,64> Cp;
  Time ct, pms;
  const Time *ctp;
  PFSSolution ts(s);
  std::fill(Cp.begin(),Cp.end(),0);
  ts.ms = 0;
  pms = 0;
  unsigned p = 1;
  n_ls *= (I.n-1); 
  while (n_ls-->0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
    ts.pi[p]=s.pi[p+1];
    ts.pi[p+1]=s.pi[p];
    ts.computeFlowtime0(I,Cp,p,I.n+1,s.ms);
    if (ts.ms<s.ms) {
      std::swap(s.pi[p],s.pi[p+1]);
      s.ms = ts.ms;
    }
    ts.pi[p]=s.pi[p];
    ts.pi[p+1]=s.pi[p+1];
    if (p+1<I.n) {
      ct = 0;
      ctp = &I.p[ts.pi[p]][1];
      for(unsigned i=I.m; i>0; i--, ctp++)
        ct = Cp[i] = std::max(Cp[i],ct)+*ctp;
      pms += ct;
      ts.ms = pms;
      p++;
    } else {
      std::fill(Cp.begin(),Cp.end(),0);
      pms = 0;
      ts.ms = 0;
      p = 1;
    }
  } 
}
inline void iSwapBest(const Instance& I, PFSSolution& s, const std::chrono::system_clock::time_point& start, unsigned timelimit, unsigned n_ls=3) {
  do {
    std::pair<unsigned,unsigned> r = bestSwap(I,s);
    if (r.second<s.ms) {
      std::swap(s.pi[r.first],s.pi[r.first+1]);
      s.ms = r.second;
    } else break;
  } while (n_ls-->0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit);
}
inline void viRZ(const Instance& I, PFSSolution &s, const std::chrono::system_clock::time_point& start, unsigned timelimit, unsigned n_insert=0) {
  PFSSolution ns(s);
  unsigned ji = 1;
  unsigned touched = 0;
  if (n_insert==0)
    n_insert = 3*I.n;
  while (touched<I.n && --n_insert>0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
    unsigned j = s.pi[ji];
    unsigned jp = std::find(ns.pi.begin()+1,ns.pi.end(),j)-ns.pi.begin();
    std::pair<unsigned,unsigned> r = bestInsert(I,ns,jp);
    if (r.second<ns.ms) {
      shiftJob(ns,jp,r.first);
      ns.ms = r.second;
      touched = 0;
    } else touched++;
    ji++;
    if (ji>I.n)
      ji = 1;
  }
  s = ns;
}
inline void powerSwap(const Instance& I, PFSSolution &s, const std::chrono::system_clock::time_point& start, unsigned timelimit) {
  for(unsigned si=1; si<=I.n && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit; si++) {
    for(unsigned sj=1; sj<=I.n; sj++) {
      if (sj == si)
        continue;
      PFSSolution ts(s);
      std::swap(ts.pi[si],ts.pi[sj]);
      ts.computeFlowtime(I);
      std::pair<unsigned,unsigned> r = bestInsert(I,ts,sj);
      if (r.second<ts.ms) {
        shiftJob(ts,sj,r.first);
        s = ts;
        s.ms = r.second;
      }
    }
  }
}
inline void insertionIT1(const Instance& I, PFSSolution &n, builder& construct, Params& params) {
  unsigned kl = 1;
  unsigned ku = I.n+1;
  unsigned n_ls = params.n_ls;
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned cms = n.ms;
    std::vector<Job> pi = n.pi;
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = pi[ji];
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };  
      construct.constructpTT_IT1(I,n,I.n,kl,ku-1,remove,infinite_time);
    }
    if (n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void insertionIT2(const Instance& I, PFSSolution &n, builder& construct, Params& params) {
  unsigned kl = 1;
  unsigned ku = I.n+1;
  unsigned n_ls = params.n_ls;
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned cms = n.ms;
    std::vector<Job> pi = n.pi;
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = pi[ji];
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };  
      construct.constructpTT_IT2(I,n,I.n,kl,ku-1,remove,infinite_time);
    }
    if (n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void lsps_iIT1(const Instance& I, PFSSolution& n, builder& construct, const char obj, const unsigned kl, const unsigned ku, const unsigned dc, unsigned n_ls, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
    unsigned cms = n.ms;
    std::vector<Job> pi = n.pi;
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = pi[ji];
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      construct.constructpTT_IT1(I,n,I.n-dc,kl,ku-1,remove,infinite_time); 
    }
    if (n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void lsps_iIT2(const Instance& I, PFSSolution& n, builder& construct, const char obj, const unsigned kl, const unsigned ku, const unsigned dc, unsigned n_ls, const std::chrono::system_clock::time_point& start, const unsigned timelimit) {
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < timelimit) {
    unsigned cms = n.ms;
    std::vector<Job> pi = n.pi;
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = pi[ji];
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      construct.constructpTT_IT2(I,n,I.n-dc,kl,ku-1,remove,infinite_time); 
    }
    if (n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void ASLS_TT(const Instance& I, PFSSolution& n, Params& params) {
  boost::multi_array<unsigned,2> C(boost::extents[I.m+1][I.n+1]); 
  std::vector<unsigned> pt(I.n+1);  
  for(unsigned j = 1; j <= I.n; j++) {
    unsigned job = n.pi[j];
    unsigned Cj = 0;
    for(unsigned i = 1; i <= I.m; i++)
      C[i][j] = Cj = std::max(Cj, C[i][j-1]) + I.p[job][i];
    pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[job]));
  }
  bool improved = true;
  unsigned n_ls = params.n_ls;
  std::vector<unsigned> cpi = n.pi;
  while(improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    improved = false;
    for(unsigned k = 1; k < I.n; k++) {
      std::swap(cpi[k], cpi[k+1]);
      std::vector<unsigned> Cm(I.m+1);
      for(unsigned i = 1; i <= I.m; i++)
        Cm[i] = std::max(Cm[i-1], C[i][k-1]) + I.p[cpi[k]][i];
      unsigned tt = pt[k-1] + std::max(0, int(Cm[I.m]-I.dd[cpi[k]]));
      for(unsigned j = k+1; j <= I.n; j++) {
        unsigned job = cpi[j];
        for(unsigned i = 1; i <= I.m; i++)
          Cm[i] = std::max(Cm[i-1], Cm[i]) + I.p[job][i];
        tt += std::max(0, int(Cm[I.m]-I.dd[job]));
      }
      if(tt < n.ms) { 
        n.pi = cpi;
        n.ms = tt;
        improved = true;
        for(unsigned j = k; j <= I.n; j++) {
          unsigned job = n.pi[j];
          unsigned Cj = 0;
          for(unsigned i = 1; i <= I.m; i++)
            C[i][j] = Cj = std::max(Cj, C[i][j-1]) + I.p[job][i];
          pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[job]));
        }
      }
      else 
        std::swap(cpi[k], cpi[k+1]);
    }
  }
}
inline void BISLS_TT(const Instance& I, PFSSolution& n, Params& params) {
  boost::multi_array<unsigned,2> C(boost::extents[I.m+1][I.n+1]); 
  std::vector<unsigned> pt(I.n+1);  
  for(unsigned j = 1; j <= I.n; j++) {
    unsigned job = n.pi[j];
    unsigned Cj = 0;
    for(unsigned i = 1; i <= I.m; i++)
      C[i][j] = Cj = std::max(Cj, C[i][j-1]) + I.p[job][i];
    pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[job]));
  }
  bool improved = true;
  unsigned n_ls = params.n_ls;
  while(improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    improved = false;
    std::vector<unsigned> cpi = n.pi;
    unsigned bk = 0, bl = 0; 
    unsigned btt = n.ms;     
    for(unsigned k = 1; k < I.n; k++) {
      for(unsigned l = k+1; l <= I.n; l++) {
        std::swap(cpi[k], cpi[l]);
        std::vector<unsigned> Cm(I.m+1);
        for(unsigned i = 1; i <= I.m; i++)
          Cm[i] = std::max(Cm[i-1], C[i][k-1]) + I.p[cpi[k]][i];
        unsigned tt = pt[k-1] + std::max(0, int(Cm[I.m]-I.dd[cpi[k]]));
        for(unsigned j = k+1; j <= I.n; j++) {
          unsigned job = cpi[j];
          for(unsigned i = 1; i <= I.m; i++)
            Cm[i] = std::max(Cm[i-1], Cm[i]) + I.p[job][i];
          tt += std::max(0, int(Cm[I.m]-I.dd[job]));
        }
        if(tt < btt) { 
          btt = tt;
          bk = k;
          bl = l;
          improved = true;
        }
        std::swap(cpi[k], cpi[l]);
      }
    }
    if(improved) {
      std::swap(n.pi[bk], n.pi[bl]);
      n.ms = btt;
      for(unsigned j = bk; j <= I.n; j++) {
        unsigned job = n.pi[j];
        unsigned Cj = 0;
        for(unsigned i = 1; i <= I.m; i++)
          C[i][j] = Cj = std::max(Cj, C[i][j-1]) + I.p[job][i];
        pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[job]));
      }
    }
  }
}
template <typename Tiebreak>
inline void CH3(const Instance& I, PFSSolution &n, builder& construct, Tiebreak tiebreak, Params& params) {
  unsigned kl = 1;
  unsigned ku = I.n+1;
  unsigned n_ls = params.n_ls;
  boost::multi_array<unsigned,2> C(boost::extents[I.m+1][I.n+1]); 
  std::vector<unsigned> pt(I.n+1);  
  bool improved = true;
  while(improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned cms = n.ms;
    std::vector<Job> cpi = n.pi;
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = cpi[ji];
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      construct.constructpTT(I,n,I.n,kl,ku-1,remove,tiebreak,infinite_time);
    }
    for(unsigned j = 1; j <= I.n; j++) {
      unsigned job = n.pi[j];
      unsigned Cj = 0;
      for(unsigned i = 1; i <= I.m; i++)
        C[i][j] = Cj = std::max(Cj, C[i][j-1]) + I.p[job][i];
      pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[job]));
    }
    cpi = n.pi;
    unsigned bk = 0, bl = 0; 
    unsigned btt = n.ms;     
    bool swap_improved = false;
    for(unsigned k = 1; k < I.n; k++) {
      for(unsigned l = k+1; l <= I.n; l++) {
        std::swap(cpi[k], cpi[l]);
        std::vector<unsigned> Cm(I.m+1);
        for(unsigned i = 1; i <= I.m; i++)
          Cm[i] = std::max(Cm[i-1], C[i][k-1]) + I.p[cpi[k]][i];
        unsigned tt = pt[k-1] + std::max(0, int(Cm[I.m]-I.dd[cpi[k]]));
        for(unsigned j = k+1; j <= I.n; j++) {
          unsigned job = cpi[j];
          for(unsigned i = 1; i <= I.m; i++)
            Cm[i] = std::max(Cm[i-1], Cm[i]) + I.p[job][i];
          tt += std::max(0, int(Cm[I.m]-I.dd[job]));
        }
        if(tt < btt) { 
          btt = tt;
          bk = k;
          bl = l;
          swap_improved = true;
        }
        std::swap(cpi[k], cpi[l]);
      }
    }
    if(swap_improved) {
      std::swap(n.pi[bk], n.pi[bl]);
      n.ms = btt;
    }
    if(n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void CH3_IT1(const Instance& I, PFSSolution &n, builder& construct, Params& params) {
  unsigned kl = 1;
  unsigned ku = I.n+1;
  unsigned n_ls = params.n_ls;
  boost::multi_array<unsigned,2> C(boost::extents[I.m+1][I.n+1]); 
  std::vector<unsigned> pt(I.n+1);  
  std::vector<unsigned> tpt(I.m+1); 
  bool improved = true;
  while(improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned cms = n.ms;
    std::vector<Job> cpi = n.pi;
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = cpi[ji];
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      construct.constructpTT_IT1(I,n,I.n,kl,ku-1,remove,infinite_time);
    }
    for(unsigned j = 1; j <= I.n; j++) {
      unsigned job = n.pi[j];
      unsigned Cj = 0;
      for(unsigned i = 1; i <= I.m; i++)
        C[i][j] = Cj = std::max(Cj, C[i][j-1]) + I.p[job][i];
      pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[job]));
    }
    cpi = n.pi;
    unsigned bk = 0, bl = 0; 
    unsigned btt = n.ms;     
    bool swap_improved = false;
    for(unsigned k = 1; k < I.n; k++) {
      for(unsigned l = k+1; l <= I.n; l++) {
        std::swap(cpi[k], cpi[l]);
        std::vector<unsigned> Cm(I.m+1);
        for(unsigned i = 1; i <= I.m; i++)
          Cm[i] = std::max(Cm[i-1], C[i][k-1]) + I.p[cpi[k]][i];
        unsigned tt = pt[k-1] + std::max(0, int(Cm[I.m]-I.dd[cpi[k]]));
        for(unsigned j = k+1; j <= I.n; j++) {
          unsigned job = cpi[j];
          for(unsigned i = 1; i <= I.m; i++)
            Cm[i] = std::max(Cm[i-1], Cm[i]) + I.p[job][i];
          tt += std::max(0, int(Cm[I.m]-I.dd[job]));
        }
        if(tt < btt) { 
          btt = tt;
          bk = k;
          bl = l;
          swap_improved = true;
        }
        std::swap(cpi[k], cpi[l]);
      }
    }
    if(swap_improved) {
      std::swap(n.pi[bk], n.pi[bl]);
      n.ms = btt;
    }
    if(n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void CH3_IT2(const Instance& I, PFSSolution &n, builder& construct, Params& params) {
  unsigned kl = 1;
  unsigned ku = I.n+1;
  unsigned n_ls = params.n_ls;
  boost::multi_array<unsigned,2> C(boost::extents[I.m+1][I.n+1]); 
  std::vector<unsigned> pt(I.n+1);  
  std::vector<unsigned> tpt(I.m+1); 
  bool improved = true;
  while(improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned cms = n.ms;
    std::vector<Job> cpi = n.pi;
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = cpi[ji];
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      construct.constructpTT_IT2(I,n,I.n,kl,ku-1,remove,infinite_time);
    }
    for(unsigned j = 1; j <= I.n; j++) {
      unsigned job = n.pi[j];
      unsigned Cj = 0;
      for(unsigned i = 1; i <= I.m; i++)
        C[i][j] = Cj = std::max(Cj, C[i][j-1]) + I.p[job][i];
      pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[job]));
    }
    cpi = n.pi;
    unsigned bk = 0, bl = 0; 
    unsigned btt = n.ms;     
    bool swap_improved = false;
    for(unsigned k = 1; k < I.n; k++) {
      for(unsigned l = k+1; l <= I.n; l++) {
        std::swap(cpi[k], cpi[l]);
        std::vector<unsigned> Cm(I.m+1);
        for(unsigned i = 1; i <= I.m; i++)
          Cm[i] = std::max(Cm[i-1], C[i][k-1]) + I.p[cpi[k]][i];
        unsigned tt = pt[k-1] + std::max(0, int(Cm[I.m]-I.dd[cpi[k]]));
        for(unsigned j = k+1; j <= I.n; j++) {
          unsigned job = cpi[j];
          for(unsigned i = 1; i <= I.m; i++)
            Cm[i] = std::max(Cm[i-1], Cm[i]) + I.p[job][i];
          tt += std::max(0, int(Cm[I.m]-I.dd[job]));
        }
        if(tt < btt) { 
          btt = tt;
          bk = k;
          bl = l;
          swap_improved = true;
        }
        std::swap(cpi[k], cpi[l]);
      }
    }
    if(swap_improved) {
      std::swap(n.pi[bk], n.pi[bl]);
      n.ms = btt;
    }
    if(n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
template <typename Tiebreak>
inline void CH4(const Instance& I, PFSSolution &n, builder& construct, Tiebreak tiebreak, Params& params) {
  unsigned kl = 1;
  unsigned ku = I.n+1;
  unsigned n_ls = params.n_ls;
  boost::multi_array<unsigned,2> C(boost::extents[I.m+1][I.n+1]); 
  std::vector<unsigned> pt(I.n+1);  
  bool improved = true;
  while(improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned cms = n.ms;
    for(unsigned j = 1; j <= I.n; j++) {
      unsigned job = n.pi[j];
      unsigned Cj = 0;
      for(unsigned i = 1; i <= I.m; i++)
        C[i][j] = Cj = std::max(Cj, C[i][j-1]) + I.p[job][i];
      pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[job]));
    }
    std::vector<Job> cpi = n.pi;
    unsigned bk = 0, bl = 0; 
    unsigned btt = n.ms;     
    bool swap_improved = false;
    for(unsigned k = 1; k < I.n; k++) {
      for(unsigned l = k+1; l <= I.n; l++) {
        std::swap(cpi[k], cpi[l]);
        std::vector<unsigned> Cm(I.m+1);
        for(unsigned i = 1; i <= I.m; i++)
          Cm[i] = std::max(Cm[i-1], C[i][k-1]) + I.p[cpi[k]][i];
        unsigned tt = pt[k-1] + std::max(0, int(Cm[I.m]-I.dd[cpi[k]]));
        for(unsigned j = k+1; j <= I.n; j++) {
          unsigned job = cpi[j];
          for(unsigned i = 1; i <= I.m; i++)
            Cm[i] = std::max(Cm[i-1], Cm[i]) + I.p[job][i];
          tt += std::max(0, int(Cm[I.m]-I.dd[job]));
        }
        if(tt < btt) { 
          btt = tt;
          bk = k;
          bl = l;
          swap_improved = true;
        }
        std::swap(cpi[k], cpi[l]);
      }
    }
    if(swap_improved) {
      std::swap(n.pi[bk], n.pi[bl]);
      std::swap(cpi[bk], cpi[bl]);
      n.ms = btt;
    }
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = cpi[ji];
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      construct.constructpTT(I,n,I.n,kl,ku-1,remove,tiebreak,infinite_time);
    }
    if(n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void CH4_IT1(const Instance& I, PFSSolution &n, builder& construct, Params& params) {
  unsigned kl = 1;
  unsigned ku = I.n+1;
  unsigned n_ls = params.n_ls;
  boost::multi_array<unsigned,2> C(boost::extents[I.m+1][I.n+1]); 
  std::vector<unsigned> pt(I.n+1);  
  bool improved = true;
  while(improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned cms = n.ms;
    for(unsigned j = 1; j <= I.n; j++) {
      unsigned job = n.pi[j];
      unsigned Cj = 0;
      for(unsigned i = 1; i <= I.m; i++)
        C[i][j] = Cj = std::max(Cj, C[i][j-1]) + I.p[job][i];
      pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[job]));
    }
    std::vector<Job> cpi = n.pi;
    unsigned bk = 0, bl = 0; 
    unsigned btt = n.ms;     
    bool swap_improved = false;
    for(unsigned k = 1; k < I.n; k++) {
      for(unsigned l = k+1; l <= I.n; l++) {
        std::swap(cpi[k], cpi[l]);
        std::vector<unsigned> Cm(I.m+1);
        for(unsigned i = 1; i <= I.m; i++)
          Cm[i] = std::max(Cm[i-1], C[i][k-1]) + I.p[cpi[k]][i];
        unsigned tt = pt[k-1] + std::max(0, int(Cm[I.m]-I.dd[cpi[k]]));
        for(unsigned j = k+1; j <= I.n; j++) {
          unsigned job = cpi[j];
          for(unsigned i = 1; i <= I.m; i++)
            Cm[i] = std::max(Cm[i-1], Cm[i]) + I.p[job][i];
          tt += std::max(0, int(Cm[I.m]-I.dd[job]));
        }
        if(tt < btt) { 
          btt = tt;
          bk = k;
          bl = l;
          swap_improved = true;
        }
        std::swap(cpi[k], cpi[l]);
      }
    }
    if(swap_improved) {
      std::swap(n.pi[bk], n.pi[bl]);
      std::swap(cpi[bk], cpi[bl]);
      n.ms = btt;
    }
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = cpi[ji];
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      construct.constructpTT_IT1(I,n,I.n,kl,ku-1,remove,infinite_time);
    }
    if(n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void CH4_IT2(const Instance& I, PFSSolution &n, builder& construct, Params& params) {
  unsigned kl = 1;
  unsigned ku = I.n+1;
  unsigned n_ls = params.n_ls;
  boost::multi_array<unsigned,2> C(boost::extents[I.m+1][I.n+1]); 
  std::vector<unsigned> pt(I.n+1);  
  bool improved = true;
  while(improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    improved = false;
    unsigned cms = n.ms;
    for(unsigned j = 1; j <= I.n; j++) {
      unsigned job = n.pi[j];
      unsigned Cj = 0;
      for(unsigned i = 1; i <= I.m; i++)
        C[i][j] = Cj = std::max(Cj, C[i][j-1]) + I.p[job][i];
      pt[j] = pt[j-1] + std::max(0, int(Cj-I.dd[job]));
    }
    std::vector<Job> cpi = n.pi;
    unsigned bk = 0, bl = 0; 
    unsigned btt = n.ms;     
    bool swap_improved = false;
    for(unsigned k = 1; k < I.n; k++) {
      for(unsigned l = k+1; l <= I.n; l++) {
        std::swap(cpi[k], cpi[l]);
        std::vector<unsigned> Cm(I.m+1);
        for(unsigned i = 1; i <= I.m; i++)
          Cm[i] = std::max(Cm[i-1], C[i][k-1]) + I.p[cpi[k]][i];
        unsigned tt = pt[k-1] + std::max(0, int(Cm[I.m]-I.dd[cpi[k]]));
        for(unsigned j = k+1; j <= I.n; j++) {
          unsigned job = cpi[j];
          for(unsigned i = 1; i <= I.m; i++)
            Cm[i] = std::max(Cm[i-1], Cm[i]) + I.p[job][i];
          tt += std::max(0, int(Cm[I.m]-I.dd[job]));
        }
        if(tt < btt) { 
          btt = tt;
          bk = k;
          bl = l;
          swap_improved = true;
        }
        std::swap(cpi[k], cpi[l]);
      }
    }
    if(swap_improved) {
      std::swap(n.pi[bk], n.pi[bl]);
      std::swap(cpi[bk], cpi[bl]);
      n.ms = btt;
    }
    for(unsigned ji = ku-1; ji>=kl; ji--) {
      Job j = cpi[ji];
      auto nend = std::remove(n.pi.begin(),n.pi.end(),j);
      std::fill(nend, n.pi.end(), 0);
      std::vector<Job> remove = { 0, j };
      construct.constructpTT_IT2(I,n,I.n,kl,ku-1,remove,infinite_time);
    }
    if(n.ms < cms)
      improved = true;
  }
}
template <typename Tiebreak>
inline void CH5(const Instance& I, PFSSolution &n, builder& construct, Tiebreak tiebreak, Params& params) {
  unsigned tt;
  unsigned n_ls = params.n_ls;
  do {
    tt = n.ms;
    insertion(I, n, construct, tiebreak, params);
    BISLS_TT(I, n, params);
  } while (n.ms < tt && n_ls-- > 0  && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit);
}
inline void CH5_IT1(const Instance& I, PFSSolution &n, builder& construct, Params& params) {
  unsigned n_ls = params.n_ls;
  bool improved = true;
  while (improved && n_ls-- > 0  && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned tt = n.ms;
    insertionIT1(I, n, construct, params);
    BISLS_TT(I, n, params);
    if(n.ms < tt)
      improved = true;
    else
      improved = false;
  }
}
inline void CH5_IT2(const Instance& I, PFSSolution &n, builder& construct, Params& params) {
  unsigned n_ls = params.n_ls;
  bool improved = true;
  while (improved && n_ls-- > 0  && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned tt = n.ms;
    insertionIT2(I, n, construct, params);
    BISLS_TT(I, n, params);
    if(n.ms < tt)
      improved = true;
    else
      improved = false;
  }
}
template <typename Tiebreak>
inline void CH6(const Instance& I, PFSSolution &n, builder& construct, Tiebreak tiebreak, Params& params) {
  unsigned n_ls = params.n_ls;
  bool improved = true;
  while (improved && n_ls-- > 0  && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned tt = n.ms;
    BISLS_TT(I, n, params);
    insertion(I, n, construct, tiebreak, params);
    if(n.ms < tt)
      improved = true;
    else
      improved = false;
  }
}
inline void CH6_IT1(const Instance& I, PFSSolution &n, builder& construct, Params& params) {
  unsigned n_ls = params.n_ls;
  bool improved = true;
  while (improved && n_ls-- > 0  && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned tt = n.ms;
    BISLS_TT(I, n, params);
    insertionIT1(I, n, construct, params);
    if(n.ms < tt)
      improved = true;
    else
      improved = false;
  }
}
inline void CH6_IT2(const Instance& I, PFSSolution &n, builder& construct, Params& params) {
  unsigned n_ls = params.n_ls;
  bool improved = true;
  while (improved && n_ls-- > 0  && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit) {
    unsigned tt = n.ms;
    BISLS_TT(I, n, params);
    insertionIT2(I, n, construct, params);
    if(n.ms < tt)
      improved = true;
    else
      improved = false;
  }
}
template<typename Tiebreak>
inline void insertion_np(const Instance& I, NPFSSolution &n, Params& params, Tiebreak tiebreak) {
  unsigned n_ls = params.n_ls_np;
  using namespace std::placeholders;
  nBuilder construct(I);
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    unsigned cms = n.ms;
    for(unsigned ji = I.n; ji >= 1; ji--) {
      Job job = n.pi[1][ji];
      for(unsigned i = 1; i <= I.m; i++) {
         unsigned dj = 1;
         for(unsigned j = 1; j <= I.n; j++) {
            if (n.pi[i][j] != job) {
               n.pi[i][dj] = n.pi[i][j];
               dj++;
            }
         }
         while (dj <= I.n) {
            n.pi[i][dj] = 0;
            dj++;
         }
      }
      std::vector<Job> remove = { 0, job };
      if (params.obj == 'm'){
        construct.constructnpNEW(n, I.n, remove, tiebreak);
      }
      else if (params.obj == 'f') {
        construct.constructnpNEWCsum(n, I.n, remove, tiebreak);
      }
      else if (params.obj == 't'){
        construct.constructnpNEWTT(n, I.n, remove, tiebreak);
      }
    }
    if (n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void CPASLS_Cmax_np(const Instance& I, NPFSSolution& n, Params& params) {
  unsigned n_ls = params.n_ls_np;
  boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
  boost::multi_array<Time,2> C(boost::extents[I.m+1][I.n+1]); 
  boost::multi_array<Time,2> h(boost::extents[I.m+1][I.n+1]), t(boost::extents[I.m+2][I.n+2]);
  boost::multi_array<unsigned,2> mi(boost::extents[I.m+1][I.n+2]);
  for (unsigned i = 1; i <= I.m; i++) {
    Time c = 0;
    for (unsigned j = 1; j <= I.n; j++) {
      c = std::max(c, n.C[i-1][n.pi[i][j]]) + I.p[n.pi[i][j]][i];
      n.C[i][n.pi[i][j]] = c; 
    }
  }
  bool improved = true;
  while (improved && n_ls-- > 0  && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = 1; j <= I.n; j++) {
        npi[i][j] = n.pi[i][j];
        invpi[i][npi[i][j]] = j;
        C[i][j] = n.C[i][j];
      }
    }
    for (unsigned i = I.m-1; i >= 1; i--) { 
      mi[i][I.n+1] = infinite_time;
      for (unsigned j = I.n; j >= 1; j--) {
        mi[i][j] = std::min(mi[i][j+1], invpi[i+1][npi[i][j]]);
      }
    }
    for(unsigned i = 1; i <= I.m; i++)
      for(unsigned j = 1; j <= I.n; j++)
        h[i][j] = std::max(h[i-1][invpi[i-1][npi[i][j]]], h[i][j-1]) + I.p[npi[i][j]][i];
    for(unsigned i = I.m; i >= 1; i--)
      for(unsigned j = I.n; j >= 1; j--)
        t[i][j] = std::max(t[i+1][invpi[i+1][npi[i][j]]], t[i][j+1]) + I.p[npi[i][j]][i];
    std::vector<std::pair<unsigned,unsigned>> r;
    std::vector<char> duplicate(I.n+1,0);
    unsigned cmax = C[I.m][npi[I.m][I.n]];
    for(unsigned i = 1; i <= I.m; i++)
      for(unsigned j = 1; j < I.n; j++)
        if(!duplicate[j] && (h[i][j] + t[i+1][j] == cmax || h[i][j+1] + t[i+1][j+1] == cmax)) {
          r.push_back(std::make_pair(j,j+1));
          duplicate[j] = 1;
        }
    boost::multi_array<unsigned,2> bs = npi; 
    boost::multi_array<unsigned,2> bC = C;   
    Time bo = n.ms;                          
    Time po = n.ms;                          
    for(unsigned k = 0; k < r.size(); k++) {
      assert(r[k].first != r[k].second);
      unsigned pos;
      std::swap(npi[1][r[k].first], npi[1][r[k].second]);
      std::swap(invpi[1][r[k].first], invpi[1][r[k].second]);
      std::swap(npi[2][r[k].first], npi[2][r[k].second]);
      std::swap(invpi[2][r[k].first], invpi[2][r[k].second]);
      pos = r[k].first;
      for(unsigned ii = 1; ii <= I.m; ii++) {
        Time ct = C[ii][npi[ii][pos-1]];
        for(unsigned ji = pos; ji <= I.n; ji++) {
          unsigned job = npi[ii][ji];
          ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
          C[ii][job] = ct;
        }
        pos = mi[ii][pos];
      }
      if (C[I.m][npi[I.m][I.n]] < bo) {
        bo = C[I.m][npi[I.m][I.n]];
        bs = npi;
        bC = C;
      }
      for (unsigned i = 3; i <= I.m-2; i++) {
        std::swap(npi[i][r[k].first], npi[i][r[k].second]);
        std::swap(invpi[i][r[k].first], invpi[i][r[k].second]);
        pos = r[k].first;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        if (C[I.m][npi[I.m][I.n]] < bo) {
          bo = C[I.m][npi[I.m][I.n]];
          bs = npi;
          bC = C;
        }
      }
      std::swap(npi[I.m-1][r[k].first], npi[I.m-1][r[k].second]);
      std::swap(invpi[I.m-1][r[k].first], invpi[I.m-1][r[k].second]);
      std::swap(npi[I.m][r[k].first], npi[I.m][r[k].second]);
      std::swap(invpi[I.m][r[k].first], invpi[I.m][r[k].second]);
      pos = r[k].first;
      for (unsigned ii = I.m-1; ii <= I.m; ii++) {
        Time ct = C[ii][npi[ii][pos-1]];
        for(unsigned ji = pos; ji <= I.n; ji++) {
          unsigned job = npi[ii][ji];
          ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
          C[ii][job] = ct;
        }
        pos = mi[ii][pos];
      }
      if (C[I.m][npi[I.m][I.n]] <= bo) {
        bo = C[I.m][npi[I.m][I.n]];
        bs = npi;
        bC = C;
      }
      std::swap(npi[1][r[k].first], npi[1][r[k].second]);
      std::swap(invpi[1][r[k].first], invpi[1][r[k].second]);
      std::swap(npi[2][r[k].first], npi[2][r[k].second]);
      std::swap(invpi[2][r[k].first], invpi[2][r[k].second]);
      pos = r[k].first;
      for (unsigned ii = 1; ii <= I.m; ii++) {
        Time ct = C[ii][npi[ii][pos-1]];
        for(unsigned ji = pos; ji <= I.n; ji++) {
          unsigned job = npi[ii][ji];
          ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
          C[ii][job] = ct;
        }
        pos = mi[ii][pos];
      }
      if (C[I.m][npi[I.m][I.n]] < bo) {
        bo = C[I.m][npi[I.m][I.n]];
        bs = npi;
        bC = C;
      }
      for (unsigned i = 3; i <= I.m-2; i++) {
        std::swap(npi[i][r[k].first], npi[i][r[k].second]);
        std::swap(invpi[i][r[k].first], invpi[i][r[k].second]);
        pos = r[k].first;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        if (C[I.m][npi[I.m][I.n]] < bo) {
          bo = C[I.m][npi[I.m][I.n]];
          bs = npi;
          bC = C;
        }
      }
      std::swap(npi[I.m-1][r[k].first], npi[I.m-1][r[k].second]);
      std::swap(invpi[I.m-1][r[k].first], invpi[I.m-1][r[k].second]);
      std::swap(npi[I.m][r[k].first], npi[I.m][r[k].second]);
      std::swap(invpi[I.m][r[k].first], invpi[I.m][r[k].second]);
      pos = r[k].first;
      for (unsigned ii = I.m-1; ii <= I.m; ii++) {
        Time ct = C[ii][npi[ii][pos-1]];
        for(unsigned ji = pos; ji <= I.n; ji++) {
          unsigned job = npi[ii][ji];
          ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
          C[ii][job] = ct;
        }
        pos = mi[ii][pos];
      }
    }
    if (bo < po) {
      n.pi = bs;
      n.C = bC;
      n.ms = bo;
      improved = true;
    }
    else {
      improved = false;
    }
  }
}
inline void ASLS_Cmax_np(const Instance& I, NPFSSolution& n, Params& params) {
  unsigned n_ls = params.n_ls_np;
  boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
  boost::multi_array<Time,2> C(boost::extents[I.m+1][I.n+1]); 
  boost::multi_array<Time,2> h(boost::extents[I.m+1][I.n+1]), t(boost::extents[I.m+2][I.n+2]);
  boost::multi_array<unsigned,2> mi(boost::extents[I.m+1][I.n+2]);
  for (unsigned i = 1; i <= I.m; i++) {
    Time c = 0;
    for (unsigned j = 1; j <= I.n; j++) {
      c = std::max(c, n.C[i-1][n.pi[i][j]]) + I.p[n.pi[i][j]][i];
      n.C[i][n.pi[i][j]] = c; 
    }
  }
  bool improved = true;
  while (improved && n_ls-- > 0  && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = 1; j <= I.n; j++) {
        npi[i][j] = n.pi[i][j];
        invpi[i][npi[i][j]] = j;
        C[i][j] = n.C[i][j];
      }
    }
    for (unsigned i = I.m-1; i >= 1; i--) { 
      mi[i][I.n+1] = infinite_time;
      for (unsigned j = I.n; j >= 1; j--) {
        mi[i][j] = std::min(mi[i][j+1], invpi[i+1][npi[i][j]]);
      }
    }
    for(unsigned i = 1; i <= I.m; i++)
      for(unsigned j = 1; j <= I.n; j++)
        h[i][j] = std::max(h[i-1][invpi[i-1][npi[i][j]]], h[i][j-1]) + I.p[npi[i][j]][i];
    for(unsigned i = I.m; i >= 1; i--)
      for(unsigned j = I.n; j >= 1; j--)
        t[i][j] = std::max(t[i+1][invpi[i+1][npi[i][j]]], t[i][j+1]) + I.p[npi[i][j]][i];
    boost::multi_array<unsigned,2> bs = npi; 
    boost::multi_array<unsigned,2> bC = C;   
    Time bo = n.ms;                          
    Time po = n.ms;                          
    for (unsigned k = 1; k < I.n; k++) {
      unsigned pos;
      std::swap(npi[1][k], npi[1][k+1]);
      std::swap(invpi[1][k], invpi[1][k+1]);
      std::swap(npi[2][k], npi[2][k+1]);
      std::swap(invpi[2][k], invpi[2][k+1]);
      pos = k;
      for (unsigned ii = 1; ii <= I.m; ii++) {
        Time ct = C[ii][npi[ii][pos-1]];
        for(unsigned ji = pos; ji <= I.n; ji++) {
          unsigned job = npi[ii][ji];
          ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
          C[ii][job] = ct;
        }
        pos = mi[ii][pos];
      }
      if (C[I.m][npi[I.m][I.n]] < bo) {
        bo = C[I.m][npi[I.m][I.n]];
        bs = npi;
        bC = C;
      }
      for (unsigned i = 3; i <= I.m-2; i++) {
        std::swap(npi[i][k], npi[i][k+1]);
        std::swap(invpi[i][k], invpi[i][k+1]);
        pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        if (C[I.m][npi[I.m][I.n]] < bo) {
          bo = C[I.m][npi[I.m][I.n]];
          bs = npi;
          bC = C;
        }
      }
      std::swap(npi[I.m-1][k], npi[I.m-1][k+1]);
      std::swap(invpi[I.m-1][k], invpi[I.m-1][k+1]);
      std::swap(npi[I.m][k], npi[I.m][k+1]);
      std::swap(invpi[I.m][k], invpi[I.m][k+1]);
      pos = k;
      for (unsigned ii = I.m-1; ii <= I.m; ii++) {
        Time ct = C[ii][npi[ii][pos-1]];
        for(unsigned ji = pos; ji <= I.n; ji++) {
          unsigned job = npi[ii][ji];
          ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
          C[ii][job] = ct;
        }
        pos = mi[ii][pos];
      }
      if (C[I.m][npi[I.m][I.n]] <= bo) {
        bo = C[I.m][npi[I.m][I.n]];
        bs = npi;
        bC = C;
      }
      std::swap(npi[1][k], npi[1][k+1]);
      std::swap(invpi[1][k], invpi[1][k+1]);
      std::swap(npi[2][k], npi[2][k+1]);
      std::swap(invpi[2][k], invpi[2][k+1]);
      pos = k;
      for (unsigned ii = 1; ii <= I.m; ii++) {
        Time ct = C[ii][npi[ii][pos-1]];
        for(unsigned ji = pos; ji <= I.n; ji++) {
          unsigned job = npi[ii][ji];
          ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
          C[ii][job] = ct;
        }
        pos = mi[ii][pos];
      }
      if (C[I.m][npi[I.m][I.n]] < bo) {
        bo = C[I.m][npi[I.m][I.n]];
        bs = npi;
        bC = C;
      }
      for (unsigned i = 3; i <= I.m-2; i++) {
        std::swap(npi[i][k], npi[i][k+1]);
        std::swap(invpi[i][k], invpi[i][k+1]);
        pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        if (C[I.m][npi[I.m][I.n]] < bo) {
          bo = C[I.m][npi[I.m][I.n]];
          bs = npi;
          bC = C;
        }
      }
      std::swap(npi[I.m-1][k], npi[I.m-1][k+1]);
      std::swap(invpi[I.m-1][k], invpi[I.m-1][k+1]);
      std::swap(npi[I.m][k], npi[I.m][k+1]);
      std::swap(invpi[I.m][k], invpi[I.m][k+1]);
      pos = k;
      for (unsigned ii = I.m-1; ii <= I.m; ii++) {
        Time ct = C[ii][npi[ii][pos-1]];
        for(unsigned ji = pos; ji <= I.n; ji++) {
          unsigned job = npi[ii][ji];
          ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
          C[ii][job] = ct;
        }
        pos = mi[ii][pos];
      }
    }
    if (bo < po) {
      n.pi = bs;
      n.C = bC;
      n.ms = bo;
      improved = true;
    }
    else {
      improved = false;
    }
  }
}
template <typename Tiebreak>
inline void Pc_np(const Instance& I, NPFSSolution& n, nBuilder& construct, Tiebreak tiebreak, std::mt19937& rng, Params& params) {
  unsigned n_ls = params.n_ls_np;
  boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
  boost::multi_array<Time,2> h(boost::extents[I.m+1][I.n+1]), t(boost::extents[I.m+2][I.n+2]);
  bool improved = true;
  while (improved && n_ls-- > 0  && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = 1; j <= I.n; j++) {
        npi[i][j] = n.pi[i][j];
        invpi[i][npi[i][j]] = j;
      }
    }
    for(unsigned i = 1; i <= I.m; i++)
      for(unsigned j = 1; j <= I.n; j++)
        h[i][j] = std::max(h[i-1][invpi[i-1][npi[i][j]]], h[i][j-1]) + I.p[npi[i][j]][i];
    for(unsigned i = I.m; i >= 1; i--)
      for(unsigned j = I.n; j >= 1; j--)
        t[i][j] = std::max(t[i+1][invpi[i+1][npi[i][j]]], t[i][j+1]) + I.p[npi[i][j]][i];
    std::vector<std::pair<unsigned,unsigned>> r;
    std::vector<char> duplicate(I.n+1,0);
    unsigned cmax = h[I.m][I.n];
    for(unsigned i = 1; i < I.m; i++)
      for(unsigned j = 1; j < I.n; j++)
        if(!duplicate[j] && (h[i][j] + t[i+1][j] == cmax || h[i][j+1] + t[i+1][j+1] == cmax)) {
          r.push_back(std::make_pair(j,j+1));
          duplicate[j] = 1;
        }
    std::shuffle(r.begin(),r.end(),rng);
    for(unsigned k = 0; k < r.size(); k++) {
      assert(r[k].first != r[k].second);
      NPFSSolution ns = n;
      Job job1 = r[k].first;
      Job job2 = r[k].second;
      for(unsigned i = 1; i <= I.m; i++) {
         unsigned dj = 1;
         for(unsigned j = 1; j <= I.n; j++) {
            if (ns.pi[i][j] != job1) {
               ns.pi[i][dj] = ns.pi[i][j];
               dj++;
            }
         }
         while (dj <= I.n) {
            ns.pi[i][dj] = 0;
            dj++;
         }
         dj = 1;
         for(unsigned j = 1; j <= I.n; j++) {
            if (ns.pi[i][j] != job2) {
               ns.pi[i][dj] = ns.pi[i][j];
               dj++;
            }
         }
         while (dj <= I.n) {
            ns.pi[i][dj] = 0;
            dj++;
         }
      }
      std::vector<Job> remove = { 0, job1, job2 };
      construct.constructnpNEW(ns, I.n-1, remove, tiebreak);
      if (ns.ms <= n.ms) {
        n = ns;
        improved = true;
      }
      else {
        improved = false;
      }
    }
  }
}
inline void AASLS(const Instance& I, NPFSSolution& n, Params& params) {
  unsigned n_ls = params.n_ls_np;
  boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
  boost::multi_array<Time,2> C(boost::extents[I.m+1][I.n+1]);
  boost::multi_array<unsigned,2> mi(boost::extents[I.m+1][I.n+2]);
  for (unsigned i = 1; i <= I.m; i++) {
    Time c = 0;
    for (unsigned j = 1; j <= I.n; j++) {
      c = std::max(c, n.C[i-1][n.pi[i][j]]) + I.p[n.pi[i][j]][i];
      n.C[i][n.pi[i][j]] = c; 
    }
  }
  for (unsigned i = 1; i <= I.m; i++) {
    npi[i][0] = 0;
  }
  for (unsigned j = 1; j <= I.n; j++) {
    C[0][j] = 0;
  }
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = 1; j <= I.n; j++) {
        npi[i][j] = n.pi[i][j];
        invpi[i][npi[i][j]] = j;
        C[i][j] = n.C[i][j];
      }
    }
    boost::multi_array<unsigned,2> bs = npi; 
    boost::multi_array<unsigned,2> bC = C;   
    Time bf = n.ms;                          
    Time pf = n.ms;                          
    for (unsigned i = I.m-1; i >= 1; i--) { 
      mi[i][I.n+1] = infinite_time;
      for (unsigned j = I.n; j >= 1; j--) {
        mi[i][j] = std::min(mi[i][j+1], invpi[i+1][npi[i][j]]);
      }
    }
    for (unsigned j = 0; j < I.n+2; j++) { 
      mi[I.m][j] = 0;
    }
    for (unsigned k = 1; k < I.n; k++) {
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned csum = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          csum += C[I.m][j];
        }
        if (csum < bf) {
          bf = csum;
          bs = npi;
          bC = C;
        }
      }
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned csum = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          csum += C[I.m][j];
        }
        if (csum < bf) {
          bf = csum;
          bs = npi;
          bC = C;
        }
      }
    }
    if (bf < pf) {
      improved = true;
      n.pi = bs;
      n.C = bC;
      n.ms = bf;
    }
    else {
      improved = false;
    }
  }
}
inline void AASLS_r(const Instance& I, NPFSSolution& n, Params& params) {
  unsigned n_ls = params.n_ls_np;
  boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
  boost::multi_array<Time,2> C(boost::extents[I.m+1][I.n+1]);
  boost::multi_array<unsigned,2> mi(boost::extents[I.m+1][I.n+2]);
  for (unsigned i = 1; i <= I.m; i++) {
    Time c = 0;
    for (unsigned j = 1; j <= I.n; j++) {
      c = std::max(c, n.C[i-1][n.pi[i][j]]) + I.p[n.pi[i][j]][i];
      n.C[i][n.pi[i][j]] = c; 
    }
  }
  for (unsigned i = 1; i <= I.m; i++) {
    npi[i][0] = 0;
  }
  for (unsigned j = 1; j <= I.n; j++) {
    C[0][j] = 0;
  }
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = 1; j <= I.n; j++) {
        npi[i][j] = n.pi[i][j];
        invpi[i][npi[i][j]] = j;
        C[i][j] = n.C[i][j];
      }
    }
    boost::multi_array<unsigned,2> bs = npi; 
    boost::multi_array<unsigned,2> bC = C;   
    Time bf = n.ms;                          
    Time pf = n.ms;                          
    for (unsigned i = I.m-1; i >= 1; i--) { 
      mi[i][I.n+1] = infinite_time;
      for (unsigned j = I.n; j >= 1; j--) {
        mi[i][j] = std::min(mi[i][j+1], invpi[i+1][npi[i][j]]);
      }
    }
    for (unsigned j = 0; j < I.n+2; j++) { 
      mi[I.m][j] = 0;
    }
    for (unsigned k = I.n-1; k >= 1; k--) {
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned csum = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          csum += C[I.m][j];
        }
        if (csum < bf) {
          bf = csum;
          bs = npi;
          bC = C;
        }
      }
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned csum = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          csum += C[I.m][j];
        }
        if (csum < bf) {
          bf = csum;
          bs = npi;
          bC = C;
        }
      }
    }
    if (bf < pf) {
      improved = true;
      n.pi = bs;
      n.C = bC;
      n.ms = bf;
    }
    else {
      improved = false;
    }
  }
}
inline void AASLS_G8(const Instance& I, NPFSSolution& n, Params& params) {
  unsigned n_ls = params.n_ls_np;
  boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
  boost::multi_array<Time,2> C(boost::extents[I.m+1][I.n+1]);
  boost::multi_array<unsigned,2> mi(boost::extents[I.m+1][I.n+2]);
  for (unsigned i = 1; i <= I.m; i++) {
    Time c = 0;
    for (unsigned j = 1; j <= I.n; j++) {
      c = std::max(c, n.C[i-1][n.pi[i][j]]) + I.p[n.pi[i][j]][i];
      n.C[i][n.pi[i][j]] = c; 
    }
  }
  for (unsigned i = 1; i <= I.m; i++) {
    npi[i][0] = 0;
  }
  for (unsigned j = 1; j <= I.n; j++) {
    C[0][j] = 0;
  }
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = 1; j <= I.n; j++) {
        npi[i][j] = n.pi[i][j];
        invpi[i][npi[i][j]] = j;
        C[i][j] = n.C[i][j];
      }
    }
    boost::multi_array<unsigned,2> bs = npi; 
    boost::multi_array<unsigned,2> bC = C;   
    Time bf = n.ms;                          
    Time pf = n.ms;                          
    for (unsigned i = I.m-1; i >= 1; i--) { 
      mi[i][I.n+1] = infinite_time;
      for (unsigned j = I.n; j >= 1; j--) {
        mi[i][j] = std::min(mi[i][j+1], invpi[i+1][npi[i][j]]);
      }
    }
    for (unsigned j = 0; j < I.n+2; j++) { 
      mi[I.m][j] = 0;
    }
    for (unsigned k = 1; k < I.n; k+=2) {
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned csum = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          csum += C[I.m][j];
        }
        if (csum < bf) {
          bf = csum;
          bs = npi;
          bC = C;
        }
      }
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned csum = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          csum += C[I.m][j];
        }
        if (csum < bf) {
          bf = csum;
          bs = npi;
          bC = C;
        }
      }
    }
    for (unsigned k = 2; k < I.n-1; k+=2) {
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned csum = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          csum += C[I.m][j];
        }
        if (csum < bf) {
          bf = csum;
          bs = npi;
          bC = C;
        }
      }
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned csum = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          csum += C[I.m][j];
        }
        if (csum < bf) {
          bf = csum;
          bs = npi;
          bC = C;
        }
      }
    }
    if (bf < pf) {
      improved = true;
      n.pi = bs;
      n.C = bC;
      n.ms = bf;
    }
    else {
      improved = false;
    }
  }
}
inline void ARNASLS(const Instance& I, NPFSSolution& n, Params& params) {
  unsigned n_ls = params.n_ls_np;
  boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
  boost::multi_array<Time,2> C(boost::extents[I.m+1][I.n+1]);
  boost::multi_array<unsigned,2> mi(boost::extents[I.m+1][I.n+2]);
  boost::multi_array<bool,2> crit(boost::extents[I.m+1][I.n+1]);
  boost::multi_array<unsigned,2> bs(boost::extents[I.m+1][I.n+1]); 
  boost::multi_array<unsigned,2> bC(boost::extents[I.m+1][I.n+1]); 
  for (unsigned i = 1; i <= I.m; i++) {
    Time c = 0;
    for (unsigned j = 1; j <= I.n; j++) {
      c = std::max(c, n.C[i-1][n.pi[i][j]]) + I.p[n.pi[i][j]][i];
      n.C[i][n.pi[i][j]] = c; 
    }
  }
  for (unsigned i = 1; i <= I.m; i++) {
    npi[i][0] = 0;
  }
  for (unsigned j = 1; j <= I.n; j++) {
    C[0][j] = 0;
  }
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = 1; j <= I.n; j++) {
        npi[i][j] = n.pi[i][j];
        invpi[i][npi[i][j]] = j;
        C[i][j] = n.C[i][j];
      }
    }
    bs = npi;       
    bC = C;         
    Time bf = n.ms; 
    Time pf = n.ms; 
    for (unsigned i = I.m-1; i >= 1; i--) { 
      mi[i][I.n+1] = infinite_time;
      for (unsigned j = I.n; j >= 1; j--) {
        mi[i][j] = std::min(mi[i][j+1], invpi[i+1][npi[i][j]]);
      }
    }
    for (unsigned j = 0; j < I.n+2; j++) { 
      mi[I.m][j] = 0;
    }
    for (unsigned l = 4; l <= I.n; l++) {
      crit[I.m][n.pi[I.m][l]] = true;
      for(unsigned i = I.m; i >= 1; i--) {
        for(unsigned j = l; j >= 1; j--) {
          unsigned job = npi[i][j];
          crit[i][npi[i][j-1]] = false;
          crit[i-1][job] = false;
          if (crit[i][job] && C[i][job]-I.p[job][i] == C[i][npi[i][j-1]])
            crit[i][npi[i][j-1]] = true;
          if (crit[i][job] && C[i][job]-I.p[job][i] == C[i-1][job])
            crit[i-1][job] = true;
        }
      }
      std::vector<std::pair<unsigned,unsigned>> r;
      crit[1][0] = true; 
      unsigned col = 1;
      for (unsigned i = 1; i <= I.m; i++) {
        for (unsigned j = col; j < I.n-1; j++) {
          if (crit[i][npi[i][j]] && crit[i][npi[i][j+1]] && !crit[i][npi[i][j-1]]) { 
            r.push_back(std::make_pair(j, j+1));
          }
          else if (crit[i][npi[i][j]] && crit[i][npi[i][j+1]] && !crit[i][npi[i][j+2]]) { 
            r.push_back(std::make_pair(j, j+1));
            col = j;
            break;
          }
        }
      }
      for (unsigned k = 0; k < r.size(); k++) {
        unsigned fj = r[k].first;
        unsigned sj = r[k].second;
        for (unsigned i = 1; i <= I.m; i++) {
          std::swap(invpi[i][npi[i][fj]], invpi[i][npi[i][sj]]);
          std::swap(npi[i][fj], npi[i][sj]);
          unsigned pos = fj;
          for (unsigned ii = i; ii <= I.m; ii++) {
            Time ct = C[ii][npi[ii][pos-1]];
            for(unsigned ji = pos; ji <= I.n; ji++) {
              unsigned job = npi[ii][ji];
              ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
              C[ii][job] = ct;
            }
            pos = mi[ii][pos];
          }
          unsigned csum = 0;
          for (unsigned j = 1; j <= I.n; j++) {
            csum += C[I.m][j];
          }
          if (csum < bf) {
            bf = csum;
            bs = npi;
            bC = C;
          }
        }
        for (unsigned i = 1; i <= I.m; i++) {
          std::swap(invpi[i][npi[i][fj]], invpi[i][npi[i][sj]]);
          std::swap(npi[i][fj], npi[i][sj]);
          unsigned pos = fj;
          for (unsigned ii = i; ii <= I.m; ii++) {
            Time ct = C[ii][npi[ii][pos-1]];
            for(unsigned ji = pos; ji <= I.n; ji++) {
              unsigned job = npi[ii][ji];
              ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
              C[ii][job] = ct;
            }
            pos = mi[ii][pos];
          }
          unsigned csum = 0;
          for (unsigned j = 1; j <= I.n; j++) {
            csum += C[I.m][j];
          }
          if (csum < bf) {
            bf = csum;
            bs = npi;
            bC = C;
          }
        }
      }
      if (bf < pf) {
        improved = true;
        n.pi = bs;
        n.C = bC;
        n.ms = bf;
      }
      else {
        improved = false;
      }
    }
  }
}
inline void AASLS_TT(const Instance& I, NPFSSolution& n, Params& params) {
  unsigned n_ls = params.n_ls_np;
  boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
  boost::multi_array<Time,2> C(boost::extents[I.m+1][I.n+1]);
  boost::multi_array<unsigned,2> mi(boost::extents[I.m+1][I.n+2]);
  for (unsigned i = 1; i <= I.m; i++) {
    Time c = 0;
    for (unsigned j = 1; j <= I.n; j++) {
      c = std::max(c, n.C[i-1][n.pi[i][j]]) + I.p[n.pi[i][j]][i];
      n.C[i][n.pi[i][j]] = c; 
    }
  }
  for (unsigned i = 1; i <= I.m; i++) {
    npi[i][0] = 0;
  }
  for (unsigned j = 1; j <= I.n; j++) {
    C[0][j] = 0;
  }
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = 1; j <= I.n; j++) {
        npi[i][j] = n.pi[i][j];
        invpi[i][npi[i][j]] = j;
        C[i][j] = n.C[i][j];
      }
    }
    boost::multi_array<unsigned,2> bs = npi; 
    boost::multi_array<unsigned,2> bC = C;   
    Time bf = n.ms;                          
    Time pf = n.ms;                          
    for (unsigned i = I.m-1; i >= 1; i--) { 
      mi[i][I.n+1] = infinite_time;
      for (unsigned j = I.n; j >= 1; j--) {
        mi[i][j] = std::min(mi[i][j+1], invpi[i+1][npi[i][j]]);
      }
    }
    for (unsigned j = 0; j < I.n+2; j++) { 
      mi[I.m][j] = 0;
    }
    for (unsigned k = 1; k < I.n; k++) {
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned tt = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          int tard = C[I.m][npi[I.m][j]] - I.dd[npi[I.m][j]];
          tt += std::max(0, tard);
        }
        if (tt < bf) {
          bf = tt;
          bs = npi;
          bC = C;
        }
      }
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned tt = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          int tard = C[I.m][npi[I.m][j]] - I.dd[npi[I.m][j]];
          tt += std::max(0, tard);
        }
        if (tt < bf) {
          bf = tt;
          bs = npi;
          bC = C;
        }
      }
    }
    if (bf < pf) {
      improved = true;
      n.pi = bs;
      n.C = bC;
      n.ms = bf;
    }
    else {
      improved = false;
    }
  }
}
inline void AASLS_r_TT(const Instance& I, NPFSSolution& n, Params& params) {
  unsigned n_ls = params.n_ls_np;
  boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
  boost::multi_array<Time,2> C(boost::extents[I.m+1][I.n+1]);
  boost::multi_array<unsigned,2> mi(boost::extents[I.m+1][I.n+2]);
  for (unsigned i = 1; i <= I.m; i++) {
    Time c = 0;
    for (unsigned j = 1; j <= I.n; j++) {
      c = std::max(c, n.C[i-1][n.pi[i][j]]) + I.p[n.pi[i][j]][i];
      n.C[i][n.pi[i][j]] = c; 
    }
  }
  for (unsigned i = 1; i <= I.m; i++) {
    npi[i][0] = 0;
  }
  for (unsigned j = 1; j <= I.n; j++) {
    C[0][j] = 0;
  }
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = 1; j <= I.n; j++) {
        npi[i][j] = n.pi[i][j];
        invpi[i][npi[i][j]] = j;
        C[i][j] = n.C[i][j];
      }
    }
    boost::multi_array<unsigned,2> bs = npi; 
    boost::multi_array<unsigned,2> bC = C;   
    Time bf = n.ms;                          
    Time pf = n.ms;                          
    for (unsigned i = I.m-1; i >= 1; i--) { 
      mi[i][I.n+1] = infinite_time;
      for (unsigned j = I.n; j >= 1; j--) {
        mi[i][j] = std::min(mi[i][j+1], invpi[i+1][npi[i][j]]);
      }
    }
    for (unsigned j = 0; j < I.n+2; j++) { 
      mi[I.m][j] = 0;
    }
    for (unsigned k = I.n-1; k >= 1; k--) {
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned tt = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          int tard = C[I.m][npi[I.m][j]] - I.dd[npi[I.m][j]];
          tt += std::max(0, tard);
        }
        if (tt < bf) {
          bf = tt;
          bs = npi;
          bC = C;
        }
      }
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned tt = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          int tard = C[I.m][npi[I.m][j]] - I.dd[npi[I.m][j]];
          tt += std::max(0, tard);
        }
        if (tt < bf) {
          bf = tt;
          bs = npi;
          bC = C;
        }
      }
    }
    if (bf < pf) {
      improved = true;
      n.pi = bs;
      n.C = bC;
      n.ms = bf;
    }
    else {
      improved = false;
    }
  }
}
inline void AASLS_G8_TT(const Instance& I, NPFSSolution& n, Params& params) {
  unsigned n_ls = params.n_ls_np;
  boost::multi_array<unsigned,2> npi(boost::extents[I.m+1][I.n+1]), invpi(boost::extents[I.m+2][I.n+1]);
  boost::multi_array<Time,2> C(boost::extents[I.m+1][I.n+1]);
  boost::multi_array<unsigned,2> mi(boost::extents[I.m+1][I.n+2]);
  for (unsigned i = 1; i <= I.m; i++) {
    Time c = 0;
    for (unsigned j = 1; j <= I.n; j++) {
      c = std::max(c, n.C[i-1][n.pi[i][j]]) + I.p[n.pi[i][j]][i];
      n.C[i][n.pi[i][j]] = c; 
    }
  }
  for (unsigned i = 1; i <= I.m; i++) {
    npi[i][0] = 0;
  }
  for (unsigned j = 1; j <= I.n; j++) {
    C[0][j] = 0;
  }
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    for (unsigned i = 1; i <= I.m; i++) {
      for (unsigned j = 1; j <= I.n; j++) {
        npi[i][j] = n.pi[i][j];
        invpi[i][npi[i][j]] = j;
        C[i][j] = n.C[i][j];
      }
    }
    boost::multi_array<unsigned,2> bs = npi; 
    boost::multi_array<unsigned,2> bC = C;   
    Time bf = n.ms;                          
    Time pf = n.ms;                          
    for (unsigned i = I.m-1; i >= 1; i--) { 
      mi[i][I.n+1] = infinite_time;
      for (unsigned j = I.n; j >= 1; j--) {
        mi[i][j] = std::min(mi[i][j+1], invpi[i+1][npi[i][j]]);
      }
    }
    for (unsigned j = 0; j < I.n+2; j++) { 
      mi[I.m][j] = 0;
    }
    for (unsigned k = 1; k < I.n; k+=2) {
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned tt = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          int tard = C[I.m][npi[I.m][j]] - I.dd[npi[I.m][j]];
          tt += std::max(0, tard);
        }
        if (tt < bf) {
          bf = tt;
          bs = npi;
          bC = C;
        }
      }
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned tt = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          int tard = C[I.m][npi[I.m][j]] - I.dd[npi[I.m][j]];
          tt += std::max(0, tard);
        }
        if (tt < bf) {
          bf = tt;
          bs = npi;
          bC = C;
        }
      }
    }
    for (unsigned k = 2; k < I.n-1; k+=2) {
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned tt = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          int tard = C[I.m][npi[I.m][j]] - I.dd[npi[I.m][j]];
          tt += std::max(0, tard);
        }
        if (tt < bf) {
          bf = tt;
          bs = npi;
          bC = C;
        }
      }
      for (unsigned i = 1; i <= I.m; i++) {
        std::swap(invpi[i][npi[i][k]], invpi[i][npi[i][k+1]]);
        std::swap(npi[i][k], npi[i][k+1]);
        unsigned pos = k;
        for (unsigned ii = i; ii <= I.m; ii++) {
          Time ct = C[ii][npi[ii][pos-1]];
          for(unsigned ji = pos; ji <= I.n; ji++) {
            unsigned job = npi[ii][ji];
            ct = std::max(ct, C[ii-1][job]) + I.p[job][ii];
            C[ii][job] = ct;
          }
          pos = mi[ii][pos];
        }
        unsigned tt = 0;
        for (unsigned j = 1; j <= I.n; j++) {
          int tard = C[I.m][npi[I.m][j]] - I.dd[npi[I.m][j]];
          tt += std::max(0, tard);
        }
        if (tt < bf) {
          bf = tt;
          bs = npi;
          bC = C;
        }
      }
    }
    if (bf < pf) {
      improved = true;
      n.pi = bs;
      n.C = bC;
      n.ms = bf;
    }
    else {
      improved = false;
    }
  }
}
inline void insertionTT_IT1_np(const Instance& I, NPFSSolution &n, Params& params) {
  unsigned n_ls = params.n_ls_np;
  using namespace std::placeholders;
  nBuilder construct(I);
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    unsigned cms = n.ms;
    for(unsigned ji = I.n; ji >= 1; ji--) {
      Job job = n.pi[1][ji];
      for(unsigned i = 1; i <= I.m; i++) {
         unsigned dj = 1;
         for(unsigned j = 1; j <= I.n; j++) {
            if (n.pi[i][j] != job) {
               n.pi[i][dj] = n.pi[i][j];
               dj++;
            }
         }
         while (dj <= I.n) {
            n.pi[i][dj] = 0;
            dj++;
         }
      }
      std::vector<Job> remove = { 0, job };
      construct.constructnpNEWTT_IT1(n, I.n, remove);
    }
    if (n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
inline void insertionTT_IT2_np(const Instance& I, NPFSSolution &n, Params& params) {
  unsigned n_ls = params.n_ls_np;
  using namespace std::placeholders;
  nBuilder construct(I);
  bool improved = true;
  while (improved && n_ls-- > 0 && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-params.start).count() < params.timelimit_np) {
    unsigned cms = n.ms;
    for(unsigned ji = I.n; ji >= 1; ji--) {
      Job job = n.pi[1][ji];
      for(unsigned i = 1; i <= I.m; i++) {
         unsigned dj = 1;
         for(unsigned j = 1; j <= I.n; j++) {
            if (n.pi[i][j] != job) {
               n.pi[i][dj] = n.pi[i][j];
               dj++;
            }
         }
         while (dj <= I.n) {
            n.pi[i][dj] = 0;
            dj++;
         }
      }
      std::vector<Job> remove = { 0, job };
      construct.constructnpNEWTT_IT2(n, I.n, remove);
    }
    if (n.ms < cms)
      improved = true;
    else
      improved = false;
  }
}
