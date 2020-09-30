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
#include <numeric>
#include <array>
#include <utility>
#include <boost/multi_array.hpp>
#include "units.hpp"
#include "instance.hpp"
struct Solution {
  Solution() : n(0), m(0), ms(infinite_time) {}
  Solution(const Instance& I) : n(I.n), m(I.m), ms(infinite_time) {}
  unsigned n,m;
  Time ms;
  Solution(Solution&& other) {
    this->swap(other);
  }
  Solution(const Solution& other) {
    n = other.n;
    m = other.m;
    ms = other.ms;
  }
  Solution& operator=(Solution other) {
    this->swap(other);
    return *this;
  }
  void swap(Solution& other) {
    using std::swap;
    swap(n,other.n);
    swap(m,other.m);
    swap(ms,other.ms);
  }
};
struct PFSSolution : public Solution {
  typedef Solution Base;
  std::vector<Job> pi;
  PFSSolution(const Instance& I) : Base(I), pi(I.n+1,0) {
    std::iota(pi.begin(),pi.end(),0);
  }
  PFSSolution(const Instance& I, const std::vector<Job>& _pi) : Base(I), pi(_pi) {}
  PFSSolution(PFSSolution&& other) {
    this->swap(other);
  }
  PFSSolution(const PFSSolution& other) : Base(other) {
    pi = other.pi;
  }
  PFSSolution& operator=(PFSSolution other) {
    this->swap(other);
    return *this;
  }
  void read(std::istream& in, bool zeroBased = true) {
    for(unsigned j=1; j<=n; j++) {
      in >> pi[j];
      if (zeroBased)
        pi[j]++;
      if (pi[j] < 1 || pi[j]>n) {
        std::cerr << "Error reading solution: invalid input: job " << pi[j] << " does not exist in a " << n << "-job instance." << std::endl;
        exit(1);
      }
    }
  }
  void swap(PFSSolution& other) {
    Base::swap(other);
    using std::swap;
    pi.swap(other.pi);
  }
  void computeMakespan(const Instance& I, std::vector<Time>& C) {
    assert(I.n == n && I.m == m && C.size()==m+1);
    fill(C.begin(),C.end(),0);
    const unsigned pn = pi.size()-1;
    for(unsigned j=1; j<=pn; j++)
      for(unsigned i=1; i<=m; i++)
        if (C[i]>C[i-1])
          C[i]+=I.p[pi[j]][i];
        else
          C[i]=C[i-1]+I.p[pi[j]][i];
    ms=C[m];
  }
  void computeMakespan(const Instance& I, std::vector<Time>& C, const unsigned end) {
    assert(I.n == n && I.m == m && C.size()==m+1);
    fill(C.begin(),C.end(),0);
    const unsigned pn = end-1;
    for(unsigned j=1; j<=pn; j++)
      for(unsigned i=1; i<=m; i++)
        if (C[i]>C[i-1])
          C[i]+=I.p[pi[j]][i];
        else
          C[i]=C[i-1]+I.p[pi[j]][i];
    ms=C[m];
  }
  void computeMakespan(const Instance& I) {
    std::vector<Time> C(I.m+1,0);
    computeMakespan(I,C);
    ms=C[m];
  }
  void computeMakespan(const Instance& I, const unsigned end) {
    std::vector<Time> C(I.m+1,0);
    computeMakespan(I,C,end);
    ms=C[m];
  }
  void reverse() {
    using std::swap;
    for(unsigned j=1; j<=n/2; j++)
      swap(pi[j],pi[n-j+1]);
  }
#define PROPFRST                       ct += *ctp++; *Cp-- = ct
#define PROPONCE if (*Cp>ct) ct = *Cp; ct += *ctp++; *Cp-- = ct
#define PROPLAST if (*Cp>ct) ct = *Cp; ct += *ctp  ; *Cp   = ct
  void computeFlowtime0(const Instance& I, std::array<Time,64> C, unsigned begin, unsigned end, unsigned ub = infinite_time) {
    const Time *ctp;
    Time ct;
    auto Cp = C.begin();
    int left = end-begin;
    if (m==5) {
      while (left-->0) {
        ctp = &I.p[pi[begin]][1];
        Cp = &C[5];
        ct = *Cp;
        PROPFRST; PROPONCE; PROPONCE; PROPONCE; PROPLAST;
        ms+=ct;
        begin++;
      }
    } else if (m==10) {
      while (left-->0) {
        ctp = &I.p[pi[begin]][1];
        Cp = &C[10];
        ct = *Cp;
        PROPFRST; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPLAST;
        ms+=ct;
        begin++;
      }
    } else if (m==20) {
      while (left-->0) {
        ctp = &I.p[pi[begin]][1];
        Cp = &C[20];
        ct = *Cp;
        PROPFRST; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE;
        PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPLAST;
        ms+=ct;
        begin++;
      }
    } else {
      while (left-->0) {
        ctp = &I.p[pi[begin]][1];
        Cp = &C[m];
        ct = *Cp;
        for(unsigned i=m; i>0; i--) {
          PROPONCE;
        }
        ms+=ct;
        begin++;
      }
    }
  }
  void computeTotalTardiness(const Instance& I, std::array<Time,64> C, unsigned begin, unsigned end, unsigned ub = infinite_time) {
    const Time *ctp;
    Time ct;
    auto Cp = C.begin();
    int left = end-begin;
    if (m==5) {
      while (left-->0) {
        ctp = &I.p[pi[begin]][1];
        Cp = &C[5];
        ct = *Cp;
        PROPFRST; PROPONCE; PROPONCE; PROPONCE; PROPLAST;
        int tard = int(ct-I.dd[pi[begin]]);
        ms+=std::max(0,tard);
        begin++;
      }
    } else if (m==10) {
      while (left-->0) {
        ctp = &I.p[pi[begin]][1];
        Cp = &C[10];
        ct = *Cp;
        PROPFRST; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPLAST;
        int tard = int(ct-I.dd[pi[begin]]);
        ms+=std::max(0,tard);
        begin++;
      }
    } else if (m==20) {
      while (left-->0) {
        ctp = &I.p[pi[begin]][1];
        Cp = &C[20];
        ct = *Cp;
        PROPFRST; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE;
        PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPONCE; PROPLAST;
        int tard = int(ct-I.dd[pi[begin]]);
        ms+=std::max(0,tard);
        begin++;
      }
    } else {
      while (left-->0) {
        ctp = &I.p[pi[begin]][1];
        Cp = &C[m];
        ct = *Cp;
        for(unsigned i=m; i>0; i--) {
          PROPONCE;
        }
        int tard = int(ct-I.dd[pi[begin]]);
        ms+=std::max(0,tard);
        begin++;
      }
    }
  }
#undef PROPFRST
#undef PROPONCE
#define PROPFRST                       ct += *(ctp + *cpp++); *Cp++ = ct
#define PROPONCE if (*Cp>ct) ct = *Cp; ct += *(ctp + *cpp++); *Cp++ = ct
  void computeFlowtime1(const Instance& I, std::array<Time,64> C0, unsigned begin, unsigned end, unsigned ub = infinite_time) {
    const Time *ctp;
    std::vector<Job>::const_iterator cpp;
    Time ct;
    int left;
    std::array<Time,512> C;
    std::fill(C.begin()+begin,C.begin()+end,0);
    auto Cp = C.begin();
    Time ms0 = ms;
    cpp = pi.begin()+begin;
    ctp = &I.q[1][0];
    left = end-begin;
    Cp = &C[begin];
    ct = C0[m];
    while (left>19) {
      PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct;
      PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct;
      PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct;
      PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct; PROPFRST; ms += ct;
      left -= 20;
    }
    while (left-->0) {
      PROPFRST; ms += ct;
    }
    if (ms>ub)
      return;
    for(unsigned i=m-1; i>0; i--) {
      ms = ms0;
      cpp = pi.begin()+begin;
      ctp = &I.q[m-i+1][0];
      left = end-begin;
      Cp = &C[begin];
      ct = C0[i];
      while (left>19) {
       PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct;
       PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct;
       PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct;
       PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct; PROPONCE; ms += ct;
       left -= 20;
     }
     while (left-->0) {
       PROPONCE;
       ms += ct;
     }
     if (ms>ub)
       return;
    }
    ms = ms0;
    for(unsigned j=begin; j<end; j++)
    ms += C[j];
  }
  void computeFlowtime(const Instance& I, unsigned end) {
    assert(I.n == n && I.m == m);
    std::array<Time,64> C;
    std::fill(C.begin(),C.end(),0);
    ms = 0;
    computeFlowtime0(I,C,1,end);
  }
  void computeFlowtime(const Instance& I) {
    computeFlowtime(I,pi.size());
  }
  void computeTotalTardiness(const Instance& I, unsigned end) {
    assert(I.n == n && I.m == m);
    std::array<Time,64> C;
    std::fill(C.begin(),C.end(),0);
    ms = 0;
    computeTotalTardiness(I,C,1,end);
  }
  void computeTotalTardiness(const Instance& I) {
    computeTotalTardiness(I,pi.size());
  }
  void write(std::ostream& o) {
    copy(pi.begin()+1,pi.end(),std::ostream_iterator<Job>(o," "));
    o << std::endl;
  }
  void removeJob(const unsigned job) {
    auto nend = pi.end();
    nend = std::remove(pi.begin()+1, nend, job);
    std::fill(nend, pi.end(), 0);
  }
};
inline std::ostream& operator<<(std::ostream& o, const PFSSolution& s) {
  o << "<";
  std::copy(s.pi.begin()+1,s.pi.end(),std::ostream_iterator<Job>(o," "));
  return o << ";ms=" << s.ms << ">";
}
struct NPFSSolution : public Solution {
  typedef Solution Base;
  boost::multi_array<unsigned,2> pi;
  boost::multi_array<Time,2>     C;
  NPFSSolution(const Instance& I) : Base(I) {
    pi.resize(boost::extents[m+1][n+1]);
    C.resize(boost::extents[m+1][n+1]);
    for(unsigned i=1; i<=m; i++)
      for(unsigned j=1; j<=n; j++)
        pi[i][j]=j;
   }
  NPFSSolution(std::istream& in, const Instance& I) : Base(I) {
    using std::cout;
    pi.resize(boost::extents[m+1][n+1]);
    C.resize(boost::extents[m+1][n+1]);
    for(unsigned i=1; i<=m; i++) {
      cout << "Machine " << i << ": ";
      for(unsigned j=1; j<=n; j++) {
        in >> pi[i][j];
        if (in.eof()) {
          std::cerr << "Error: Short input" << std::endl;
          exit(1);
        }
        pi[i][j] = ((pi[i][j]-1)/m)+1;
        cout << pi[i][j] << " ";
      }
      cout << std::endl;
    }
  }
  NPFSSolution(const Instance& I, const boost::multi_array<unsigned,2>& opi) : Base(I) {
    pi.resize(boost::extents[m+1][n+1]);
    C.resize(boost::extents[m+1][n+1]);
    for(unsigned i=1; i<=m; i++)
      for(unsigned j=1; j<=n; j++)
        pi[i][j] = opi[i][j];
  }
  void read(std::istream& in, bool zeroBased = true) {
    for(unsigned i=1; i<=m; i++) {
      for(unsigned j=1; j<=n; j++) {
        in >> pi[i][j];
        if (in.eof()) {
          std::cerr << "Error: Short input" << std::endl;
          exit(1);
        }
        if (zeroBased)
          pi[i][j]++;
        if (pi[i][j] < 1 || pi[i][j]>n) {
          std::cerr << "Error: invalid input: job " << pi[i][j] << " does not exists in a " << n << "-job instance." << std::endl;
          exit(1);
        }
      }
    }
  }
  NPFSSolution(NPFSSolution&& other) {
    this->swap(other);
  }
  NPFSSolution(const NPFSSolution& other) : Base(other) {
    pi.resize(boost::extents[other.pi.shape()[0]][other.pi.shape()[1]]);
    pi = other.pi;
    C.resize(boost::extents[other.C.shape()[0]][other.C.shape()[1]]);
    C = other.C;
  }
  NPFSSolution(const Instance& I, const PFSSolution& other) : Base(other) {
    pi.resize(boost::extents[m+1][n+1]);
    for(unsigned i=1; i<=m; i++)
      for(unsigned j=1; j<=n; j++)
        pi[i][j]=other.pi[j];
     C.resize(boost::extents[m+1][n+1]);
   }
   NPFSSolution& operator=(NPFSSolution other) {
    this->swap(other);
    return *this;
  }
  void reverse() {
    using std::swap;
    for(unsigned i=1; i<=m/2; i++)
      for(unsigned j=1; j<=n; j++)
        swap(pi[i][j],pi[m-i+1][j]);
    for(unsigned i=1; i<=m; i++)
      for(unsigned j=1; j<=n/2; j++)
        swap(pi[i][j],pi[i][n-j+1]);
  }
  void swap(NPFSSolution& other) {
    Base::swap(other);
    using std::swap;
    pi.resize(boost::extents[other.pi.shape()[0]][other.pi.shape()[1]]);
    swap(pi,other.pi);
    C.resize(boost::extents[other.C.shape()[0]][other.C.shape()[1]]);
    swap(C,other.C);
  }
  void computeMakespan(const Instance& I) {
    assert(I.n == n && I.m == m);
    for(unsigned j=0; j<=n; j++)
      C[0][j]=0;
    computeMakespan(I, 1);
  }
  void computeMakespan(const Instance& I, const unsigned mi) {
    assert(I.n == n && I.m == m);
    for(unsigned i = mi; i <= m; i++) {
      Time c = 0;
      for(unsigned j = 1; j <= n; j++) {
        c = std::max(c, C[i-1][pi[i][j]]) + I.p[pi[i][j]][i];
        C[i][pi[i][j]] = c;
      }
    }
    ms = C[m][pi[m][n]];
  }
  void computeFlowtime(const Instance& I) {
    assert(I.n == n && I.m == m);
    for(unsigned j=0; j<=n; j++)
      C[0][j]=0;
    computeFlowtime(I, 1);
  }
  void computeFlowtime(const Instance& I, const unsigned mi) {
    assert(I.n == n && I.m == m);
    for(unsigned i = mi; i <= m; i++) {
      Time c = 0;
      for(unsigned j = 1; j <= n; j++) {
        c = std::max(c, C[i-1][pi[i][j]]) + I.p[pi[i][j]][i];
        C[i][pi[i][j]] = c;
      }
    }
    ms = 0;
    for(unsigned j = 1; j <= n; j++)
      ms+=C[m][j];
  }
  void computeTotalTardiness(const Instance& I) {
    assert(I.n == n && I.m == m);
    for(unsigned j=0; j<=n; j++)
      C[0][j]=0;
    computeTotalTardiness(I, 1);
  }
  void computeTotalTardiness(const Instance& I, const unsigned mi) {
    assert(I.n == n && I.m == m);
    for(unsigned i = mi; i <= m; i++) {
      Time c = 0;
      for(unsigned j = 1; j <= n; j++) {
        c = std::max(c, C[i-1][pi[i][j]]) + I.p[pi[i][j]][i];
        C[i][pi[i][j]] = c;
      }
    }
    ms = 0;
    for(unsigned j = 1; j <= n; j++) {
      int tard = C[m][j]-I.dd[pi[m][j]];
      ms+=std::max(0,tard);
    }
  }
  boost::multi_array<bool,2> criticalPath(const Instance& I) {
    boost::multi_array<bool,2> c;
    c.resize(boost::extents[m+1][n+1]);
    c[m][pi[m][n]]=true;
    for(unsigned i=m; i>=1; i--) {
      for(unsigned j=n; j>=1; j--) {
        unsigned tj = pi[i][j];
        if (c[i][tj] && C[i][tj]-I.p[tj][i] == C[i][pi[i][j-1]])
          c[i][pi[i][j-1]]=true;
        if (c[i][tj] && C[i][tj]-I.p[tj][i] == C[i-1][tj])
          c[i-1][tj]=true;
      }
    }
    return c;
  }
  void write(std::ostream& o) {
    for(unsigned i=1; i<=m; i++) {
      for(unsigned j=1; j<=n; j++)
       o << pi[i][j] << " ";
     o << std::endl;
   }
  }
  void removeJob(const unsigned job) {
    for(unsigned i = 1; i <= m; i++) {
      unsigned dj = 1;
      for(unsigned j = 1; j <= n; j++) {
        if (pi[i][j] != job) {
          pi[i][dj] = pi[i][j];
          dj++;
        }
      }
      while (dj <= n) {
        pi[i][dj] = 0;
        dj++;
      }
    }
  }
};
