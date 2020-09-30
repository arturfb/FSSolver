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
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <iomanip>
using namespace std;
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace boost;
#include "constructive.hpp"
#include "threemachine.hpp"
unsigned a = 9, b = 3, c = 7;
float a_tt = 0.0, b_tt = 0.15, c_tt = 1.25, e_tt = 4.0;
struct beamnode : m3::psolution {
  typedef m3::psolution Base;
  vector<Time> C;         
  vector<Time> Sp;        
  float SIT;              
  float SCT;              
  beamnode(const Instance& I) : Base(I.n), C(I.m+1,0), Sp(I.m+1,0), SIT(0), SCT(0) {
    for(unsigned i=1; i<=I.m; i++)
      for(unsigned j=1; j<=I.n; j++)
        Sp[i] += I.p[j][i];
  }
  beamnode(beamnode&& other) : Base(other) {
    this->swap(other);
  }
  beamnode(const beamnode& other) : Base(other) {
    C = other.C;
    Sp  = other.Sp;
    SIT = other.SIT;
    SCT = other.SCT;
  }
  beamnode& operator=(beamnode other) {
    this->swap(other);
    return *this;
  }
  void swap(beamnode& other) {
    using std::swap;
    Base::swap(other);
    C.swap(other.C);
    swap(Sp,other.Sp);
    swap(SIT,other.SIT);
    swap(SCT,other.SCT);
  }
  void addi(unsigned ji, const Instance& I) {
    Job j = pi[ji];
    float IT = 0.0;
    for(unsigned i=1; i<=I.m; i++)
      if (C[i-1]<=C[i])
        C[i] += I.p[j][i];
      else {
        IT += float(I.m*(C[i-1]-C[i]))/float(i-1+numFixed()*(float(I.m-i+1)/float(I.n-2)));
        C[i] = C[i-1] + I.p[j][i];
      }
    for(unsigned i=1; i<=I.m; i++)
      Sp[i] -= I.p[j][i];
    float Ca = 0;
    if (numFree()>0)
      for(unsigned i=1; i<=I.m; i++)
        Ca=std::max(Ca,float(C[i]))+float(Sp[i])/float(numFree());
    SCT += C[I.m] + Ca;
    SIT = (1.0f-float(b)/float(I.n))*(SIT+IT*(I.n-numFixed()-2));
    Base::addi(ji);
  }
  float Bvalue(unsigned ji, const Instance& I) {
    Job j = pi[ji];
    float IT = 0.0;
    Time Cj = 0;
    for(unsigned i=1; i<=I.m; i++)
      if (Cj<=C[i])
        Cj = C[i]+I.p[j][i];
      else {
        IT += float(I.m*(Cj-C[i]))/float(i-1+numFixed()*(float(I.m-i+1)/float(I.n-2)));
        Cj += I.p[j][i];
      }
    return a*SCT+SIT+c*Cj+IT*(I.n-numFixed()-2);
  }
  friend std::ostream& operator<<(std::ostream& o, const beamnode& n);
};
std::ostream& operator<<(std::ostream& o, const beamnode& n) {
  o << n.numFixed() << "(";
  if (n.full())
    copy(n.pi.begin()+1,n.pi.end(),std::ostream_iterator<Job>(o,","));
  else {
    copy(n.pi.begin()+1,n.pi.begin()+n.n1+1,std::ostream_iterator<Job>(o,","));
    o << "|";
    copy(n.pi.begin()+n.n1+1,n.pi.begin()+n.n2,std::ostream_iterator<Job>(o,","));
    o << "|";
    copy(n.pi.begin()+n.n2,n.pi.end(),std::ostream_iterator<Job>(o,","));
  }
  return o << "] " << n.n1;
}
void beam_search(const Instance& I, PFSSolution& cs, unsigned w) {
  typedef tuple<float,unsigned,unsigned> beamext;
  vector<beamext> beamcand; 
  vector<beamnode> beam(w,I), newbeam(w,I);
  for(unsigned ji=1; ji<=I.n; ji++) {
    Job j = ji; 
    float IT = 0.0;
    Time Cj = 0;
    for(unsigned i=1; i<=I.m; i++) {
      Cj += I.p[j][i];
      if (i<I.m)
        IT += float(Cj)/i; 
    }
    IT *= I.m;
    float xi = float(I.n-2)/4.0*IT+float(Cj);
    beamext be = make_tuple(xi,Cj,ji);
    auto ip = lower_bound(beamcand.begin(),beamcand.end(),be);
    beamcand.insert(ip,be);
    if (beamcand.size()>w)
      beamcand.pop_back();
  }
  assert(beamcand.size()==w);
  for(unsigned i=0; i<w; i++)
    beam[i].addi(get<2>(beamcand[i]),I);
  for(unsigned k=2; k<=I.n; k++) {
    beamcand.clear();
    for(unsigned i=0; i<w; i++) {
      for(unsigned ji=beam[i].n1+1; ji<beam[i].n2; ji++) {
        float B = beam[i].Bvalue(ji,I);
        if (beamcand.size()<w || B<get<0>(beamcand.back())) {
          beamext be = make_tuple(B,i,ji);
          auto ip = lower_bound(beamcand.begin(),beamcand.end(),be);
          beamcand.insert(ip,be);
          if (beamcand.size()>w)
            beamcand.pop_back();
        }
      }
    }
    assert(beamcand.size()==w);
    for(unsigned i=0; i<w; i++) {
      newbeam[i]=beam[get<1>(beamcand[i])];
      newbeam[i].addi(get<2>(beamcand[i]),I);
    }
    beam.swap(newbeam);
  }
  cs.ms=infinite_time;
  for(unsigned i=0; i<w; i++) {
    PFSSolution ns(I);
    ns.pi = beam[i].pi;
    ns.computeFlowtime(I);
    if (ns.ms<cs.ms)
      cs = ns;
  }
}
struct beamnode_tt : m3::psolution {
  typedef m3::psolution Base;
  vector<Time> C; 
  unsigned TT; 
  unsigned TE; 
  float    TI; 
  unsigned TW; 
  beamnode_tt(const Instance& I) : Base(I.n), C(I.m+1,0), TT(0), TE(0), TI(0.0), TW(0) {
  }
  beamnode_tt(beamnode_tt&& other) : Base(other) {
    this->swap(other);
  }
  beamnode_tt(const beamnode_tt& other) : Base(other) {
    C = other.C;
    TT = other.TT;
    TE = other.TE;
    TI = other.TI;
    TW = other.TW;
  }
  beamnode_tt& operator=(beamnode_tt other) {
    this->swap(other);
    return *this;
  }
  void swap(beamnode_tt& other) {
    using std::swap;
    Base::swap(other);
    C.swap(other.C);
    swap(TT,other.TT);
    swap(TE,other.TE);
    swap(TI,other.TI);
    swap(TW,other.TW);
  }
  void addi(unsigned ji, const Instance& I) {
    const Job j = pi[ji];
    const unsigned K = numFixed()+1;
    unsigned Cj = 0;
    float it = 0.0;
    for(unsigned i = 1; i <= I.m; i++) {
      if(i > 1)
        it += float(I.m*std::max(0, int(Cj-C[i]))) / float(i-1+(K-1)*float(I.m-i+1)/float(I.n-2));
      Cj = std::max(Cj, C[i]) + I.p[j][i];
      C[i] = Cj;
    }
    if(K>1) {
    TT += std::max(0, int(Cj-I.dd[j]));
    TE += std::max(0, int(I.dd[j]-Cj));
    TI += it;
    }
    Base::addi(ji);
    TW = 0;
    for(unsigned k = n1+1; k < n2; k++) {
      Cj = 0;
      for(unsigned i = 1; i <= I.m; i++)
        Cj = std::max(Cj,C[i]) + I.p[pi[k]][i];
      TW += std::max(0,int(Cj-I.dd[pi[k]]));
    }
  }
  float Gvalue(unsigned ji, const Instance& I) {
    const Job j = pi[ji];
    const unsigned K = numFixed()+1;
    unsigned Cj = 0; 
    float it = 0.0;
    for(unsigned i = 1; i <= I.m; i++) {
      if(i > 1) 
        it += float(I.m*std::max(0, int(Cj-C[i]))) / float(i-1+(K-1)*float(I.m-i+1)/float(I.n-2));
      Cj = std::max(Cj, C[i]) + I.p[j][i];
    }
    float F1 = TI*(float(I.n-K-1)/float(I.n));
    float F2 = a_tt*TE*(float(2*I.n-K-1)/float(2*I.n));
    float F3 = b_tt*TT*(float(K-1+I.n)/float(2*I.n));
    float L1 = float(I.n-K-1)*it;
    float L2 = c_tt*std::max(0, int(I.dd[j]-Cj));
    float L3 = (e_tt/float(I.n-K+1))*TW;
    float F = F1 + F2 + F3;
    float L = L1 + L2 + L3;
    return F+L;
  }
  friend std::ostream& operator<<(std::ostream& o, const beamnode_tt& n);
};
void beam_search_tt(const Instance& I, PFSSolution& cs, unsigned w) {
  typedef tuple<float,unsigned,unsigned> beamext;
  vector<beamext> beamcand; 
  vector<beamnode_tt> beam(w,I), newbeam(w,I);
  float bej = std::numeric_limits<float>::max(), bw = std::numeric_limits<float>::max();
  unsigned bji = 1;
  for(unsigned ji = 1; ji <= I.n; ji++) {
    Job j = ji; 
    Time tpt = 0;    
    float wit = 0.0; 
    for(unsigned i = 1; i <= I.m; i++) {
      if(i > 1)
        wit += float(tpt*I.m)/(i-1);
      tpt += I.p[j][i];
    }
    wit *= float(I.n-2)/4.0;
    float ej = tpt + wit;
    if(ej < bej || (ej == bej && wit < bw)) { 
      bej = ej;
      bw = wit;
      bji = ji;
    }
  }
  beam[0].addi(bji, I);
  for(unsigned ji = beam[0].n1+1; ji < beam[0].n2; ji++) {
    float G = beam[0].Gvalue(ji, I);
    if (beamcand.size() < w || G < get<0>(beamcand.back())) {
      beamext be = make_tuple(G, 0, ji);
      auto ip = lower_bound(beamcand.begin(), beamcand.end(), be);
      beamcand.insert(ip, be);
      if (beamcand.size() > w)
        beamcand.pop_back();
    }
  }
  assert(beamcand.size() == w);
  for(unsigned i = 0; i < w; i++) {
    newbeam[i] = beam[get<1>(beamcand[i])];
    newbeam[i].addi(get<2>(beamcand[i]),I);
  }
  beam.swap(newbeam);
  for(unsigned k = 3; k <= I.n-1; k++) {
    beamcand.clear();
    for(unsigned i = 0; i < w; i++) {
      for(unsigned ji = beam[i].n1+1; ji < beam[i].n2; ji++) {
        float G = beam[i].Gvalue(ji, I);
        if (beamcand.size() < w || G < get<0>(beamcand.back())) {
          beamext be = make_tuple(G, i, ji);
          auto ip = lower_bound(beamcand.begin(), beamcand.end(), be);
          beamcand.insert(ip, be);
          if (beamcand.size() > w)
            beamcand.pop_back();
        }
      }
    }
    assert(beamcand.size()==w);
    for(unsigned i = 0; i < w; i++) {
      newbeam[i] = beam[get<1>(beamcand[i])];
      newbeam[i].addi(get<2>(beamcand[i]),I);
    }
    beam.swap(newbeam);
  }
  cs.ms = infinite_time;
  for(unsigned i = 0; i < w; i++) {
    PFSSolution ns(I);
    ns.pi = beam[i].pi;
    ns.computeTotalTardiness(I);
    if (ns.ms < cs.ms)
      cs = ns;
  }
}

