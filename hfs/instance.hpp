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
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <boost/multi_array.hpp>
#include "units.hpp"
inline  Time getTime(std::istream& in) {
  std::string token;
  in >> token;
  if (token[0]=='i' || token[0]=='I')
    return infinite_time;
  std::stringstream num(token);
  Time t;
  num >> t;
  return t;
}
struct Instance {
  unsigned n; 
  unsigned m; 
  boost::multi_array<Time,2> p; 
  boost::multi_array<Time,2> q; 
  boost::multi_array<Time,3> l; 
  std::vector<unsigned> dd;     
  Instance(unsigned _n=0, unsigned _m=0) : n(_n), m(_m), p(boost::extents[n+1][m+2]), q(boost::extents[m+1][n+1]), l(boost::extents[n+1][m+1][m+1]), dd(n+1) {}
  Instance(std::istream& in) {
    read(in);
  }
  Instance(const Instance& I, unsigned i1, unsigned i2) : n(I.n), m(i2-i1+1), p(boost::extents[n+1][m+2]), q(boost::extents[m+1][n+1]), l(boost::extents[n+1][m+1][m+1]) {
    for(unsigned j=1; j<=n; j++)
      for(unsigned i=i1; i<=i2; i++)
	      p[j][i-i1+1] = q[i-i1+1][j] = I.p[j][i];
    for(unsigned j=1; j<=n; j++)
      for(unsigned i=i1; i<=i2; i++)
	      for(unsigned k=i1; k<=i2; k++)
	        l[j][i-i1+1][k-i1+1] = I.l[j][i][k];
  }
  void read(std::istream& in) {
    in >> n >> m;
    read_extended(in);
  }
  void read_extended(std::istream& in) {
    p.resize(boost::extents[n+1][m+2]);
    q.resize(boost::extents[m+1][n+1]);
    for(unsigned j=1; j<=n; j++)
      for(unsigned i=1; i<=m; i++) {
	      unsigned dummy;
	      in >> dummy;
	      assert(dummy==i-1 || dummy==i); 
	      p[j][i] = q[i][j] = getTime(in);
      }
    for(unsigned i=0; i<=m; i++)
      q[i][0]=0;
    for(unsigned j=0; j<=n; j++)
      p[j][0]=p[j][m+1]=0;
    std::fill(&p[0][0],&p[0][m+1],0);
    l.resize(boost::extents[n+1][m+1][m+1]);
    for(unsigned j=1; j<=n; j++)
      for(unsigned i=1; i<m; i++)
	      for(unsigned k=i+1; k<=m; k++) {
	        l[j][i][k]=0;
	        for(unsigned o=i+1; o<k; o++)
	          l[j][i][k]+=p[j][o];
	      }
  }
  void read_duedates(std::istream& in) {
    in >> n >> m;
    assert(n>0);
    assert(n<1000000);
    dd.resize(n+1); 
    unsigned dummy;
    std::string dummyStr;
    int dummyInt;
    p.resize(boost::extents[n+1][m+2]);
    q.resize(boost::extents[m+1][n+1]);
    for(unsigned j=1; j<=n; j++) {
      for(unsigned i=1; i<=m; i++) {
        in >> dummy;
        assert(dummy==i-1);
        p[j][i] = q[i][j] = getTime(in);
      }
    }
    for(unsigned i=0; i<=m; i++)
      q[i][0]=0;
    for(unsigned j=0; j<=n; j++)
      p[j][0]=p[j][m+1]=0;
    std::fill(&p[0][0],&p[0][m+1],0);
    l.resize(boost::extents[n+1][m+1][m+1]);
    for(unsigned j=1; j<=n; j++)
      for(unsigned i=1; i<m; i++)
        for(unsigned k=i+1; k<=m; k++) {
          l[j][i][k]=0;
          for(unsigned o=i+1; o<k; o++)
            l[j][i][k]+=p[j][o];
        }
    in >> dummyStr;
    assert(dummyStr.compare("Reldue") == 0);
    for(unsigned j=1; j<=n; j++) {
      in >> dummyInt;
      assert(dummyInt == -1);
      in >> dd[j];
      assert(dd[j]>=0);
      in >> dummyInt;
      assert(dummyInt == -1);
      in >> dummyInt;
      assert(dummyInt == 1);
    }
  }
  void reverse() {
    using std::swap;
    for(unsigned i=1; i<=m/2; i++)
      for(unsigned j=1; j<=n; j++) {
	      swap(p[j][i],p[j][m-i+1]);
	      swap(q[i][j],q[m-i+1][j]);
      }
  }
  Time totalTime(unsigned j=0) const {
    Time T = 0;
    if (j==0)
      for(unsigned i=1; i<=m; i++)
	      for(unsigned j=1; j<=n; j++)
	        T += p[j][i];
    else
      for(unsigned i=1; i<=m; i++)
	      T += p[j][i];
    return T;
  }
  unsigned LB_Cmax() {
    std::vector<unsigned> b(m+1);
    std::vector<unsigned> a(m+1);
    std::vector<unsigned> Tm(m+1);
    std::vector<unsigned> tpt(n+1);
    for(unsigned i = 1; i <= m; i++) {
      for(unsigned j = 1; j <= n; j++) {
        Tm[i] += p[j][i];
        tpt[j] += p[j][i];
      }
      if(i < m)
        b[i+1] = *std::min_element(tpt.begin()+1,tpt.end());
    }
    std::vector<unsigned> _tpt = tpt;
    for(unsigned i = 1; i <= m; i++) {
      for(unsigned j = 1; j <= n; j++)
        _tpt[j] -= p[j][i];
      if(i < m)
        a[i] = *std::min_element(_tpt.begin()+1,_tpt.end());
    }
    std::vector<unsigned> v(m+1);
    for(unsigned i = 1; i <= m; i++)
      v[i] = b[i] + Tm[i] + a[i];
    return std::max(*std::max_element(v.begin(),v.end()), *std::max_element(tpt.begin(),tpt.end()));
  }
};
