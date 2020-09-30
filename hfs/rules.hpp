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
#include <random>
#include <boost/accumulators/accumulators.hpp> 
#include <boost/accumulators/statistics/mean.hpp> 
#include <boost/accumulators/statistics/variance.hpp> 
#include <boost/accumulators/statistics/stats.hpp>
#include "instance.hpp"
inline std::vector<Job> getTotalTimeOrder(const Instance& I) {
  std::vector<Time> T(I.n+1,0);
  for(unsigned j = 1; j <= I.n; j++)
    for(unsigned i = 1; i <= I.m; i++)
      T[j] += I.p[j][i];
  std::vector<Job> pi(I.n+1);
  std::iota(pi.begin(), pi.end(),0);
  sort(pi.begin()+1, pi.end(), [&T](Job i, Job j) {
      return T[i]>T[j] || (T[i]==T[j] && i<j); 
    });
  return pi;
}
void sortFF(std::vector<unsigned>& pi, std::vector<unsigned>& T, const unsigned left, const unsigned right) {
  unsigned i = left; 
  unsigned j = right;
  const unsigned mid = T[(left+right)/2];
  while(i <= j) {
    while(T[i] > mid && i < right)
      i++;
    while(T[j] < mid && j > left)
      j--;
    if(i <= j) {
      std::swap(T[i], T[j]);
      std::swap(pi[i], pi[j]);
      i++;
      j--;
    }
  }
  if(left < j)
    sortFF(pi, T, left, j);
  if(right > i)
    sortFF(pi, T, i, right);
}
inline std::vector<Job> getTotalTimeOrderFF(const Instance& I) {
  std::vector<Time> T(I.n+1,0);
  for(unsigned j = 1; j <= I.n; j++)
    for(unsigned i = 1; i <= I.m; i++)
      T[j] += I.p[j][i];
  std::vector<Job> pi(I.n+1);
  std::iota(pi.begin(), pi.end(),0);
  sortFF(pi, T, 1, I.n);
  return pi;
}
inline std::vector<Job> getCorrectedMeanTimeOrder(const Instance& I) {
  std::vector<boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::mean,boost::accumulators::tag::variance> > > acc(I.n+1);
  for(unsigned j=1; j<=I.n; j++)
    for(unsigned i=1; i<=I.m; i++)
      acc[j](I.p[j][i]);
  std::vector<double> T(I.n+1,0);
  for(unsigned j=1; j<=I.n; j++)
    T[j]=boost::accumulators::mean(acc[j])-boost::accumulators::variance(acc[j]);
  std::vector<Job> pi(I.n+1);
  std::iota(pi.begin(),pi.end(),0);
  sort(pi.begin()+1, pi.end(), [&T](Job i, Job j) {
      return
	T[i]>T[j] ||
	(T[i]==T[j] && i<j);
    });
  return pi;
}
inline std::vector<Job> getKK1(const Instance& I, std::vector<char>& first) {
  const unsigned c = (I.m-1)*(I.m-2)/2;
  const unsigned n = I.n;
  const unsigned m = I.m;
  std::vector<Time> a(n+1), b(n+1), T(n+1), mT(n+1);
  for(unsigned j=1; j<=n; j++)
    for(unsigned i=1; i<=m; i++) {
      a[j] += (c+m-i)*I.p[j][i];
      b[j] += (c+i-1)*I.p[j][i];
      T[j] += I.p[j][i];
      mT[j] = std::max(mT[j],I.p[j][i]);
    }
  std::vector<Job> pi(n+1);
  std::iota(pi.begin(),pi.end(),0);
  std::sort(pi.begin()+1, pi.end(), [&a,&b,&T,&mT](const Job& j1, const Job& j2) {
      return
	std::min(a[j1],b[j1])>std::min(a[j2],b[j2]) ||
	(std::min(a[j1],b[j1])==std::min(a[j2],b[j2]) && j1<j2);
    });
  for(unsigned j=1; j<=n; j++)
    first[j] = a[j] <= b[j];
  return pi;
}
inline std::vector<Job> getKK2(const Instance& I, std::vector<char>& first) {
  const unsigned n = I.n;
  const unsigned m = I.m;
  std::vector<int> T(n+1,0), mT(n+1,0);
  std::vector<int> U(n+1,0);
  unsigned s = m/2;
  unsigned t = (m+1)/2;
  for(unsigned j=1; j<=n; j++) {
    for(unsigned i=1; i<=m; i++) {
      T[j] += I.p[j][i];
      mT[j] = std::max(mT[j],int(I.p[j][i]));
    }
    T[j] *= 4*s-3;
    for(unsigned h=1; h<=s; h++)
      U[j] += int(4*h-3)*(int(I.p[j][s+1-h])-int(I.p[j][t+h]));
  }
  std::vector<Job> pi(n+1);
  std::iota(pi.begin(),pi.end(),0);
  std::sort(pi.begin()+1, pi.end(), [&U,&T,&mT](const Job& a, const Job& b) {
      return
	 std::min(T[a]+U[a],T[a]-U[a])> std::min(T[b]+U[b],T[b]-U[b]) ||
	  (std::min(T[a]+U[a],T[a]-U[a])==std::min(T[b]+U[b],T[b]-U[b]) && T[a] > T[b]) || 
	  (std::min(T[a]+U[a],T[a]-U[a])==std::min(T[b]+U[b],T[b]-U[b]) && T[a] == T[b] && mT[a] > mT[b]);
    });
  for(unsigned j=1; j<=n; j++)
    first[j] = (U[j] <= 0);
  return pi;
}
inline std::vector<int> getCurvature(const Instance& I) {
  const unsigned n = I.n;
  const unsigned m = I.m;
  std::vector<int> c(n+1,0);
  for(unsigned j=1; j<=n; j++) {
    std::vector<Time> C(m+1,0);
    for(unsigned i=1; i<=m; i++)
      C[i] = C[i-1]+I.p[j][i];
    for(unsigned i=1; i<=m; i++) {
      double alpha = double(i-1)/double(m-1);
      Time t = (1-alpha)*C[1]+alpha*C[m];
      c[j] += int(C[i])-int(t);
    }
  }
  return c;
}
inline std::vector<double> getAvgDev(const Instance& I) {
  std::vector<boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance> > > acc(I.n+1);
  for(unsigned j=1; j<=I.n; j++)
    for(unsigned i=1; i<=I.m; i++)
      acc[j](I.p[j][i]);
  std::vector<double> T(I.n+1,0);
  for(unsigned j=1; j<=I.n; j++)
    T[j]=boost::accumulators::mean(acc[j]);
  return T;
}
inline std::vector<Job> getEDD(const Instance& I) {
  std::vector<Time> D(I.n+1,0);
  for(unsigned j = 1; j <= I.n; j++)
    D[j] = I.dd[j];
  std::vector<Job> pi(I.n+1);
  std::iota(pi.begin(), pi.end(),0);
  sort(pi.begin()+1, pi.end(), [&D](Job i, Job j) {
      return D[i] < D[j] || (D[i] == D[j] && i < j);
    });
  return pi;
}
inline unsigned firstposition(const unsigned job, const unsigned position) {
  return position;
}
inline unsigned lastposition(const unsigned job, const unsigned position) {
  return std::numeric_limits<unsigned>::max()-position;
}
inline unsigned randomtiebreak(const unsigned job, const unsigned position) {
  return lrand48()%1000;
}
inline unsigned usebool(const unsigned job, const unsigned position, const std::vector<char>& first) {
  if (first[job])
    return position;
  else
    return std::numeric_limits<unsigned>::max()-position;
}
