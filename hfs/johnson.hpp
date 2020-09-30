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
/**
 * \file johnson.hpp
 *   \author Marcus Ritt <mrpritt@inf.ufrgs.br>
 *   \version $Id: johnson.hpp 6288 2016-01-04 01:28:28Z ritt $
 *   \date Time-stamp: <2020-09-29 20:15:35 ritt>
 *
 * Several variants of Johnson's algorithm for F2||C_max.
 */
#pragma once
#include <algorithm>
#include <vector>
#include "instance.hpp"
#include "units.hpp"
template <class ForwardIterator>
inline void johnson(const Instance& I, ForwardIterator first, ForwardIterator last, unsigned mA, unsigned mB) {
  auto mid = partition(first,last,[&I,&mA,&mB](const Job& i) { return I.p[i][mA]<I.p[i][mB];});
  std::sort(first,mid,[&I,&mA,&mB](const Job& i, const Job& j) { return I.p[i][mA]<I.p[j][mA];});
  std::sort(mid,last,[&I,&mA,&mB](const Job& i, const Job& j) { return I.p[i][mB]>I.p[j][mB];});
}
inline void johnson(const Instance& I, std::vector<Job>& pi, unsigned mA, unsigned mB) {
  johnson(I,pi.begin()+1,pi.end(),mA,mB);
}
inline std::vector<Job> johnson(const Instance& I, unsigned mA, unsigned mB) {
  std::vector<Job> pi(I.n+1);
  std::iota(pi.begin(),pi.end(),0);
  johnson(I,pi,mA,mB);
  return pi;
}
template <class ForwardIterator>
inline void johnson(ForwardIterator first, ForwardIterator last, const std::vector<Time>& A, const std::vector<Time>& B) {
  assert(A.size()==B.size() && size_t(distance(first,last))==A.size());
  auto mid = partition(first,last,[&A,&B](const Job& i) { return A[i]<B[i];});
  std::sort(first,mid,[&A,&B](const Job& i, const Job& j) { return A[i]<A[j];});
  std::sort(mid,last,[&A,&B](const Job& i, const Job& j) { return B[i]>B[j];});
}
inline void johnson(std::vector<Job>& pi, const std::vector<Time>& A, const std::vector<Time>& B) {
  assert(A.size()==B.size() && pi.size()==A.size());
  johnson(pi.begin()+1,pi.end(),A,B);
}
inline std::vector<Job> johnson(const std::vector<Time>& A, const std::vector<Time>& B) {
  assert(A.size()==B.size());
  std::vector<Job> pi(A.size());
  std::iota(pi.begin(),pi.end(),0);
  johnson(pi,A,B);
  return pi;
}
Time twomachine(const Instance& I, std::vector<Job>& pi, unsigned mA, unsigned mB);
template <class ForwardIterator>
Time twomachine(const Instance& I, ForwardIterator first, ForwardIterator last, unsigned mA, unsigned mB, Time h2=0, Time t1=0);
inline Time twomachine(const Instance& I, unsigned mA, unsigned mB) {
  std::vector<Job> pi(I.n+1,0);
  std::iota(pi.begin(),pi.end(),0);
  return twomachine(I,pi,mA,mB);
}
