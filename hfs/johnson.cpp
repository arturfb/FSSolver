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
 * \file johnson.cpp
 *   \author Marcus Ritt <mrpritt@inf.ufrgs.br> 
 *   \version $Id: emacs 6018 2015-10-14 18:39:39Z ritt $
 *   \date Time-stamp: <2020-09-29 20:16:01 ritt>
 */
#include "johnson.hpp"
Time twomachine(const Instance& I, std::vector<Job>& pi, unsigned mA, unsigned mB) {
  return twomachine(I,pi.begin()+1,pi.end(),mA,mB);
}
template <class ForwardIterator>
Time twomachine(const Instance& I, ForwardIterator first, ForwardIterator last, unsigned mA, unsigned mB, Time h2, Time q1) {
  assert(mA <= mB);
  if (mA == mB) {
    Time t = 0;
    while (first != last) {
      t += I.p[*first][mA];
      first++;
    }
    return t;
  } else if (mA + 1 == mB) {
    johnson(I,first,last,mA,mB);
    unsigned C1 = 0, C2 = 0;
    while (first != last) {
      C1 += I.p[*first][mA];
      C2 = std::max(C2,C1)+I.p[*first][mB];
      first++;
    }
    return C2;
  } else {
    auto mid = partition(first,last,[&I,&mA,&mB](const Job& i) { return I.p[i][mA]<I.p[i][mB];});
    std::sort(first,mid,[&I,&mA,&mB](const Job& i, const Job& j) { return I.p[i][mA]+I.l[i][mA][mB]<I.p[j][mA]+I.l[j][mA][mB];});
    std::sort(mid,last,[&I,&mA,&mB](const Job& i, const Job& j) { return I.p[i][mB]+I.l[i][mA][mB]>I.p[j][mB]+I.l[j][mA][mB];});
    unsigned C1 = 0, C2 = h2; 
    while (first != last) {
      C1 += I.p[*first][mA];
      C2 = std::max(C2,C1+I.l[*first][mA][mB])+I.p[*first][mB];
      first++;
    }
    C1 += q1;
    C2 = std::max(C2,C1);
    return C2;
  }
}
