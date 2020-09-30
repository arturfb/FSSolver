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
 * \file threemachine.hpp
 *   \author Marcus Ritt <mrpritt@inf.ufrgs.br>
 *   \version $Id: threemachine.hpp 6316 2016-01-10 19:19:38Z ritt $
 *   \date Time-stamp: <2020-09-29 20:12:02 ritt>
 */
#pragma once
#include <vector>
#include <list>
#include <iostream>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <boost/intrusive/list.hpp>
#include "instance.hpp"
#include "solution.hpp"
#include "units.hpp"
#include "johnson.hpp"
namespace m3 {
  struct psolution {
    std::vector<Job> pi;
    unsigned n1,n2;
    psolution(unsigned n, unsigned k) : pi(k+1), n1(0), n2(k+1) {
      assert(k<=n);
      iota(pi.begin(),pi.end(),0);
    }
    psolution(unsigned n) : psolution(n,n) {}
    psolution(psolution&& other) {
      this->swap(other);
    }
    psolution(const psolution& other) {
      pi = other.pi;
      n1 = other.n1;
      n2 = other.n2;
    }
    psolution& operator=(psolution other) {
      this->swap(other);
      return *this;
    }
    void swap(psolution& other) {
      using std::swap;
      pi.swap(other.pi);
      swap(n1,other.n1);
      swap(n2,other.n2);
    }
    void add(Job j) {
      n1++;
      for(unsigned ji=n1+1; ji<n2; ji++)
	    if (pi[ji]==j) {
	      std::swap(pi[n1],pi[ji]);
	    return;
	}
    }
    void addi(unsigned ji) {
      n1++;
      std::swap(pi[n1],pi[ji]);
    }
    void dda(Job j) {
      n2--;
      for(unsigned ji=n1+1; ji<n2; ji++)
	    if (pi[ji]==j) {
	      std::swap(pi[n2],pi[ji]);
	    return;
	}
    }
    void ddai(unsigned ji) {
      n2--;
      std::swap(pi[n2],pi[ji]);
    }
    bool empty() const {
      if (n1>0)
	      return false;
      if (n2<pi.size())
	      return false;
      return true;
    }
    unsigned totalJobs() const {
      return pi.size()-1;
    }
    bool full() const {
      return n1+1 >= n2;
    }
    bool determined() const {
      return n1+1 >= n2;
    }
    unsigned numFixed() const {
      return n1+pi.size()-n2;
    }
    unsigned numFree() const {
      return n2-n1-1;
    }
    bool free(Job j) {
      return std::find(pi.begin()+n1+1,pi.begin()+n2,j) != pi.begin()+n2;
    }
  };
}
