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
#include <algorithm>
#include <fstream>
#include <iterator>
#include <limits>
#include <utility>

#include "hfs/instance.hpp"
#include "hfs/solution.hpp"
#include "hfs/ls.hpp"
#include "hfs/parameters.hpp"

#include "Components.hpp"

namespace search {

inline PFSSolution iga(Instance& I, SGEDecider& d, Params& params, const std::chrono::system_clock::time_point& start, std::mt19937& rng, unsigned lb = 0, bool verbose=false) {
   using namespace std::placeholders;
   unsigned steps = 0, sbest = 0, saccepted = 0;
   builder construct(I);
   PFSSolution cs = initial_sol(I, d, params, rng, construct, start);
   std::vector<unsigned> past(params.l, cs.ms);
   if (verbose) std::cout << cs.ms << " " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start).count() << " ";
   if (d.dv_local_search > 0)
      local_search(I, d, params, cs, rng, construct, steps, start);
   PFSSolution bs = cs;
   PFSSolution n(I);
   std::chrono::system_clock::time_point tbest = std::chrono::system_clock::now();
   do {
      n = cs;
      iga_step(I, d, params, n, rng, construct, start);
      local_search(I, d, params, n, rng, construct, steps, start);
      if (n.ms < bs.ms) {
         sbest = steps;
         bs = n;
         tbest = std::chrono::system_clock::now();
      }
      saccepted += acceptance(n, cs, bs.ms, d, params, past, steps % params.l);
      past[steps % params.l] = cs.ms;
      steps++;
   } while (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < params.timelimit && bs.ms>lb);
   if (verbose) {
      std::cout << steps << " " << saccepted << " " << bs.ms << " " << sbest << " " << std::chrono::duration_cast<std::chrono::seconds>(tbest-start).count() << " " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() << std::endl;
   } else if (!params.np) {
      std::cout << bs.ms << std::endl;
   }
   return bs;
}

inline NPFSSolution iga_np(Instance& I, SGEDecider& d, Params& params, const std::chrono::system_clock::time_point& start, std::mt19937& rng, unsigned lb = 0, bool verbose=false) {
   unsigned steps = 0, sbest = 0, saccepted = 0;
   nBuilder construct(I);
   PFSSolution is = iga(I, d, params, start, rng);
   if (verbose) std::cout << is.ms << " " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start).count() << " ";
   NPFSSolution cs(I, is);
   NPFSSolution bs = cs;
   if (d.dv_local_search_np > 0) 
      local_search_np(I, d, params, cs, rng, construct, steps, start);
   std::vector<unsigned> past(params.l_np, cs.ms);
   std::chrono::system_clock::time_point tbest = std::chrono::system_clock::now();
   do {
      NPFSSolution n = cs;
      iga_step_np(I, d, params, n, rng, construct, start);
      local_search_np(I, d, params, n, rng, construct, steps, start);
      if (n.ms < bs.ms) {
         sbest = steps;
         bs = n;
         tbest = std::chrono::system_clock::now();
      }
      saccepted += acceptance_np(n, cs, bs.ms, d, params, past, steps % params.l_np);
      past[steps % params.l_np] = cs.ms;
      steps++;
   } while (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() < params.timelimit_np && bs.ms > lb);
   if (verbose) {
      std::cout << steps << " " << saccepted << " " << bs.ms << " " << sbest << " " << std::chrono::duration_cast<std::chrono::seconds>(tbest-start).count() << " " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start).count() << std::endl;
   } else {
      std::cout << bs.ms << std::endl;
   }
   return bs;
}
} 
