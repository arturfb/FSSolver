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
#include "Search.hpp"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <limits>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace boost;

int main(int argc, char *argv[]) {
   std::vector<std::string> iname(1);      
   std::vector<std::string> ename(1);
   std::vector<std::string> tname(1);
   std::vector<unsigned> dv(19);     

   Params params;

   po::options_description desc("Options");
   desc.add_options()
      ("help,h", "show help")
      ("obj", po::value<char>(&params.obj)->default_value('m'), "objective: m for makespan, f for total completion time, or t for total tardiness")
      ("instance,i", po::value<std::vector<std::string>>(&iname), "path to instance file")
      ("iga", "run iga")
      ("seed", po::value<unsigned>(), "random seed")
      ("timelimit,t", po::value<unsigned>(&params.timelimit)->default_value(10), "time limit (s)")
      ("timefactor", po::value<unsigned>(&params.timefactor)->default_value(30), "time per shop operation (ms)")
      ("both", po::value<bool>(&params.both)->default_value(false), "use standard and reverse instance for initial heuristic (0,1)")
      ("dc", po::value<unsigned>(&params.dc)->default_value(8), "number of jobs to destruct/reconstruct")
      ("dcf", po::value<double>(), "fraction of jobs destructed/reconstructed in IGA")
      ("n_ls", po::value<unsigned>(&params.n_ls)->default_value(std::numeric_limits<unsigned>::max()), "number of times local search is repeated")
      ("alpha", po::value<double>(&params.alpha)->default_value(0.5), "temperature multiplier")
      ("l", po::value<unsigned>(&params.l)->default_value(100), "l for late acceptance")
      ("rrtd", po::value<unsigned>(&params.rrtd)->default_value(10), "Record-to-Record Travel deviation")
      ("thres", po::value<double>(&params.thres)->default_value(0.005), "Threshold acceptance")
      ("alternate", po::value<unsigned>(&params.alternate)->default_value(2), "Interval at which the second local search is performed")
      ("x", po::value<unsigned>(), "value x for LR (Liu & Reeves, 2001)")
      ("p", po::value<unsigned>(&params.p)->default_value(3), "perturbation intensity for pool ILS")
      ("pf", po::value<double>(), "perturbation factor for ILS (p=I.n*pf)")
      ("ps", po::value<unsigned>(&params.pool_size)->default_value(5), "pool size")
      ("w", po::value<unsigned>(), "beam width")
      ("k", po::value<unsigned>(&params.k)->default_value(1), "factor k defines beam width n/k")
      ("rcon", po::value<unsigned>(&params.rcon)->default_value(25), "number of repeated constructions for initial solution")
      ("verbose", "verbose output")
      ("np", "use non-permutation schedules")
      ("npf", po::value<double>(), "percentage of global time limit to spend to exploring NP schedules")
      ("dc_np", po::value<unsigned>(&params.dc_np)->default_value(8), "number of jobs to destruct/reconstruct")
      ("alpha_np", po::value<double>(&params.alpha_np)->default_value(0.5), "temperature multiplier")
      ("n_ls_np", po::value<unsigned>(&params.n_ls_np)->default_value(std::numeric_limits<unsigned>::max()), "number of times local search is repeated")
      ("l_np", po::value<unsigned>(&params.l_np)->default_value(100), "l for late acceptance")
      ("rrtd_np", po::value<unsigned>(&params.rrtd_np)->default_value(10), "Record-to-Record Travel deviation")
      ("thres_np", po::value<double>(&params.thres_np)->default_value(0.005), "Threshold acceptance")
      ("alternate_np", po::value<unsigned>(&params.alternate_np)->default_value(4), "Interval at which the second local search is performed")
      ("verbose_np", "verbose output")
      ("dv01", po::value<unsigned>(&dv[0])->default_value(0), "Decider value: initial solution")
      ("dv02", po::value<unsigned>(&dv[1])->default_value(0), "Decider value: priority order for initial solution")
      ("dv03", po::value<unsigned>(&dv[2])->default_value(0), "Decider value: tiebreaker for initial solution")
      ("dv04", po::value<unsigned>(&dv[3])->default_value(0), "Decider value: tiebreaker for perturbation")
      ("dv05", po::value<unsigned>(&dv[4])->default_value(0), "Decider value: tiebreaker for local search procedure 1")
      ("dv06", po::value<unsigned>(&dv[5])->default_value(0), "Decider value: tiebreaker for local search procedure 2")
      ("dv07", po::value<unsigned>(&dv[6])->default_value(0), "Decider value: perturbation")
      ("dv08", po::value<unsigned>(&dv[7])->default_value(0), "Decider value: number of local search methods")
      ("dv09", po::value<unsigned>(&dv[8])->default_value(0), "Decider value: local search procedure 1")
      ("dv10", po::value<unsigned>(&dv[9])->default_value(0), "Decider value: local search procedure 2")
      ("dv11", po::value<unsigned>(&dv[10])->default_value(0), "Decider value: acceptance criterion")
      ("dv12", po::value<unsigned>(&dv[11])->default_value(0), "Decider value: tiebreaker for perturbation (NP)")
      ("dv13", po::value<unsigned>(&dv[12])->default_value(0), "Decider value: tiebreaker for local search procedure 1 (NP)")
      ("dv14", po::value<unsigned>(&dv[13])->default_value(0), "Decider value: tiebreaker for local search procedure 2 (NP)")
      ("dv15", po::value<unsigned>(&dv[14])->default_value(0), "Decider value: number of local search methods (NP)")
      ("dv16", po::value<unsigned>(&dv[15])->default_value(0), "Decider value: local search procedure 1 (NP)")
      ("dv17", po::value<unsigned>(&dv[16])->default_value(0), "Decider value: local search procedure 2 (NP)")
      ("dv18", po::value<unsigned>(&dv[17])->default_value(0), "Decider value: acceptance criterion (NP)")
      ("dv19", po::value<unsigned>(&dv[18])->default_value(0), "Decider value: perturbation (NP)");

   po::variables_map vm;
   po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
   po::notify(vm);

   if (vm.count("help")) {
      std::cout << desc << std::endl;
      return 0;
   }

   std::mt19937 rng;
   unsigned seed;
   if (!vm.count("seed")) {
      seed = time(0);
      std::ifstream f("/dev/urandom");
      if (f.good()) {
         f.read((char *)(&seed), sizeof(unsigned int));
      }
   } else {
      seed = vm["seed"].as<unsigned>();
   }
   srand48(seed);
   rng.seed(seed);

   if (params.obj != 'm') {
      if (params.obj == 'f' || params.obj == 't') {
         params.both = false;
      }
      else {
         std::cerr << "Unknown objective!" << std::endl;
         return 2;
      }
   }
   if (vm.count("verbose"))
      params.verbose = 1;
   else
      params.verbose = 0;

   if (vm.count("verbose_np"))
      params.verbose_np = 1;
   else
      params.verbose_np = 0;


   if (vm.count("np"))
      params.np = 1;
   else
      params.np = 0;

   Instance I;
   if (!vm.count("instance")) {
      I.read(std::cin);
   }
   else {
      std::filebuf fb;
      if (fb.open(iname[0], std::ios::in)) {
         std::istream is(&fb);
         if (params.obj != 't')
            I.read(is);
         else
            I.read_duedates(is);
         fb.close();
      }
   }


   if (vm.count("dcf")) {
      params.dc = int(ceil(I.n*vm["dcf"].as<double>()));
   }
   if (vm.count("pf")) {
      params.p = int(ceil(I.n*vm["pf"].as<double>()));
   }
   if (vm.count("timefactor")) {
      params.timelimit = (I.n*I.m*params.timefactor+999)/1000;
   }
   if (vm.count("npf")) {
      params.timelimit_np = params.timelimit;
      params.timelimit = unsigned(ceil(params.timelimit * (1.0 - vm["npf"].as<double>())));
   }
   else if (vm.count("np")) {
      params.timelimit_np = params.timelimit;
      params.timelimit = unsigned(ceil(params.timelimit * 0.5));      
   }
   if (vm.count("k")) {
      params.w = I.n / vm["k"].as<unsigned>();
   }
   if (vm.count("w")) {
      params.w = vm["w"].as<unsigned>();
   }
   if (vm.count("x")) {
      params.x = vm["x"].as<unsigned>();
   } 
   else {
      params.x = unsigned(ceil((I.n+I.m-1)/I.m));
   }
   
   if (dv[10] == 0 || (vm.count("np") && dv[17] == 0)) {
      if(params.obj != 't') {
         const double pavg = double(I.totalTime()) / (I.n*I.m);
         if (params.obj == 'm')
            params.T = params.alpha*pavg/10;
         else
            params.T = params.alpha*pavg/10*double(I.n);
      }
      else {
         unsigned LB = I.LB_Cmax();
         double sum = 0;
         for(unsigned j = 1; j <= I.n; j++)
            sum += LB - I.dd[j];
         params.T = sum/(10*I.n);
      }
      
      if (vm.count("np"))
         params.T_np = params.T;
   }



   SGEDecider d(I, dv, params.obj, params.np);
   std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
   params.start = start;

   if (!vm.count("iga")) {
      std::cerr << "No algorithm specified! Running iga..." << std::endl;
   }

   if(!vm.count("np")) {
      PFSSolution s(I);
      s = search::iga(I, d, params, start, rng, 0, params.verbose);
   }
   else {
      
      NPFSSolution s_np(I);
      s_np = search::iga_np(I, d, params, start, rng, 0, params.verbose_np);
   }
 
   return 0;  
}
