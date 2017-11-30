
/*

This program was created at:  Tue Sep  8 21:05:23 2015
This program was created by:  Zev N. Kronenberg


Contact: zev.kronenber@gmail.com

Organization: Unviersity of Utah
    School of Medicine
    Salt Lake City, Utah


The MIT License (MIT)

Copyright (c) <2015> <Zev N. Kronenberg>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*/

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <algorithm>
#include "split.h"
#include "gpatInfo.hpp"

struct options{
  std::string file;
  double cut;
}globalOpts;

using namespace std;

static const char *optString = "hf:s:";

void printHelp(void){

  cerr << endl << endl;
  cerr << "INFO: help" << endl;
  cerr << "INFO: description:" << endl;
  cerr << "      Creates genomic segments (bed file) for regions with high wcFst  " << endl;

  cerr << "Output : 8 columns :                 "    << endl;
  cerr << "     1. Seqid                        "    << endl;
  cerr << "     2. Start (zero based)           "    << endl;
  cerr << "     3. End   (zero based)           "    << endl;
  cerr << "     4. Average iHS                  "    << endl;
  cerr << "     5. Average high Fst (iHS > -s)  "    << endl;
  cerr << "     6. N iHS values in segment      "    << endl;
  cerr << "     7. N high iHS values in segment "    << endl;
  cerr << "     8. Segment length               "    << endl;

  cerr << "INFO: usage:  segmentFst -s 2 -f iHS.normalized.output.txt " << endl;
  cerr << endl;
  cerr << "INFO: required: -f            -- Output from normalizeIHS     "   << endl;
  cerr << "INFO: optional: -s            -- High absolute iHS cutoff [2] "    << endl;

  cerr << endl;

  printVersion();
}

//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
    {
    int opt = 0;
    globalOpts.file = "NA";
    opt = getopt(argc, argv, optString);
    while(opt != -1){
      switch(opt){
      case 's':
	{
	  string op = optarg;    
	  globalOpts.cut = atof(op.c_str());
	  break;
	}
      case 'h':
	{
	  printHelp();
	  exit(1);
	  break;
	}
	
      case 'f':
	{
	  globalOpts.file = optarg;
	  break;
	}
      case '?':
	{
	  break;
	}
      }	
	opt = getopt( argc, argv, optString ); 
   }
return 1;
}
//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  :

 Function does   :

 Function returns:

*/

bool growWindow(vector<double> & values, 
		int * begin            , 
		int * end              , 
		int * nhigh            , 
		int * nlow             ,
		double * hSum          ,
		double * lSum           ){

  *begin -= 1;
  *end   += 1;

  if(*begin < 0){
    return false;
  }
  if(*end >= values.size()){
    return false;
  }

  for(int index = *begin; index <= *end; index++){
    if(values[index] > globalOpts.cut){
      *nhigh += 1;
      *hSum += values[index];
    }
    else{
      *nlow += 1;
      *lSum += values[index];
    }
  }
  if((*nhigh)*2 < (*nlow)){
    return false;
  }
  return true;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : takes the sorted Fst scores, positions, and seqids

 Function does   : segments and prints bed

 Function returns:

*/


void process(vector<int> & pos, vector<double> & value, vector<string> & seqid)
{
  

  // i is the index of the outter loop/aka variant sites.
  // always start the seed with 9 SNPs the window grows to 10 in "growWindow"
  for(int i = 9; i < value.size()-9; i++){

    int begin = i -9;
    int end   = i +9;
    
    int nHigh = 0;
    int nLow  = 0;
    
    double HighSum = 0;
    double LowSum  = 0;

    bool anyGroth = false;
    
    while(growWindow(value, &begin, &end, 
		     &nHigh, &nLow, &HighSum, &LowSum)){
      anyGroth = true;
    }
    // the current window had an extention
    if(anyGroth){
    // reset the index beyond current window
      i = end + 1;

      if(begin < 0){
	begin = 0;
      }
      if(end >= value.size()){
	end = (value.size() - 1);
      }

      double avgFstHigh = HighSum / double(nHigh);
      double avgFst     = (HighSum + LowSum) / (double(nHigh)+double(nLow));



      cout << seqid[begin]   << "\t"
           << pos[begin]  -1 << "\t"
           << pos[end]    -1 << "\t"
           << avgFst         << "\t"
           << avgFstHigh     << "\t"
           << nHigh + nLow   << "\t"
           << nHigh          << "\t"
           << (pos[end] - pos[begin])
           << endl;      

    } 
  }
}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
  globalOpts.cut = 2;
  int parse = parseOpts(argc, argv);

  string last;
  int lastPos = 0;

  vector<string> seqid;
  vector<int>      pos;
  vector<double> value;

  if(globalOpts.file.empty()){
    printHelp();
    exit(1);
  }

  string line;
  ifstream myfile (globalOpts.file);
  if (myfile.is_open())
    {
      while ( getline (myfile,line) )
	{
	  vector<string> lineDat = split(line, '\t');
	  if(lineDat.size() != 9){
	    cerr << "FATAL: not valid normalized iHS file" << endl;
	    exit(1);
	  }
	  if(last.compare(lineDat[0]) != 0){
	    last = lineDat[0];
	    process(pos, value, seqid);
	    pos.clear();
	    value.clear();
	    seqid.clear();
	    lastPos = 0;
	  }
	  else{
	    if(atoi(lineDat[1].c_str()) < lastPos){
	      cerr << "FATAL: normalized iHS input must be sorted by position." << endl;
	      exit(1);
	    }
	    lastPos = atoi(lineDat[1].c_str());
	    seqid.push_back(lineDat[0]);
	    pos.push_back(atoi(lineDat[1].c_str()));
	    value.push_back(abs(atof(lineDat[6].c_str())));	  
	  }
	}

      std::cerr << "INFO: about to segment: " << pos.size() << " scores." << std::endl;
      process(pos, value, seqid);
      myfile.close();
    }
  else{
    cerr << "FATAL: could not open file: " << globalOpts.file << endl;
    exit(1);
  }
  return 0;
}
