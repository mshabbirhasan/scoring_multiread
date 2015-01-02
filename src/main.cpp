// ==========================================================================
//                             scoring_multiread
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Saima Sultana Tithi <saima5@vt.edu>
// ==========================================================================

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/modifier.h>
#include <sstream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <math.h>
#include <map>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;

typedef multimap<string, string> AmbiguousMap;
typedef AmbiguousMap::iterator mapAmbIter;

typedef multimap<string, string> UniqueMap;
typedef UniqueMap::iterator mapUniIter;

char getRevComp(char c) 
{
    switch(c) {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'C':
            return 'G';
        case 'G':
            return 'C';    
        default:
            return c;
    }
}

double convertPhred33(char c) 
{
    int x = c - 33;
    return pow(10.0, (-x/10.0));
}

vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

void prior(string seq, string genomicSeq, double phredArray[], double priorArr[])
{
	double pAC = 0.00025;
	double pAT = 0.00025;
	double pAG = 0.0005;
	double pCA = 0.00025;
	double pCT = 0.0005;
	double pCG = 0.00025;
	double pTC = 0.0005;
	double pTA = 0.00025;
	double pTG = 0.00025;
	double pGA = 0.0005;
	double pGT = 0.00025;
	double pGC = 0.00025;
	double pSNP = 0.001;
	double CG = 0.85;
	double CH = 0.02;
    
	for (int i = 0; i < seq.length(); i++)
    {
		priorArr[i] = 0.0;
	}
    
    for (int i = 0; i < seq.length(); i++)
    {
        if (phredArray[i] > 0)
        {
            if (genomicSeq[i] == 'A')
            {
                if (seq[i] == 'A')
                    priorArr[i] = 1 - pSNP;
                else if (seq[i] == 'C')
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pAC * CG;
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pAC * CG;
                    else
                        priorArr[i] = pAC * CH;                 
                }
                else if (seq[i] == 'G')
                    priorArr[i] = pAG;
                else
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pAT + pAC * (1-CG);
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pAT + pAC * (1-CG);
                    else
                        priorArr[i] = pAT + pAC * (1-CH);
                }
            }
            else if (genomicSeq[i] == 'C')
            {
                if (seq[i] == 'A')
                    priorArr[i] = pCA;
                else if (seq[i] == 'C')
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = (1-pSNP) * CG;
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = (1-pSNP) * CG;
                    else
                        priorArr[i] = (1-pSNP) * CH;                
                }
                else if (seq[i] == 'G')
                    priorArr[i] = pCG;
                else
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pCT + (1-pSNP) * (1-CG);
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pCT + (1-pSNP) * (1-CG);
                    else
                        priorArr[i] = pCT + (1-pSNP) * (1-CH);
                }
            }
            else if (genomicSeq[i] == 'G')
            {
                if (seq[i] == 'A')
                    priorArr[i] = pGA;
                else if (seq[i] == 'C')
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pGC * CG;
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pGC * CG;
                    else
                        priorArr[i] = pGC * CH;                
                }
                else if (seq[i] == 'G')
                    priorArr[i] = 1 - pSNP;
                else
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pGT + pGC * (1-CG);
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pGT + pGC * (1-CG);
                    else
                        priorArr[i] = pGT + pGC * (1-CH);
                }
            }
            else //T at genomic location
            {
                if (seq[i] == 'A')
                    priorArr[i] = pTA;
                else if (seq[i] == 'C')
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pTC * CG;
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pTC * CG;
                    else
                        priorArr[i] = pTC * CH;                
                }
                else if (seq[i] == 'G')
                    priorArr[i] = pTG;
                else
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = (1-pSNP) + pTC * (1-CG);
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = (1-pSNP) + pTC * (1-CG);
                    else
                        priorArr[i] = (1-pSNP) + pTC * (1-CH);
                }
            }
        }
        else
		{
            priorArr[i] = 1;
		}
    }
    
}

double sum(double observedProb[], int size)
{
    double sum = 0.0;
    for(int i = 0; i < size; i++)
    {
        sum += observedProb[i];
    }
    return sum;
}

int main(int argc, char const ** argv)
{
	if (argc != 4)
    {
        std::cerr << "USAGE: ./main FILE.fa ambiguous_read_file unique_overlap_read_file \n";
        return 1;
    } 
	
	// Try to load index and create on the fly if necessary.
	seqan::FaiIndex faiIndex;
	if (seqan::read(faiIndex, argv[1]) != 0)
	{
		if (build(faiIndex, argv[1]) != 0)
		{
		    std::cerr << "ERROR: Index could not be loaded or built.\n";
		    return 1;
		}
		if (write(faiIndex) != 0)  // Name is stored from when reading.
		{
		    std::cerr << "ERROR: Index could not be written do disk.\n";
		    return 1;
		}
	}

	AmbiguousMap myAmbMap;
    
    string line;
    ifstream myfile(argv[2]);
    if (myfile.is_open()) {
        while (getline(myfile, line)) {
            string::size_type pos = line.find_first_of('\t');
            string key = line.substr(0, pos);
            string value = line.substr(pos);
            value.erase(0, value.find_first_not_of('\t'));
            myAmbMap.insert(make_pair(key, value));
        }
        myfile.close();
    }

	UniqueMap myUniMap;
    
    ifstream myfileUni(argv[3]);
    if (myfileUni.is_open()) {
        while (getline(myfileUni, line)) {
            vector<string> tokens = split(line, '\t');
            string key = tokens.at(3) + "_" + tokens.at(0) + "_" + tokens.at(1);
            string value = tokens.at(4) + "\t" + tokens.at(6) + "\t" + tokens.at(10) 
                    + "\t" + tokens.at(11) + "\t" + tokens.at(12) + "\t" + tokens.at(13) 
                    + "\t" + tokens.at(14) + "\t" + tokens.at(15) + "\t" + tokens.at(16);
            myUniMap.insert(make_pair(key, value));
        }
        myfileUni.close();
    }

	ofstream outputfile;
  	outputfile.open ("Result.txt");

	ofstream outputfile1;
  	outputfile1.open ("Result_with_highest_probability.txt");

	string header = "#ID\tCHROM\tPOS\tProbability(999 means probability cannot be calculated)\n";
	string header1 = "#ID\tCHROM\tPOS\tHighest_Probability\n";

	outputfile << header;
	outputfile1 << header1;

    mapAmbIter m_it, s_it;
    mapUniIter mUni_it, sUni_it;

    for (m_it = myAmbMap.begin();  m_it != myAmbMap.end();  m_it = s_it)
    {
        string theKey = (*m_it).first;

        pair<mapAmbIter, mapAmbIter> keyRange = myAmbMap.equal_range(theKey);

        // Iterate over all map elements with key == theKey
		vector<double> likelihoodArr;
		string output = "";
        for (s_it = keyRange.first;  s_it != keyRange.second;  ++s_it)
        {
            string theValue = (*s_it).second;
            vector<string> tokens = split(theValue, '\t');
            unsigned flag = atoi(tokens.at(0).c_str());
            string seq = tokens.at(8);
            string phred33 = tokens.at(9);

			double multiSeqErr[phred33.length()];
    
            for (int i = 0; i < phred33.length(); i++)
            {
               multiSeqErr[i] = convertPhred33(phred33[i]); 
            }

			//get the genome seq using seqan library
            string genomicSeq;			

			// Translate sequence name to index.
			unsigned idx = 0;
			if (!getIdByName(faiIndex, tokens.at(1), idx))
			{
				std::cerr << "ERROR: Index does not know about sequence " << tokens.at(1) << "\n";
				return 1;
			}
			//unsigned seqLength = sequenceLength(faiIndex, idx);
			//cout << seqLength << "\n";
			
			// Load characters of chromosome
			seqan::CharString seqChrPrefix;
			unsigned startPos = atoi(tokens.at(2).c_str()) - 2;
			unsigned endPos = startPos + seq.length() + 2;

			if (readRegion(seqChrPrefix, faiIndex, idx, startPos, endPos) != 0)
			{
				std::cerr << "ERROR: Could not load reference sequence.\n";
				return 1;
			}
			//cout << "Seq:" << seqChrPrefix << "\n";
						
			std::stringstream stream;
			stream<<seqChrPrefix;

			genomicSeq = stream.str();
			//std::cout << genomicSeq.length() << "\n";

			if(flag == 16 || flag == 272) {
				string revGenomicSeq(genomicSeq);
				for (int i = 0; i < genomicSeq.length(); i++)
                {
                   revGenomicSeq[genomicSeq.length() - 1 -i] = getRevComp(genomicSeq[i]); 
                }
				genomicSeq = revGenomicSeq;
			}
			genomicSeq = genomicSeq.substr(1);
			
			double priorArr[seq.length()];
            prior(seq, genomicSeq, multiSeqErr, priorArr);		

            if(flag == 16 || flag == 272) {   
    			double revPriorArr[seq.length()];
				for (int i = 0; i < seq.length(); i++)
                {
                   revPriorArr[seq.length() - 1 -i] = priorArr[i]; 
                }
				for (int i = 0; i < seq.length(); i++)
                {
                   priorArr[i] = revPriorArr[i]; 
                }

				double revMultiSeqErr[phred33.length()];
				for (int i = 0; i < phred33.length(); i++)
                {
                   revMultiSeqErr[phred33.length() - 1 -i] = multiSeqErr[i]; 
                }
				for (int i = 0; i < phred33.length(); i++)
                {
                   multiSeqErr[i] = revMultiSeqErr[i]; 
                }

                for (int i = 0; i < seq.length(); i++)
                {
                   seq[tokens.at(8).length() - 1 -i] = getRevComp(tokens.at(8)[i]); 
                }
            }
			
            double observed[seq.length()];
			//initialize observed array
			for(int a = 0; a < seq.length(); a++)
			{
				observed[a] = 0.0;
			}

            double posterior[seq.length()];
			//initialize posterior array
			for(int a = 0; a < seq.length(); a++)
			{
				posterior[a] = 0.0;
			}
            
            string theUniKey = theKey + "_" + tokens.at(1) + "_" + tokens.at(2);
            if (myUniMap.count(theUniKey) != 0)
            {
                pair<mapUniIter, mapUniIter> keyRangeUni = myUniMap.equal_range(theUniKey);
                
                for (int i = 0; i < seq.length(); i++) 
                {
                    if (multiSeqErr[i] > 0)
                    {
                        double observedProb[seq.length()];
						//initialize observedProb array
						for(int a = 0; a < seq.length(); a++)
						{
							observedProb[a] = 0.0;
						}
	
                        int baseTotal = 0;
                        
                        // Iterate over all map elements with key == theUniKey
                        for (sUni_it = keyRangeUni.first;  sUni_it != keyRangeUni.second;  ++sUni_it)
                        {
                            string theUniValue = (*sUni_it).second;
                            vector<string> uniTokens = split(theUniValue, '\t');
                            
                            double uniSeqErr[uniTokens.at(8).length()];
    
                            for (unsigned x = 0; x < uniTokens.at(8).length(); x++)
                            {
                               uniSeqErr[x] = convertPhred33(uniTokens.at(8)[x]); 
                            }
                            unsigned ambPos = atoi(tokens.at(2).c_str());
                            unsigned uniPos = atoi(uniTokens.at(1).c_str());
                            if(uniPos <= (ambPos+i) && (uniPos+uniTokens.at(7).length()-1) >= (ambPos+i))
                            {
                                if(uniTokens.at(7).at(ambPos+i-uniPos) == seq[i])
                                {
                                    observedProb[baseTotal] = 1 - uniSeqErr[ambPos+i-uniPos] 
                                            - multiSeqErr[i] + (uniSeqErr[ambPos+i-uniPos] * multiSeqErr[i]);
                                }
                                else
                                {
                                    observedProb[baseTotal] = uniSeqErr[ambPos+i-uniPos] 
                                            + multiSeqErr[i] - (uniSeqErr[ambPos+i-uniPos] * multiSeqErr[i]);
                                }
                                baseTotal++;
                            }
                        }
                        if (baseTotal != 0)
                        {
                            observed[i] = sum(observedProb, seq.length()) / baseTotal;
                        }
                        else
                        {
                            observed[i] = 0.999999999999999;
                        }
						double value = (priorArr[i] * observed[i]) / (priorArr[i] * observed[i] + (1-priorArr[i]) * (1-observed[i]));
						//std::cout << std::setprecision(15) << value << "\n";
						
						if (value > 0)
							posterior[i] = log10(value);
                    }
                    else
                    {
                        posterior[i] = 0.0;
                    }
                }//end for
				double likelihood = sum(posterior, seq.length());
				output = theKey + "\t" + tokens.at(1) + "_" + tokens.at(2) + "\t";
				outputfile << output;
				outputfile << std::setprecision(15) << likelihood << "\n";
				likelihoodArr.push_back(likelihood);
            }
            else
            {
				output = theKey + "\t" + tokens.at(1) + "_" + tokens.at(2) + "\t";
				outputfile << output;
				outputfile << 999 << "\n";
            }
		
        }//end for
		if(!likelihoodArr.empty()) {
			vector<double>::const_iterator it;
			it = std::max_element(likelihoodArr.begin(), likelihoodArr.end());
			outputfile1 << output;
			outputfile1 << std::setprecision(15) << *it << "\n";	
		}
    }//end for

	outputfile.close(); 
	outputfile1.close();
	cout << "Output is written in Result.txt and Result_with_highest_probability.txt\n"; 

    return 0;
}
