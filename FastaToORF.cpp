// FastaToORFs
// Adam Jones (adamjones@uidaho.edu)
// Modified 4/23/2020
// v1.1

// This program takes a fasta file and converts each sequence into its longest ORF
// no arguments = interactive mode

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>


std::string translate(std::string nucleotides)
{
	size_t i;
	std::string translation;
	std::string codon;
	std::string aa;
	
	for (i = 0; i < nucleotides.length(); i = i + 3)
	{
		aa = "X";
		codon = "XXX";
		if (i + 2 < nucleotides.length())
		{
			codon[0] = nucleotides[i];
			codon[1] = nucleotides[i + 1];
			codon[2] = nucleotides[i + 2];
		}

		if (codon == "TTT") aa = "F";
		if (codon == "TTC") aa = "F";
		if (codon == "TTA") aa = "L";
		if (codon == "TTG") aa = "L";
		if (codon == "TCT") aa = "S";
		if (codon == "TCC") aa = "S";
		if (codon == "TCA") aa = "S";
		if (codon == "TCG") aa = "S";
		if (codon == "TAT") aa = "Y";
		if (codon == "TAC") aa = "Y";
		if (codon == "TAA") aa = "*"; // I am using * for stop
		if (codon == "TAG") aa = "*";
		if (codon == "TGT") aa = "C";
		if (codon == "TGC") aa = "C";
		if (codon == "TGA") aa = "*";
		if (codon == "TGG") aa = "W";

		if (codon == "CTT") aa = "L";
		if (codon == "CTC") aa = "L";
		if (codon == "CTA") aa = "L";
		if (codon == "CTG") aa = "L";
		if (codon == "CCT") aa = "P";
		if (codon == "CCC") aa = "P";
		if (codon == "CCA") aa = "P";
		if (codon == "CCG") aa = "P";
		if (codon == "CAT") aa = "H";
		if (codon == "CAC") aa = "H";
		if (codon == "CAA") aa = "Q";
		if (codon == "CAG") aa = "Q";
		if (codon == "CGT") aa = "R";
		if (codon == "CGC") aa = "R";
		if (codon == "CGA") aa = "R";
		if (codon == "CGG") aa = "R";

		if (codon == "ATT") aa = "I";
		if (codon == "ATC") aa = "I";
		if (codon == "ATA") aa = "I";
		if (codon == "ATG") aa = "M";
		if (codon == "ACT") aa = "T";
		if (codon == "ACC") aa = "T";
		if (codon == "ACA") aa = "T";
		if (codon == "ACG") aa = "T";
		if (codon == "AAT") aa = "N";
		if (codon == "AAC") aa = "N";
		if (codon == "AAA") aa = "K";
		if (codon == "AAG") aa = "K";
		if (codon == "AGT") aa = "S";
		if (codon == "AGC") aa = "S";
		if (codon == "AGA") aa = "R";
		if (codon == "AGG") aa = "R";

		if (codon == "GTT") aa = "V";
		if (codon == "GTC") aa = "V";
		if (codon == "GTA") aa = "V";
		if (codon == "GTG") aa = "V";
		if (codon == "GCT") aa = "A";
		if (codon == "GCC") aa = "A";
		if (codon == "GCA") aa = "A";
		if (codon == "GCG") aa = "A";
		if (codon == "GAT") aa = "D";
		if (codon == "GAC") aa = "D";
		if (codon == "GAA") aa = "E";
		if (codon == "GAG") aa = "E";
		if (codon == "GGT") aa = "G";
		if (codon == "GGC") aa = "G";
		if (codon == "GGA") aa = "G";
		if (codon == "GGG") aa = "G";

		if (i == 0)
			translation = aa;
		else
			translation = translation + aa;

	}
	return translation;
}

class fastq_read {
public:
	std::string header;
	std::string sequence;
	std::string qual_header;
	std::string quality;
	int ORFstart;
	int ORFend;
	int ORFlength;
	bool ORFforward;
	std::string reversecomplement;
	std::string ORFseq;
	std::vector<std::string> all_Orfs;
	std::vector<int> all_Orfs_starts;
	std::vector<int> all_Orfs_ends;
	std::vector<int> all_Orfs_lengths;
	bool all_Orfs_reverse;
	int all_Orfs_readingframe;
	bool all_Orfs_ATG_start;
	int all_Orfs_counter;

	void AllORFs(bool Reverse, bool ATG_required, int ReadingFrame, int Min_ORF_Length) // The function returns the length of the longest ORF in the current reading frame
	{
		int ii;

		all_Orfs.clear();
		all_Orfs_starts.clear();
		all_Orfs_ends.clear();
		all_Orfs_reverse = Reverse;
		all_Orfs_readingframe = ReadingFrame;
		all_Orfs_ATG_start = ATG_required;
		all_Orfs_counter = 0;

		int tempORFstart, tempORFend;
		std::string teststring;
		if (Reverse)
			teststring = reversecomplement.substr(ReadingFrame,reversecomplement.length()-ReadingFrame);
		else
			teststring = sequence.substr(ReadingFrame,sequence.length()-ReadingFrame);

		// Find all stop codons in the test string
		// Also keep in mind that positions in the
		// test string are offset from the original string by
		// the reading frame.

		std::string codon("XXX");
		int seqlength = static_cast<int>(teststring.length());

		// Find the locations of all stop codons

		int NumberStops = 0;
		std::vector<int> StopList;
		int NumberATGs = 0;
		std::vector<int> ATGlist;
		StopList.clear();
		ATGlist.clear();

		for (ii = 0; ii < seqlength; ii = ii+3)
		{
			codon = "XXX";
			codon[0] = teststring[ii];
			if (ii+1 < seqlength)
				codon[1] = teststring[ii+1];
			if (ii+2 < seqlength)
				codon[2] = teststring[ii+2];

			if (codon == "TAG" || codon == "TAA" || codon == "TGA" || codon == "NNN" ||
				codon == "NNA" || codon == "NNG" || codon == "NAN" || codon == "NGN" ||
				codon == "TNN" || codon == "NAA" || codon == "NAG" || codon == "NGA" ||
				codon == "TNG" || codon == "TNA" || codon == "TAN" || codon == "TGN")
			{
				StopList.push_back(ii);
				NumberStops++;
			}

			if (codon == "ATG")
			{
				ATGlist.push_back(ii);
				NumberATGs++;
			}
		}

		if (ATG_required)  // find all transcripts starting with ATG longer than Min_ORF_length
		{
			all_Orfs_counter = 0;
			int last_stop = -1;
			bool same_as_last;
			int counter;
			bool stop_after_atg;
			for (ii = 0; ii < NumberATGs; ii++)
			{
				if (NumberStops == 0)
				{
					stop_after_atg = false;
				}
				else
				{
					if (StopList[NumberStops-1] > ATGlist[ii])
						stop_after_atg = true;
					else
						stop_after_atg = false;
				}

				if (stop_after_atg)  // There is at least one stop after the ATG
				{
					counter = 0;
					while (StopList[counter] <= ATGlist[ii])
						counter++;
					tempORFstart = ATGlist[ii];
					if (StopList[counter]+2 < seqlength)
						tempORFend = StopList[counter]+2; // this is the index of the last base of the stop in the tempsequence
					else
						tempORFend = seqlength - 1;
					
					if (tempORFend == last_stop)
						same_as_last = true;
					else
						same_as_last = false;
					last_stop = tempORFend;
				}
				else  // There are no stops after this ATG
				{
					tempORFstart = ATGlist[ii];
					tempORFend = seqlength-1;
					if (tempORFend == last_stop)
						same_as_last = true;
					else
						same_as_last = false;
					last_stop = tempORFend;
				}

				// Add temp ORF to the list of ORFs if it's long enough
				if (tempORFend-tempORFstart+1 >= Min_ORF_Length && !same_as_last)
				{
					all_Orfs_starts.push_back(tempORFstart+ReadingFrame);
					if (tempORFend+ReadingFrame < static_cast<int>(sequence.length()) && tempORFend+ReadingFrame < static_cast<int>(reversecomplement.length()))
					{
						all_Orfs_ends.push_back(tempORFend+ReadingFrame);
					}
					else
					{
						if (sequence.length()-1 < reversecomplement.length()-1)
							all_Orfs_ends.push_back(static_cast<int>(sequence.length())-1);
						else
							all_Orfs_ends.push_back(static_cast<int>(reversecomplement.length())-1);
					}

					all_Orfs.push_back(teststring.substr(tempORFstart,tempORFend-tempORFstart+1));
					all_Orfs_counter++;
				}
				
			}  // end of ii
		} // end of ATG required
		else
		{
			// Here we're finding all ORFs, whether or not there is an ATG, so it's basically just the spans between stop codons.
			all_Orfs_counter = 0;

			// If there are no stop codons, then the whole sequence is the reading frame.

			if (NumberStops == 0)
			{
				if (seqlength > Min_ORF_Length)
				{
					all_Orfs_starts.push_back(ReadingFrame); // It starts with the first base of this reading frame
					all_Orfs_ends.push_back(static_cast<int>(sequence.length())-1);  // It goes to the very end of the sequence
					all_Orfs.push_back(teststring);
					all_Orfs_counter++;
				}
			}
			else 
			{
				// if there are any stop codons, then the first ORF is from the start to the first stop.
				tempORFstart = 0;
				tempORFend = StopList[0] + 2;
				if (tempORFend-tempORFstart+1 > Min_ORF_Length && StopList[0] > 0)
				{
					all_Orfs_starts.push_back(tempORFstart+ReadingFrame);
					all_Orfs_ends.push_back(tempORFend+ReadingFrame);
					all_Orfs.push_back(teststring.substr(tempORFstart,tempORFend-tempORFstart+1));
					all_Orfs_counter++;
				}

				// Next comes the reading frames between stop codons

				for (ii = 1; ii < NumberStops; ii++)
				{
					tempORFstart = StopList[ii-1] + 3;  // the codon after the previous stop is the start
					tempORFend = StopList[ii] + 2;  // The end is the third position in the current stop

					if (tempORFend-tempORFstart+1 > Min_ORF_Length)
					{
						all_Orfs_starts.push_back(tempORFstart+ReadingFrame);
						all_Orfs_ends.push_back(tempORFend+ReadingFrame);
						all_Orfs.push_back(teststring.substr(tempORFstart,tempORFend-tempORFstart+1));
						all_Orfs_counter++;
					}
				}

				// And finally the one from the last stop codon to the end of the sequence

				tempORFstart = StopList[NumberStops-1] + 3;
				tempORFend = seqlength - 1;
				if (tempORFend - tempORFstart + 1 > Min_ORF_Length)
				{
					all_Orfs_starts.push_back(tempORFstart+ReadingFrame);
					all_Orfs_ends.push_back(static_cast<int>(sequence.length())-1);  // It goes to the very end
					all_Orfs.push_back(teststring.substr(tempORFstart,tempORFend-tempORFstart+1));
					all_Orfs_counter++;
				}

			}
		}

		// Calculate lengths of all retained ORFs
		all_Orfs_lengths.clear();
		for (ii = 0; ii < all_Orfs_counter; ii++)
		{
			all_Orfs_lengths.push_back(all_Orfs_ends[ii]-all_Orfs_starts[ii]+1);
		}
	}
};


int main( int argc, char* argv[] )
{
	int i, j, k;

	bool headerflag;
	std::string currentsequence;
	int numbertranscripts;
	int underscoreloc;
	bool nextlocus;
	std::string ORFoutFile;
	std::string AAoutFile;
	int min_orf_length;
	int retained_orf_counter;
	bool require_ATG;
	bool help_mode;
	bool interactive_mode;
	std::string query;
	std::string min_orf_length_string;
	std::string requireATG_string;
	std::string remove_Ns_string, remove_stop_string;
	bool remove_Ns;
	bool remove_stop;

	char *line = new char[1000000];

	std::string InFileName;

	help_mode = false;
	interactive_mode = false;
	if (argc <= 1)
	{
		std::cout << "\n(i)nteractive mode or (h)elp?\n";
		std::cin >> query;
		if (query[0] == 'h' || query[0] == 'H')
			help_mode = true;
		else
			interactive_mode = true;
	}

	// Set defaults
	min_orf_length = 0;
	min_orf_length_string = "0";
	InFileName = "input.fasta";
	ORFoutFile = "out_ORFs";
	require_ATG = false;
	requireATG_string = "n";
	remove_Ns = false;
	remove_Ns_string = "n";
	remove_stop = false;
	remove_stop_string = "n";


	// parse arguments

	std::string tempstring1;
	std::string tempstring2;

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h" || tempstring1 == "--help")
			help_mode = true;
		for (i = 1; i < argc - 1; i++)
		{
			tempstring1 = argv[i];
			tempstring2 = argv[i+1];
			if (tempstring1 == "-i")
				InFileName = tempstring2;
			if (tempstring1 == "-o")
				ORFoutFile = tempstring2;
			if (tempstring1 == "-m")
				min_orf_length_string = tempstring2;
			if (tempstring1 == "-a")
				requireATG_string = tempstring2;
			if (tempstring1 == "-r")
				remove_Ns_string = tempstring2;
			if (tempstring1 == "-s")
				remove_stop_string = tempstring2;
		}
	}

	if (help_mode)
	{
		std::cout << "\n\nFastaToORF, version 1.1\n";
		std::cout << "Adam Jones -- adamjones@uidaho.edu\n\n";
		std::cout << "This program takes a fasta file of DNA sequences\n";
		std::cout << "and converts each read into its longest ORF.\n\n";
		std::cout << "-i:\tinput filename [default: input.fasta]\n";
		std::cout << "-o:\toutput filename for ORFs [default: out_ORFs]\n";
		std::cout << "-m:\tminimum ORF length (in nucleotides) required to keep ORF [default: 0]\n";
		std::cout << "-a:\trequire ATG at the beginning of each ORF (y or n) [default: n]\n";
		std::cout << "-r:\tremove DNA sequences containing ambiguous bases (y or n) [default: n]\n";
		std::cout << "-s:\tremove stop codons (y or n) [default: n]\n";
		std::cout << "no arguments:\tinteractive mode\n";
		std::cout << "\nenter any character to exit\n";

		std::cin >> i;
		return 0;
	}

	if (interactive_mode)
	{
		std::cout << "Input file:\n";
		std::cin >> InFileName;

		std::cout << "\nOutput filename:\n";
		std::cin >> ORFoutFile;
			
		std::cout << "\nMin ORF length:\n";
		std::cin >> min_orf_length_string;

		std::cout << "\nRequire ATG start (y or n)?:\n";
		std::cin >> requireATG_string;

		std::cout << "\nRemove sequences with ambiguous bases (y or n)?:\n";
		std::cin >> remove_Ns_string;

		std::cout << "\nRemove stop codons (y or n)?:\n";
		std::cin >> remove_stop_string;

	}

	if (requireATG_string[0] == 'y' || requireATG_string[0] == 'Y')
		require_ATG = true;
	else
		require_ATG = false;

	if (remove_Ns_string[0] == 'n' || remove_Ns_string[0] == 'N')
		remove_Ns = false;
	else
		remove_Ns = true;

	if (remove_stop_string[0] == 'n' || remove_stop_string[0] == 'N')
		remove_stop = false;
	else
		remove_stop = true;

	std::stringstream tempss(min_orf_length_string);
	tempss >> min_orf_length;

	if (interactive_mode)
	{
		std::cout << "\n\nParameters:\n";
		std::cout << "\nIn File:\t" << InFileName;
		std::cout << "\nORF Outfile:\t" << ORFoutFile;
		std::cout << "\nMin ORF Len:\t" << min_orf_length;
		std::cout << "\nRequire ATG:\t";
		if (require_ATG)
			std::cout << "yes";
		else
			std::cout << "no";
		std::cout << "\nRemove Ambiguous Reads?:\t";
		if (remove_Ns)
			std::cout << "yes";
		else
			std::cout << "no";

		std::cout << "\nRemove Stop Codons?:\t";
		if (remove_stop)
			std::cout << "yes";
		else
			std::cout << "no";

		std::cout << "\n\nContinue (y/n)?\n";
		std::cin >> query;
		if (query[0] == 'n' || query[0] == 'N')
			return 0;
	}

	AAoutFile = ORFoutFile + "_aa.fasta";
	ORFoutFile = ORFoutFile + "_nucl.fasta";

	char* infile_name = new char[InFileName.length()+2];
	char* orffile_name = new char[ORFoutFile.length()+2];
	char* aafile_name = new char[AAoutFile.length() + 2];

	for (i = 0; i < static_cast<int>(InFileName.length()); i++)
		infile_name[i] = InFileName[i];
	infile_name[i] = '\0';

	for (i = 0; i < static_cast<int>(ORFoutFile.length()); i++)
		orffile_name[i] = ORFoutFile[i];
	orffile_name[i] = '\0';

	for (i = 0; i < static_cast<int>(AAoutFile.length()); i++)
		aafile_name[i] = AAoutFile[i];
	aafile_name[i] = '\0';

	std::ifstream infile;
	infile.open(infile_name);

	std::ofstream orfoutfile;
	orfoutfile.open(orffile_name);

	std::ofstream aaoutfile;
	aaoutfile.open(aafile_name);

	bool ambiguous_read;
	retained_orf_counter = 0;
	int loci_done = 0;
	while (!infile.eof())
	{
		fastq_read currentread;
		currentread.header = "NULL";
		currentread.sequence = "NULL";
		currentread.qual_header = "NULL";
		currentread.quality = "NULL";
		if (!infile.eof())
		{
			infile.getline(line, 999999);
			currentread.header = line;
		}

		if (!infile.eof())
		{
			infile.getline(line, 999999);
			currentread.sequence = line;
		}

		ambiguous_read = false;
		int original_seq_length;
		if (currentread.sequence.length() > 5 && currentread.header[0] == '>')
		{
				
			std::string tempseq;

			// Make sure there are no weird characters in the sequence
			k = static_cast<int>(currentread.sequence.length());
			original_seq_length = k;
			bool knownbase;
			for (j = 0; j < k; j++)
			{
				knownbase = false;
				if (currentread.sequence[j] == 'A' || currentread.sequence[j] == 'a')
				{
					tempseq = tempseq + 'A';
					knownbase = true;
				}
				if (currentread.sequence[j] == 'G' || currentread.sequence[j] == 'g')
				{
					tempseq = tempseq + 'G';
					knownbase = true;
				}
				if (currentread.sequence[j] == 'T' || currentread.sequence[j] == 't')
				{
					tempseq = tempseq + 'T';
					knownbase = true;
				}
				if (currentread.sequence[j] == 'C' || currentread.sequence[j] == 'c')
				{
					tempseq = tempseq + 'C';
					knownbase = true;
				}
				if (!knownbase)
				{
					tempseq = tempseq + 'N';
					ambiguous_read = true;
				}
			}
			currentread.sequence = tempseq;

			// Reverse Complement of Sequence
			currentread.reversecomplement.clear();
			k = static_cast<int>(currentread.sequence.length());
			if (k < original_seq_length)
				ambiguous_read = true;
			for (j = 0; j < k; j++)
			{
				if (currentread.sequence[j] == 'A' || currentread.sequence[j] == 'a')
					tempseq = tempseq + 'T';
				if (currentread.sequence[j] == 'G' || currentread.sequence[j] == 'g')
					tempseq = tempseq + 'C';
				if (currentread.sequence[j] == 'T' || currentread.sequence[j] == 't')
					tempseq = tempseq + 'A';
				if (currentread.sequence[j] == 'C' || currentread.sequence[j] == 'c')
					tempseq = tempseq + 'G';
				if (currentread.sequence[j] == 'N' || currentread.sequence[j] == 'n')
					tempseq = tempseq + 'N';
			}

			for (j = 0; j < k; j++)
			{
				currentread.reversecomplement = currentread.reversecomplement + tempseq[k - j - 1];
			}

			// Figure out the longest open reading frame 

			int longestORF[6];     // one for each reading frame.  0-2: forward, 3-5: reverse.
			int longestORFstart[6];
			int longestORFend[6];

			int seqlngth;

			int verylongestorf;
			int verylongestorfstart;
			int verylongestorfend;
			bool verylongestorfforward;

			for (j = 0; j < 6; j++)  // go through all the reading frames and get the longest orf from each one
			{
				if (j < 3)
					currentread.AllORFs(false, require_ATG, j, 6);
				else
					currentread.AllORFs(true, require_ATG, j - 3, 6);
				longestORF[j] = -1;
				longestORFstart[j] = -1;
				longestORFend[j] = -1;
				for (k = 0; k < currentread.all_Orfs_counter; k++)
				{
					if (currentread.all_Orfs_lengths[k] > longestORF[j])
					{
						longestORF[j] = currentread.all_Orfs_lengths[k];
						longestORFstart[j] = currentread.all_Orfs_starts[k];
						longestORFend[j] = currentread.all_Orfs_ends[k];
					}
				}
			} // j

			// Figure out the longest ORF for the current contig in any reading frame

			verylongestorf = -1;
			verylongestorfstart = -1;
			verylongestorfend = -1;
			verylongestorfforward = true;
			for (j = 0; j < 6; j++)
			{
				if (longestORF[j] > verylongestorf)
				{
					verylongestorf = longestORF[j];
					verylongestorfstart = longestORFstart[j];
					verylongestorfend = longestORFend[j];
					if (j < 3)
						verylongestorfforward = true;
					else
						verylongestorfforward = false;
				}
			}

			currentread.ORFstart = verylongestorfstart;
			currentread.ORFend = verylongestorfend;
			currentread.ORFlength = verylongestorf;
			currentread.ORFforward = verylongestorfforward;

			currentread.ORFseq.clear();
			seqlngth = static_cast<int>(currentread.sequence.length());
			if (currentread.ORFlength > 0) // Make sure there is at least one orf
			{
				if (currentread.ORFforward)
				{
					for (j = 0; j < currentread.ORFlength; j++)
					{
						k = j + currentread.ORFstart;
						if (k < seqlngth)
							currentread.ORFseq = currentread.ORFseq + currentread.sequence[k];
					}
				}
				else
				{
					seqlngth = static_cast<int>(currentread.reversecomplement.length());
					for (j = 0; j < currentread.ORFlength; j++)
					{
						k = j + currentread.ORFstart;
						if (k < seqlngth)
							currentread.ORFseq = currentread.ORFseq + currentread.reversecomplement[k];
					}
				}
			}

			// Write the longest ORF to the output file

			std::string temp_aa_seq;
			if (!ambiguous_read || !remove_Ns)
			{
				if (static_cast<int>(currentread.ORFseq.length()) >= min_orf_length)
				{
					currentread.header[0] = '>';
					if (currentread.ORFforward)
						orfoutfile << currentread.header << "_for_" << currentread.ORFstart + 1 << "_" << currentread.ORFend + 1 << "\n";
					else
						orfoutfile << currentread.header << "_rev_" << currentread.ORFstart + 1 << "_" << currentread.ORFend + 1 << "\n";
					orfoutfile << currentread.ORFseq << "\n";

					temp_aa_seq = translate(currentread.ORFseq);

					if (remove_stop) // check to see if the last aa is a stop and remove it if parameter set
					{
						if (temp_aa_seq.back() == '*')
						{
							temp_aa_seq = temp_aa_seq.substr(0, temp_aa_seq.length() - 1);
						}
					}

					if (currentread.ORFforward)
						aaoutfile << currentread.header << "_for_" << currentread.ORFstart + 1 << "_" << currentread.ORFend + 1 << "\n";
					else
						aaoutfile << currentread.header << "_rev_" << currentread.ORFstart + 1 << "_" << currentread.ORFend + 1 << "\n";
					aaoutfile << temp_aa_seq << "\n";

					retained_orf_counter++;
				}
			}

				loci_done++;
			
		}
			
		

	} // end of while
	

	infile.close();
	orfoutfile.close();


	delete[] line;
	delete[] infile_name;
	delete[] orffile_name;

	std::cout << "\n\nDone!\n\n";
	std::cout << retained_orf_counter << " reads retained.\n";
	std::cout << loci_done << " total reads processed.\n";

	if (interactive_mode)
	{
		std::cout << "\nEnter any integer to exit!\n";
		std::cin >> i;
	}
	return 0;

}