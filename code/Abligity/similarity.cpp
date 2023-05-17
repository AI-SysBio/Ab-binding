#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/serialization/map.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/program_options.hpp>

#include <iomanip>
#include <list>
#include <algorithm>
#include <numeric>
#include <map>
#include <unordered_map>
#include <cmath>
#include <utility>

using namespace std;

void save(map<string,double> const& obj, string const& fname) {
    std::ofstream ofs(fname, std::ios::binary);
    {
        boost::iostreams::filtering_ostreambuf fos;

        // push the ofstream and the compressor
        fos.push(boost::iostreams::zlib_compressor(boost::iostreams::zlib::best_compression));
        fos.push(ofs);

        // start the archive on the filtering buffer:
        boost::archive::binary_oarchive bo(fos);
        bo << obj;
    }
}

map<string,double> load(string const& fname) {
	map<string,double> obj;
	
	std::ifstream ifs(fname, std::ios::binary);
	try{
		//boost::iostreams::filtering_istreambuf fis;
		boost::iostreams::filtering_streambuf<boost::iostreams::input> fis;
		
		// push the ifstream and the decompressor
		fis.push(boost::iostreams::zlib_decompressor());
		fis.push(ifs);

		// start the archive on the filtering buffer:
		boost::archive::binary_iarchive bi(fis);
		
		bi >> obj;
	}	catch(const std::exception &exc){
		cerr << "Error in reading " + fname << endl;
		cerr << exc.what() << endl;
	}
		return obj;
}


bool fexists(const string& filename) {
  ifstream ifile(filename.c_str());
  return (bool)ifile;
}

double get_similarity(map<string,double> map1,map<string,double> map2,double alpha,double beta){
	double a=0.0,b=0.0,commonab=0.0,similarity=0.0;
	for(auto it1 = map1.begin(); it1!= map1.end(); it1++) {
		auto it2 = map2.find(it1->first);
		if (it2 != map2.end()) {
			if(it1->second < 0.0 || it2->second < 0.0) {cout << it1->first << " "<< it2->second << endl;continue;}
			commonab += min(it1->second,it2->second);
			if(alpha != 0.0) a += max(0.0,it1->second-it2->second);
			if(beta != 0.0) b += max(0.0,it2->second-it1->second);
		}
		else if(alpha != 0.0) {
			a += it1->second;
		}
	}
	if(beta != 0.0) {
		for(auto it2 = map2.begin(); it2!=map2.end();it2++){
			if(map1.find(it2->first) == map1.end()) b+= it2->second;
		}
	}
	//cout <<input1<<"\t"<<input2<<"\t"<<a<<"\t"<<b<<"\t"<<commonab<<"\t"<<commonab/(alpha*a+beta*b+commonab)<<endl;
	similarity = commonab/(alpha*a+beta*b+commonab);
	return similarity;
}
// Format results
void format_results(string output,double** similaritymatrix, int firstdim, int seconddim)
{
	std::ofstream outputfile;
	outputfile.open (output);
	for (int h=0; h<firstdim; ++h){
		for (int w=0; w<seconddim; ++w){
			outputfile << similaritymatrix[h][w] << ",";
		}
		outputfile << "\n";
	}
	outputfile.close();
}


int main(int argc, char** argv){
	boost::program_options::options_description desc("Options");
	desc.add_options()
		("help,h","Help screen")
		("input1,i",boost::program_options::value<string>(),"Input map 1")
		("input2,j",boost::program_options::value<string>(),"Input map 2")
		("list,l",boost::program_options::value<string>(),"Text file with a list of .ht")
		("output,o",boost::program_options::value<string>(),"Output Similarity matrix")
		("alpha,a",boost::program_options::value<double>()->default_value(0.5),"alpha")
		("beta,b",boost::program_options::value<double>()->default_value(0.5),"beta")
		("cache,c",boost::program_options::bool_switch()->default_value(false),"cache .ht");
	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc),vm);
	boost::program_options::notify(vm);
	if(vm.count("help")) {cout << desc << endl; return 0;}
	
	// Parsing options
   double alpha = vm["alpha"].as<double>(), beta = vm["beta"].as<double>();
   bool cache = vm["cache"].as<bool>(); 
   string input1,input2,output="",file_list; int filelistlength;
   
   // Algorithm variables
   double similarity = 0.0; double** similaritymatrix;
   map<string,double> map1, map2; vector<string> files; string fileloc,fileloc1,fileloc2;
   vector<map<string,double>> htcaches;
   int mode = 0; // Determines which mode to run in: 1 = pair; 2 = pairwise; 3 = one-against-all
   
   if (vm.count("output")){ output = vm["output"].as<string>(); }
	if (vm.count("input1")){
		input1 = vm["input1"].as<string>();
		assert(fexists(input1));
		mode += 1;
		map1 = load(input1);
	} 
	if (vm.count("input2")){
		input2 = vm["input2"].as<string>();
		assert(fexists(input2));
		map2 = load(input2);
	} 
	if (vm.count("list")){
		file_list = vm["list"].as<string>();
		assert(fexists(file_list));
		mode += 2;
		// Read list of files
		ifstream file_list_handle(file_list);
		while (file_list_handle.good()) {
			getline(file_list_handle,fileloc); // last item
			if(fileloc.size() < 2 || fileloc.substr(fileloc.size()-2) != "ht"){continue;}
			assert(fexists(fileloc));
			files.push_back(fileloc);
		}
		filelistlength = files.size();
		// Read all heatmaps - potentially using up a massive amount of RAM
		if (cache){
		cout << "Cached ";
		for (auto const& file : files){
			map2 = load(file);
			htcaches.push_back(map2);
		}}
	}
	// Execute according to the mode

	switch (mode) {
		case 1: {
			cout << "Pair comparison\n";
			similarity = get_similarity(map1,map2,alpha,beta);
			cout << similarity << endl;
			break;}
		case 2: {cout << "Pairwise comparison\n";
			if (output==""){ cout << "Must have output destination\n"; return 1;}
			similaritymatrix = new double*[filelistlength];
			for (int i=0; i<filelistlength; i++){
				similaritymatrix[i] = new double[filelistlength];
			}
			
			// Similarity matrix
			for (int fileid1=0; fileid1<filelistlength; fileid1++){
				if (cache) {map1 = htcaches[fileid1];}
				else {
					fileloc1 = files[fileid1];
					map1 = load(fileloc1);
				}
				for (int fileid2=fileid1+1; fileid2<filelistlength; fileid2++){
					if (cache) {map2 = htcaches[fileid2];}
					else {
						fileloc2 = files[fileid2];
						map2 = load(fileloc2);
					}
					
					similarity = get_similarity(map1,map2,alpha,beta);
					similaritymatrix[fileid1][fileid2] = similarity;
					similaritymatrix[fileid2][fileid1] = similarity;
			}}
			format_results(output,similaritymatrix,filelistlength,filelistlength);
			break;}
		case 3: {cout << "One-vs-all comparison\n";
			if (output==""){ cout << "Must have output destination\n"; return 1;}
			similaritymatrix = new double*[1];
			similaritymatrix[0] = new double[filelistlength];
			
			// Similarity matrix
			for (int fileid=0; fileid<filelistlength; fileid++){
				if (cache) {map2 = htcaches[fileid];}
				else {
					fileloc = files[fileid];
					map2 = load(fileloc);
				}
				similarity = get_similarity(map1,map2,alpha,beta);
				similaritymatrix[0][fileid] = similarity;
			}
			format_results(output,similaritymatrix,1,filelistlength);
			break;}
		default: {cout << "Invalid input combinations";
			return 1;} 
    }
    

	
	return 0;
	
}
