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
#include <locale>   
#include <cassert>
#include <regex>

using namespace std;
double** distance_matrix;
map<string,string> vtype;
list<string> allowed_residuegroups = {"5","7a","7b"};
// Triangulate
struct Pip{
	string pymol, chain; double resnum;
	vector<char> atomtypes;
	int id;
	double xcoord,ycoord,zcoord;
	Pip(string pipid,string pippymol,string pipatomtypes,string pipxcoord,string pipycoord,string pipzcoord){
		id = stoi(pipid);
		pymol = pippymol;
		copy(pipatomtypes.begin(),pipatomtypes.end(),back_inserter(atomtypes));
		xcoord = stod(pipxcoord);		
		ycoord = stod(pipycoord);		
		zcoord = stod(pipzcoord);
		smatch sm; 
		regex resnumre("(\\w+)/(.*)/\\w+");
		regex_match (pymol,sm,resnumre);
		chain = sm[1].str();
		resnum = stoi(sm[2].str());
	}
	~Pip(){}
};


bool sortpiplist(const Pip& a, const Pip& b) {
	return (a.chain<b.chain) || ((a.chain==b.chain) && (a.resnum < b.resnum));
}

// Look up table for vertex type
/*
[A] Aliphatic: A,G,I,L,M,V
[B] Polar: C,N,P,Q,S,T
[C] Aromatic: F,W,Y
[D] Negatively charged: D,E
[E] Positively charged: K,R,H

*/
const map<string,string> vtype5 = {
	{ "AA","a" },{ "AB","b" },{ "AC","c" },{ "AD","d" },{ "AE","e" },
					 { "BB","f" },{ "BC","g" },{ "BD","h" },{ "BE","i" },
					 				  { "CC","j" },{ "CD","k" },{ "CE","l" },
					 				  				   { "DD","m" },{ "DE","n" },
					 				  				  				    { "EE","o" }
};
/*
[A]cidic: D,E
[B]asic: H, K, R
[H]ydrophobic: I,L,M,P,V
A[M]ine: N, Q
[N]ucleophilic: C, S, T
A[R]omatic: F,W,Y
[S]mall: A, G
*/
const map<string,string> vtype7a = {
	{ "AA","a" },{ "AB","b" },{ "AH","c" },{ "AM","d" },{ "AN","e" },{ "AR","f" },{ "AS","g" },
					 { "BB","h" },{ "BH","i" },{ "BM","j" },{ "BN","k" },{ "BR","l" },{ "BS","m" },
					 				  { "HH","n" },{ "HM","o" },{ "HN","p" },{ "HR","q" },{ "HS","r" },
					 				  				   { "MM","s" },{ "MN","t" },{ "MR","u" },{ "MS","v" },
					 				  				  				    { "NN","w" },{ "NR","x" },{ "NS","y" },
					 				  				  				  				     { "RR","z" },{ "RS","+" },
					 				  				  				  				  				      { "SS","-" }
};
/*
[A] Aliphatic: A, G, I, L, P, V
[B] Aromatic: F, W, Y
[C] Sulfur: C, M
[D] Hydroxyl: S, T
[E] Basic: H, K, R
[F] Acidic: D, E
[G] Amide: N, Q
*/
const map<string,string> vtype7b = {
	{ "AA","a" },{ "AB","b" },{ "AC","c" },{ "AD","d" },{ "AE","e" },{ "AF","f" },{ "AG","g" },
					 { "BB","h" },{ "BC","i" },{ "BD","j" },{ "BE","k" },{ "BF","l" },{ "BG","m" },
					 				  { "CC","n" },{ "CD","o" },{ "CE","p" },{ "CF","q" },{ "CG","r" },
					 				  				   { "DD","s" },{ "DE","t" },{ "DF","u" },{ "DG","v" },
					 				  				  				    { "EE","w" },{ "EF","x" },{ "EG","y" },
					 				  				  				  				     { "FF","z" },{ "FG","+" },
					 				  				  				  				  				      { "GG","-" }
};


// Sort pips by atom types
struct piptri{
	char atomtype;
	Pip pip;
	piptri(int atomtype_, Pip pip_) : atomtype(atomtype_), pip(pip_) {}
	bool operator<(const piptri& p) const { return atomtype < p.atomtype;}
};

struct sortpipfunc{
	inline bool operator() (const piptri* atpip1, const piptri* atpip2)
    {
        return (*atpip1 < *atpip2);
    }
};

double get_dist(Pip& pip1, Pip& pip2) {
	return sqrt(pow(pip1.xcoord-pip2.xcoord,2)+pow(pip1.ycoord-pip2.ycoord,2)+pow(pip1.zcoord-pip2.zcoord,2));
};

void calc_distance_matrix(vector<Pip> Pipvec){
	int length = Pipvec.size();
	// initialize double arrays
	distance_matrix = new double*[length];
	for (int i=0; i<=length-1;i++){
		distance_matrix[i] = new double[length];
	}
	int id1, id2; double distance;
	for(int i=0; i<=length-1;i++){
		Pip pip1 = Pipvec[i]; 
		id1 = pip1.id;
		for(int j=i+1; j<=length-1;j++) {
			Pip pip2 = Pipvec[j]; 
			id2 = pip2.id;
			distance = get_dist(pip1,pip2);
			if (distance < 1.5){distance = -1;}
			distance_matrix[id1][id2] = distance;
			distance_matrix[id2][id1] = distance;
	}}
}

bool checkdist(double ds) {return ( ds>=1.5 );}

void diffpip(Pip& pip1, Pip& pip2, double result[3]) {
	result[0] = pip2.xcoord - pip1.xcoord;
	result[1] = pip2.ycoord - pip1.ycoord;
	result[2] = pip2.zcoord - pip1.zcoord;
}
double dotproduct(double coord1[3], double coord2[3]) {
	double product = 0.0;
	for(int i=0; i<3; i++) product += coord1[i] * coord2[i];
	return product;
}
void crossproduct(double coord1[3], double coord2[3], double result[3]) {
	result[0] = coord1[1]*coord2[2] - coord1[2]*coord2[1];
	result[1] = coord1[0]*coord2[2] - coord1[2]*coord2[0];
	result[2] = coord1[0]*coord2[1] - coord1[1]*coord2[0];
}

string get_distbin(vector<double> distbins, double binsize) {
	string distbin;
	for(double eachdistbin : distbins) {distbin += "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[int(round(eachdistbin/binsize))-1];}
	return distbin;
}


map<string,double> eachtri(vector<Pip> Pipvec, bool gaussian, double height, double binsize) {
	map<string,double> hm; string hk; char v;
	
	// Iterate through atom types in these pips
	for(int i_=0; i_<=Pipvec[0].atomtypes.size()-1;i_++){
	for(int j_=0; j_<=Pipvec[1].atomtypes.size()-1;j_++){
	for(int k_=0; k_<=Pipvec[2].atomtypes.size()-1;k_++){
		
		vector<piptri*> piptrivec;
		piptrivec.push_back(new piptri(Pipvec[0].atomtypes[i_],Pipvec[0]));
		piptrivec.push_back(new piptri(Pipvec[1].atomtypes[j_],Pipvec[1]));
		piptrivec.push_back(new piptri(Pipvec[2].atomtypes[k_],Pipvec[2]));
		sort(piptrivec.begin(),piptrivec.end(),sortpipfunc());
		
		// Generate 3 edges for the triangle

		string vts; double dist;
		vector<double> diststores, point_heights; vector<string> distbins;
		tuple<string,double> vdiststore;
		vector< tuple<string,double> > vdiststores;
		// This looks ugly
		for(int i__=0; i__<=piptrivec.size()-2;i__++)
		for(int j__=i__+1; j__<=piptrivec.size()-1;j__++) {
			//dist = get_dist(piptrivec[i__]->pip,piptrivec[j__]->pip);
			dist = distance_matrix[piptrivec[i__]->pip.id][piptrivec[j__]->pip.id];
			if (dist < 0) {return hm;}
			string vt = string() + piptrivec[i__]->atomtype + piptrivec[j__]->atomtype;
			cout << vt << endl;
			assert(vtype.find(vt)!=vtype.end());
			vdiststore = make_tuple(vtype.find(vt)->second,dist);
			vdiststores.push_back(vdiststore);
		}
		// Sort vertex types and edge types
		sort(vdiststores.begin(),vdiststores.end());
		for(int vdiststorei = 0; vdiststorei < vdiststores.size(); vdiststorei++ ){
			vts += get<0>(vdiststores[vdiststorei]);
			diststores.push_back(get<1>(vdiststores[vdiststorei]));
		}
		// Add the main one
		point_heights.push_back(1);
		vector<double> mainbin{double(diststores[0]),double(diststores[1]),double(diststores[2])};

		distbins.push_back(get_distbin(mainbin, binsize));
		// Then its neighbours
		if(gaussian) {
			vector< vector <double> > bounds(3);
			double dist, x, point_height;
			// Reinitialize for comparison later
			vector<double> maindist{round(diststores[0]/binsize),round(diststores[1]/binsize),round(diststores[2]/binsize)};
			double low, high;
			for (auto dsi = 0; dsi < diststores.size(); dsi++){
				double ds = diststores[dsi];
				low = round(ds/binsize)-1;
				high = round(ds/binsize)+1;
				
				if ( checkdist(ds-binsize)) bounds[dsi].push_back(low);
				if ( checkdist(ds)) bounds[dsi].push_back(round(ds/binsize));
				if ( checkdist(ds+binsize)) bounds[dsi].push_back(high);
				
			}
			for(auto a: bounds[0]){
			for(auto b: bounds[1]){
			for(auto c: bounds[2]){
				vector<double> gaussbins{a,b,c};
				if (gaussbins == maindist) {continue;}
				point_heights.push_back(0.5*height); // Step
				distbins.push_back(get_distbin(gaussbins, 1.0));// Because the bin has already been divided by binsize
			}}}
			bounds.clear();
		}

		for (int i=0; i < distbins.size(); i++){
			hk = vts + distbins[i];
			hm[hk] += point_heights[i];
		}
		distbins.clear();
	}}}
	return hm;
	
}

map<string,double> removebinsize(map<string,double>map1, double upper, double lower) {
	map<string,double> map1p;
	for(auto it1=map1.begin();it1!=map1.end();it1++) {
		if (it1->second > lower && it1->second < upper)	map1p[it1->first] = it1->second;
	}
	return map1p;
}

// IO
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

	std::ifstream ifs(fname, std::ios::binary);
	{
		//boost::iostreams::filtering_istreambuf fis;
		boost::iostreams::filtering_streambuf<boost::iostreams::input> fis;
		
		// push the ifstream and the decompressor
		fis.push(boost::iostreams::zlib_decompressor());
		fis.push(ifs);

		// start the archive on the filtering buffer:
		boost::archive::binary_iarchive bi(fis);
		map<string,double> obj;
		bi >> obj;
		return obj;
	}
}
bool fexists(const string& filename) {
  ifstream ifile(filename.c_str());
  return (bool)ifile;
}

bool is_numeric(string rownum){
	for (auto n:rownum){
		if(isdigit(n)) return false;
	}
	return true;
}

int main(int argc, char** argv){

	boost::program_options::options_description desc("Options");
	desc.add_options()
		("help,h","Help screen")
		("input,i",boost::program_options::value<string>(),"Input PIP list")
		("output,o",boost::program_options::value<string>(),"Output hashtable")
		("verbose,v",boost::program_options::value<bool>()->default_value(false),"Verbose")
		("gaussian,g",boost::program_options::value<bool>()->default_value(false),"Introduce gaussian counts to neighbours")
		("rmbin,r",boost::program_options::value<bool>()->default_value(false),"If remove bins with counts outside of defined range (define the range with minbin and maxbin options")
		("minbin,n",boost::program_options::value<double>()->default_value(0),"If remove bins with count less than input")
		("maxbin,m",boost::program_options::value<double>()->default_value(10000),"If remove bins with count more than input")
		("binsize,b",boost::program_options::value<double>()->default_value(1.0),"Grid/Bin size")
		("residuegroup,x",boost::program_options::value<string>()->default_value("7b"),"Residue Grouping to be used");
	
	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc),vm);
	boost::program_options::notify(vm);
	if(vm.count("help")) {cout << desc << endl; return 0;}
	
	bool verbose = vm["verbose"].as<bool>();
	bool gaussian = vm["gaussian"].as<bool>();
	bool rmbin = vm["rmbin"].as<bool>();
	double minbin = vm["minbin"].as<double>();
	double maxbin = vm["maxbin"].as<double>();
	double binsize = vm["binsize"].as<double>();
	string inputpip(vm["input"].as<string>());
	string residuegroup(vm["residuegroup"].as<string>());
	
	assert(fexists(inputpip));
	assert((std::find(allowed_residuegroups.begin(), allowed_residuegroups.end(), residuegroup) != allowed_residuegroups.end()));
	// Residue group
	if(residuegroup=="5"){ vtype = vtype5; }
	else if(residuegroup=="7a") { vtype = vtype7a; }
	else if(residuegroup=="7b") { vtype = vtype7b; }
	else {cout << "Error in residue group selection\n"; return 1;}

	ifstream pipfile(inputpip);
	map<string,double> hm, localhm; // hash table 
	// IO
	string outputht = vm["output"].as<string>();
	// Var
	double height = 1.0;
	string rownum,pippymol,pipatomtypes,pipxcoord,pipycoord,pipzcoord;
	list<Pip> Piplist; vector<Pip> TriNode;
	
	while (pipfile.good()) {
		getline(pipfile,rownum,',');
		getline(pipfile,pippymol,',');
		getline(pipfile,pipatomtypes,',');
		getline(pipfile,pipxcoord,',');
		getline(pipfile,pipycoord,',');
		getline(pipfile,pipzcoord); // last item
		if(is_numeric(rownum)){continue;}
		Pip piplocal(rownum,pippymol,pipatomtypes,pipxcoord,pipycoord,pipzcoord);
		Piplist.push_back(piplocal);
	}
	
	// List to Vector
	vector<Pip> Pipvec{begin(Piplist),end(Piplist)};
	sort(Pipvec.begin(),Pipvec.end(),sortpiplist);

	// Iterate through combinations
	int N = Pipvec.size();
	assert(N>=3);
	
	// Calculate distance matrix
	calc_distance_matrix(Pipvec);

	// Triangulate
	
	for(int i=0; i<=N-3;i++){ 
	for(int j=i+1; j<=N-2;j++){
	for(int k=j+1; k<=N-1;k++){
		TriNode.push_back(Pipvec[i]);
		TriNode.push_back(Pipvec[j]);
		TriNode.push_back(Pipvec[k]);
		localhm = eachtri(TriNode, gaussian, height, binsize);
		for (const auto& kv: localhm){
			hm[kv.first] += kv.second;
		}
		localhm.clear();
		TriNode.clear();
	}}}
	if (rmbin) {hm = removebinsize(hm,maxbin,minbin);}
	// MapIO
	save(hm,outputht);
	
	if (verbose){
		cout << "Hash table is written to "<< outputht << endl;
		cout<< "End"<<endl;
	}

	return 0;
	
}

