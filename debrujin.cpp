#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<unordered_map>
#include<unordered_set>
using namespace std;

#define N 8500
#define L_s 100 // length of short read
#define K 29	// length of k-mers

vector<string> short_read(2 * N);

// node in de brujin graph
class Node{
public:
	Node() {}
	Node(string s) {k_mer = s;}
	string k_mer;
	unordered_set<string> ingree;
	unordered_set<string> outgree;
};

unordered_map<string, Node> map;
vector<string> out_data;

void read_file(){
	ifstream short_1("data1\\short_1.fasta");
	ifstream short_2("data1\\short_2.fasta");
	int i = 0;
	string tmp;
	while(getline(short_1, tmp)){
		if(i & 1){
			short_read[(i >> 1)] = tmp;
		}
		i++;
	}
	while(getline(short_2, tmp)){
		if(i & 1){
			short_read[(i >> 1)] = tmp;
		}
		i++;
	}
}


void build_brujin(){
	for(string& s : short_read){
		vector<string> kmers(L_s - K + 1);
		for(int i = 0; i <= L_s - K; i++){
			kmers[i] = s.substr(i, K);
		}
		for(int i = 0; i < L_s - K; i++){
			string s_1 = kmers[i];
			string s_2 = kmers[i + 1];
			
			if(map.count(s_1)){
				Node& t = map[s_1];
				if(!t.outgree.count(s_2)){
					t.outgree.insert(s_2);
				}
			}else{
				map[s_1] = Node(s_1);
				map[s_1].outgree.insert(s_2);
			}
			
			if(map.count(s_2)){
				Node& t = map[s_2];
				if(!t.ingree.count(s_1)){
					t.ingree.insert(s_1);
				}
			}else{
				map[s_2] = Node(s_2);
				map[s_2].ingree.insert(s_1);
			}
		}
	}
}

void maximal_nonbranch_paths(){
	for(pair<const string, Node>& p : map){
		Node node = p.second;
		// node is not 1-in-1-out
		if(node.ingree.size() != 1 || node.outgree.size() != 1){
			for(string out : node.outgree){
				// find a maximal_nonbranch_path started from this node
				string path = p.first;
				string tmp_s = out;
				
				while(1){
					path += tmp_s[K - 1];
					Node& tmp_node = map[tmp_s];
					if(tmp_node.ingree.size() != 1 || tmp_node.outgree.size() != 1) break;
					tmp_s = (*tmp_node.outgree.begin());
				}
				if(path.length() > L_s)
					out_data.push_back(path);
			}
		}
	}
}

void write_file(){
	ofstream writer("contig.fasta", ios::out);
	int i = 0;
	for(string& s : out_data){
		writer << ">contig_" << i << endl;
		writer << s << endl; 
		i++;
	}
}

int main(){
	read_file();
	build_brujin();
	maximal_nonbranch_paths();
	write_file();
	cout << out_data.size() << endl;
	return 0;	
} 
