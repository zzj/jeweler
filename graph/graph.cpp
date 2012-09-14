#include "graph.hpp"
#include "../constants.hpp"
Graph::~Graph(){
	for (auto i=nodes.begin();i!=nodes.end();i++){
		delete (*i);
	}
}

ExonNode * Graph::find_exon_node(int start, int end, int origin){
	auto i=nodes.begin();
	for (; i!=nodes.end();i++){
		if ((*i)->start==start &&
			(*i)->end==end &&
			(*i)->origin==origin){
			return (*i);
		}
		
	}
	return GRAPH_EXON_NOT_FOUND;
}

ExonNode * Graph::add_exon_node(int start, int end, int origin, 
								set<JewelerAlignment *> &reads,
								set< JewelerAlignment *> &allele_reads){
	ExonNode * ret;
	if ((ret=find_exon_node(start, end, origin))==GRAPH_EXON_NOT_FOUND){
		ret=new ExonNode(start, end, origin);
	}
	ret->insert_allele_reads(allele_reads);
	ret->insert_reads(reads);
	nodes.insert(ret);
	return ret;
}

void Graph::add_edge(ExonNode * first, ExonNode *second, int num_reads){
	if (first->origin==EXON_NO_INFO || second->origin==EXON_NO_INFO ||
		(first->origin==second->origin)
		){
		if (first->out_nodes_checklist.find(second)==first->out_nodes_checklist.end()){
			first->out_nodes.push_back(second);
			first->out_edge_weight.push_back(num_reads);
			first->out_nodes_checklist.insert(second);
		}
		if (second->in_nodes_checklist.find(first)==second->in_nodes_checklist.end()){
			second->in_nodes.push_back(first);
			second->in_edge_weight.push_back(num_reads);
			second->in_nodes_checklist.insert(first);
		}
	}
	else {
		fprintf(stdout, "Invalid Edge!\n");
		exit(0);
	}
}


void Graph::add_edge(ExonNode * first, ExonNode *second, vector<JewelerAlignment *> reads){
	if (first->origin==EXON_NO_INFO || second->origin==EXON_NO_INFO ||
		(first->origin==second->origin)
		){
		for (size_t i=0;i<reads.size();i++){
			// TODO: finish this function
		}
	}
	else {
		fprintf(stdout, "Invalid Edge!\n");
		exit(0);
	}
}

void Graph::dump_graph(FILE *foutput){
	fprintf(foutput,"digraph G{\n");
	vector<ExonNode *> queue;
	set<ExonNode *> checklist;
	int top=0;
	for (auto i=nodes.begin();i!=nodes.end();i++){
		string color="black";
		if ((*i)->origin!=EXON_NO_INFO){
			if ((*i)->origin==EXON_MATERNAL){
				color="pink";
			}
			else {
				color="lightblue";
			}
		}
		fprintf(foutput,"\t\"%s\" [color=%s];\n",(*i)->detach().c_str(),color.c_str());
	}
	queue=get_starting_nodes();
	
	while(queue.size() != top){
		ExonNode *current=queue[top];
		top++;
		string current_node_name=current->detach();

		for (auto j=current->out_nodes.begin();
			 j!=current->out_nodes.end();
			 j++){
			string next_node_name=(*j)->detach();
			fprintf(foutput,"\t\"%s\" -> \"%s\";\n", 
					current_node_name.c_str(), next_node_name.c_str());
			if (checklist.find((*j))==checklist.end()){
				queue.push_back(*j);
				checklist.insert(*j);
			}
		}
	}
	fprintf(foutput,"}\n");

}

vector<ExonNode *> Graph::get_starting_nodes(){
	vector<ExonNode *> ret;
	auto i=nodes.begin();
	for (;i!=nodes.end();i++){
		if ((*i)->in_nodes.size()==0){
			ret.push_back(*i);
		}
	}
	return ret;
}

int Graph::get_all_paths(vector<Path> & records){
	records.clear();
	vector<ExonNode *> starting_nodes;
	starting_nodes=get_starting_nodes();
	vector<ExonNode *> path;
	for (auto i =starting_nodes.begin(); i!=starting_nodes.end(); i++){
		traverse_graph((*i),  path, records, (*i)->origin);
	}
	return 0;
}
	
int Graph::traverse_graph( ExonNode * current_node, vector<ExonNode *> &path, vector<Path> &records, int origin){
	
	if (origin!=EXON_NO_INFO) {
		if (current_node->origin!=EXON_NO_INFO &&
			current_node->origin!=origin){
			return 0;
		}
	}

	if (current_node->origin!=EXON_NO_INFO) origin=current_node->origin;

	path.push_back(current_node);
	
	if (current_node->out_nodes.size()==0){
		records.push_back(Path(path));
	}
	else {
		for (auto i=current_node->out_nodes.begin(); 
			 i!=current_node->out_nodes.end();
			 i++){
			traverse_graph((*i), path, records, origin);
		}
	}
	path.pop_back();
	return 0;
}
