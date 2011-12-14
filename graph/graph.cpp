#include "graph.hpp"

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

ExonNode * Graph::add_exon_node(int start, int end, int origin){
	ExonNode * ret;
	if ((ret=find_exon_node(start, end, origin))!=GRAPH_EXON_NOT_FOUND){
		return ret;
	}
	ret=new ExonNode(start, end, origin);
	nodes.insert(ret);
	return ret;
}

int Graph::add_edge(ExonNode * first, ExonNode *second, int num_reads){
	if (first->origin==EXON_NO_INFO || second->origin==EXON_NO_INFO ||
		(first->origin==second->origin)
		){
		first->out_nodes.push_back(second);
		first->out_edge_weight.push_back(num_reads);
		second->in_nodes.push_back(first);
		second->in_edge_weight.push_back(num_reads);
	}
	else {
		fprintf(stdout, "Invalid Edge!\n");
		exit(0);
	}
}


int Graph::add_edge(ExonNode * first, ExonNode *second, vector<BamAlignment *> reads){
	if (first->origin==EXON_NO_INFO || second->origin==EXON_NO_INFO ||
		(first->origin==second->origin)
		){
		for (int i=0;i<reads.size();i++){
			// TODO: finish this function
		}
	}
	else {
		fprintf(stdout, "Invalid Edge!\n");
		exit(0);
	}
}


int Graph::dump_graph(FILE *foutput){
	fprintf(foutput,"digraph G{\n");
	vector<ExonNode *> queue;
	set<ExonNode *> checklist;
	int top=0;
	queue=get_starting_nodes();
	while(queue.size()!=top){
		ExonNode *current=queue[top];
		top++;
		string current_node_name=current->detach();
		for (auto j=current->out_nodes.begin();
			 j!=current->out_nodes.end();
			 j++){
			string next_node_name=(*j)->detach();
			fprintf(foutput,"\t%s -> %s\n", 
					current_node_name.c_str(), next_node_name.c_str());
			if (checklist.find((*j))!=checklist.end()){
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
