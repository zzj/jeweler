
#include "sewing_machine.hpp"
#include "../proto/jeweler.pb.h"
#include "../zleveldb.hpp"
#include <memory>
#include <boost/filesystem.hpp>
using boost::filesystem::path;
using boost::filesystem::create_directory;
using boost::filesystem::remove_all;
using boost::filesystem::exists;


locator::locator(const JewelerAlignment &al, const RefVector &rv) {
	this->data.set_chr(rv[al.RefID].RefName);
	this->data.set_position(al.Position);
	this->data.set_cigar_string(get_cigar_string(al));
    int edit_distance;
	al.GetTag("NM", edit_distance);
    this->data.set_distance(edit_distance);
    this->data.set_is_first(al.IsFirstMate());
}

void locator::dump(Jeweler::SewingMachineData::Locator *new_data) {
    new_data->CopyFrom(this->data);
    return ;
}

bool locator::is_first() {
    return this->data.is_first();
}
bool compare_locator::operator()(const locator *a, const locator *b) {
	if (a->data.chr() != b->data.chr()) 
        return a->data.chr() < b->data.chr();
	else return a->data.position()<b->data.position();
}

SewingMachine::SewingMachine() {
}


SewingMachine::SewingMachine(BamReader &reader, string &db_folder) {
	this->initialize(reader, db_folder);
}


void SewingMachine::initialize(BamReader &reader, string &db_folder) {
    this->db_folder = db_folder;
	BamInfo::initialize(reader);
    this->load_zdb(db_folder);
    this->zdb->clear();
	//this->num_multiple_reads_per_chr.resize(references.size(),0);
}

void SewingMachine::add_read(const string &name,
                             shared_ptr<Jeweler::SewingMachineData> smd) {
    if (this->is_multiple_alignment(smd)) {
        this->multiple_alignment_names.insert(name);
    }
}

void SewingMachine::load_zdb(string &db_file){
    this->zdb = new ZLevelDB(db_file);
    this->db_folder = db_file;
    if (!exists(this->snapshot_file())) {
        this->snapshot();
    }
    else {
        this->load_snapshot();
    }
}

string SewingMachine::snapshot_file() {
    return this->db_folder + "/sewing_machine.data";
}

void SewingMachine::snapshot() {
    this->zdb->foreach<Jeweler::SewingMachineData, SewingMachine >(this, &SewingMachine::add_read);
    fstream out(this->snapshot_file(),
                ios::out | ios::binary | ios::trunc);
    for (auto i = this->multiple_alignment_names.begin();
         i != this->multiple_alignment_names.end();
         i ++) {
        unique_ptr<Jeweler::String> s(new Jeweler::String);
        s->set_data(*i);
        write_protobuf_data(&out, s.get());
    }
    this->initialized_size = this->multiple_alignment_names.size();
}

void SewingMachine::load_snapshot() {
    fstream input(this->snapshot_file(),
                  ios::in | ios::binary);
    unique_ptr<Jeweler::String> s(new Jeweler::String);
    while (load_protobuf_data(&input, s.get()) != -1) {
        this->multiple_alignment_names.insert(s->data());
    }
}

SewingMachine::~SewingMachine() {
    if (this->initialized_size != this->multiple_alignment_names.size())
        this->snapshot();
    delete this->zdb;
}


string SewingMachine::get_reference_sequence(FastaReference &fr, JewelerAlignment &al) {
	string ret;
	// TODO : read from reference genome to get the mismatch position
	// or other useful information.
	return ret;
}


int SewingMachine::add_alignment(const JewelerAlignment &al) {
    unique_ptr<locator>l(new locator(al, this->references));
    shared_ptr<Jeweler::SewingMachineData> smd = \
        this->zdb->get<Jeweler::SewingMachineData>(al.Name);
    if (smd.get() == NULL) {
        smd = shared_ptr<Jeweler::SewingMachineData>(new Jeweler::SewingMachineData);
    }
    l->dump(smd->add_locator());
    this->zdb->set<Jeweler::SewingMachineData>(al.Name, smd.get());
    this->add_read(al.Name, smd);
  	return 0;
}


int SewingMachine::build_alignment_connection_map_core() {
	// auto i = seqs.begin();
	// size_t j,k;
	// alignment_connection_map.clear();
	// for ( ; i != seqs.end(); i++) {
	// 	if (i->second.size() < 2) continue;
	// 	for (j = 0; j < i->second.size(); j++) {
	// 		// statistics
	// 		num_multiple_reads_per_chr[i->second[j]->ref_id]++;

	// 		// map
	// 		for (k = 0;k < i->second.size(); k++) {
	// 			if (j != k)
	// 				alignment_connection_map[i->second[j]].push_back(i->second[k]);
	// 		}
	// 	}
	// }
	return 0;
}


int SewingMachine::output_alignment_connection_map_core(FILE * file) {
	// auto i=alignment_connection_map.begin();
	// for (; i!=alignment_connection_map.end(); i++) {
	// 	size_t j=0;
	// 	for (; j!=i->second.size(); j++) {
	// 		fprintf(file, "%s\t%d\t%s\t%d\n",
	// 				references[i->first->ref_id].RefName.c_str(),
	// 				i->first->position,
	// 				references[i->second[j]->ref_id].RefName.c_str(),
	// 				i->second[j]->position);
	// 	}
	// }
	return 0;
}


int SewingMachine::output_alignment_connection_map(FILE * file) {
	build_alignment_connection_map_core();
	output_alignment_connection_map_core(file);
	return 0;
}


bool SewingMachine::is_multiple_alignment(const string &name) {
    return multiple_alignment_names.find(name) != multiple_alignment_names.end();
}

bool SewingMachine::is_multiple_alignment(const shared_ptr<Jeweler::SewingMachineData> smd) {
    bool has_first = false, has_second = false;
    for (int i = 0; i < smd->locator_size(); i++) {
        if (smd->locator(i).is_first()) has_first = true;
        else has_second =true;
    }
    if ((has_first && !has_second) || (!has_first && has_second))
        return smd->locator_size() >= 2;
    return smd->locator_size() >= 3;
}
