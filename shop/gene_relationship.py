import constants

class GeneRelationship:
    def get_num_consistent_mismatches(self, sg, origin, gene):
        target_ret = set()
        origin_ret = set()
        target_locations = sg.mismatch_analyzer.get(gene.name, set())
        origin_locations = sg.mismatch_analyzer.get(origin, set())
        for c in gene.coverage:
            for k in c.target_position:
                if k in target_locations:
                    target_ret.add(k)
            if c.position in origin_locations:
                origin_ret.add(c.position)
        return (len(origin_ret), len(target_ret))

    def __init__(self, data, gene, sg):
        '''Data is actually a protobuf data type, which is
        BraceletData, gene is the SubStructure: transcript, sg is a shared_graph class

        '''
        self.origin_name = data.name
        self.target_name = gene.name
        self.origin_chr = data.chr
        self.target_chr = gene.chr
        self.origin_start = data.start_position
        self.target_start = gene.start_position
        self.origin_end = data.end_position
        self.target_end = gene.end_position
        self.origin_coverage_shared_rate = gene.origin_coverage_shared_rate
        self.origin_region_shared_rate = gene.origin_region_shared_rate
        self.target_coverage_shared_rate = gene.target_coverage_shared_rate
        self.target_region_shared_rate = gene.target_region_shared_rate
        self.origin_num_isoforms = len(gene.origin_num_exon)
        self.target_num_isoforms = len(gene.target_num_exon)
        self.origin_num_reads = data.num_read
        self.target_num_reads = gene.num_read
        self.origin_num_exons = sum(gene.origin_num_exon) / len(gene.origin_num_exon)
        self.target_num_exons = sum(gene.target_num_exon) / len(gene.target_num_exon)
        target_ncm, origin_ncm = \
                    self.get_num_consistent_mismatches(sg, data.name, gene)
        self.origin_num_consistent_mismatches = origin_ncm
        self.target_num_consistent_mismatches = target_ncm
        self.origin_matched_score = sg.cuffcompare_result.get(data.name).matched_score
        self.target_matched_score = sg.cuffcompare_result.get(gene.name).matched_score
        self.origin_gene_name = sg.cuffcompare_result.get(data.name).gene_name
        self.target_gene_name = sg.cuffcompare_result.get(gene.name).gene_name
        self.origin_transcript_name = sg.cuffcompare_result.get(data.name).transcript_name
        self.target_transcript_name = sg.cuffcompare_result.get(gene.name).transcript_name
        self.relationshiop_type = sg.gene2gene.get((self.origin_name, self.target_name),
                                                   constants.UU)

    @property
    def X(self):
        return [self.origin_num_reads,
                self.target_num_reads,
                self.origin_num_reads * 1.0 /
                self.target_num_reads,
                self.origin_num_reads >
                self.target_num_reads,
                self.target_coverage_shared_rate,
                self.origin_coverage_shared_rate,
                self.target_coverage_shared_rate -
                self.origin_coverage_shared_rate,
                self.target_coverage_shared_rate * 1.0 /
                self.origin_coverage_shared_rate,
                self.target_coverage_shared_rate >
                self.origin_coverage_shared_rate,
                self.target_region_shared_rate,
                self.origin_region_shared_rate,
                self.target_region_shared_rate -
                self.origin_region_shared_rate,
                self.target_region_shared_rate >
                self.origin_region_shared_rate,
                self.target_num_isoforms == 1,
                self.origin_num_isoforms == 1,
                self.origin_num_exons == 1,
                self.target_num_exons == 1,
                self.target_num_consistent_mismatches,
                self.origin_num_consistent_mismatches,
                self.target_num_consistent_mismatches -
                self.origin_num_consistent_mismatches,
        ]

    def Y(self, correct_gene_set=None):
        if correct_gene_set is not None:
            return 1 if self.origin_gene_name in correct_gene_set and self.target_gene_name not in correct_gene_set  else 0
        return 1 if self.relationshiop_type == constants.GP else 0

    def used_for_training(self, correct_gene_set=None, black_list=None):
        if correct_gene_set is not None:
            ## if two scripts are too close, should be removed because
            ## they might be the same script
            if self.origin_chr == self.target_chr:
                return False
            if self.origin_gene_name not in correct_gene_set and self.target_gene_name not in correct_gene_set:
                return False
            if not black_list:
                return True
            return self.origin_name not in black_list

        return self.relationshiop_type in (constants.GP,
                                           constants.PG,
                                           constants.GG)
