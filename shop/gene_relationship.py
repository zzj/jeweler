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
        self.origin_coverage_shared_rate = gene.origin_coverage_shared_rate
        self.origin_region_shared_rate = gene.origin_region_shared_rate
        self.target_coverage_shared_rate = gene.target_coverage_shared_rate
        self.target_region_shared_rate = gene.target_region_shared_rate
        self.origin_num_isoforms = len(gene.origin_num_exon)
        self.target_num_isoforms = len(gene.target_num_exon)
        self.origin_num_exons = sum(gene.origin_num_exon) / len(gene.origin_num_exon)
        self.target_num_exons = sum(gene.target_num_exon) / len(gene.target_num_exon)
        target_ncm, origin_ncm = \
                    self.get_num_consistent_mismatches(sg, data.name, gene)
        self.origin_num_consistent_mismatches = origin_ncm
        self.target_num_consistent_mismatches = target_ncm
        self.relationshiop_type = sg.gene2gene.get((self.origin_name, self.target_name),
                                                   constants.UU)

    @property
    def X(self):
        return [self.target_coverage_shared_rate - self.origin_coverage_shared_rate,
                self.origin_region_shared_rate,
                self.target_region_shared_rate,
                self.origin_num_isoforms == 1,
                self.target_num_isoforms == 1,
                self.origin_num_exons == 1,
                self.target_num_exons == 1,
                self.origin_num_consistent_mismatches,
                self.target_num_consistent_mismatches,]

    @property
    def X_1(self):
        return [self.origin_coverage_shared_rate -
                self.target_coverage_shared_rate,
                self.origin_region_shared_rate -
                self.target_region_shared_rate,
                self.origin_num_isoforms == 1,
                self.target_num_isoforms == 1,
                self.origin_num_exons == 1,
                self.target_num_exons == 1,
                self.origin_num_consistent_mismatches -
                self.target_num_consistent_mismatches,]

    @property
    def Y(self):
        return 1 if self.relationshiop_type == constants.GP else 0

    @property
    def used_for_training(self):
        return self.relationshiop_type in (constants.GP,
                                           constants.PG,
                                           constants.PP,
                                           constants.GG)
