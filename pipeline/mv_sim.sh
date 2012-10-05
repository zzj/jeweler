mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0000_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0000_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.0> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.0.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.0.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0000_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0000_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0000_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.0> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.0.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.0.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0000_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0000_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0000_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.0> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.0.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.0.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0000_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0001_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0001_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.1> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.1.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.1.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0001_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0001_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0001_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.1> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.1.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.1.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0001_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0001_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0001_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.1> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.1.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.1.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0001_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0002_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0002_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.2> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.2.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.2.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0002_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0002_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0002_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.2> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.2.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.2.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0002_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0002_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0002_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.2> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.2.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.2.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0002_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0003_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0003_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.3> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.3.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.3.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0003_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0003_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0003_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.3> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.3.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.3.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0003_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0003_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0003_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.3> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.3.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.3.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0003_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0004_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0004_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.4> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.4.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.4.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0004_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0004_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0004_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.4> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.4.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.4.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0004_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0004_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0004_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.4> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.4.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.4.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0004_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0005_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0005_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.5> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.5.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.5.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0005_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0005_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0005_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.5> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.5.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.5.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0005_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0005_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0005_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.5> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.5.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.5.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0005_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0006_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0006_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.6> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.6.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.6.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0006_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0006_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0006_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.6> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.6.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.6.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0006_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0006_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0006_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.6> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.6.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.6.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0006_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0007_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0007_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.7> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.7.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.7.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0007_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0007_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0007_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.7> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.7.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.7.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0007_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0007_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0007_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.7> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.7.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.7.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0007_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0008_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0008_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.8> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.8.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.8.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0008_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0008_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0008_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.8> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.8.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.8.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0008_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0008_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0008_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.8> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.8.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.8.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0008_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0009_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0009_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.9> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.9.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.9.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0009_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0009_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0009_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.9> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.9.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.9.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0009_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0009_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0009_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.9> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.9.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.9.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0009_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0010_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0010_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.10> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.10.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.10.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0010_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0010_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0010_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.10> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.10.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.10.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0010_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0010_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0010_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.10> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.10.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.10.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0010_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0011_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0011_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.11> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.11.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.11.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0011_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0011_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0011_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.11> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.11.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.11.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0011_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0011_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0011_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.11> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.11.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.11.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0011_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0012_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0012_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.12> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.12.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.12.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0012_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0012_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0012_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.12> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.12.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.12.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0012_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0012_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0012_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.12> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.12.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.12.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0012_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0013_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0013_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.13> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.13.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.13.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0013_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0013_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0013_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.13> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.13.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.13.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0013_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0013_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0013_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.13> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.13.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.13.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0013_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0014_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0014_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.14> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.14.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.14.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0014_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0014_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0014_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.14> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.14.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.14.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0014_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0014_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0014_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.14> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.14.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.14.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0014_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0015_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0015_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.15> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.15.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.15.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0015_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0015_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0015_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.15> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.15.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.15.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0015_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0015_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0015_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.15> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.15.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.15.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0015_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0016_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0016_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.16> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.16.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.16.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0016_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0016_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0016_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.16> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.16.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.16.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0016_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0016_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0016_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.16> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.16.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.16.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0016_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0017_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0017_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.17> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.17.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.17.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0017_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0017_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0017_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.17> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.17.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.17.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0017_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0017_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0017_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.17> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.17.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.17.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0017_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0018_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0018_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.18> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.18.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.18.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0018_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0018_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0018_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.18> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.18.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.18.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0018_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0018_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0018_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.18> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.18.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.18.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0018_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//HH_0019_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0019_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.19> /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.19.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/WSB/real/abundance.txt.19.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//HH_0019_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//GG_0019_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0019_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.19> /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.19.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/PWK/real/abundance.txt.19.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//GG_0019_F.abundance
mv /lustre/scr/z/z/zzj/jeweler/data/simulation_bam//FF_0019_F/accepted_hits.bam /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0019_F_merged.bam
sort /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.19> /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.19.sorted
join -2 2 /lustre/scr/z/z/zzj/RNAseqSim/output/CAST/real/abundance.txt.19.sorted gid-tid-mapping.txt > /lustre/scr/z/z/zzj/jeweler/data/simulation_bam_merge//FF_0019_F.abundance
