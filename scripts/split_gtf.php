<?php


// argument format
// php split_get.php alias gtf_file_from_cufflinks left_side_ref_sequence right_side_ref_sequence bam_file

// the script will create a folder named by the given alias under folder result/cuffsequence


@mkdir('result');
@mkdir('result/cuffsequence')
$cufffolder="result/cuffsequence/".$argv[1]."/";
$output_info=$cufffolder.$argv[1].".info";


$transcriptfile=$argv[2]
$foutput=fopen($output_info,"w+");
$father=$argv[3];
$mother=$argv[4];
$bamfile=$argv[5];
mkdir($cufffolder);

$fd=file($transcriptfile, FILE_IGNORE_NEW_LINES);
$lastid="";
$tranfile=NULL;
$summary=array();
$nodiff=array();
$diff=array();

for ($i=1;$i<=31594;$i++){
	$summary["CUFF.".$i]=0;
}
$id=0;
$chr="";
$min=-1;
$max=-1;
foreach ($fd as $line){
	$data=split("\t", $line);

	if ($data[2]=='transcript'){
		$newdata=split(";",$data[8]);
		$geneid=split("\"",$newdata[0]);
		$geneid=$geneid[1];
		
		@$summary[$geneid]++;
		if ($lastid!=$geneid){
			
			if ($tranfile!=NULL) {
				fclose($tranfile);
				printf($id."\n");
				$id++;
				$fa_seq=$genefolder.$lastid.'.'.$father.".fa";
				$ma_seq=$genefolder.$lastid.'.'.$mother.".fa";
				$fa_map=$genefolder.$lastid.'.'.$father.".map";
				$ma_map=$genefolder.$lastid.'.'.$mother.".map";
				$read_seq=$genefolder.$lastid.".seq.fasta";
				system("gffread -w ".$fa_seq." -g $father ".$genefolder.$lastid);
				system("gffread -w ".$ma_seq." -g $mother ".$genefolder.$lastid);
				system('samtools view '.$bamfile.' '.$chr.':'.$min.'-'.$max.' |awk \'{OFS="\\t"; print ">"$1"\\n"$10}\' - > '.$read_seq);
				$result=array();
				exec('diff '.$fa_seq." ".$ma_seq,$result);
				system("./blat ".$fa_seq." ".$read_seq." -t=dna -q=dna ".$fa_map);
				system("./blat ".$ma_seq." ".$read_seq." -t=dna -q=dna ".$ma_map);
				if (count($result)<2){
					$nodiff[]=$lastid;
				}
				else {
					$diff[]=$lastid;
				}
				fprintf($foutput,$lastid."\t".$genefolder."\t".$genefolder.$lastid."\t".$fa_seq."\t".$ma_seq."\t".$read_seq."\t".$fa_map."\t".$ma_map."\n");
			}
			$genefolder=$cufffolder.$geneid."/";
			if (!file_exists($genefolder))
				mkdir($genefolder);
			$tranfile=fopen($genefolder.$geneid,"w+");
			$lastid=$geneid;
			$chr=$data[0];
			$min=$data[3];
			$max=$data[4];
		}
		else{
			$chr=$data[0];
			if($min<0){
				$min=$data[3];
			}
			else {
				if ($data[3]<$min){
					$min=$data[3];
				}
			}

			if($max<0){
				$max=$data[4];
			}
			else {
				if ($data[4]<$max){
					$max=$data[4];
				}
			}

		}
	}
	fprintf($tranfile,$line."\n");
}

//asort($summary);
//print_r($summary);

?>