<?php
function cmp($a,$b) {
	if ($a[15]!=$b[15])
		return $a[15]>$b[15];
	else 
		return $a[9]<$b[9];
}

$full=file($argv[1], FILE_IGNORE_NEW_LINES);

$database=array();
$lastid="";
$output="";
$total_id="";
foreach ($full as $line){
	$data=split("\t",$line);
	if ($data[11]!=0 || $data[10]!=$data[12]) 
		continue;
	if ($data[1]>3)
		continue;

	$database[]=$data;
}

usort($database,"cmp");

$cuffid=0;
$transcriptid=0;
foreach ($database as $data){

	$id=$data[9];
	$chr=$data[13];
	$start=$data[15];
	$end=$data[16];
	$srand=$data[8];
	if (strstr($id,"seq",true)==strstr($lastid,"seq",true)) {
		$transcriptid++;
	}
	else {
		$cuffid++;
		$transcriptid=1;
	}
	$lastid=$id;

	$output.=$chr."\t"."Cufflinks"."\t"."transcript"."\t".
		$start."\t".$end."\t"."1000"."\t".$srand."\t"."."."\t";
	
	$meta=sprintf('gene_id "CUFF.%d"; transcript_id "CUFF.%d.%d"; FPKM "7.6807933123"; frac "1.000000"; conf_lo "6.891722"; conf_hi "8.469865"; cov "14.775373";',$cuffid,$cuffid,$transcriptid);

	$output.=$meta."\n";
	// remove last ","
	$start_exons=split(",",substr($data[20],0,-1));
	$exon_size=split(",",substr($data[18],0,-1));
	for($i=0;$i<count($start_exons);$i++){
		$output.=$chr."\t"."Cufflinks"."\t"."exon"."\t".$start_exons[$i]."\t".($start_exons[$i]+$exon_size[$i])."\t"."1000"."\t".$srand."\t"."."."\t";
		$meta=sprintf('gene_id "CUFF.%d"; transcript_id "CUFF.%d.%d"; exon_number "%d"; FPKM "7.6807933123"; frac "1.000000"; conf_lo "6.891722"; conf_hi "8.469865"; cov "14.775373";',$cuffid,$cuffid,$transcriptid,$i+1);

		$output.=$meta."\n";
	}
}

$fd=fopen($argv[2],"w+");
fprintf($fd,$output);
fclose($fd);

?>