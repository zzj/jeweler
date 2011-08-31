<?php

$fdr=fopen($argv[1],"r");
$fdw=fopen($argv[2],"w+");

while(($buffer=fgets($fdr,100000))!=false){

	$buffer=trim($buffer);
	if ($buffer==">chrM") fclose($fdw);
	if (strlen($buffer) % 70==0 && strlen($buffer)>70){
		for ($i=0;$i<strlen($buffer);$i=$i+70){
			fprintf($fdw,substr($buffer,$i,70)."\n");
		}
	}
	else {
		if (strlen($buffer)%70!=0){
			print(strlen($buffer)."\n");
			print($buffer."\n");
		}
		fprintf($fdw,$buffer."\n");
	}

}

?>

