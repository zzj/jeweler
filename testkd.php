<?php

$f=file("command_list");
mkdir ('temp');
$id=1;
foreach ($f as $l){
	//system("bsub -C 0 -q week -M 16 -n4 -R \"span[hosts=1]\"  -o temp/test".$id." ".$l);	
	system("bsub -q  week -M 4 -C 0 -o temp/test".$id." ".$l);
	$id++;
}
?>