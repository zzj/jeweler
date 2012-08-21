<?php

$f=file("command_list");
mkdir ('temp');
$id=1;
foreach ($f as $l){
	system("bsub -C 0 -q hour -M 4 -n 1 -R \"span[hosts=1]\"  -o temp/test".$id." ".$l);	
	//system("bsub -q day -M 96 -C 0 -o temp/test".$id." ".$l);
	$id++;
}
?>