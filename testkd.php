<?php

$f=file("command_list");
mkdir ('temp');
$id=1;
foreach ($f as $l){
	//system("bsub -q bigmem  -n2 -R \"span[hosts=1]\"  -o temp/test".$id." ".$l);	
	system("bsub -q day -C 0 -o temp/test".$id." ".$l);
	$id++;
}
?>