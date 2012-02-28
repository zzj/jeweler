<?php

$f=file("command_list");
mkdir ('temp');
$id=1;
foreach ($f as $l){
	//system("bsub -q day  -n4 -R \"span[hosts=1]\"  -o temp/test".$id." ".$l);	
	system("bsub -q bigmem -o temp/test".$id." ".$l);
	$id++;
}
?>