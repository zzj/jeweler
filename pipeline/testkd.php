<?php

$f=file("command_list");
mkdir ('temp');
$id=1;
foreach ($f as $l){
	system("bsub -q day -M 8 -n 6 -R \"span[hosts=1]\"  -o temp/test".$id." ".$l);
	//system("bsub -q day  -o temp/test".$id." \"".$l."\"");
	$id++;
}
?>