<?php

$f=file("command_list");
mkdir ('temp');
$id=1;
foreach ($f as $l){
	system("bsub -q day -M 24 -n 1 -R \"span[hosts=1]\"  -o temp/test".$id." ".$l);
	//system("bsub -q day  -o temp/test".$id." \"".$l."\"");
	$id++;
}
?>
