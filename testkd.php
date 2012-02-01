<?php

$f=file("command_list");
mkdir ('temp');
$id=1;
foreach ($f as $l){
	system("bsub -q week -R mem96 -o temp/test".$id." ".$l);
	$id++;
}
?>