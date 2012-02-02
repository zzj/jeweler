<?php
sleep(3600);
$f=file("command_list");
mkdir ('temp');
$id=1;
foreach ($f as $l){
	system("bsub -q bigmem -o temp/test".$id." ".$l);
	$id++;
}
?>