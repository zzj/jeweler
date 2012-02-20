for i in {1..100} 
do 
echo  ~/bin/bin/R CMD BATCH --no-save --no-restore \'--args idx=$i\' active_study.R 
done