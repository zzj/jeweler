for i in {1..200} 
do 
echo  ~/bin/bin/R CMD BATCH --no-save --no-restore \'--args idx=$i\' plot_cuffcompare.R done
done
