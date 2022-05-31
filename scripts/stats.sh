index=ref/combined_pave_hpv.fa.fai
inref=ref/combined_pave_hpv.fa

## split bam
while read files;
do
    outdir=$files
    inbam=$files.sorted.bam
    samtools index $inbam
    mkdir $outdir || true
    while read line;
    do
        chr=$(echo $line | awk '{print $1}')
        samtools view -b $inbam $chr > $outdir/$chr.bam
        samtools index $outdir/$chr.bam
        samtools depth -a $outdir/$chr.bam > $outdir/$chr.perbase_depth
        grep $chr $inref | cut -d' ' -f 2- > $outdir/$chr.species

#remove empty bam files
        depthfile=$(wc -l $outdir/$chr.perbase_depth | cut -d' ' -f1)
        if [ "$depthfile" -eq "0" ]; then
            rm $outdir/$chr*
        fi
    done <$index
done < ../files.txt

##get stats
for bam in $outdir/*.bam;
do
#echo $chr
#mean depth
samtools depth -a $bam | awk '{{c++;s+=$3}}END{{print s/c}}' > $bam.mean_depth
#breadth of coverage for bases with coverage > mincov
#mincov=0
samtools depth -a $bam | awk '{{c++; if($3>0) total+=1}}END{{print (total/c)*100}}' > $bam.0.breadth
samtools depth -a $bam | awk '{{c++; if($3>=5) total+=1}}END{{print (total/c)*100}}' > $bam.5.breadth
samtools depth -a $bam | awk '{{c++; if($3>=10) total+=1}}END{{print (total/c)*100}}' > $bam.10.breadth
samtools depth -a $bam | awk '{{c++; if($3>=20) total+=1}}END{{print (total/c)*100}}' > $bam.20.breadth
samtools depth -a $bam | awk '{{c++; if($3>=50) total+=1}}END{{print (total/c)*100}}' > $bam.50.breadth
samtools depth -a $bam | awk '{{c++; if($3>=100) total+=1}}END{{print (total/c)*100}}' > $bam.100.breadth


#nanostat for pairwise ID
NanoStat --bam $bam > $bam.nanostats
grep 'Average percent identity' $bam.nanostats | awk '{{print $4}}' > $bam.pwid

done

Rscript --vanilla scripts/collect_stats.R $outdir $outdir.pave.stats.txt
