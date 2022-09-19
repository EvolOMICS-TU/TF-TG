#!/bin/bash
#SBATCH --job-name=%j_apt
#SBATCH --ntasks=40
#SBATCH --nodelist=compute1

##############################################################################################################
################## Defining the path variables ########################################################
##############################################################################################################

#replace "meme-5.1.1" with correct version like "meme-5.4.1"
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.1.1:$PATH

#enter the name of the project
proj="trial_ashutosh"

#Enter the path to the the list of DEGs
list=~/Project/Arabidopsis/Data_list/genelist

#Enter the path to the list of Transcription Factors (TFs)
pwms=~/Project/Arabidopsis/PWMs/pwm

#enter the path to the file containing TFs from cis-bp
tf=~/Project/Arabidopsis/PWMs/TF_Information.txt

#Enter the path to the retrieved upstream sequences
seq=~/Project/Arabidopsis/Data_list/sequences

#Define the path to the folder containing the PWMs
mkdir ~/$proj/
mkdir ~/$proj/Sequences/
mkdir ~/$proj/Sequences/background/
mkdir ~/$proj/Sequences/meme_motifs/
mkdir ~/$proj/Results/
mkdir ~/$proj/extra/
mkdir ~/$proj/TFs/

#Retrieve the upstream sequences in 
sequence=~/$proj/Sequences
pseudog=~/$proj/extra
result=~/$proj/Results
mscan=~/$proj/TFs
echo "Paths stored"

cat $seq > $sequence/ref-seq

################################################################################################################
#################### Editing the genelist for any windows-based error/duplications #############################
################################################################################################################
:<<'END'
tr -d '\15\32' < $list > genelist_w1

tr '[:lower:]' '[:upper:]' < genelist_w1 > genelist_w2

sort -u genelist_w2 > $list

rm genelist_w1

rm genelist_w2 

#cp $list/genelist $mscan/genelist
END
###############################################################################################################
################### Retrieving the start and end coordinates of genes #########################################
###############################################################################################################

#The subsequent step requires RSAT installation

#retrieve-seq  -org Arabidopsis_thaliana.TAIR10.42 -feattype gene -type upstream -format fasta -label id,name -from -1000 -to -1 -noorf -i $list -o $sequence/ref-seq

grep -v "WARNING" $sequence/ref-seq > $sequence/fasta
cat $sequence/fasta > $sequence/ref-seq
rm $sequence/fasta 

echo "Upstream regions of the target genes extracted"

cd $pseudog
grep -i ">" $sequence/ref-seq | sed "s/[|>:;]/ /g" | awk -F "\t" '{print $1,$16,$17,$18,$1}' | awk -F "\t" '{if($2<$3) print}' > $pseudog/Gene_Coord
sed -i 's/D/+/g' $pseudog/Gene_Coord
sed -i 's/R/-/g' $pseudog/Gene_Coord
echo "Gene Coordinates of the upstream regions extracted"

echo "Files with upstream coordinates transferred to Pseudogenomes"

echo "Extracting the list of TFs from the DEGs"

rm $mscan/tflist
grep -if $list $tf | awk '{if($9 != "N") print $4,$6,$7,$9}' >> $mscan/tflist
#| grep -i "HORVU"

:<<'END'
while read -r line
do
#echo $line
	grep $line $tf | awk '{if($9 != "N") print $4,$6,$7,$9}' >> $mscan/tflist
done < $list
END
w=($(wc $mscan/tflist))

if [[ $w -gt 0 ]]
then

	echo "List of transcription factors extracted"
	
fi

###########################################################################################################################
####################### Compilation of all motifs in a file in TRANSFAC format ############################################
###########################################################################################################################


rm $sequence/collected

#cd $mscan

while read -r line
do

	a=($line)
	motif="${a[0]}.txt"
	w=($(wc $pwms/$motif))
	if [[ $w -gt 1 ]]
	then

		echo "AC ${a[1]}" >> $sequence/collected
		echo "XX" >> $sequence/collected
		echo "ID ${a[2]}" >> $sequence/collected
		echo "XX" >> $sequence/collected
		cat $pwms/$motif >> $sequence/collected
		echo "XX" >> $sequence/collected
	fi
	
done < $mscan/tflist

echo "//" >> $sequence/collected

sed -i "s/Pos/PO/g" $sequence/collected

echo "Concatenation of motifs done"

cd $sequence
ls | grep "sequences" | grep -v "list" > $sequence/list_of_sequences

##################################################################################################################################
################### Scanning of sequences by single compiled motif file using FIMO ###############################################
##################################################################################################################################

fasta-get-markov -m 1 -dna $sequence/ref-seq $sequence/background/ref-seq
transfac2meme -bg $sequence/background/ref-seq -use_acc $sequence/collected > $sequence/meme_motifs/ref-seq
fimo -o $result/ref-seq -bfile $sequence/background/ref-seq $sequence/meme_motifs/ref-seq $sequence/ref-seq

#rm $mscan/collected

###################################################################################################################################
###################################### Removing duplicated rows from result #######################################################
###################################################################################################################################




