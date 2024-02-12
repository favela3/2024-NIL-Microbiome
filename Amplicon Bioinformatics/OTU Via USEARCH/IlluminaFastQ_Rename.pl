#!/usr/bin/perl
# This script is to rename all the reads in a Illumina Miseq output as following format:
# SampleName_ReadCount;



####1. Import sample name from a file;
# All the samples names will be saved in an Array called "SampleNames";
# The number of samples will be saved in "$SampleSize";
print "Please input the file name that contains all the sample names: \n";
open (INFILE, <STDIN>);
my @SampleNames=<INFILE>;
close(INFILE);
chomp @SampleNames;
my $SampleSize=@SampleNames;


####2. For each sample names, search for the corresponding sample file;
my $i=0;
my $suffix=".extendedFrags";
my $ReadCount=0;

for ($i=0; $i<$SampleSize; $i++){
    my $infilename="@SampleNames[$i]${suffix}.fastq";
    my $outfilename="@SampleNames[$i].renamed.fastq";
    open (INFILE, "<$infilename");
    open (OUTFILE,">$outfilename");
    $ReadCount=0;
    while(my $line=<INFILE>){
        if ($line =~ /^\@A00419/){
            $ReadCount++;
            $line="@@SampleNames[$i]_Seq$ReadCount\n";        }
        print OUTFILE $line;
    }
    close (INFILE);
    close (OUTFILE);
    print("Finished renaming sequences for sample @SampleNames[$i]\n");
}
