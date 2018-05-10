#!/usr/bin/perl

print STDERR "
#####################################################
# insertCSAH2dat.pl -l<csah_list> [-s] < sp_datfile # 
#                    > sp_datfile_with_csah         #
#                                                   #
# -s causes to write CSAH-containing dat entries    #
#    only                                           #
#                                                   #
# Version 1.2                           07.06.2011. #
#####################################################

";

use Getopt::Std;

getopts('l:svh');

# reading CSAH table
open (LF,"$opt_l") || die "Cannot open CSAH list file $opt_l\n";

while(<LF>){

    chomp;
    unless ($_=~/^\#/){
	#($ac,$name,$con,$s4c,$ftc)=split(/\t/,$_);
         ($ac,$con,$s4c,$ftc)=split(/\t/,$_);
	 $CSAH{$ac}.=":$con";
    }
}
close(LF);

#Reading Swiss-Prot dat file
$lastft="";
$writeentry=1;
while(<>){
    chomp;
    if ($_=~/^ID/){
	    $idline=$_;
	    $writeentry=0;
	    $CSAHtowrite="";
    }
    # Have to correctly deal with multi-line AC enries: writing id only once and storing CSAH
    if ($_=~/^AC/){
	$name=$_;$name=~s/^AC +([A-Z0-9]+)\;.*$/$1/;
	if ($CSAH{$name} =~ /[0-9]/){$writeentry=1;$CSAHtowrite=$CSAH{$name};print "$idline\n";}
	#elsif ($opt_s == 1){$writeentry=0}
	elsif ($opt_s == 0){$writeentry=1}
	print STDERR "AC: $name CSAH: $CSAHtowrite\n" if $opt_v;
	#print "$idline\n" if $writeentry;
	$lastft="";
    }

    elsif($_=~/^FT /){
	if (($_!~/ COILED /) && ($lastft=~/ COILED /)){
            print STDERR "Writing CSAH record for $name after COILED lines\n" if $opt_v;
	    &writecsahlines;
	    $csahwritten{$name}=1;
	}
	$lastft=$_;
    }
    #if (($_=~/^SQ/) && ($CSAH{$name} =~/[0-9]/) && (!$csahwritten{$name})){
    if ($_=~/^SQ/){
	print STDERR "Writing CSAH record for $name before sequence\n" if $opt_v;
	&writecsahlines unless $csahwritten{$name};
	$csahwritten{$name}=1;
    }


    print "$_\n" if $writeentry;

}

sub writecsahlines{
    print STDERR "Writing CSAH lines for entry $name\n" if $opt_v; 
    $CSAHtowrite=~s/^\://;
    #$CSAH{$name}=~s/^\://;
    foreach $csah (split(/\:/,$CSAHtowrite)){
        ($cs,$ce)=split(/ \- /,$csah);
	#print  "FT   COILED      195    255       Potential.\n";
	printf ("FT   CSAH     %6d %6d       Potential.\n",$cs,$ce);
    }
}#_ sub writecsahlines

