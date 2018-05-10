#!/usr/bin/perl

print STDERR "
################################################################
# csahdetect.pl --mode=C[onsensus] | S[can4csah] | F[t_charge] #
#               --scan4csah=<scan4csah.pl command>             #
#               --ft_charge=<ft_charge.pl command>             #
#               --s4coutfile=<scan4csah output file to read>   #
#               --ftcoutfile=<ft_charge output file to read>   #
#               --helicalP=<P-value for nonhelicity filtering> #
#               --ftcmhP=<helicalP scaling in segment merging> # 
#               --infasta=<input FASTA file>                   #
#               --maskedfasta=<output masked FASTA file>       #
#               --csahfasta=<output FASTA file with CSAHs>     #
#               > output_CSAH_table                            #
#                                                              #
# Invokes scan4csah and/or ft_charge to detect Charged Single  #
# Alpha-Helices in the input FASTA file.                       #
# For detailed usage & credits, invoke with --help.            #
#                                                              #
# Version 2.0, 21 October 2015.                                #
################################################################

";

$opt_F="";$opt_O="";
$|=1;
use Getopt::Long;
&InitCFscores;

$CSAHDIR="/home/gaspari/.CSAH";

$SCAN4CSAH="$CSAHDIR/scan4csah.pl --paramfile $CSAHDIR/scan4csah_evdtable.txt";
$FTCHARGE ="$CSAHDIR/ft_charge.pl -e $CSAHDIR/ftcharge_evdtable.txt";
$MODE="Consensus";
$FASTA="";
$MINCONSLEN=30;
$HELICALP=0.5;
$FTCMHP=0.5;
$S4COUTFILE="";
$FTCOUTFILE="";
$FASTAROW=60; # row length in output FASTA files

#Getting options
GetOptions('help|?' => \$help,
           'mode=s' => \$MODE,
           'verbose' => \$VERBOSE,
           'scan4csah=s' => \$SCAN4CSAH,
           'ft_charge=s' => \$FTCHARGE,
           's4coutfile=s' => \$S4COUTFILE,
           'ftcoutfile=s' => \$FTCOUTFILE,
           'minconslen=i' => \$MINCONSLEN,
           'helicalP=f' => \$HELICALP,
           'ftcmhP=f' => \$FTCMHP,
           'infasta=s' => \$FASTA,
           'maskedfasta=s' => \$XFASTA,
           'csahfasta=s' => \$CFASTA
          );


&Usage if $help;

exit if ($FASTA =~ /^$/);

if ($MODE !~ /^[CSF]/i){
  die "Invalid mode, please choose from 'Consensus/Scan4csah/Ftcharge'. Invoke with --help for more info\n";
}

$scan4csah="scan4csah.out";
$ftcharge="ftcharge.out";



# Invoking scan4csah & ft_charge
if ($MODE =~ /^[CS]/i){
 if ($S4COUTFILE =~ /^$/){
   $S4COUTFILE=$FASTA;$S4COUTFILE=~s/\.fst/\.s4c/;
   $S4CLOG=$FASTA;$S4CLOG=~s/\.fst/_s4c.log/;
   print STDERR ("Invoking: 'perl $SCAN4CSAH < $FASTA 2> $S4CLOG > $S4COUTFILE'\n") if $VERBOSE; 
   system ("perl $SCAN4CSAH < $FASTA 2> $S4CLOG > $S4COUTFILE "); 
 }
 else{
   print STDERR ("Will read precomputed SCAN4CSAH output from $S4COUTFILE\n") if $VERBOSE; 
 }
 # For consensus mode, only need to invoke ft_charge.pl for sequences with scan4csah hits.
 # Thus, a filtered fasta file can be genaretad and used further based on scan4csah results.
 if ($MODE =~ /^C/i){
   open (S4,"$S4COUTFILE") || die "cannot open input file $S4COUTFILE\n";
   while(<S4>){
      chomp;
      if ($_ !~ /^\#/){
	$seqid=$_;$seqid=~s/ .*$//;
	$seqid=processinputseqid($seqid);
	$KEEPFORFTC{$seqid}=1;
      }
   }
   close(S4);
   $FASTA_S4C_FILTERED=$FASTA."_s4cfiltered";
   open (FF,"$FASTA") || die "cannot open FASTA input file $FASTA\n";
   open (F4,">$FASTA_S4C_FILTERED") || die "cannot create filteted FASTA file $FASTA_S4C_FILTERED\n";
   $writetof4=0;
   while(<FF>){
       chomp;
       if ($_=~/^>/){
	$seqid=$_;$seqid=~s/^>//;
	$seqid=processinputseqid($seqid);
	if ($KEEPFORFTC{$seqid}){$writetof4=1}
	else{$writetof4=0}
       }
       print F4 "$_\n" if $writetof4;
   }
   close(F4);
   close(FF);
   # Hits & everything will only be found in the filtered file in the consensus mode
   $FASTA=$FASTA_S4C_FILTERED;
 }#_consensus mode
}#If S4C to be used

if ($MODE =~ /^[CF]/i){
 if ($FTCOUTFILE =~ /^$/){
   $FTCOUTFILE=$FASTA;$FTCOUTFILE=~s/\.fst/\.ftc/;
   $FTCLOG=$FASTA;$FTCLOG=~s/\.fst/_ftc.log/;
   print STDERR ("Invoking: 'perl $FTCHARGE < $FASTA 2> $FTCLOG > $FTCOUTFILE'\n")  if $VERBOSE;
   system ("perl $FTCHARGE < $FASTA 2> $FTCLOG > $FTCOUTFILE");  
 }
 else{
   print STDERR ("Will read precomputed FT_CHARGE output from $FTCOUTFILE\n") if $VERBOSE; 
 }
}

# Reading scan4csah output file
if ($MODE =~ /^[CS]/){
 open (S4,"$S4COUTFILE") || die "cannot open input file $S4COUTFILE\n";
 while(<S4>){
    chomp;
    if ($_ !~ /^\#/){
	$seqid=$_;$seqid=~s/ .*$//;
	$seqid=processinputseqid($seqid);
        print STDERR "#S4C $seqid\n" if $VERBOSE;
        $all=(split(' ',$_))[17];
	$SCANHITSIN{$seqid}=$all;
    }# scan4csah hit
 }
 close(S4);
}

# Reading ft_charge output file
if ($MODE =~ /^[CF]/){
 open (FT,"$FTCOUTFILE") || die "cannot open $FTCOUTFILE\n";
 #processing ft_charge output: merging regions (generating nonoverlapping predictions)
 %ftn=();%FTS=();%FTE=();
 %TEMPS=();%TEMPE=();$TEMPN=0;%TEMPQ=();
 $prevseqid="none";
 while(<FT>){
    chomp;
    #next if ($_!~/^\>/);
    next if ($_!~/[A-Z]+$/);
    $hit=$_;$hit=~/ .*$/;
    $hit=~s/\>//;$hit=~s/ .*$//;
    $seqid=$hit;
    $seqid=processinputseqid($seqid);
    if (($prevseqid ne "none") && ($seqid ne $prevseqid) && ($TEMPN > 0)){
      &ftc_merge($prevseqid);
    }
    print STDERR "FTC SEQID $seqid\n" if $VERBOSE;
    $range=$hit;#$range=~s/^.*\_//;
    $range=~s/^[a-zA-Z0-9\_\|]+\_([0-9]+\-[0-9]+)$/$1/;
    ($ftstart,$ftend)=split(/\-/,$range);
    $TEMPS{$TEMPN}=$ftstart;$TEMPE{$TEMPN}=$ftend;
    # Storing sequences to be able to introduce another helicity filter
    $sequence=$_;$sequence=~s/^.* ([A-Z]+)$/$1/;
    $TEMPQ{$TEMPN}=$sequence;
    #print STDERR "$_ $TEMPN > $sequence | $TEMPQ{$TEMPN}\n";
    $TEMPN++;
    $prevseqid=$seqid;
 }#_while FT
 # last sequence
 if ($TEMPN > 0){
    &ftc_merge($seqid);
 }
}
close(FT);

# finding consensus
$tcsn=0;
$csthwritten=0;
# CONSENSUS MODE
if ($MODE =~ /^C/){
 foreach $seqid (keys %ftn){
    if ($SCANHITSIN{$seqid} ne ""){
      $csn=0;
      @CS=();@CE=();@SSE=();@FSE=();@ACCESSION=();
      for ($ftx=0; $ftx <= $ftn{$seqid}; $ftx++){
      $ftstart=$FTS{$seqid}[$ftx];$ftend=$FTE{$seqid}[$ftx];
      foreach $scsegment (split(/\,/,$SCANHITSIN{$seqid})){
	    $scsegment=~s/^.*\://;
            ($scstart,$scend)=split(/\-/,$scsegment);
	    print STDERR "#SEQID $seqid > sc: $scstart - $scend | ft: $ftstart - $ftend\n" if $VERBOSE;
            $consstart=-1;
            if ($FTS{$seqid}[$ftx] >= $scstart){
               # case 1: SC segment contains FT segment
               if ($FTE{$seqid}[$ftx] <= $scend){
                 $consstart=$FTS{$seqid}[$ftx];
                 $consend  =$FTE{$seqid}[$ftx];
               }
               # case 2: SC segment contains start but not end of FT segment
               elsif ($FTS{$seqid}[$ftx] <= $scend){
                 $consstart=$FTS{$seqid}[$ftx];
                 $consend  =$scend;
               }
            }
            else { # FT segment starts before SC segment
               # case 3: FT segment contains SC segment
               if ($FTE{$seqid}[$ftx] >= $scend){
                 $consstart=$scstart;
                 $consend  =$scend;
               }
               # case 4: FT segment contains start but not end of SC segment
               elsif ($FTE{$seqid}[$ftx] >= $scstart){
                 $consstart=$scstart;
                 $consend  =$FTE{$seqid}[$ftx];
               }
            }
           # Overlap found
            if ($consstart > -1){
                # New CSAH if over minimum length
                if (($consend-$consstart+1) >= $MINCONSLEN){
		    $CFscoreP=gethelicalEVDP_FASTA($FASTA,$seqid,$consstart,$consend);
		    if ($CFscoreP < $HELICALP){
			$csn++;
			$tcsn++;
			$ACCESSION[$csn]=$seqid;
			$CS[$csn]=$consstart;$CE[$csn]=$consend;
			$SSE[$csn]=$scstart."-".$scend;
			$FSE[$csn]=$FTS{$seqid}[$ftx]."-".$FTE{$seqid}[$ftx];
                        # Marking CSAH positions
                        $CONSENSUS{$seqid}.=":"."$consstart-$consend";
		    }
                }
            }#_if consstart > -1
	}
    }
   &writecsahtable;
  }
 }
}#_consensus mode


# SCAN4CSAH MODE
elsif ($MODE =~ /^S/){
 foreach $seqid (keys %SCANHITSIN){
      $csn=0;
      @CS=();@CE=();@SSE=();@FSE=();@ACCESSION=();
      foreach $scsegment (split(/\,/,$SCANHITSIN{$seqid})){
            $consstart=-1;
	    $scsegment=~s/^.*\://;
            ($scstart,$scend)=split(/\-/,$scsegment);
            $consstart=$scstart;$consend=$scend;
            if ($consstart > -1){
                # New CSAH if over minimum length
                if (($consend-$consstart+1) >= $MINCONSLEN){
		    $CFscoreP=gethelicalEVDP_FASTA($FASTA,$seqid,$consstart,$consend);
		    if ($CFscoreP < $HELICALP){
			$csn++;
			$tcsn++;
			$ACCESSION[$csn]=$seqid;
			$CS[$csn]=$consstart;$CE[$csn]=$consend;
			$SSE[$csn]=$scstart."-".$scend;
			$FSE[$csn]="-";
                        # Marking CSAH positions
                        $CONSENSUS{$seqid}.=":"."$consstart-$consend";
		    }
                }
            }
	}
   &writecsahtable;
  }
}#_scan4csah mode

# FT_CHARGE MODE
elsif ($MODE =~ /^F/){
 foreach $seqid (keys %ftn){
      $csn=0;
      @CS=();@CE=();@SSE=();@FSE=();@ACCESSION=();
      for ($ftx=0; $ftx <= $ftn{$seqid}; $ftx++){
      $consstart=-1;
      $ftstart=$FTS{$seqid}[$ftx];$ftend=$FTE{$seqid}[$ftx];
            $consstart=$ftstart;$consend=$ftend;
            if ($consstart > -1){
                # New CSAH if over minumum length
                if (($consend-$consstart+1) >= $MINCONSLEN){
		    $CFscoreP=gethelicalEVDP_FASTA($FASTA,$seqid,$consstart,$consend);
		    if ($CFscoreP < $HELICALP){
                        #print STDERR "CSAH! $consstart - $consend \n" if $VERBOSE; 
			$csn++;
			$tcsn++;
			$ACCESSION[$csn]=$seqid;
			$CS[$csn]=$consstart;$CE[$csn]=$consend;
			$SSE[$csn]="-";
			$FSE[$csn]=$ftstart."-".$ftend;
                        $CONSENSUS{$seqid}.=":"."$consstart-$consend";
		    }
                }
            }
    }
    &writecsahtable;
 }
}#_ft_charge mode



# Writing masked FASTA file if hits are found
if (($tcsn) && (($XFASTA !~ /^$/) ||  ($CFASTA !~ /^$/))){

 open (FST,"$FASTA") || die "cannot open $fstfile";
 if ($XFASTA !~ /^$/){
     open (XFT,">$XFASTA") || die "cannot create masked output FASTA file $XFASTA\n";
     print STDERR "Writing masked FASTA file...\n" if $verbose;
 }
 if ($CFASTA !~ /^$/){
     open (CFT,">$CFASTA") || die "cannot create output CSAH fasta file $CFASTA\n";
     print STDERR "Writing CSAH FASTA file...\n" if $verbose;
 }

 # reading and processing fasta input file

 $seq="";
 while(<FST>){
    chomp;
    if ($_ =~ />/){
	if ($seq ne ""){
	    &writeseq;
	    $seq="";
	}
        $header=$_;
	$seqid=$_;
	$seqid=~s/^>//;
	$seqid=processinputseqid($seqid);
        #$seqid=(split('\|',$_))[1];
    }
    else{
        $seq.=$_;
    }

 }
 # last but not least
 if ($seq ne ""){
   &writeseq;
 }
 close(XFT);

}#_if fst output required

# Writing summary
if ($tcsn){
 $tseq=keys %CONSENSUS;
 print STDERR "Found $tcsn CSAHs in $tseq proteins.\n";
 print STDERR "Masked FASTA file written to $XFASTA.\n" if ($XFASTA !~ /^$/);
 print STDERR "CSAH FASTA file written to $CFASTA.\n" if ($CFASTA !~ /^$/);
}
else{
 print STDERR "No CSAHs found.\n";
 print STDERR "No masked FASTA file written.\n" if ($XFASTA !~ /^$/);
 print STDERR "No CSAH FASTA file written.\n" if ($CFASTA !~ /^$/);
}


exit;

############################################
# S U B R O U T I N E S
############################################

#writing CSAH table to STDOUT
sub writecsahtable{
   if ($csn > 0){
        print "#Accession\tConsensus\tscan4csah\tft_charge\n" if (!$csthwritten);
        $csthwritten=1;
        for ($x=1; $x <=$csn; $x++){
            print "$ACCESSION[$x]\t$CS[$x] - $CE[$x]\t$SSE[$x]\t$FSE[$x]\n";
        }#_for x
    }
}#_sub writecsahtable
############################################

# Writing masked and/or CSAH sequence
sub writeseq{
    my $maskedseq=$seq;
    my $cstart,$cend,$clen,$cx;
    # write only if CSAH found
    if ($CONSENSUS{$seqid} ne ""){
	$CONSENSUS{$seqid}=~s/^://;
	@ranges=split(/\:/,$CONSENSUS{$seqid});
        foreach $range (@ranges){
            ($cstart,$cend)=split(/\-/,$range);
	    $clen=$cend-$cstart+1;
            substr($maskedseq,$cstart-1,$clen)="x"x$clen;
	    if ($CFASTA !~ /^$/){
		$csahseq=substr($seq,$cstart-1,$clen);
		print CFT ">$seqid $range\n";
		$cx=0;
		while($cx < $clen){
		    print CFT substr($csahseq,$cx,$FASTAROW)."\n";
		    $cx+=$FASTAROW;
		}
	    }#_if CFASTA
        }

	if ($XFASTA !~ /^$/){
	    print XFT "$header\n";
	    $cx=0;
	    while($cx < length($seq)){
		print XFT substr($maskedseq,$cx,$FASTAROW)."\n";
		$cx+=$FASTAROW;
	    }
	}
    }# if consensus found
}#_sub writeseq

#############################xx
# A routine to handle ft-charge output with multiple window sizes
# merges segments
sub ftc_merge{
$my_seqid=shift;
print STDERR "ftc_merge invoked, my_seqid: $my_seqid | seqid: $seqid | prevseqid: $prevseqid\n" if $opt_v;
$ftn{$my_seqid}=-1;
foreach $TEMPN (sort {$TEMPS{$a} <=> $TEMPS{$b}} keys %TEMPS){
    if (helicalEVDP($TEMPQ{$TEMPN}) < ($HELICALP*$FTCMHP)){
	if ($ftn{$my_seqid} < 0){$ftn{$my_seqid}=0;$FTS{$my_seqid}[0]=$TEMPS{$TEMPN};$FTE{$my_seqid}[0]=$TEMPE{$TEMPN}}
	if (($TEMPS{$TEMPN} <= (1+$FTE{$my_seqid}[$ftn{$my_seqid}])) && ($TEMPE{$TEMPN} > $FTE{$my_seqid}[$ftn{$my_seqid}])){
	    $FTE{$my_seqid}[$ftn{$my_seqid}]=$TEMPE{$TEMPN};
	}
	elsif ($TEMPS{$TEMPN} > $FTE{$my_seqid}[$ftn{$my_seqid}]){
	    $ftn{$my_seqid}++;
	    $FTS{$my_seqid}[$ftn{$my_seqid}]=$TEMPS{$TEMPN};$FTE{$my_seqid}[$ftn{$my_seqid}]=$TEMPE{$TEMPN};
	}
    }
    #else{print "Huha\n";}
    #print "$HELICALP | $FTS{$my_seqid}[$ftn{$my_seqid}] - $FTE{$my_seqid}[$ftn{$my_seqid}]\n"; 
}#_foreach

%TEMPS=();%TEMPE=();$TEMPN=0;%TEMPQ=();

}#sub ftc_merge

############################
sub processinputseqid{
    my $inputseqid=shift;
    $inputseqid=~s/^[a-z]+\|//;
    if ($inputseqid=~/^[^\|]+\|/){
      # should work for UniProt & NCBI ids...
      $inputseqid=~s/^([^\|]+)\|.*$/$1/;
    }
    else{
      $inputseqid=substr($inputseqid,0,15);
    }
    return($inputseqid);
}#_sub processinputseqid;

############################
sub InitCFscores{
#Chou-Fasman helix propensity scores
# TiBS June 1977, pp. 128-131, Table I.
%CFscore=("E",1.51,
          "M",1.45,
          "A",1.42,
          "L",1.21,
          "K",1.16,
          "F",1.13,
          "Q",1.11,
          "W",1.08,
          "I",1.08,
          "V",1.06,
          "D",1.01,
          "H",1.00,
          "R",0.98,
          "T",0.83,
          "S",0.77,
          "C",0.70,
          "Y",0.69,
          "N",0.67,
          "P",0.57,
          "G",0.57
    )
}

############################
# function to return helicap P for a given sequence provided as a string
sub helicalEVDP{
    my $sq=shift;
    my $cfs=0;
    for (my $cfx=0; $cfx < length($sq); $cfx++){
	$cfs+=$CFscore{substr($sq,$cfx,1)};
    }
    $cfs/=length($sq);

    my $loc=1.060; my $scale=0.072; my$shape=-0.299;
    my $z = ($cfs - $loc) / $scale;
    my $P_VALUE = 1 - exp(-(1+$shape*$z)**(-1/$shape));
    print STDERR "SQ: $sq: $P_VALUE\n" if $VERBOSE;
    if ($P_VALUE =~ /nan/i){$P_VALUE=0}
    return $P_VALUE;
}

# function to return helicap P for a given sequence from a FASTA seq
sub gethelicalEVDP_FASTA{
    my $FASTAfile=shift;
    my $sid=shift;
    my $sqstart=shift;
    my $sqend=shift;

    open (FH,"$FASTAfile") || die "Cannot open FASTA file $FASTAfile\n";
    my $readseq=0;
    my $fullprotseq="";
    while(<FH>){
	chomp;
	if ($_=~/^>/){
	    $id=$_;$id=~s/^>//;$id=processinputseqid($id);
	    if ($id eq $sid){$readseq=1}
	    else{$readseq=0}
	}
	elsif ($readseq){
	    $fullprotseq.=$_;
	}
    }
    close(FH);

    my $csahsq=substr($fullprotseq,$sqstart-1,$sqend-$sqstart+1);
    return helicalEVDP($csahsq);


}

############################
sub Usage{

print "

CSAHdetect package, version 2.0

================================================================================
CREDITS

csahdetect.pl is a wrapper to invoke scan4csah and/or ft_charge and
return a uniformly formatted output.

When using this program, scan4csah or ft_charge, 
kindly cite one of the following references:

    * Dániel Dudola, Gábor Tóth, László Nyitray, Zoltán Gáspári:
      Consensus prediction of charged single α-helices with CSAHserver
      in preparation 

    * Zoltán Gáspári, Dániel Süveges, András Perczel, László Nyitray, Gábor 
      Tóth:
      Charged single alpha-helices in proteomes revealed by a consensus 
      prediction approach.
      Biochem. Biophys. Acta - Proteins and Proteomics (2012) 1824:637-646.

    * Dániel Süveges, Zoltán Gáspári, Gábor Tóth, László Nyitray:
      Charged single α-helix: a versatile protein structural motif
      Proteins (2009) 74:905-916.

================================================================================
INSTALLATION

The only specific requirement is that FFT.pm should be avalable for ft_charge
to be able to run. FFT.pm is obtainable from CPAN (http://www.cpan.org).
However, the install script can invke 'sudo cpan Math::FFT' for a 
straightforward install, if you chose so.

Installation steps:

- Unpack CSAHdetect.zip in a suitable diretory
- Type 'perl INSTALL.PL' from that direcotry.
- The install script will
  - check the availability of Math::FFT and offer to install it using cpan 
    (you can skip this but you will have to install Math::FFT manually 
    later to use the FT_CHARGE method)
  - ask you where to put the scan4csah.pl and ft_charge.pl executables 
    and the EVD parameter files
  - ask you where to put the csahdetect.pl executable
  - will put the executables into their respective locations and apply 
    the given settings.

If everything went OK, you can now invoke csahdetect.pl --help for usage.

The notes below are for users who might want to use the scan4csah.pl 
and ft_charge.pl programs separately for some purpose.
It is advised that both scan4csah.pl and ft_charge.pl are copied to a 
location where they can be easily executable (e.g. /usr/local/bin).
It is also advised that their respective EVD parameter files are moved to
a standard location and their default location is changed by editing
the default values of the corresponding variables in the scripts:

[scan4csah.pl, line 206]: \#my  = \"/home/szpari/csahserver/scan4csah_evdtable.txt\";
                (remove the \"\#\" to make the change active, by defult, scan4csah.pl
                           is able to run without an external parameter file)


[ft_charge.pl, line 46]: $opt_e=\"/home/szpari/csahserver/ftcharge_evdtable.txt\";

================================================================================
USAGE

Running csahdetect.pl

Invoking csahdetect.pl without any options will give you an overview of 
options and input/output files.

By default, you can just type
csahdetect.pl --infasta=<input_fasta_file>

The algorithms used are selected with the '--mode' option:
 C[onsensus] (default): invoke both programs
 S[can4csah] : invoke scan4csah
 F[t_charge] : invoke ft_charge

(only the first character is meaningful).

As the FT_CHARGE method is much slower than SCAN4CSAH, in the consensus 
mode SCAN4CSAH is invoked first and then a fasta file is generated
containing only those sequences where CSAH segments were predicted by SCAN4CSAH.
FT_CHARGE is then invoked on this reduced sequence set to save time.
(If you are not happy with this, you can get a full FT_CHARGE-based prediction  
by invoking 'csahdetect.pl --mode=F' or running ft_charge.pl separately.)

By default, the programs are invoked without any special options, meaning they
use their defaults (they use their default evd files). 
If you use non-standard locations or want to change any
of their respective options to non-default values, you might either specify 
these on the command line or change the default values of the relevant
variables in this script.

The user might also chose to use precomputed SCAN4CSAH and/or FT_CHARGE 
outputs, in this case the output files specified with --s4coutfile and
--ftcoutfile will be read in and SCAN4CSAH and/or FT_CHARGE will not be 
invoked. This is useful e.g. to extract consensus from runs of the two
algorithms on different sequence sets. 

The minimum length of consensus CSAHs is set to 30 (the default minumum length
in scan4csah as it is shorter than the default window size of 32-64 in 
ft_charge), you can change it using the --minconslen option.

Examples:

 - using command-line options:

 > csahdetect.pl --scan4csah='/home/pompom/bin/scan4csah.pl --minlen=40' 
   --ft_charge='/programs/ft_charge.pl -p 0.1' --infasta=...
 
 - rewriting defaults in csahdetect.pl:
 [...]
 \$SCAN4CSAH=\"/home/pompom/scan4csah.pl\";
 \$FTCHARGE=\"/programs/ft_charge.pl\";
 [...]

If option --maskedfasta is given, a masked FASTA file (with residues in CSAH 
regions masked as 'x') will also be written.

For more information on the usage of scan4csah.pl and ft_charge.pl kindly refer
to their own description (invoke 'scan4csah.pl --help' or 'ft_charge.pl' -h)

================================================================================
DISCLAIMER AMD RIGHTS

This program, as well as scan4csah.pl and ft_charge.pl are provided on an 
\"as is\" basis in the hope that they will prove meaningful.
The authors are not responsible for any damage or malfunction caused
by the usage of these programs.
All programs in this csahdetect package are distributed freely and the user is 
free to make any modifications provided (s)he 1) acknowledges the use of the 
programs by citing the references above and 2) distributes the modified
versions as free software.
 
================================================================================
";

exit(1);

}#_sub Usage

