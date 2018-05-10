#!/usr/bin/perl -w
#
# scan4csah.pl 2.3.2                       (C) Gabor Toth, 2005-2010
#
# Scans protein sequence(s) for charged amino acid patterns and
# putative salt bridges indicative of stable alpha-helix (CSAH domain).
#
#[ HISTORY:
#
#  v2.3.2 (27-Jul-2010)
#    Bug: when the first CSAH started with residue 1, it was not detected.
#    Fixed.
# 
#  v2.3.1 (13-Jan-2010)
#    Bug: the start and end positions of the best CSAH are still
#    not correct. Fixed.
#
#  v2.3 (12-Jan-2010)
#    Bug: CSAH postitions are incorrectly reported
#    (start and end positions are less than the true values by one
#    due to 0-based position numbering).
#    Fixed. 
#
#  v2.2 (11-Jan-2010)
#    Correcting HTML output.
#
#  v2.1 (22-Dec-2009)
#    Adjusting script output to web server (quiet run etc.).
#
#  v2.0 (15-Dec-2009)
#    Modified from 'scan4sah2.pl' v1.0.1.  New name: 'scan4csah.pl'!
#    The program can use pre-computed score distribution parameters
#    and fitted EVD parameters to calculate Z-score and P-value for
#    each score.
#
#
#  v1.0.1 (21-Sep-2007)
#    The description of real column 17 was missing from the POD
#    documentation. Fixed.
#
#  v1.0 (20-Sep-2007)
#    The behaviour using option --swissprot is modified.
#    Sequence ID in new UniProtKB/SwissProt fasta-format files
#    follows the format <acc_num>|<seq_id>.
#
#  v0.9 (1-Jan-2007)
#    Bug: default max. rel. score ($maxrelscore) was zero, and this
#    erroneously prevented printing of the real relative score if it
#    was negative. Fixed.
#
#  v0.8 (13-Oct-2006)
#    Made conformant with scan4sah3.pl v0.3:
#    Score all triad combinations:
#    (i,i+4,i+8), (i,i+4,i+7), (i,i+3,i+7), (i,i+3,i+6);
#    in both orientations: (a-b-a), (b-a-b)
#
#    Default parameters have been set according to the grid search result
#    of 4-Jan-2006 ('gridsearch7.csh').
#
#  v0.7 (7-Dec-2005)
#    Hits with regions having relative score greater than threshold
#    were reported but max. rel. score could be much less then threshold
#    because the longest region was selected as best.
#    This discrepancy is corrected now: 
#    the best region is selected by max. rel. score instead of max length.
#    Score and max. rel. score are written for this region.
#    Note that specifying a minimum length limit is important not to be
#    mislead by strong but very short regions and to obtain meaningful results.
#
#  v0.6 (27-Nov-2005)
#    Negative scrores are introduced for destabilizing effects:
#    (a) opposite (+-) charges at (i,i+1) and (i,i+2) positions;
#    (b) identical (++ or --) charges at (i,i+3) and (i,i+4) positions.
#
#  Script renamed to 'scan4sah2.pl' (24-Nov-2005)
#
#  v0.52 (09-Sep-2005)
#    Integer scores are forced in command-line options.
#
#  v0.51 (09-Sep-2005)
#    Default score threshold has been increased to 4.
#    Default values are corrected in POD.
#
#  v0.5 (06-Sep-2005)
#    Scoring system modified: each score must be an integer.
#    Previous scores have been multiplied by a factor of 4.
#    Since the basis of the relative score is the length, the relative 
#    score threshold also has to be increased accordingly
#    (i.e. multiplied by a factor of 4).
#
#  v0.4 (28-Aug-2005)
#    Scores for each series of consecutive dyads are calculated separately.
#    The score of the best (longest?) series and the total score are reported
#    in the output.
#
#  v0.32 (25-Aug-2005)
#    Extension of option --html. Output table in HTML format with links to
#    individual HTML files of each annotated sequences.
#
#  v0.31 (24-Aug-2005)
#    New option: HTML output of annotated sequences.
#
#  v0.3 (23-Aug-2005)
#    Absolute (and relative) score is calculated only from those series
#    of consecutive dyad-forming residues that form a series longer than
#    the minimum length ($minlen) specified for a putative CSAH.
#
#    Removed: histogram plot of scores. Moved to a separate script: 'plot4sah.pl'.
#
#  v0.25 (19-Aug-2005)
#    New relative score is calculated only from those series of consecutive
#    dyad-forming residues that form a series longer than the minimum length
#    ($minlen) specified for a putative CSAH.
#    Rel. score: number of dyad-forming residues divided by the total length
#    of the series of dyads formed by such residues.
#
#    Old relative score was: total score divided by the number of all charged
#    residues.
#
#  v0.2 (19-Aug-2005)
#    Added: histogram plot of scores.
#
#  v0.11 (15-Aug-2005)
#    Documentation: output format described.
#
#  v0.1 (15-Aug-2005)
#    Added: length constraint.
#    Tabular output.
#
#  v0.02 (14-Aug-2005)
#
#  v0.01 (14-Aug-2005):
#    Modified from scan4charge.pl v0.41.
#]

$VERSION = '2.3.1';

$SERVER_HOST = 'himalia';  # valid values: 'emboss' or 'himalia' or '' (empty string)

############################
# Default values:          #
############################
# Scores:                  #
our $dyad4a = 1;           #  i,i+4 dyad (-+) (integer score!)
our $dyad3a = 9;           #  i,i+3 dyad (-+) (integer score!)
our $dyad4b = 0;           #  i,i+4 dyad (+-) (integer score!)
our $dyad3b = 3;           #  i,i+3 dyad (+-) (integer score!)
our $triada = 7;           #  i,i+4,i+8 or i,i+3|4,i+7 or i,i+3,i+6 (-+-) triad (integer score!)
our $triadb = 0;           #  i,i+4,i+8 or i,i+3|4,i+7 or i,i+3,i+6 (+-+) triad (integer score!)
our $dyad4p = 3;           #  penalty for i,i+4 dyad (++ or --) (pos. integer!)
our $dyad3p = 0;           #  penalty for i,i+3 dyad (++ or --) (pos. integer!)
our $dyad2p = 2;           #  penalty for i,i+2 dyad (+- or -+) (pos. integer!)
our $dyad1p = 0;           #  penalty for i,i+1 dyad (+- or -+) (pos. integer!)
############################
our $gap = 5;              #  max. gap between two consecutive dyad-forming
                           #   residues
our $minlen = 30;          #  minimum length of a candidate region with
                           #   consecutive dyads (30?)
#our $threshold = 2.5;     #  score threshold (relative)
our $threshold = 0;        #  score threshold (relative)
                           #
our %MATRIX = (            #  list of charged amino acids with their charges
	       'D' => -1,  #
	       'E' => -1,  #
	       'K' => 1,   #
	       'R' => 1,   #
#	       'H' => 0.5  #  H is not needed now, and
                           #   in any case, do not use floating point score!
	       );          #
our $SWISS = 0;            #  0: do not write accession number;
                           #  1: swissprot input file,
                           #   write acc. num. into output table
our $PVAL = 0;             #  0: do not calculate P-value; 1: calculate P-value
                           #   (set --paramfile to enable!)
our $VERBOSE = 0;          #  0: do not write negative results into output
our $QUIET = 0;            #  0: normal run; 1: quite run (less messages to stderr)
                           #
our $HTML = 0;             #  do not write HTML output
if ($SERVER_HOST eq 'emboss') {        #
    our $bdir = "/srv/www/htdocs";     # base directory for output
    our $htmldir = "/csah/temp/html";  # HTML directory within DocumentRoot
}                                      #
elsif ($SERVER_HOST eq 'himalia') {    #
    our $bdir = "/srv/www/himalia/htdocs"; # base directory for output
    our $htmldir = "/temp/html";       # HTML directory within DocumentRoot
}                                      #
else {                                 #
    our $bdir = ".";                   # base directory for output
    our $htmldir = "/html";            # HTML directory within DocumentRoot
}                                      #
########################################


$\ = "\n";                      # output record separator for print
$, = " ";                       # output field separator for print

use Getopt::Long;
use Pod::Usage;

&Init;
&Read_sequences;

# End of main program.


# Subroutines:
#--------------------
sub Init {
    our %PARAMS = ();  # Distributon and EVD parameters
    my $man = 0;
    my $help = 0;
    my $matrix = 0;
    my @AminoAcids = ('A','B','C','D','E','F','G','H','I','K','L','M',
                      'N','P','Q','R','S','T','U','V','W','Y','Z','X');
    my $aa = '';
    #my $paramfile = "/home/szpari/csahserver/scan4csah_evdtable.txt";
#    our @Scores = ();

    GetOptions('help|?' => \$help,
	       'man' => \$man,
	       'verbose' => \$VERBOSE,
	       'quiet' => \$QUIET,
	       'matrix=s' => \$matrix,
	       'gap=i' => \$gap,
	       'minlen=i' => \$minlen,
	       'threshold=f' => \$threshold,
	       'dyad4acidic=i' => \$dyad4a,
	       'dyad4basic=i' => \$dyad4b,
	       'dyad3acidic=i' => \$dyad3a,
	       'dyad3basic=i' => \$dyad3b,
	       'triadacidic=i' => \$triada,
	       'triadbasic=i' => \$triadb,
	       'dyad4pen=i' => \$dyad4p,
	       'dyad3pen=i' => \$dyad3p,
	       'dyad2pen=i' => \$dyad2p,
	       'dyad1pen=i' => \$dyad1p,
	       'html=s' => \$html,
	       'swissprot' => \$SWISS,
	       'paramfile=s' => \$paramfile,
	       );
#	       ) or pod2usage(2);

    &Message unless $QUIET;

    pod2usage(1) if $help;
#    pod2usage(-exitstatus => 0, -verbose => 2) if $man;

    if ($gap < 4) {
	die "Gap (option --gap) must be greater than 3!\n";
    }

    if ($html) {
        $HTML = 1;
        $html .= '.html';
        open(HTML,">$html") or die "Cannot create HTML file $html: $!\n";
#	our $oldfh = select HTML;
	$| = 1;

	print HTML "<HTML>\n<HEAD><TITLE>scan4csah.pl output</TITLE></HEAD>\n<BODY>\n";
	print HTML "<FONT FACE=\"Courier New, Courier\">";
	mkdir "${bdir}${htmldir}";
    }

    if ($matrix) {
	my $score = 0;
	undef(%MATRIX);

	open (MATRIX,"$matrix") || die "Cannot open file $matrix (option --matrix): $!\n";

	while (<MATRIX>) {
	    chomp;
	    next if (/^$/ || /^\#/);
	    ($aa,$score,undef) = split;
	    $MATRIX{$aa} = $score;
	}
	close(MATRIX);
    }

    foreach $aa (@AminoAcids) {
	$MATRIX{$aa} = 0 unless defined($MATRIX{$aa});
    }

    unless ($QUIET) {
	# to decrease size of error log file:
	print STDERR "Min. length: $minlen";
	print STDERR "Max. gap: $gap";
	print STDERR "D4a: $dyad4a  D4b: $dyad4b  D3a: $dyad3a  D3b: $dyad3b  Ta: $triada  Tb: $triadb";
	print STDERR "P4: $dyad4p  P3: $dyad3p  P2: $dyad2p  P1: $dyad1p";
	print STDERR "Score threshold: $threshold";
    }

    if ($paramfile) {
	print STDERR "Parameter table is read from file \"$paramfile\"." unless $QUIET;
	&Read_parameter_table;
	$PVAL = 1;                # calculate P-value
    }

    # header of output table
    if ($HTML) {
#	my @C1 = ('Seq. ID','Acc. num.','+-','aa','positive','negative','dyad4a','dyad4b','dyad3a','dyad3b','triad4a','triad4b','max. score','max. rel. score','score','rel. score','number of CSAHs','best CSAH','all CSAHs');
#	my @C2 = ('Seq. ID','+-','aa','positive','negative','dyad4a','dyad4b','dyad3a','dyad3b','triad4a','triad4b','max. score','max. rel. score','score','rel. score','number of CSAHs','best CSAH','all CSAHs');
	my @C1 = ('Sequence<br>ID','Accession<br>number','amino<br>acid<br>residues','basic<br>residues<br>(positive)','acidic<br>residues<br>(negative)','dyad4<br>-+','dyad4<br>+-','dyad3<br>-+','dyad3<br>+-','triad<br>-+-','triad<br>+-+','score<br>of<br>best<br>CSAH','relative<br>score<br>of best<br>CSAH','total<br>score','relative<br>score','number<br>of<br>CSAHs','best CSAH<br>[length: start-end]','all CSAHs<br>[length: start-end]');
	my @C2 = ('Sequence<br>ID','amino<br>acid<br>residues','basic<br>residues<br>(positive)','acidic<br>residues<br>(negative)','dyad4<br>-+','dyad4<br>+-','dyad3<br>-+','dyad3<br>+-','triad<br>-+-','triad<br>+-+','score<br>of<br>best<br>CSAH','relative<br>score<br>of best<br>CSAH','total<br>score','relative<br>score','number<br>of<br>CSAHs','best CSAH<br>[length: start-end]','all CSAHs<br>[length: start-end]');
#	print HTML "<PRE>";
	if ($SWISS) {
#	    printf HTML "<P>#%-11s %-5s +- %4s %4s %4s %3s %3s %3s %3s %3s %3s %5s %6s %5s %6s %3s %-12s %s",'seqID','acc','aa','pos','neg','d4a','d4b','d3a','d3b','t4a','t4b','m.sc.','m.r.sc','score','rel.sc','num','best_SAH','all_SAHs';
	    printf HTML "<TABLE CELLSPACING=0 BORDER=1><TR BGCOLOR=\"lightgreen\">";
	    foreach (@C1) {
		printf HTML "<TH>%s</TH>", $_;
	    }
#	    printf HTML "</TR>";
	}
	else {
#	    printf HTML "<P>#%-17s +- %4s %4s %4s %3s %3s %3s %3s %3s %3s %5s %6s %5s %6s %3s %-12s %s",'seqID','aa','pos','neg','d4a','d4b','d3a','d3b','t4a','t4b','m.sc.','m.r.sc','score','rel.sc','num','best_SAH','all_SAHs';
	    printf HTML "<TABLE CELLSPACING=0 BORDER=1><TR BGCOLOR=\"lightgreen\">";
	    foreach (@C2) {
		printf HTML "<TH>%s</TH>", $_;
	    }
#	    printf HTML "</TR>";
	}
	if ($PVAL) {
#	    printf HTML "   %7s %7s\n",'Z-score','P-value';
	    printf HTML "<TH>%s</TH><TH>%s</TH></TR>\n",'Z-score','P-value';
	}
	else {
	    printf HTML "</TR>\n";
	}
    }
    # write text output even if HTML output is requested:
    if ($SWISS) {
	printf "#%-11s %-5s +- %4s %4s %4s %3s %3s %3s %3s %3s %3s %5s %6s %5s %6s %3s %-12s %s",'seqID','acc','aa','pos','neg','d4a','d4b','d3a','d3b','t4a','t4b','m.sc.','m.r.sc','score','rel.sc','num','best_SAH','all_SAHs';
    }
    else {
	printf "#%-17s +- %4s %4s %4s %3s %3s %3s %3s %3s %3s %5s %6s %5s %6s %3s %-12s %s",'seqID','aa','pos','neg','d4a','d4b','d3a','d3b','t4a','t4b','m.sc.','m.r.sc','score','rel.sc','num','best_SAH','all_SAHs';
    }
    if ($PVAL) {
	printf "   %7s %7s\n",'Z-score','P-value';
    }
    else {
	printf "\n";
    }

} #_sub Init


#--------------------
sub Read_parameter_table {
    my $i = 0;
    my $len = my $gap = 0;
    my $mean = my $sd = my $loc = my $scale = my $shape = 0;
    my $key = '';

    %PARAMS = ();

    open (PARAM,"$paramfile") || die "Cannot open file $paramfile (option --paramfile): $!\n";
    while (<PARAM>) {
	chomp;
	next if (/^$/ || /^\#/);
	$i++;
	($len,$gap,$mean,$sd,$loc,$scale,$shape,undef) = split;
	$key = "${len}:${gap}";
	$PARAMS{$key}->{'mean'}  = $mean;
	$PARAMS{$key}->{'sd'}    = $sd;
	$PARAMS{$key}->{'loc'}   = $loc;
	$PARAMS{$key}->{'scale'} = $scale;
	$PARAMS{$key}->{'shape'} = $shape;
    }
    close(PARAM);

    print STDERR "$i lines have been read from parameter file \"$paramfile\"." unless $QUIET;

} #_sub Read_parameter_table


#--------------------
sub Read_sequences {
#    @Time = localtime(time);
#    printf STDERR "Processing of sequences started at %.2d:%.2d:%.2d\n\n", @Time[2,1,0];

    my $firstline = 1;
    my $SEQ = 0;
    local $locus = '';
    local $acc = '';
    local $def = '';
    local $l = 0;
    local $sequence = '';

  SEQ_READ:
    while (<>) {
	chomp;
	s/[\r\t ]+$//g;		# remove whitespaces at the end of line
	next if (/^$/);         # ignore blank lines

	if ($firstline) {     # check first line of input file for correct format
	    if ($_ !~ /^>/) {
		die "\nWrong format, check input file!\n";
	    }
	    $firstline = 0;
	}
	if (/^>/) {
	    if ($SEQ) {
		&Process_sequence;
		$SEQ = 0;
	    }
	    $SEQ = 0;
	    ($locus,$def) = split(' ',$_,2);
	    $locus =~ s/^>//;
	    $locus =~ s/sp\|//;
	    if ($SWISS) {
		($acc,$locus) = split(/\|/,$locus,2);
	    }
	    next;
	}
	else {
	    $sequence .= $_;         # read the whole sequence into one string
	    $SEQ = 1;
	}
    }   # label SEQ_READ

    if ($SEQ) {
	&Process_sequence;
    }

    print STDERR ("\nA total of $l sequences processed.\n") unless $QUIET;

    if ($HTML) {
#	print HTML "</PRE>";
	print HTML "</TABLE>";
	print HTML "\n</FONT>";
	print HTML "</BODY></HTML>";
	close(HTML);
#	select $oldfh;
    }


#    @Time = localtime(time);
#    printf STDERR ("Processing finished at %.2d:%.2d:%.2d\n", @Time[2,1,0]);

} #_sub Read_sequences


#--------------------
sub Process_sequence {
    my $i = 0;
    my $k = 0;
    local @Seq = ();
    my $aa = '';
    local $lastresidue = 0;
    local @SeqCharge = ();
    local @Charged = ();
    local @InDyad = ();
    local @InRegion = ();
    my $charged = 0;
    my $positive = 0;
    my $negative = 0;
    local @Dyad = ();
    local @S_Dyad = ();

    local $dyads_4_acidic = 0;
    local $dyads_4_basic = 0;
    local $dyads_3_acidic = 0;
    local $dyads_3_basic = 0;
    local $dyads_2_acidic = 0;
    local $dyads_2_basic = 0;
    local $dyads_1_acidic = 0;
    local $dyads_1_basic = 0;
    local $triads_4_acidic = 0;
    local $triads_4_basic = 0;

    local %Dyads_4_acidic = ();
    local %Dyads_4_basic = ();
    local %Dyads_3_acidic = ();
    local %Dyads_3_basic = ();
    local %Dyads_2_acidic = ();
    local %Dyads_2_basic = ();
    local %Dyads_1_acidic = ();
    local %Dyads_1_basic = ();
    local %Triads_4_acidic = ();
    local %Triads_4_basic = ();

    local $bad_dyads_4_acidic = 0;
    local $bad_dyads_4_basic = 0;
    local $bad_dyads_3_acidic = 0;
    local $bad_dyads_3_basic = 0;
    local %BadDyads_4_acidic = ();
    local %BadDyads_4_basic = ();
    local %BadDyads_3_acidic = ();
    local %BadDyads_3_basic = ();

    my %seen = ();
    local @Runs = ();
    local @CorrRuns = ();
    local %DYADRUN = ();
    my $runs = my $runs4h = '';
    my $sah_num = 0;
    my %was = ();
    my $plusminus = '';
    my %SCORE = ();
    my %RELSCORE = ();
    my $score = 0;
    my $relscore = 0;
    my $totalscore = 0;
    my $totrelscore = 0;
    local $total_len = 0;
    my $maxscore = 0;
    my $maxrelscore = -1000; # score can be negative because of the penalties
                    # this is set to allow the first real value to replace it
    my $zsore = my $pvalue = 0;
    my $key = '';
    my $THIS = 0;
    my $LEN = 0;
    my $seqhtml = '';

    $l++;
    $sequence =~ tr/1 \t\n//d;  # remove 1, whitespaces and newlines 
                                #  from string
#    my $init_len = length($sequence);
    @Seq = split(//,$sequence);
    undef($sequence);


    $lastresidue = $#Seq;                    # last index in array @Seq
#    my $final_len = $lastresidue + 1;

    for ($i = 0; $i <= $lastresidue; $i++) {   # step along sequence
	$aa = $Seq[$i];
	if (defined($MATRIX{$aa})) {
	    $SeqCharge[$i] = $MATRIX{$aa};   # fill array with aa charge for each position
	}
	else {
	    $SeqCharge[$i] = 0;
	}
	push @Charged, $i unless ($SeqCharge[$i] == 0);
	$positive++ if ($SeqCharge[$i] > 0);
	$negative++ if ($SeqCharge[$i] < 0);

	$InDyad[$i] = 0;  # is this position in a dyad to be reported?
                          # (fill the array with starting values: 0=no)
	$InRegion[$i] = 0;  # is this position in a putative CSAH region to be reported?
                            # (fill the array with starting values: 0=no)
    }

    foreach $i (@Charged) {
	if ($i <= $lastresidue-4 && $SeqCharge[$i]*$SeqCharge[$i+4] < 0) {
	    push @Dyad, $i unless $seen{$i}++;
	    push @Dyad, ($i+4) unless $seen{$i+4}++;
	} 
	if ($i <= $lastresidue-3 && $SeqCharge[$i]*$SeqCharge[$i+3] < 0) {
	    push @Dyad, $i unless $seen{$i}++;
	    push @Dyad, ($i+3) unless $seen{$i+3}++;
	} 
    }

    my $prev = -1;
    my $beg = -1;
    my $end = -1;
    my $len = 0;
    my $pos = '';
    my $r = my $cr = '';
    my $region = '';
    my @DyadRun = ();
    my $maxlen = 0;
    local $maxregion = '';
    @seen = ();

#    @S_Dyad = sort by_number @Dyad;
    @S_Dyad = sort {$a <=> $b} @Dyad;
    foreach $i (@S_Dyad) {
	if ($prev >= 0) {
	    if ($i - $prev <= $gap) {
		if ($beg < 0) {
		    $beg = $prev;
#		    push @DyadRun, $beg unless $seen{$beg}++;
		    push @DyadRun, $beg;
		}
		$end = $i;
		push @DyadRun, $end;
	    }
	    else {
		$len = $end - $beg + 1;
		if ($len >= $minlen) {
		    $LEN = 1;
#		    $total_len += $len;
		    $r = "$len:$beg-$end";                # 0-based position numbering
		    $cr = $len.":".($beg+1)."-".($end+1);  # corrected (1-based) numbering
		    push @Runs, $r;              # array of all series of dyads longer than min_length
		    push @CorrRuns, $cr;         # array of all series of dyads (corrected positions)
		    $DYADRUN{$r} = [ @DyadRun ]; # arrays of dyad-forming positions for each series
		    if ($len > $maxlen) {        # find longest series of dyads
			$maxlen = $len;
#			$maxregion = $r;
		    }

#		    print STDERR "$locus $len:$beg-$end";  # DEBUG
		}

		@DyadRun = ();
		@seen = ();
		$end = 0;
		$len = 0;

		$beg = $i;
		push @DyadRun, $beg;
	    }
	} # if ($prev >= 0)

	$prev = $i;

    } # foreach $i (@S_Dyad)

#    if ($beg && $end) {          # unprocessed region left
    if (($beg >= 0) && $end) {          # unprocessed region left, bug fix for v2.3.2
	$len = $end - $beg + 1;
	if ($len >= $minlen) {
	    $LEN = 1;
#	    $total_len += $len;
	    $r = "$len:$beg-$end";                # 0-based position numbering
	    $cr = $len.":".($beg+1)."-".($end+1);  # corrected (1-based) numbering
	    push @Runs, $r;              # array of all series of dyads longer than min_length
	    push @CorrRuns, $cr;         # array of all series of dyads (corrected positions)
	    $DYADRUN{$r} = [ @DyadRun ];
	    if ($len > $maxlen) {     # find longest series of dyads
		$maxlen = $len;
#		$maxregion = $r;
	    }

#	    print STDERR "$locus $len $beg-$end";  # DEBUG
	}
	$beg = $end = $len = 0;
    }

    if (@Runs) {
	foreach $r (@Runs) {
	    foreach $i (@{ $DYADRUN{$r} }) {
		$InDyad[$i] = 1;    # is this position in a dyad to be reported? (1=yes)
	    }
	}
#	$runs = join(',',@Runs);      # 0-based position numbering
	$runs = join(',',@CorrRuns);  # corrected (1-based) position numbering
	&Find_dyads;
	$sah_num = $#Runs + 1;
    }
    else {
	$runs = '-';
	$sah_num = 0;
    } # if (@Runs)

    foreach $r (@Runs) {
	($len,$pos) = split(':',$r);
	($beg,$end) = split('-',$pos);
	$cr = $len.":".($beg+1)."-".($end+1);  # corrected (1-based) numbering

	# Positive scores:

	# opposite charges at i,i+4 and i,i+3 positions
#	$score = $dyad4 * ($#{ $Dyads_4_acidic{$r} } + 1 + $#{ $Dyads_4_basic{$r} } + 1);
#	$score += $dyad3 * ($#{ $Dyads_3_acidic{$r} } + 1 + $#{ $Dyads_3_basic{$r} } + 1);

	# modified in v0.8
	# orientation scores replaced by scoring the two orientations of i,i+4 and i,i+3 dyads separately
	$score = $dyad4a * ($#{ $Dyads_4_acidic{$r} } + 1);
	$score += $dyad3a * ($#{ $Dyads_3_acidic{$r} } + 1);
	$score += $dyad4b * ($#{ $Dyads_4_basic{$r} } + 1);
	$score += $dyad3b * ($#{ $Dyads_3_basic{$r} } + 1);

	# acidic (negatively charged) residue first (N-term)
#	$score += $acidic4first * ($#{ $Dyads_4_acidic{$r} } + 1);
#	$score += $acidic3first * ($#{ $Dyads_3_acidic{$r} } + 1);

	# alternating charges at i,i+4,i+8 positions
#	$score += $triad * ($#{ $Triads_4_acidic{$r} } + 1 + $#{ $Triads_4_basic{$r} } + 1);
	# alternating charges at i,i+4,i+8 or i,i+3|4,i+7 or i,i+3,i+6 positions
	$score += $triada * ($#{ $Triads_4_acidic{$r} } + 1);
	$score += $triadb * ($#{ $Triads_4_basic{$r} } + 1);

	# Penalties (negative scores):

	# identical charges at i,i+4 and i,i+3 positions
	$score -= $dyad4p * ($#{ $BadDyads_4_acidic{$r} } + 1 + $#{ $BadDyads_4_basic{$r} } + 1);
	$score -= $dyad3p * ($#{ $BadDyads_3_acidic{$r} } + 1 + $#{ $BadDyads_3_basic{$r} } + 1);

	# opposite charges at i,i+2 and i,i+1 positions
	$score -= $dyad2p * ($#{ $Dyads_2_acidic{$r} } + 1 + $#{ $Dyads_2_basic{$r} } + 1);
	$score -= $dyad1p * ($#{ $Dyads_1_acidic{$r} } + 1 + $#{ $Dyads_1_basic{$r} } + 1);


	$SCORE{$r} = $score;
	$relscore = round($score/$len,3);
	$RELSCORE{$r} = $relscore;
	if ($RELSCORE{$r} >= $threshold) {
	    $THIS = 1;
	}
	$totalscore += $score;
	$total_len += $len;

#	if ($r eq $maxregion) {   # modified in v0.7
#	    $maxscore = $score;
#	    $maxrelscore = $relscore;
#	}

#	if ($score > $maxscore) {        # find the region with maximum score
	if ($relscore > $maxrelscore) {  # find the region with maximum relative score
#	    $maxregion = $r;             # 0-based position numbering
	    $maxregion = $cr;            # corrected (1-based) numbering
	    $maxscore = $score;
	    $maxrelscore = $relscore;
	}

    } # foreach $r (@Runs)

    if ($total_len > 0) {
	$totrelscore = round($totalscore/$total_len,3);
	if ($totrelscore >= $threshold) {
	    $THIS = 1;
	}
    }
    else {
	$totrelscore = 0;
    }

    $key = "${minlen}:${gap}";
    if ($PVAL) {
	$zscore = ($maxrelscore - $PARAMS{$key}->{'mean'})/$PARAMS{$key}->{'sd'};
	$pvalue = &P_value($zscore,$PARAMS{$key}->{'loc'},$PARAMS{$key}->{'scale'},$PARAMS{$key}->{'shape'});
    }

#    push @Scores, $score;

    if ($THIS && $LEN) {
	$plusminus = '+';
    }
    elsif ($VERBOSE) {
	$plusminus = '-';
    }
    if ($plusminus && $LEN) {
	if ($HTML) {
	    $runs4h = $runs;
	    $runs4h =~ s/,/,<br>/g;
	    $runs4h =~ s/-/&ndash;/g;
	    if ($SWISS) {
#		@Data = ($acc,($lastresidue+1),$positive,$negative,$dyads_4_acidic,$dyads_4_basic,$dyads_3_acidic,$dyads_3_basic,$triads_4_acidic,$triads_4_basic,$maxscore,$maxrelscore,$totalscore,$totrelscore,$sah_num,$maxregion,$runs4h);
#		printf HTML "<P><A HREF=\"$htmldir/%s.html\">%-12s</A> %-6s %1s %4d %4d %4d %3d %3d %3d %3d %3d %3d %5d %6.3f %5d %6.3f  %2d %-12s %s",$locus,$locus,$acc,$plusminus,($lastresidue+1),$positive,$negative,$dyads_4_acidic,$dyads_4_basic,$dyads_3_acidic,$dyads_3_basic,$triads_4_acidic,$triads_4_basic,$maxscore,$maxrelscore,$totalscore,$totrelscore,$sah_num,$maxregion,$runs;
#		printf HTML "<TR><TD><A HREF=\"$htmldir/%s.html\">%s</A></TD><TD>%s</TD><TD>%s</TD><TD>%4d</TD><TD>%4d</TD><TD>%4d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%5d</TD><TD>%6.3f</TD><TD>%5d</TD><TD>%6.3f</TD><TD>%2d</TD><TD>%s</TD><TD>%s</TD>",$locus,$locus,$acc,$plusminus,($lastresidue+1),$positive,$negative,$dyads_4_acidic,$dyads_4_basic,$dyads_3_acidic,$dyads_3_basic,$triads_4_acidic,$triads_4_basic,$maxscore,$maxrelscore,$totalscore,$totrelscore,$sah_num,$maxregion,$runs4h;
#		printf HTML "<TR><TD><A HREF=\"$htmldir/%s.html\">%s</A></TD><TD>%s</TD><TD>%4d</TD><TD>%4d</TD><TD>%4d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%5d</TD><TD>%6.3f</TD><TD>%5d</TD><TD>%6.3f</TD><TD>%2d</TD><TD>%s</TD><TD>%s</TD>",$locus,$locus,$acc,($lastresidue+1),$positive,$negative,$dyads_4_acidic,$dyads_4_basic,$dyads_3_acidic,$dyads_3_basic,$triads_4_acidic,$triads_4_basic,$maxscore,$maxrelscore,$totalscore,$totrelscore,$sah_num,$maxregion,$runs4h;
		printf HTML "<TR><TD align=\"center\"><A HREF=\"$htmldir/%s.html\">%s</A></TD>",$locus,$locus;
		printf HTML "<TD align=\"center\">%s</TD><TD align=\"center\">%4d</TD><TD align=\"center\">%4d</TD><TD align=\"center\">%4d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%5d</TD><TD align=\"center\">%6.3f</TD><TD align=\"center\">%5d</TD><TD align=\"center\">%6.3f</TD><TD align=\"center\">%2d</TD><TD align=\"center\"><nobr>%s</nobr></TD><TD align=\"center\"><nobr>%s</nobr></TD>",$acc,($lastresidue+1),$positive,$negative,$dyads_4_acidic,$dyads_4_basic,$dyads_3_acidic,$dyads_3_basic,$triads_4_acidic,$triads_4_basic,$maxscore,$maxrelscore,$totalscore,$totrelscore,$sah_num,$maxregion,$runs4h;
	    }
	    else {
#		printf HTML "<P><A HREF=\"$htmldir/%s.html\">%-19s</A> %1s %4d %4d %4d %3d %3d %3d %3d %3d %3d %5d %6.3f %5d %6.3f  %2d %-12s %s",$locus,$locus,$plusminus,($lastresidue+1),$positive,$negative,$dyads_4_acidic,$dyads_4_basic,$dyads_3_acidic,$dyads_3_basic,$triads_4_acidic,$triads_4_basic,$maxscore,$maxrelscore,$totalscore,$totrelscore,$sah_num,$maxregion,$runs;
#		printf HTML "<TR><TD><A HREF=\"$htmldir/%s.html\">%s</A></TD><TD>%s</TD><TD>%4d</TD><TD>%4d</TD><TD>%4d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%5d</TD><TD>%6.3f</TD><TD>%5d</TD><TD>%6.3f</TD><TD>%2d</TD><TD>%s</TD><TD>%s</TD>",$locus,$locus,$plusminus,($lastresidue+1),$positive,$negative,$dyads_4_acidic,$dyads_4_basic,$dyads_3_acidic,$dyads_3_basic,$triads_4_acidic,$triads_4_basic,$maxscore,$maxrelscore,$totalscore,$totrelscore,$sah_num,$maxregion,$runs4h;
#		printf HTML "<TR><TD><A HREF=\"$htmldir/%s.html\">%s</A></TD><TD>%4d</TD><TD>%4d</TD><TD>%4d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%3d</TD><TD>%5d</TD><TD>%6.3f</TD><TD>%5d</TD><TD>%6.3f</TD><TD>%2d</TD><TD>%s</TD><TD>%s</TD>",$locus,$locus,($lastresidue+1),$positive,$negative,$dyads_4_acidic,$dyads_4_basic,$dyads_3_acidic,$dyads_3_basic,$triads_4_acidic,$triads_4_basic,$maxscore,$maxrelscore,$totalscore,$totrelscore,$sah_num,$maxregion,$runs4h;
		printf HTML "<TR><TD align=\"center\"><A HREF=\"$htmldir/%s.html\">%s</A></TD>",$locus,$locus;
		printf HTML "<TD align=\"center\">%4d</TD><TD align=\"center\">%4d</TD><TD align=\"center\">%4d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%3d</TD><TD align=\"center\">%5d</TD><TD align=\"center\">%6.3f</TD><TD align=\"center\">%5d</TD><TD align=\"center\">%6.3f</TD><TD align=\"center\">%2d</TD><TD align=\"center\"><nobr>%s</nobr></TD><TD align=\"center\"><nobr>%s</nobr></TD>",($lastresidue+1),$positive,$negative,$dyads_4_acidic,$dyads_4_basic,$dyads_3_acidic,$dyads_3_basic,$triads_4_acidic,$triads_4_basic,$maxscore,$maxrelscore,$totalscore,$totrelscore,$sah_num,$maxregion,$runs4h;
	    }
	    if ($PVAL) {
#		printf HTML " %7.3f %7.3f\n",$zscore,$pvalue;
		printf HTML "<TD align=\"center\">%7.3f</TD><TD align=\"center\">%7.3f</TD></TR>\n",$zscore,$pvalue;
	    }
	    else {
#		printf HTML "\n";
		printf HTML "</TR>\n";
	    }
	}
        # write text output even if HTML output is requested:
	if ($SWISS) {
	    printf "%-12s %-6s %1s %4d %4d %4d %3d %3d %3d %3d %3d %3d %5d %6.3f %5d %6.3f  %2d %-12s %s",$locus,$acc,$plusminus,($lastresidue+1),$positive,$negative,$dyads_4_acidic,$dyads_4_basic,$dyads_3_acidic,$dyads_3_basic,$triads_4_acidic,$triads_4_basic,$maxscore,$maxrelscore,$totalscore,$totrelscore,$sah_num,$maxregion,$runs;
	}
	else {
	    printf "%-19s %1s %4d %4d %4d %3d %3d %3d %3d %3d %3d %5d %6.3f %5d %6.3f  %2d %-12s %s",$locus,$plusminus,($lastresidue+1),$positive,$negative,$dyads_4_acidic,$dyads_4_basic,$dyads_3_acidic,$dyads_3_basic,$triads_4_acidic,$triads_4_basic,$maxscore,$maxrelscore,$totalscore,$totrelscore,$sah_num,$maxregion,$runs;
	}
	if ($PVAL) {
	    printf " %7.3f %7.3f\n",$zscore,$pvalue;
	}
	else {
	    printf "\n";
	}
    }

    if ($HTML && $plusminus && $LEN) {
	$seqhtml = $bdir . $htmldir . "/" . $locus . ".html";
	open(SEQHTML,">$seqhtml");
	print SEQHTML "<HTML>\n<HEAD><TITLE>$locus - annotated sequence</TITLE></HEAD>\n<BODY>\n";
	print SEQHTML "<FONT FACE=\"Courier New, Courier\">";
	if (defined($def)) {
	    if ($SWISS) {
		print SEQHTML "<B>&gt;$locus $acc</B> $def";
	    }
	    else {
		print SEQHTML "<B>&gt;$locus</B> $def";
	    }
	}
	else {
	    if ($SWISS) {
		print SEQHTML "<B>&gt;$locus $acc</B>";
	    }
	    else {
		print SEQHTML "<B>&gt;$locus</B>";
	    }
	}
	for ($i = 0; $i <= $lastresidue; $i++) {   # step along sequence
	    if ($i % 60 == 0) {
		print SEQHTML "<BR>";
	    }
	    if ($InRegion[$i]) {
		printf SEQHTML "<B>";
		if ($InDyad[$i]) { printf SEQHTML "<U>" }
		if ($SeqCharge[$i] > 0) {
		    printf SEQHTML "<FONT COLOR=\"#0000FF\">";
		}
		if ($SeqCharge[$i] < 0) {
		    printf SEQHTML "<FONT COLOR=\"#FF0000\">";
		}
	    }
	    printf SEQHTML $Seq[$i];
	    if ($InRegion[$i]) {
		if ($SeqCharge[$i] > 0 || $SeqCharge[$i] < 0 ) {
		    printf SEQHTML "</FONT>";
		}
		if ($InDyad[$i]) { printf SEQHTML "</U>" }
		printf SEQHTML "</B>";
	    }
	}
	print SEQHTML "\n";
	print SEQHTML "</FONT>";
	print SEQHTML "\n</BODY></HTML>";
	close(SEQHTML);
    }

} #_sub Process_sequence


sub Find_dyads {

    my $i = 0;
    my $k = 0;
    my $beg = 0;
    my $end = 0;
    my $len = 0;
    my $r = '';
    my %seen = ();

    foreach $r (@Runs) {
	(undef,$region) = split(':',$r);
	($beg,$end) = split('-',$region);
	for ($i = $beg; $i <= $end; $i++) {
	    $InRegion[$i] = 1; # is this position in a putative CSAH region to be reported? (1=yes)
	}

#	print STDERR "$locus   $r   ".join(',',@{ $DYADRUN{$r} });    # DEBUG
	foreach $i (@{ $DYADRUN{$r} }) {
	    if ($SeqCharge[$i] > 0) {                       # Positive charge first (N-term)
		if ($i <= $end-4 && $SeqCharge[$i+4] < 0) {
                    $dyads_4_basic++;
		    push @{ $Dyads_4_basic{$r} },$i;
		}
		if ($i <= $end-3 && $SeqCharge[$i+3] < 0) {
                    $dyads_3_basic++;
		    push @{ $Dyads_3_basic{$r} },$i;
		}
		if ($i <= $end-8 && $SeqCharge[$i+4] < 0 && $SeqCharge[$i+8] > 0) {
                    $triads_4_basic++;
		    push @{ $Triads_4_basic{$r} },$i;
		}
		if ($i <= $end-7 && $SeqCharge[$i+4] < 0 && $SeqCharge[$i+7] > 0) {
                    $triads_4_basic++;
		    push @{ $Triads_4_basic{$r} },$i;
		}
		if ($i <= $end-7 && $SeqCharge[$i+3] < 0 && $SeqCharge[$i+7] > 0) {
                    $triads_4_basic++;
		    push @{ $Triads_4_basic{$r} },$i;
		}
		if ($i <= $end-6 && $SeqCharge[$i+3] < 0 && $SeqCharge[$i+6] > 0) {
                    $triads_4_basic++;
		    push @{ $Triads_4_basic{$r} },$i;
		}
		if ($i <= $end-2 && $SeqCharge[$i+2] < 0) { # +- at i,i+2
                    $dyads_2_basic++;
		    push @{ $Dyads_2_basic{$r} },$i;
		}
		if ($i <= $end-1 && $SeqCharge[$i+1] < 0) { # +- at i,i+1 (adjacent residues)
                    $dyads_1_basic++;
		    push @{ $Dyads_1_basic{$r} },$i;
		}
		if ($i <= $end-4 && $SeqCharge[$i+4] > 0) { # ++ at i,i+4
                    $bad_dyads_4_basic++;
		    push @{ $BadDyads_4_basic{$r} },$i;
		}
		if ($i <= $end-3 && $SeqCharge[$i+3] > 0) { # ++ at i,i+3
                    $bad_dyads_3_basic++;
		    push @{ $BadDyads_3_basic{$r} },$i;
		}
	    }
	    if ($SeqCharge[$i] < 0) {                       # Negative charge first (N-term)
		if ($i <= $end-4 && $SeqCharge[$i+4] > 0) {
                    $dyads_4_acidic++;
		    push @{ $Dyads_4_acidic{$r} },$i;
		}
		if ($i <= $end-3 && $SeqCharge[$i+3] > 0) {
                    $dyads_3_acidic++;
		    push @{ $Dyads_3_acidic{$r} },$i;
		}
		if ($i <= $end-8 && $SeqCharge[$i+4] > 0 && $SeqCharge[$i+8] < 0) {
                    $triads_4_acidic++;
		    push @{ $Triads_4_acidic{$r} },$i;
		}
		if ($i <= $end-7 && $SeqCharge[$i+4] > 0 && $SeqCharge[$i+7] < 0) {
                    $triads_4_acidic++;
		    push @{ $Triads_4_acidic{$r} },$i;
		}
		if ($i <= $end-7 && $SeqCharge[$i+3] > 0 && $SeqCharge[$i+7] < 0) {
                    $triads_4_acidic++;
		    push @{ $Triads_4_acidic{$r} },$i;
		}
		if ($i <= $end-6 && $SeqCharge[$i+3] > 0 && $SeqCharge[$i+6] < 0) {
                    $triads_4_acidic++;
		    push @{ $Triads_4_acidic{$r} },$i;
		}
		if ($i <= $end-2 && $SeqCharge[$i+2] > 0) { # -+ at i,i+2
                    $dyads_2_acidic++;
		    push @{ $Dyads_2_acidic{$r} },$i;
		}
		if ($i <= $end-1 && $SeqCharge[$i+1] > 0) { # -+ at i,i+1 (adjacent residues)
                    $dyads_1_acidic++;
		    push @{ $Dyads_1_acidic{$r} },$i;
		}
		if ($i <= $end-4 && $SeqCharge[$i+4] < 0) { # -- at i,i+4
                    $bad_dyads_4_acidic++;
		    push @{ $BadDyads_4_acidic{$r} },$i;
		}
		if ($i <= $end-3 && $SeqCharge[$i+3] < 0) { # -- at i,i+3
                    $bad_dyads_3_acidic++;
		    push @{ $BadDyads_3_acidic{$r} },$i;
		}
	    }
	}
    }

} #_sub Find_dyads

#--------------------
sub P_value {

# Usage:
#   $z_score = ($score - $mean)/$sd;
#   $p_value = &P_value($z_score,$loc,$scale,$shape);

    my $x = $_[0];               # Input parameter
    my $loc = $_[1];             # location parameter of GEV distribution
    my $scale = $_[2];           # scale parameter of GEV distribution
    my $shape = $_[3];           # shape parameter of GEV distribution

    my $z = 0;
    my $pvalue = 0;

    $z = ($x - $loc) / $scale;
    $pvalue = 1 - exp(-(1+$shape*$z)**(-1/$shape));

    return $pvalue;

} #_sub P_value

#--------------------
sub round {
    # Why do we need this function?
    # It seems that e.g. 'sprintf "%4.2f",$fl' incorrectly
    #  rounds 0.125 to 0.12 (should be 0.13).
    #
    # Use: round(0.123456,4) to get 0.1235

    my $numbertoround = $_[0];        # Input number to be rounded
    my $roundtodigits = $_[1];        # The number of decimal places kept
    my $roundednum = 0;
    my $factor = 0;

    $factor = 10**$roundtodigits;
    $roundednum = int($numbertoround*$factor+0.5)/$factor;

    return $roundednum;
} #_sub round

#--------------------
sub by_number {
    $a <=> $b;
}

#--------------------
sub Message {
    print STDERR "
scan4csah.pl $VERSION

Scans protein sequences for charged amino acid patterns indicative
of stable alpha-helix (CSAH domain).

Usage:  scan4csah.pl [options] < infile > outfile

";
} #_sub Message

__END__


=head1 NAME

scan4csah.pl - Scans protein sequence(s) for charged amino acid patterns
               indicative of stable alpha-helix (CSAH domain)

=head1 SYNOPSIS

scan4csah.pl [options] < infile > outfile

Options:
  --help
  --man
  --verbose
  --matrix=input_matrix_file
  --gap=max_gap
  --minlen=min_length_of_CSAH
  --threshold=score_threshold
  --dyad4acidic=dyad4a_score
  --dyad3acidic=dyad3a_score
  --dyad4basic=dyad4b_score
  --dyad3basic=dyad3b_score
  --triadacidic=triada_score
  --triadbasic=triadb_score
  --dyad4pen=dyad4_penalty
  --dyad3pen=dyad3_penalty
  --dyad2pen=dyad2_penalty
  --dyad1pen=dyad1_penalty
  --html=html_output_file
  --swissprot
  --paramfile=input_parameter_table_file

=head1 OPTIONS

=over 8

=item B<--help>

Prints a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=item B<--verbose>

Verbose output (write lines with negative result, too)

=item B<--matrix> I<input_matrix_file>

A 2-column matrix of residue scores (format: AMINO_ACID CHARGE);
if not specified, the following default values are used:
D,E = -1, R,K = +1

=item B<--gap> I<max_gap>

Maximum gap between two consecutive residues that are members of any dyad (default: 5, minimum: 4)

=item B<--minlen> I<min_length_of_CSAH>

Minimum length of a candidate region with consecutive dyads to be accepted as putative CSAH domain

=item B<--threshold> I<score_threshold>

Minimum relative score threshold to report a sequence entry
(absolute score per number of charged residues)
(default: 2.5)

=item B<--dyad4acidic> I<dyad4a_score>

Score given for an i,i+4 (-+) dyad of oppositely charged residues

=item B<--dyad4basic> I<dyad4b_score>

Score given for an i,i+4 (+-) dyad of oppositely charged residues

=item B<--dyad3acidic> I<dyad3a_score>

Score given for an i,i+3 (-+) dyad of oppositely charged residues

=item B<--dyad3basic> I<dyad3b_score>

Score given for an i,i+3 (+-) dyad of oppositely charged residues

=item B<--triadacidic> I<triada_score>

Extra score given for an i,i+3|4,i+6|7|8 triad (-+-)

=item B<--triadbasic> I<triadb_score>

Extra score given for an i,i+3|4,i+6|7|8 triad (+-+)

=item B<--dyad4pen> I<dyad4_penalty>

Penalty given for an i,i+4 dyad of identically charged residues (default: 2)

=item B<--dyad3pen> I<dyad3_penalty>

Penalty given for an i,i+3 dyad of identically charged residues (default: 2)

=item B<--dyad2pen> I<dyad2_penalty>

Penalty given for an i,i+2 dyad of oppositely charged residues (default: 2)

=item B<--dyad1pen> I<dyad1_penalty>

Penalty given for an i,i+1 dyad of oppositely charged residues (default: 2)

=item B<--html> I<html_output_file>

Base name for the HTML-format output file (".html" will be added to this name). This output includes the output table, otherwise written to standard output, with links to individual HTML files of annotated sequences. Charged amino acid residues in putative CSAHs are emphasized in the latter files: B<red> -> negatively charged, B<blue> -> positively charged, B<underlined> -> dyad-forming residue, B<bold> -> region of putative CSAH. HTML sequence files are written into a subdirectory named "html". This directory will be created and its content will be overwritten without asking for confirmation. HTML output is not written by default.

=item B<--swissprot>

Input sequence file has SwissProt format sequence headers with accession number as second field after sequence ID. Write protein accession number into output table. Accession number is not written by default.

=item B<--paramfile> I<input_parameter_table_file>

Input table with precomputed parameters of max. rel. score distribution and precomputed parameters of EVD distribution fitted upon the Z-score distribution. Precomputed values are normally calculated by running 'scan4csah.pl' and 'plot4csah.pl' on all protein sequences in Swiss-Prot. (Columns in the table: 1. min. CSAH length, 2. max. gap size within a CSAH, 3. mean of all max. rel. scores, 4. standard deviation of all max. rel. scores, 5. EVD location parameter, 6. EVD scale parameter, 7. EVD shape parameter.)


=back

=head1 INPUT

Input file: FASTA format protein sequence file

=head1 OUTPUT

Output file (written to standard output): table with the following columns:

  1. seqID   : sequence ID
  2. acc     : accession number  (only if --swissprot is set)
  3. +-      : + or -  (+: putative CSAH domain-containing protein;
                        -: all other, written only if option --verbose is given)
  4. aa      : length of sequence (number of amino acid residues)
  5. pos     : number of positively charged residues (according to score matrix)
  6. neg     : number of negatively charged residues (according to score matrix)
  7. d4a     : number of [i,i+4] dyads (N-term residue is acidic)
  8. d4b     : number of [i,i+4] dyads (N-term residue is basic)
  9. d3a     : number of [i,i+3] dyads (N-term residue is acidic)
 10. d3b     : number of [i,i+3] dyads (N-term residue is basic)
 11. t4a     : number of [i,i+4,i+8] triads (N-term residue is acidic)
 12. t4b     : number of [i,i+4,i+8] triads (N-term residue is basic)
 13. m.sc.   : maximum score (absolute score of longest putative CSAH)
 14. m.r.sc. : maximum relative score (relative score of highest-scoring putative CSAH)
               (absolute score per length of CSAH region)
 15. score   : combined absolute score of all putative CSAH regions
 16. rel.sc. : combined relative score of all putative CSAH regions
               (combined absolute score per total lengths of putative CSAH regions)
 17. num     : number of putative CSAH regions
 18. best SAH: region of best (longest) putative CSAH
               (format: LENGTH1:START1-END1)
 19. all SAHs: region(s) of putative CSAH(s)
               (format: LENGTH1:START1-END1,LENGTH2:START2-END2 ...)

 20. Z-score : Z-score (m.r.sc.) calculated using precomputed distribution parameters
 21. P-value : P-value (m.r.sc.) calculated using precomputed EVD parameters
               (the last 2 columns are written only if --paramfile is set)


=head1 DESCRIPTION

B<scan4csah.pl> will scan protein sequences for charged amino acid patterns
indicative of stable alpha-helix (SAH domain).

=head1 VERSION

2.3.1

=head1 DATE

13 January 2010

=head1 AUTHOR

Gabor Toth  < tothg[at]abc[dot]hu >

=cut
