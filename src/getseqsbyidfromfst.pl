#!/usr/bin/perl

print STDERR"
# getgenefromfst.pl -f<file_with_IDs> < fst_file > fst_file_IDsonly

";

use Getopt::Std;
$opt_f="";
getopts('f:hv') || die "Optional error\n";
exit if $opt_h;
open (IF,"$opt_f") || die ("Cannot open inout file $opt_f\n");
while(<IF>){
 chomp;
 $getgene{$_}=1;
}
close(IF);

$write=0;
while(<>){
    chop;
    if ($_ =~ /^>/){
        $write=0;
        $id=$_;
        $id=~s/^.*\>[a-z]+\|([0-9A-Z\-]+)\|.*$/$1/;
        #$id=~s/^.*\>([0-9A-Z\-]+)\|.*$/$1/;
        if ($getgene{$id}){
          $write=1;
        }
#        foreach $getid (keys %getgene){
#          if ($_ =~ /\W$getid\W/){
#            {$write=1}
#          }
#        }
    }
    print "$_\n" if $write;
}


