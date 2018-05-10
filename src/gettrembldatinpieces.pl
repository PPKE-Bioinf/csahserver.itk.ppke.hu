#!/usr/bin/perl

$linenum=0;
open (T,">test.ids") || die "Cannot create test.ids\n";

while(<>){

 print T "$_";
 $linenum++;
 if ($linenum == 5000){
   close(T);
   system "./retrieve_uniprot_dat_for_id_list.pl test.ids";
   $linenum=0;
   open (T,">test.ids") || die "Cannot create test.ids\n";
 }

}
if($linenum > 0){
 system "./retrieve_uniprot_dat_for_id_list.pl test.ids";
}
close(T);
