#!usr/bin/perl
use strict;
use Getopt::Long; 
use Bio::SeqIO;

my($USAGE) = "Incomplete number of arguments in $0\n Arguments: [1]Input alignment (FASTA format), [2]Output file (FASTA) \n\n";

my $num_args = $#ARGV + 1;
if ($num_args != 2) {
print $USAGE;
	exit;
	}


my $infile = $ARGV[0];
my$outfile = $ARGV[1];

my @sequences;
my $codon;

my $nuc  = Bio::SeqIO->new(-file => $infile , '-format' => 'fasta');
my $out = Bio::SeqIO->new(-file => ">$outfile" , '-format' => 'fasta');


while (my $nucseq = $nuc->next_seq) {
		push(@sequences, $nucseq);
}

foreach my $nucseq (@sequences){
	my $nuc_id=$nucseq->id();
	my $nuc_str=uc($nucseq->seq);
		
		my $codon;
		
  for (my $i=0;$i<(length($nuc_str)-2);$i+=3) {
	$codon=substr($nuc_str,$i,3);
	if ($codon =~ /(((U|T)A(A|G|R))|((T|U)GA))/){
		print "stop codon (",$codon,") detected in ",$nuc_id,"\n";
		substr($nuc_str, $i, 3) = "---";
}
}

my $newseq = Bio::Seq->new(-seq => "$nuc_str",                          
                         -display_id => $nuc_id);
  $out->write_seq($newseq); 
}
