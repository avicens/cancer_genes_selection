# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5 
use strict;
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

my($USAGE) = "Incomplete number of arguments in $0\n Arguments: [1]Protein_stable_ID, [2]Ensembl_database \n\n";

my $num_args = $#ARGV + 1;
if ($num_args != 2) {
print $USAGE;
	exit;
	}

# Read in the Protein ID and specific Ensembl database from the argument on the command line.
my($protid) = $ARGV[0];
my($spdb) = $ARGV[1];

my $confFile = "/home/uvi/be/avs/tools/biomart-perl/conf/martURLLocation.xml";
my $action='cached';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

		
	$query->setDataset("${spdb}_gene_ensembl");
	$query->addFilter("ensembl_peptide_id", [$protid]);
	$query->addAttribute("ensembl_gene_id");
	$query->addAttribute("ensembl_transcript_id");
	$query->addAttribute("coding");
	$query->addAttribute("ensembl_peptide_id");
	$query->addAttribute("external_gene_name");

$query->formatter("FASTA");

my $query_runner = BioMart::QueryRunner->new();

############################## GET RESULTS ##########################

$query_runner->execute($query);
$query_runner->printHeader();
$query_runner->printResults();
$query_runner->printFooter();
#####################################################################
