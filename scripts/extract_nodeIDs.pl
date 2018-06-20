#!/usr/bin/perl
# Extract nodeIDs information from branch-site analysis

use strict;
use warnings;

#Declare usage for standard input
my($USAGE) = "Incomplete number of arguments in $0\n Arguments: [1]Node-IDs file\n";

my $num_args = $#ARGV + 1;
if ($num_args != 1) {
print $USAGE;
	exit;
	}

# declare and initialize variables
my @annotation = (  );
my @ids = (  );
my $filename = $ARGV[0];

parse_ids(\@annotation, \@ids, $filename);

# Print the IDs annotation to check if we got it okay.
print @ids;
exit;

###SUBROUTINES###

# parse_ids
#
# parse IDs annotation from Branch information file

sub parse_ids {

    my($annotation, $ids, $filename) = @_;

       
    # declare and initialize variables
    my $in_id = 0; 
    my @BranchFile = (  );
    
    # Get the Branch-site output into an array from a file
    @BranchFile = get_file_data($filename);
    
    # Extract all the sequence lines
    foreach my $line (@BranchFile) {

        if( $line =~ /Descendant/ ) { # If $line is begin of next record SUMMARY BY MODEL,
            last; #break out of the foreach loop.
        } elsif( $in_id) { # If we know we're in the LRT information,
            push(@ids, $line); # add the current line to @lrt.
        } elsif ( $line =~ /Leaf/ ) { # If $line begins the LRT information,
            $in_id = 1; # set the $in_lrt flag.
        } else{ # Otherwise
            push( @$annotation, $line); # add the current line to @annotation.
        }
    }
    # remove only-whitespace lines 
    s{^\s+$|\s{5}|\s{4}}{}g foreach @ids;
}

# get_file_data
#
# A subroutine to get data from a file given its filename

sub get_file_data {

    my($filename) = @_;

    use strict;
    use warnings;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}
