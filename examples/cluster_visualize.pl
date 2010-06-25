#!/usr/bin/perl -w

use lib '../blib/lib','../blib/arch';

use Algorithm::KMeans qw(kmeans visualize cluster_data_generator);


# Mask:

# The mask tells the module which columns of the data file
# are are to be used for clustering, which columns are to be
# ignored and which column contains the symbolic ID tag for
# a data point.  If the ID is in column 1 and you are
# clustering 3D data, the mast would be "I111".  Note the
# first character in the mask in this case is `I' for "Id".
# If, on the other hand, you wanted to ignore the first data
# coordinate for clustering, the mask would be "I011".  The
# symbolic ID can be in any column --- you just have to
# place the character `I' at the right place:

my $mask = "I111";
my $datafile = "mynewdatafile.dat";

kmeans( datafile => $datafile,
        mask     =>  $mask,
        terminal_output => 1,  
        K => 0 );

visualize( datafile => $datafile,
           mask     =>  $mask );

