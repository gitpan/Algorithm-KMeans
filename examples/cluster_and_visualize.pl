#!/usr/bin/perl -w

#use lib '../blib/lib', '../blib/arch';

use strict;
use Algorithm::KMeans;

#my $datafile = "mydatafile2.dat";
my $datafile = "mydatafile3.dat";
#my $datafile = "mydatafile1.dat";
#my $datafile = "Nadeem.txt";
#my $datafile = "features_temp.dat";
#my $datafile = "fernando.dat";

# Mask:

# The mask tells the module which columns of the data file
# are are to be used for clustering, which columns are to be
# ignored and which column contains the symbolic ID tag for
# a data point.  If the ID is in column 1 and you are
# clustering 3D data, the mast would be "N111".  Note the
# first character in the mask in this case is `N' for
# "Name".  If, on the other hand, you wanted to ignore the
# first data coordinate for clustering, the mask would be
# "N011".  The symbolic ID can be in any column --- you just
# have to place the character `N' at the right place:

my $mask = "N11";       # for mydatafile3.dat
#my $mask = "N111";     # for mydatafile1.dat --- use all three data cols
#my $mask = "N011";     # for mydatafile1.dat --- use all only last two cols
#my $mask = "N100";      # for mydatafile1.dat, Nadeem.txt
#my $mask = "N10";       # for fernando.dat and features_temp.dat datafiles
my $clusterer = Algorithm::KMeans->new( datafile => $datafile,
                                        mask     => $mask,
                                        K        => 2,
                                        cluster_seeding => 'smart',
                                        terminal_output => 1,
#                                        write_clusters_to_files => 1,
                                        debug => 0,
    );

$clusterer->read_data_from_file();

# If you want to access the clusters in your own script:
my ($clusters, $cluster_centers) = $clusterer->kmeans();

# Once you have the clusters in your own top-level script,
# you can now examine the contents of the clusters by the
# following sort of code:
foreach my $cluster (@$clusters) {
    print "Cluster:   @$cluster\n\n"
}


#$clusterer->get_initial_cluster_centers_1_40(3);


# VISUALIZATION:

# Visualization mask:

# In most cases, you would not change the value of the mask
# between clustering and visualization.  But, if you are
# clustering multi-dimensional data and you wish to
# visualize the projection of of the data on each plane
# separately, you can do so by changing the value of the
# visualization mask.  The number of on bits in the
# visualization must not exceed the number of on bits in the
# original data mask.

#my $visualization_mask = "111";   # for mydatafile1.dat with all 3 data cols
my $visualization_mask = "11";   
#my $visualization_mask = "1";  #for fernando.dat, features_temp.dat, Nadeem.txt

$clusterer->visualize_clusters($visualization_mask);


