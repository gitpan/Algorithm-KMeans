#!/usr/bin/perl -w

#use lib '../blib/lib', '../blib/arch';

use strict;
use Algorithm::KMeans;

my $datafile = "mydatafile1.dat";

my $mask = "N111";
my $clusterer = Algorithm::KMeans->new( datafile => $datafile,
                                        mask     => "N111",
                                        Kmin     => 3,
                                        Kmax     => 4,
                                        terminal_output => 1,
                                        write_clusters_to_files => 1,
    );

$clusterer->read_data_from_file();

# If you want to access the clusters in your own script:
my ($clusters, $cluster_centers) = $clusterer->kmeans();
foreach my $cluster (@$clusters) {
    print "Cluster:   @$cluster\n\n"
}

# The following call also prints out the best value for K:
my $K_best = $clusterer->get_K_best();

$clusterer->show_QoC_values();


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

my $visualization_mask = "111";
$clusterer->visualize_clusters($visualization_mask);

