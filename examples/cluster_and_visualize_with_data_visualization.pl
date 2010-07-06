#!/usr/bin/perl -w

#use lib '../blib/lib', '../blib/arch';

use strict;
use Algorithm::KMeans;

my $datafile = "mydatafile1.dat";

#  READ the comment block associated with the following call in
#
#             cluster_visualize.pl

my $mask = "N111";
my $clusterer = Algorithm::KMeans->new( datafile => $datafile,
                                        mask     => "N111",
                                        K        => 3,
                                        terminal_output => 1,
                                        do_variance_normalization => 1,
#                                        write_clusters_to_files => 1,
    );

$clusterer->read_data_from_file();
$clusterer->kmeans();


# VISUALIZATION:

# Visualization mask:

# Read the comment block in cluster_and_visualize() that is
# associated with the setting up of the visualization mask.

my $visualization_mask = "111";

# In order to see the effects of variance normalization of
# the data (each data coordinate is normalized by the
# standard-deviation along that coordinate axis), it is
# sometimes useful to see both the raw data and its
# normalized form.  The following two calls accomplish that:

$clusterer->visualize_data($visualization_mask, 'original');

$clusterer->visualize_data($visualization_mask, 'normed');


# Finally, you can visualize the clusters.  BUT NOTE THAT
# THE VISUALIZATION MASK FOR CLUSTER VISUALIZATION WILL, IN
# GENERAL, BE INDEPENDENT OF THE VISUALIZATION MASK FOR
# VIEWING THE DATA:

$clusterer->visualize_clusters($visualization_mask);

