package Algorithm::KMeans;

#---------------------------------------------------------------------------
# Copyright (c) 2010 Avinash Kak. All rights reserved.
# This program is free software.  You may modify and/or
# distribute it under the same terms as Perl itself.
# This copyright notice must remain attached to the file.
#
# Algorithm::KMeans is a pure Perl implementation for
# clustering multi-dimensional data.
#---------------------------------------------------------------------------

use 5.10.0;
use strict;
use warnings;

use Exporter;
our @ISA = qw(Exporter);

our @EXPORT_OK = qw(kmeans visualize cluster_data_generator);

our $VERSION = '1.0';

my ($terminal_output, %all_data);

# from perl docs:
my $num_regex =  '^[+-]?\ *(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$'; 

sub kmeans {
    my %args = @_;
    die "\n\n$0 requires exactly four command line arguments. " .
               "\n\nExample:\n\n " .
               "      kmeans.pl filename.dat -I0111  1  K=3 " .
               "\n\nARG 1: the name of the data file; " .
               "\n\nARG 2: a mask array expressed as a single string that " .
               "starts with the character `-' and the rest of which " .
               "consists of a single letter I, for the data record " .
               "identifier field and 0s and 1s to indicate which " .
               "fields in the data file records to ignore and which " .
               "fields to treat as numerical data for clustering; " .
               "\n\nARG 3: either 0 or 1, with 1 indicating you wish " .
               "to see the results on the terminal screen and 0 otherwise;" .
               "\n\nARG 4: a string like K=N where N is zero if want the " .
               "program to figure out the best K (the number of clusters), " .
               "or a specific integer for some fixed value for K" 
                                          unless keys %args == 4; 
    my $datafile = $args{datafile};
    my $mask     = $args{mask};
    my $K = $args{K};
    $terminal_output = $args{terminal_output};

    my @mask = split //, $mask;

    open INPUT, $datafile 
        or die "unable to open file $datafile: $!\n";
    chomp( my @raw_data = <INPUT> );
    close INPUT;

    # Transform strings into number data
    foreach my $record (@raw_data) {
        next if $record =~ /^#/;
        my @data_fields;
        my @fields = split /\s+/, $record;
        die "\nABORTED: Mask size does not correspond to row record size\n" 
            if $#fields != $#mask;
        my $record_id;
        foreach my $i (0..@fields-1) {
            if ($mask[$i] eq '0') {
                next;
            } elsif ($mask[$i] eq 'I') {
                $record_id = $fields[$i];
            } elsif ($mask[$i] eq '1') {
                push @data_fields, $fields[$i];
            } else {
                die "misformed mask for reading the data file\n";
            }
        }
        my @nums = map {/$num_regex/;$_} @data_fields;
        $all_data{ $record_id } = \@nums;
    }

    my @all_data_ids = keys %all_data;

    if ($K == 0) {
        iterate_through_K( \@all_data_ids );
    } elsif ( $K =~ /\d+/) {
        cluster_for_fixed_K( \@all_data_ids, $K );
    } else {
        die "Incorrect call syntax used.  See documentation.\n";
    }
}

sub cluster_for_fixed_K {
    my @all_data_ids = @{shift @_};
    my $K = shift;
    my $N = @all_data_ids;
    die "You need at least 8 data samples. The number of data points " .
        "must satisfy the relation N = 2xK**2 where K is the " .
        "number of clusters.  The smallest value for K is 2.\n"
        if $N <= 8;
    my $Kmax = int( sqrt( $N / 2.0 ) );
    print "Value of Kmax is: $Kmax\n";
    my $QoC;
    my @QoC_values;
    my @array_of_clusters;
    my @array_of_cluster_centers;
    my $clusters;
    my $new_clusters;
    my $cluster_centers;
    my $new_cluster_centers;
    print "Clustering for K = $K\n";
    foreach my $trial (1..20) {
        my ($new_clusters, $new_cluster_centers) = 
                              cluster_for_given_K( $K, \@all_data_ids);
        my $newQoC = cluster_quality( $new_clusters, 
                                          $new_cluster_centers );
        if ( (!defined $QoC) || ($newQoC < $QoC) ) {
            $QoC = $newQoC;
            $clusters = deep_copy( $new_clusters );
            $cluster_centers = deep_copy( $new_cluster_centers );
        } 
    }
    if ($terminal_output) {
        print "\nDisplaying final clusters for best K (= $K) :\n";
        display_clusters( $clusters );
        display_cluster_centers( $cluster_centers );
        print "QoC value: $QoC\n";
#       print "QoC values array for different K: " . 
#               "@{[ map {my $x = sprintf('%.4f', $_); $x} @QoC_values ]}\n";
    }
    write_clusters_out_to_files( $clusters );
}

# The following subroutine is the top-level routine to call
# if you want the system to figure out on its own what value
# to use for K, the number of clusters.  It will try every
# possible value of K between 2 and the maximum possible
# depending on the number of data points available. For
# example, if the number of data points is 10,000, it will
# try all possible values of K between 2 and 70. For how the
# maximum value is set for K, see the comments made under
# Description.  Note also how this function makes 20
# different tries for each value of K as a defense against
# the problem of the final result corresponding to some
# local minimum in the values of the QoC metric.  Out of
# these 20 tries for each K, it retains the clusters and the
# cluster centers for only that try that yields the smallest
# value for the QoC metric.  After estimating the "best" QoC
# values for all possible K in this manner, it then finds
# the K for which the QoC is the minimum.  This is taken to
# be the best value for K.  Finally, the output clusters are
# written out to separate files.
sub iterate_through_K {
    my @all_data_ids = @{shift @_};
    my $N = @all_data_ids;
    die "You need at least 8 data samples. The number of data points " .
        "must satisfy the relation N = 2xK**2 where K is the " .
        "number of clusters.  The smallest value for K is 2.\n"
        if $N <= 8;
    my $Kmax = int( sqrt( $N / 2.0 ) );
    print "Value of Kmax is: $Kmax\n";
    my @QoC_values;
    my @array_of_clusters;
    my @array_of_cluster_centers;
    foreach my $K (2..$Kmax) {
        my $QoC;
        my $clusters;
        my $cluster_centers;
        print "Clustering for K = $K\n";
        foreach my $trial (1..20) {
            my ($new_clusters, $new_cluster_centers) = 
                              cluster_for_given_K( $K, \@all_data_ids);
            my $newQoC = cluster_quality( $new_clusters, 
                                          $new_cluster_centers );
            if ( (!defined $QoC) || ($newQoC < $QoC) ) {
                $QoC = $newQoC;
                $clusters = deep_copy( $new_clusters );
                $cluster_centers = deep_copy( $new_cluster_centers );
            } 
        }
        push @QoC_values, $QoC;
        push @array_of_clusters, $clusters;
        push @array_of_cluster_centers, $cluster_centers;
    }
    my ($min, $max) = minmax( \@QoC_values );
    die "Unsuccessful. Try again.\n" if ($max - $min ) < 0.00001;
    my $K_minus_2_best = get_index_at_value($min, \@QoC_values );
    my $K_best = $K_minus_2_best + 2;
    if ($terminal_output) {
        print "\nDisplaying final clusters for best K (= $K_best) :\n";
        display_clusters( $array_of_clusters[$K_minus_2_best] );
        display_cluster_centers( $array_of_cluster_centers[$K_minus_2_best] );
        print "\nBest clustering achieved for K=$K_best with QoC = $min\n";
        print "QoC values array for different K starting with K=2:  @QoC_values\n";
#       print "QoC values array for different K: " . 
#               "@{[ map {my $x = sprintf('%.4f', $_); $x} @QoC_values ]}\n";
    }
    write_clusters_out_to_files( $array_of_clusters[$K_minus_2_best] );
}

# The following function returns the value of QoC for a
# given partitioning of the data into K clusters.  It
# calculates two things: the average value for the distance
# between a data point and the center of the cluster in
# which the data point resides, and the average value for
# the distances between the cluster centers.  We obviously
# want to minimize the former and maximize the latter.  All
# of the "from center" distances within each cluster are
# stored in the variable $sum_of_distances_for_one_cluster.
# When this variable, after it is divided by the number of
# data elements in the cluster, is summed over all the
# clusters, we get a value that is stored in
# $avg_dist_for_cluster.  The inter-cluster-center distances
# are stored in the variable $inter_cluster_center_dist.
sub cluster_quality {
    my $clusters = shift;
    my $cluster_centers = shift;
    my $K = @$cluster_centers;          # Number of clusters
    my $cluster_radius = 0;
    foreach my $i (0..@$clusters-1) {
        my $sum_of_distances_for_one_cluster = 0;
        foreach my $ele (@{$clusters->[$i]}) {
            $sum_of_distances_for_one_cluster += 
                distance( $ele, $cluster_centers->[$i] );
        }
       $cluster_radius += 
           $sum_of_distances_for_one_cluster / @{$clusters->[$i]};
    }
    my $inter_cluster_center_dist = 0;
    foreach my $i (0..@$cluster_centers-1) {
        foreach my $j (0..@$cluster_centers-1) {
            $inter_cluster_center_dist += 
              distance2( $cluster_centers->[$i], 
                         $cluster_centers->[$j] );
        }
    }
    my $avg_inter_cluster_center_dist = $inter_cluster_center_dist /
                    ( $K * ($K-1) / 2.0 );
    return $cluster_radius / $avg_inter_cluster_center_dist;
}

# This is the function to call if you already know what
# value you want to use for K, the number of expected
# clusters.  The purpose of this function is to do the
# initialization of the cluster centers and to carry out the
# initial assignment of the data to the clusters with the
# initial cluster centers.  The initialization consists of 3
# steps: Construct a random sequence of K integers between 0
# and N-1 where N is the number of data points to be
# clustered; 2) Call get_initial_cluster_centers() to index
# into the data array with the random integers to get a list
# of K data points that would serve as the initial cluster
# centers; and (3) Call assign_data_to_clusters_initial() to
# assign the rest of the data to each of the K clusters on
# the basis of the proximity to the cluster centers.
sub cluster_for_given_K {
    my $K = shift;
    my @all_data_ids = @{shift @_};
    my @cluster_center_indices = 
                initialize_cluster_centers( $K, scalar(@all_data_ids) );
    my $cluster_centers = get_initial_cluster_centers( 
                                              \@cluster_center_indices, 
                                              \@all_data_ids );

    my $clusters = assign_data_to_clusters_initial( $cluster_centers, 
                                              \@all_data_ids );
    my $cluster_nonexistant_flag = 0;
    foreach my $trial (0..2) {
        ($clusters, $cluster_centers) = 
                              assign_data_to_clusters( $clusters, $K );
        my $num_of_clusters_returned = @$clusters;
        foreach my $cluster (@$clusters) {
            $cluster_nonexistant_flag = 1 if ((!defined $cluster) 
                                             ||  (@$cluster == 0));
        }
        last unless $cluster_nonexistant_flag;
    }
    return ($clusters, $cluster_centers);
}

# Returns a set of K random integers.  These serve as
# indices to reach into the data array.  A data element
# whose index is one of the random numbers returned by this
# routine serves as an initial cluster center.  Note the
# quality check it runs on the list of K random integers
# constructed.  We first make sure that all K random
# integers are different.  Subsequently, we carry out a
# quality assessment of the K random integers constructed.
# This quality measure consists of the ratio of the values
# spanned by the random integers to the value of N, the
# total number of data points to be clustered.  Currently,
# if this ratio is less than 0.3, we discard the K integers
# and try again.
sub initialize_cluster_centers {
    my $K = shift;
    my $data_store_size = shift;
    my @cluster_center_indices;
    while (1) {
        foreach my $i (0..$K-1) {
            $cluster_center_indices[$i] = int rand( $data_store_size );
            next if $i == 0;
            foreach my $j (0..$i-1) {
                while ( $cluster_center_indices[$j] == 
                        $cluster_center_indices[$i] ) {
                    my $old = $cluster_center_indices[$i];
                    $cluster_center_indices[$i] = int rand($data_store_size);
                }
            }
        }
        my ($min,$max) = minmax(\@cluster_center_indices );
        my $quality = ($max - $min) / $data_store_size;
        last if $quality > 0.3;
    }
    return @cluster_center_indices;
}

# This routine merely reaches into the data array with the
# random integers, as constructed by the previous routine,
# serving as indices and fetching values corresponding to
# those indices.  The fetched data samples serve as the
# initial cluster centers.
sub get_initial_cluster_centers {
    my @cluster_center_indices = @{shift @_};
    my @all_data_ids = @{shift @_};
    my @result;
    foreach my $i (@cluster_center_indices) {    
        push @result, $all_data{$all_data_ids[$i]};        
    }
    return \@result;
}

# The purpose of this routine is to form initial clusters by
# assigning the data samples to the initial clusters formed
# by the previous routine on the basis of the best proximity
# of the data samples to the different cluster centers.
sub assign_data_to_clusters_initial {
    my @cluster_centers = @{ shift @_ };
    my @all_data_ids =  @{ shift @_ };
    my @clusters;
    foreach my $ele (@all_data_ids) {
        my $best_cluster;
        my @dist_from_clust_centers;
        foreach my $center (@cluster_centers) {
            push @dist_from_clust_centers, distance($ele, $center);
        }
        my ($min, $best_center_index) = minimum( \@dist_from_clust_centers );
        push @{$clusters[$best_center_index]}, $ele;
    }
    return \@clusters;
}    

# This is the main routine that along with the
# update_cluster_centers() routine constitute the two key
# steps of the K-Means algorithm.  In most cases, the
# infinite while() loop will terminate automatically when
# the cluster assignments of the data points remain
# unchanged. For the sake of safety, we keep track of the
# number of iterations. If this number reaches 100, we exit
# the while() loop anyway.  In most cases, this limit will
# not be reached.
sub assign_data_to_clusters {
    my $clusters = shift;
    my $K = shift;
    my $final_cluster_centers;
    my $iteration_index = 0;
    while (1) {
        my $new_clusters;
        my $assignment_changed_flag = 0;
        my $current_cluster_center_index = 0;
        my $cluster_size_zero_condition = 0;
        my $how_many = @$clusters;
        my $cluster_centers = update_cluster_centers( $clusters );
        $iteration_index++;
        foreach my $cluster (@$clusters) {
            my $current_cluster_center = 
                          $cluster_centers->[$current_cluster_center_index];
            foreach my $ele (@$cluster) {
                my @dist_from_clust_centers;
                foreach my $center (@$cluster_centers) {
                    push @dist_from_clust_centers, distance($ele, $center);
                }
                my ($min, $best_center_index) = 
                              minimum( \@dist_from_clust_centers );
                my $best_cluster_center = 
                                 $cluster_centers->[$best_center_index];
                if (vector_equal($current_cluster_center, 
                                         $best_cluster_center)){
                    push @{$new_clusters->[$current_cluster_center_index]}, 
                                  $ele;
                } else {
                    $assignment_changed_flag = 1;             
                    push @{$new_clusters->[$best_center_index]}, $ele;
                }
            }
            $current_cluster_center_index++;
        }
        # Now make sure that we still have K clusters since K is fixed:
        next if ((@$new_clusters != @$clusters) && ($iteration_index < 100));
        # Now make sure that none of the K clusters is an empty cluster:
        foreach my $newcluster (@$new_clusters) {
            $cluster_size_zero_condition = 1 if ((!defined $newcluster) 
                                             or  (@$newcluster == 0));
        }
        next if (($cluster_size_zero_condition) && ($iteration_index < 100));
        last if $iteration_index == 100;
        # Now do a deep copy of new_clusters into clusters
	$clusters = deep_copy( $new_clusters );
        last if $assignment_changed_flag == 0;
    }
    $final_cluster_centers = update_cluster_centers( $clusters );
    return ($clusters, $final_cluster_centers);
}

# After each new assignment of the data points to the
# clusters on the basis of the current values for the
# cluster centers, we call the routine shown here for
# updating the values of the cluster centers.
sub update_cluster_centers {
    my @clusters = @{ shift @_ };
    my @new_cluster_centers;
    foreach my $cluster (@clusters) {
        die "Cluster became empty --- untenable condition " .
            "for a given K.  Try again. \n" if !defined $cluster;
        my $cluster_size = @$cluster;
        die "Cluster size is zero --- untenable.\n" if $cluster_size == 0;
        my @new_cluster_center = @{add_vectors( $cluster )};
        @new_cluster_center = map {my $x = $_/$cluster_size; $x} 
                                  @new_cluster_center;
        push @new_cluster_centers, \@new_cluster_center;
    }        
    return \@new_cluster_centers;
}

##################  Cluster Visualization Code #################

#  Note that it makes sense to call visualize() only AFTER
#  you have called kmeans() The latter function deposits the
#  symbolic labels of the data points in the different
#  clusters in different files that are named ClusterX.dat
#  for values of X starting with X=0.
#
#
#  IMPORTANT:  visualize() can only be run on the same 
#              source data file on which you last ran
#              kmeans()
#
#  The visualize() implementation automatically figures out
#  whether it should do a 2D plot or a 3D plot.  If the
#  number of on bits in the mask that is supplied as one of
#  the arguments is greater than 2, it does a 3D plot for
#  the first three data coordinates.  That is, the clusters
#  will be displayed in the 3D space formed by the first
#  three data coordinates. On the other hand, if the number
#  of on bits in the mask is exactly 2, it does a 2D plot.
#  Should it happen that only one on bit is specified for
#  the mask, visualize() aborts.
#
#  The visualization code consists of reading the original
#  source data file (which is presumed to have the same
#  structure as in the call to the clustering routine) and
#  creating the hash %all_data.  Subsequently, we read each
#  ClusterX.dat file that was created by the kmeans()
#  subroutine.  Note that the cluster files contain only the
#  symbolic names for the individual records in the source
#  data file.  We therefore next reach into the %all_data
#  hash table and get the data coordinates associated with
#  each symbolic label in the ClusterX.dat file.  The
#  numerical data thus generated is then written out to a
#  temp file.  When doing so we must remember to insert TWO
#  BLANK LINES between the data blocks corresponding to the
#  different ClusterX.dat files.  This constraint is imposed
#  on us by Gnuplot when plotting data from the same file.
#  Note that we want to use different point styles for the
#  data points in different cluster files.
#
#  Subsequently, we call upon the Perl interface provided by
#  the Graphics::GnuplotIF module to plot the data clusters.
sub visualize {
    my %args = @_;
    die "\n\n$0 Needs the name of the original data file with symbolic ids for all the data points that are subject to clustering and the mask string that has the same meaning as in the call to kmeans()" unless keys %args == 2; 

    use Graphics::GnuplotIF;

    my $master_datafile = $args{datafile};
    my $mask = $args{mask};

    my @mask = split //, $mask;

    open INPUT, $master_datafile 
         or die "unable to open file $master_datafile: $!\n";
    chomp( my @raw_data = <INPUT> );
    close INPUT;

    my %all_data;
    # Transform strings into number data
    foreach my $record (@raw_data) {
        next if $record =~ /^#/;
        my @data_fields;
        my @fields = split /\s+/, $record;
        die "\nABORTED: Mask size does not correspond to row record size\n" 
            if $#fields != $#mask;
        my $record_id;
        foreach my $i (0..@fields-1) {
            if ($mask[$i] eq '0') {
                next;
            } elsif ($mask[$i] eq 'I') {
                $record_id = $fields[$i];
            } elsif ($mask[$i] eq '1') {
                push @data_fields, $fields[$i];
            } else {
                die "misformed mask for reading the data file\n";
            }
        }
        my @nums = map {/[\d.]/;$_} @data_fields;
        $all_data{ $record_id } = \@nums;
    }

    #count the number of data fields to plot in each record
    my $data_field_width = 0;
    foreach my $c (@mask) {
        $data_field_width++ if $c eq '1';
    }

    my @all_data_ids = keys %all_data;

    my @cluster_files = glob "Cluster*.dat";
    #print "cluster files: @cluster_files\n";
    my $K = @cluster_files;

    my $temp_file = "temp_" . $master_datafile;
    unlink $temp_file if -e $temp_file;

    open OUTPUT, ">$temp_file"
           or die "Unable to open a temp file in this directory: $!\n";
    foreach my $file (@cluster_files) {
        open INPUT, $file or die "Cluster file $file disappeared: $!\n";
        while (<INPUT>) {
            chomp;
            print OUTPUT "@{$all_data{$_}}";
            print OUTPUT "\n";
        }
        close INPUT;
        print OUTPUT "\n\n";
    }
    close OUTPUT;

    my $plot = Graphics::GnuplotIF->new( persist => 1 );

    $plot->gnuplot_cmd( "set noclip" );
    $plot->gnuplot_cmd( "set pointsize 2" );

    my $arg_string = "";
    if ($data_field_width > 2) {
        foreach my $i (0..$K-1) {
            my $j = $i + 1;
            $arg_string .= "\"$temp_file\" index $i using 1:2:3 notitle with points lt $j pt $j, ";
        }
    } elsif ($data_field_width == 2) {
        foreach my $i (0..$K-1) {
            my $j = $i + 1;
            $arg_string .= "\"$temp_file\" index $i using 1:2 notitle with points lt $j pt $j, ";
        }
    } elsif ($data_field_width == 1 ) {
        foreach my $i (0..$K-1) {
            my $j = $i + 1;
            $arg_string .= "\"$temp_file\" index $i using 1 notitle with points lt $j pt $j, ";
        }
    }

    $arg_string = $arg_string =~ /^(.*),[ ]+$/;
    $arg_string = $1;

    if ($data_field_width > 2) {
        $plot->gnuplot_cmd( "splot $arg_string" );
    } elsif ($data_field_width == 2) {
        $plot->gnuplot_cmd( "plot $arg_string" );
    } elsif ($data_field_width == 1) {
        die "No provision for plotting 1-D data\n";
    }
}


###########  Generating Synthetic Data for Clustering  ############

#  The data generated corresponds to a multivariate
#  distribution.  The mean and the covariance of each
#  Gaussian in the distribution are specified individually
#  in a parameter file.  See the example parameter file
#  param.txt in the examples directory.  Just edit this
#  file for your own needs.
#
#  The multivariate random numbers are generated by calling
#  the Math::Random module.  As you would expect, that
#  module will insist that the covariance matrix you
#  specify be symmetric and positive definite.
sub cluster_data_generator {
    my %args = @_;

    my $input_parameter_file = $args{input_parameter_file};
    my $output_file = $args{output_datafile};
    my $N = $args{number_data_points_per_cluster};

    open INPUT, $input_parameter_file
        || "unable to open parameter file: $!";

    my @all_params = <INPUT>;
    @all_params = grep { $_ !~ /^[ ]*#/ } @all_params;
    chomp @all_params;
    my $param_string = join ' ', @all_params;
    my @cluster_strings = split /[ ]*cluster[ ]*/, $param_string;
    @cluster_strings = grep  $_, @cluster_strings;

    my $K = @cluster_strings;
    die "Too many clusters requested" if $K > 12;
    my @point_labels = ('a'..'z');

    print "Number of Gaussians used for the synthetic data: $K\n";
    my @means;
    my @covariances;
    my $data_dimension;
    foreach my $i (0..$K-1) {
        my @num_strings = split /  /, $cluster_strings[$i];
        my @cluster_mean = map {/$num_regex/;$_} split / /, $num_strings[0];
        $data_dimension = @cluster_mean;
        push @means, \@cluster_mean;
        my @covariance_nums = map {/$num_regex/;$_} split / /, $num_strings[1];
        die "dimensionality error" if @covariance_nums != 
                                      ($data_dimension ** 2);
        
        my $cluster_covariance;
        foreach my $j (0..$data_dimension-1) {
            foreach my $k (0..$data_dimension-1) {        
                $cluster_covariance->[$j]->[$k] = 
                         $covariance_nums[$j*$data_dimension + $k];
            }
        }
        push @covariances, $cluster_covariance;
    }

    use Math::Random;           # for normal and uniform densities

    random_seed_from_phrase( 'hellojello' );

    my @data_dump;
    foreach my $i (0..$K-1) {
        my @m = @{shift @means};
        my @covar = @{shift @covariances};
        my @new_data = Math::Random::random_multivariate_normal( $N, @m, @covar );
        my $p = 0;
        my $label = $point_labels[$i];
        @new_data = map {unshift @$_, $label.$i; $i++; $_} @new_data;
        push @data_dump, @new_data;     
    }

    fisher_yates_shuffle( \@data_dump );

    open OUTPUT, ">$output_file";
    foreach my $ele (@data_dump) {
        foreach my $coord ( @$ele ) {
            print OUTPUT "$coord ";
        }
        print OUTPUT "\n";
    }
    close OUTPUT;
}

######################   Support Routines  ########################

# This routine is really not necessary in light of the new
# `~~' operator in Perl.  Will use the new operator in the
# next version.
sub vector_equal {
    my $vec1 = shift;
    my $vec2 = shift;
    die "wrong data types for distance calculation\n" if @$vec1 != @$vec2;
    foreach my $i (0..@$vec1-1){
        return 0 if $vec1->[$i] != $vec2->[$i];
    }
    return 1;
}

sub add_vectors {
    my @arr_of_ids = @{shift @_};      # array of data element names
    my @result;
    my $data_dimensionality = @{$all_data{$arr_of_ids[0]}};
    foreach my $i (0..$data_dimensionality-1) {
        $result[$i] = 0.0;
    }
    foreach my $id (@arr_of_ids) {
        my $ele = $all_data{$id};
        my $i = 0;
        foreach my $component (@$ele) {
            $result[$i] += $component;
            $i++;
        }
    }
    return \@result;
}

sub display_cluster_centers {
    my @cluster_center_arr = @{shift @_};
    my $i = 1;
    foreach my $ele (@cluster_center_arr) {
        print "Cluster center $i: " .
               "@{[ map {my $x = sprintf('%.4f', $_); $x} @$ele ]}\n";
        $i++;
    }
}

# The following routine is for computing the distance
# between a data point specified by its symbolic name in the
# master datafile and a point (such as the center of a
# cluster) expressed as a vector of coordinates:
sub distance {
    my $ele1_id = shift @_;            # symbolic name of data sample
    my @ele1 = @{$all_data{$ele1_id}};
    my @ele2 = @{shift @_};
    die "wrong data types for distance calculation\n" if @ele1 != @ele2;
    my $how_many = @ele1;
    my $squared_sum = 0;
    foreach my $i (0..$how_many-1) {
        $squared_sum += ($ele1[$i] - $ele2[$i])**2;
    }    
    my $dist = sqrt $squared_sum;
    return $dist;
}

# The following routine does the same as above but now both
# arguments are expected to be arrays of numbers:
sub distance2 {
    my @ele1 = @{shift @_};
    my @ele2 = @{shift @_};
    die "wrong data types for distance calculation\n" if @ele1 != @ele2;
    my $how_many = @ele1;
    my $squared_sum = 0;
    foreach my $i (0..$how_many-1) {
        $squared_sum += ($ele1[$i] - $ele2[$i])**2;
    }    
    return sqrt $squared_sum;
}
sub get_index_at_value {
    my $value = shift;
    my @array = @{shift @_};
    foreach my $i (0..@array-1) {
        return $i if $value == $array[$i];
    }
}

# For displaying the individual clusters on a terminal
# screen.  Each cluster is displayed through the symbolic
# names associated with the data points.
sub display_clusters {
    my @clusters = @{shift @_};
    my $i = 1;
    foreach my $cluster (@clusters) {
        @$cluster = sort @$cluster;
        my $cluster_size = @$cluster;
        print "\n\nCluster $i ($cluster_size records):\n";
        foreach my $ele (@$cluster) {
            print "  $ele";
        }
        $i++
    }
    print "\n\n";
}

sub write_clusters_out_to_files {
    my @clusters = @{shift @_};        
    unlink glob "Cluster*.dat";
    foreach my $i (1..@clusters) {
        my $filename = "Cluster" . $i . ".dat";
        print "Writing cluster $i to file $filename\n";
        open FILEHANDLE, "| sort > $filename"
            or die "Unable to pen file: $!";
        foreach my $ele (@{$clusters[$i-1]}) {        
            print FILEHANDLE "$ele\n";
        }
        close FILEHANDLE;
    }
}

# Meant only for constructing a deep copy of an array of arrays
sub deep_copy {
    my $ref_in = shift;
    my $ref_out;
    foreach my $i (0..@{$ref_in}-1) {
        foreach my $j (0..@{$ref_in->[$i]}-1) {
            $ref_out->[$i]->[$j] = $ref_in->[$i]->[$j];
        }
    }
    return $ref_out;
}

# Returns the minimum value and its positional index in an array
sub minimum {
    my $arr = shift;
    my $min;
    my $index;
    foreach my $i (0..@{$arr}-1) {
        if ( (!defined $min) || ($arr->[$i] < $min) ) {
            $index = $i;
            $min = $arr->[$i];
        }
    }
    return ($min, $index);
}

sub minmax {
    my $arr = shift;
    my $min;
    my $max;
    foreach my $i (0..@{$arr}-1) {
        if ( (!defined $min) && (!defined $max) ) {
            $min = $arr->[$i];
            $max = $arr->[$i];
        } elsif ( $arr->[$i] < $min ) {
            $min = $arr->[$i];
        } elsif ( $arr->[$i] > $max ) {
            $max = $arr->[$i];
        }
    }
    return ($min, $max);
}

# from perl docs:
sub fisher_yates_shuffle {                
    my $arr =  shift;                
    my $i = @$arr;                   
    while (--$i) {                   
        my $j = int rand( $i + 1 );  
        @$arr[$i, $j] = @$arr[$j, $i]; 
    }
}

1;

__END__

=head1 NAME

Algorithm::KMeans - Clustering multi-dimensional data with a pure-Perl implementation

=head1 SYNOPSIS

  use Algorithm::KMeans qw(kmeans visualize cluster_data_generator);

  # Set the mask to tell system which columns of a datafile to use
  # for clustering and which column contains a symbolic ID for each
  # data record.  For example, if the ID is in the first column, if
  # you want the second column to be ignored, and for 3D clustering
  # from the rest:
  my $mask = "I0111";

  # If you want the module to figure out the optimum number of clusters
  # from the data in the file supplied as $datafile:
  kmeans( datafile => $datafile,
          mask     =>  $mask,
          terminal_output => 1,
          K => 0 );

  # If you know how many clusters you want (in this case 3):
  kmeans( datafile => $datafile,
          mask     =>  $mask,
          terminal_output => 1,
          K => 3 );

  # To view the clusters formed:
  visualize( datafile =>  $datafile,
             mask     =>  $mask );

  # If you want to generate your own multivariate data for clustering,
  # you can call
  my $N = 60;              # If you want 60 data points per cluster
  cluster_data_generator( input_parameter_file => $parameter_file,
                          output_datafile => $datafile,
                          number_data_points_per_cluster => $N );

=head1 DESCRIPTION

B<Algorithm::KMeans> is a I<perl5> module for the clustering
of numerical data in multidimensional spaces.  Since the
module is entirely in Perl (in the sense that it is not a
Perl wrapper around a C library that actually does the
clustering), the code in the module can easily be modified
to experiment with several aspects of automatic clustering.
For example, one can change the criterion used to measure
the "distance" between two data points, the stopping
condition for accepting final clusters, the criterion used
for measuring the quality of the clustering achieved, etc.

A K-Means clusterer is a poor man's implementation of the EM
algorithm.  EM stands for Expectation Maximization. For the
case of Gaussian data, the results obtained with a good
K-Means implementation should match those obtained with the
EM algorithm. Clustering with K-Means takes place
iteratively and involves two steps: 1) assignment of data
samples to clusters; and 2) Recalculation of the cluster
centers.  The assignment step can be shown to be akin to the
Expectation step of the EM algorithm, and the calculation of
the cluster centers akin to the Maximization step of the EM
algorithm.

Of the two key steps of the K-Means algorithm, the
assignment step consists of assigning each data point to
that cluster from whose center the data point is the
closest.  That is, during assignment, you compute the
distance between the data point and each of the current
cluster centers.  You assign the data sample on the basis of
the minimum value of the computed distance.  The second step
consists of re-computing the cluster centers for the newly
modified clusters.

Obviously, before the two-step approach can proceed, we need
to initialize the both the cluster center values and the
clusters that can then be iteratively modified by the
two-step algorithm.  How this initialization is carried out
is very important.  The implementation here uses a random
number generator to find K random integers between 0 and N
where N is the total number of data samples that need to be
clustered and K the number of clusters you wish to form.
The K random integers are used as indices for the data
samples in the overall data array --- the data samples thus
selected are treated as seed cluster centers.  This
obviously requires a prior knowledge of K.

How to specify K is one of the most vexing issues in any
approach to clustering.  In some case, we can set K on the
basis of prior knowledge.  But, more often than not, no such
prior knowledge is available.  When the programmer does not
explicitly specify a value for K, the approach taken in the
current implementation is to try all possible values between
2 and some largest possible value that makes statistical
sense.  We then choose that value for K which yields the
best value for the QoC (Quality of Clustering) metric.  It
is generally believed that the largest value for K should
not exceed sqrt(N/2) where N is the number of data point to
be clustered.

How to set the QoC metric is obviously a critical issue unto
itself.  In the current implementation, the value of QoC is
a ratio of the average radius of the clusters and the
average distance between the cluster centers.  But note that
this is a good criterion only when the data exhibits the
same variance in all directions.  When the data variance is
different directions, but still remains the same for all
clusters, a more appropriate QoC can be formulated using
other distance metrics such as the Mahalanobis distance.

Every iterative algorithm requires a stopping criterion.
The criterion implemented here is that we stop iterations
when there is no re-assignment of the data points during the
assignment step.

Ordinarily, the output produced by a K-Means clusterer will
correspond to a local minimum for the QoC values, as opposed
to a global minimum.  The current implementation protects
against that, but only in a very small way, by trying
different randomly selected initial cluster centers and then
selecting the one that gives the best overall QoC value.

=head1 METHODS

The module provides the following methods for clustering,
for cluster visualization, and for the generation of data
for testing a clustering algorithm:

=over

=item   Algorithm::KMeans::kmeans( datafile => $data_file,
                                   mask     =>  $mask,
                                   terminal_output => 1,
                                   K => 0 );

where the keyword argument "K=>0" tells the module that you
want it to figure out the optimum number of clusters to
form.  The datafile keyword names the file that contains the
data that needs to be clustered.  The data file is expected
to contain entries in the following format

   c20  0  10.7087017086940  9.63528386251712  10.9512155258108
   c7   0  12.8025925026787  10.6126270065785  10.5228482095349
   b9   0  7.60118206283120  5.05889245193079  5.82841781759102
   ....
   ....

where the first column contains the symbolic ID tag for each
data record and the rest of the columns the numerical
information.  As to which columns are actually used for
clustering is decided by the string value of the mask.  For
example, if we wanted to cluster on the basis of the entries
in just the last three columns, the mask value would be
`I0111' where the character `I' indicates that the ID tag is
in the first column, the character '0' that the second
column is to be ignored, and '1's that the last three
columns are to be used for clustering.

The parameter `terminal_value' is boolean and determines
what you will see on the terminal screen of the window in
which you make these method calls.  If you set it to 1, you
will see different clusters as lists of the symbolic IDs and
you will also see the cluster centers, along with the QoC
(Quality of Clustering) value for the clusters displayed.
If this parameter is set to 0, you will see only minimal
information.  In either case, the clusters will be written
out to files named Cluster0.dat, Cluster1.dat, Cluster2.dat,
etc.  Before the clusters are written to these files, the
module destroys all files with such names in the directory
in which you call the module.

=item   Algorithm::KMeans::kmeans( datafile => $data_file,
                                   mask     =>  $mask,
                                   terminal_output => 1,
                                   K => $K );

for a non-zero integer value for the keyword `K'.  This has
a profound effect of the behavior of the module, as it will
now form exactly $K clusters.  


=item  Algorithm::KMeans::visualize( datafile => $datafile,
                                     mask     =>  $mask );

for visualizing the clusters formed.  The datafile is the
same as used in the call to kmeans().  This datafile is used
by the visualization script to extract the numerical
coordinates associated with the symbolic ID tags in the
cluster files.  The mask here does not have to be identical
to the one used for clustering, but must be a subset of that
mask.  This is convenient for visualizing the clusters formed
in two- or three-dimensional subspaces of the original
space, assuming that the clustering was carried out in a
space whose dimensionality is greater than 3.

=item Algorithm::KMeans::cluster_data_generator( 
                          input_parameter_file => $parameter_file_name,
                          output_datafile =>  $datafile,
                          number_data_points_per_cluster => $N );

for generating multivariate data for clustering if you wish
to play with synthetic data for clustering.  The input
parameter file contains the means and the variances for the
different Gaussians you wish to use for the synthetic data.
See the file param.txt provided in the examples directory.
It will be easiest for you to just edit this file for your
data generation needs.  In addition to the format of the
parameter file, the main constraint you need to observe in
specifying the parameters is that the dimensionality of the
covariance matrix must correspond to the dimensionality of
the mean vectors.  The multivariate random numbers are
generated by calling the Math::Random module.  As you would
expect, this module requires that the covariance matrices
you specify in your parameter file be symmetric and positive
definite.  Should the covariances in your parameter file not
obey this condition, the Math::Random module will let you
know.

=back

=head1 HOW ARE THE CLUSTERS OUTPUT?

The module dumps the clusters into files named

    Cluster0.dat
    Cluster1.dat
    Cluster2.dat
    ...
    ...

in the directory in which you execute the module.  The
number of such files will equal the number of clusters
formed.  On each run of the kmeans() function, all such
existing files in the directory are destroyed before fresh
ones created.  Each cluster file consists of the symbolic ID
tags of the data points in that cluster.

It would be trivial to modify the module so that a call to
kmeans() returns a list of hashes, each hash representing a
different cluster, and with each hash holding the symbolic
ID tags for the keys and their data point coordinates for
the values.  If there is demand for such a feature, it will
be added to the next version of this module.

=head1 REQUIRED

This module requires the following two modules:

   Math::Random
   Graphics::GnuplotIF

the former for generating the multivariate random numbers
and the latter for the visualization of the clusters.

=head1 EXAMPLES

See the examples directory in the distribution for how to
make calls to the clustering and the visualization routines.
The examples directory also includes a parameter file,
param.txt, for generating synthetic data for clustering.
Just edit this file if you would like to generate your own
multivariate data for clustering.


=head1 CAVEATS

Please note that this clustering module is not meant for
very large datafiles.  Being an all-Perl implementation, the
goal here is not the speed of execution.  On the contrary,
the goal is to make it easy to experiment with the different
facets of K-Means clustering.  If you need to process a
large data file, you'd be better off with a module like
Algorithm::Cluster.  But note that when you use a wrapper
module in which it is a C library that is actually doing the
job of clustering for you, it is more daunting to experiment
with various aspects of clustering.

Clustering usually does not work well when the data is
highly anisotropic, that is, when the data has very
different variances along its different dimensions.  This
problem becomes particularly severe when the different
clusters you expect to see in the data have non-uniform
anisotropies.  When the anisotropies are uniform, one could
improve on the implementation shown in the current module by
first normalizing the data with respect to the different
variances along the different dimensions.  One could also
try to cluster the data in a low-dimensional space formed by
a principal components analysis of the data, especially
after the variances are normalized out.  Depending on how
the current module is received, its future versions may
include those enhancements.


=head1 BUGS

Please send me email if you find any bugs.  They will
be dealt with promptly.

=head1 INSTALLATION

The usual

    perl Makefile.PL
    make
    make test
    make install

if you have root access.  If not, 

    perl Makefile.PL prefix=/some/other/directory/
    make
    make test
    make install

=head1 AUTHOR

Avinash Kak, kak@purdue.edu

If you send email, please place the string "KMeans" in your
subject line to get past my spam filter.

=head1 COPYRIGHT

This library is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

 Copyright 2010 Avinash Kak

=cut

