use Test::Simple tests => 3;

use lib '../blib/lib','../blib/arch';

use Algorithm::KMeans qw(kmeans visualize cluster_data_generator);;

my $parameter_file = "param.txt";
my $datafile = "testdata.dat";
my $mask = "I111";



# Test 1 (Data Generation):
cluster_data_generator( input_parameter_file => $parameter_file,
                        output_datafile => $datafile,
                        number_data_points_per_cluster => 30 );
open IN, $datafile;
my @data_records = <IN>;
ok( @data_records == 90,  'Data generation works' );


# Test 2 (K-Means Clustering):
kmeans( datafile => $datafile,
        mask     =>  $mask,
        terminal_output => 0,  
        K => 3 );
my @clusters = glob "Cluster*.dat";
ok( @clusters == 3,  'Clustering works' );


# Test 3 (Data Visualization)
eval {
    visualize( datafile => $datafile,
               mask     =>  $mask );
};
ok( !$@,  'Visualization works' );

unlink glob "*.dat";
unlink glob "Cluster*";
unlink glob "temp_*";
