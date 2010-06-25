#!/usr/bin/perl -w

#  There should be no need to call this script ordinarily.

#  When the Algorithm::KMeans module creates new cluster files,
#  it automatically delete all previously created such files.
#  Such files are named ClusterX.dat for X starting with X = 0.
#  The files temp_* are created by the visualization script.


unlink glob "Cluster*.dat";

unlink glob "temp_*";

