#!/usr/bin/perl

# ./edge_to_infomap.pl orig_graph.txt vertices.txt sigma thresh

$min_dist = 0.06; #0

if (defined($ARGV[3]) && ($ARGV[3] != 0)) {
    # threshold is upper limit!
    $thresh = $ARGV[3];
} else {
    # just set at infinity, so everything passes
    $thresh = 9**9**9; 
}

if (defined $ARGV[2]) {
    $sigma = $ARGV[2];
} else {
    $sigma = 1;
}

if (defined $ARGV[1]) {
    $vertices = $ARGV[1];
} else {
    $vertices = "/home/st23955/code/der_eval/REF/trn_ev08+10_4w.A.lst";
}

# hash vertices to line number...
$counter = 1;
$vtx_hash = {};
open(VTX, "$vertices") || die("Could not open file $vertices");
while(defined($line = <VTX>)) {
    chomp($line);
    $vtx_hash->{$line} = $counter;
    $counter++;
}
close(VTX);

# get and print out number of vertices
$num_vertices = $counter - 1;
printf(STDOUT "*Vertices %s\n", $num_vertices);

# print out vertex names -- useless, but a necessary format formality
for ($count = 0; $count < $num_vertices; $count++) {
    printf(STDOUT "%s \"%s\"\n", $count+1, $count+1);
}

# get and then print number of edges...
$num_edges = 0;
open(GTXT, "$ARGV[0]") || die("Could not open file $ARGV[0]");
while(defined($line = <GTXT>)) {
    chomp($line);

    # split to get vertex labels
    ($vertex1, $vertex2, $edge) = split(",", $line);
    
    # split edge description to just get weight 
    ($nil, $weight) = split(" ", $edge);
    chop($weight); # remove closing braces

    # if both vertices exist as keys in our hash
    # AND edge distance is less than threshold
    if ( exists($vtx_hash->{$vertex1}) && exists($vtx_hash->{$vertex2}) && ($weight < $thresh) ) {
        $num_edges++;
    }
}
close(GTXT);
printf(STDOUT "*Edges %s\n", $num_edges);

# go through graph AGAIN and print edges + weights
open(GTXT, "$ARGV[0]") || die("Could not open file $ARGV[0]");
while(defined($line = <GTXT>)) {
    chomp($line);

    # split to get vertex labels
    ($vertex1, $vertex2, $edge) = split(",", $line);
    
    # split edge description to just get weight 
    ($nil, $weight) = split(" ", $edge);
    chop($weight); # remove closing braces

    # if both vertices exist as keys in our hash
    # AND edge distance is less than threshold
    if ( exists($vtx_hash->{$vertex1}) && exists($vtx_hash->{$vertex2}) && ($weight < $thresh) ) {

        if ($sigma == 0) {
            printf(STDOUT "%s %s 1\n", $vtx_hash->{$vertex1}, $vtx_hash->{$vertex2});
        } else {
            printf(STDOUT "%s %s %3.5f\n", $vtx_hash->{$vertex1}, $vtx_hash->{$vertex2}, exp(-1 * ($weight - $min_dist) / $sigma));
        }
    }
}
close(GTXT);
