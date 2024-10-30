#!/bin/bash

# Input file with the collapsed data
collapsed_file="pav_collapsed_output.txt"

# Use awk to find the minimum values for each class and track the corresponding GCA
awk '
NR > 1 {  # Skip the header row
  # Initialize the min values only during the first line of data processing
  if (NR == 2) {
    min_core = $2 != "NA" ? $2 : "NA"
    min_core_gca = $1
    min_softcore = $3 != "NA" ? $3 : "NA"
    min_softcore_gca = $1
    min_dispensable = $4 != "NA" ? $4 : "NA"
    min_dispensable_gca = $1
    min_private = $5 != "NA" ? $5 : "NA"
    min_private_gca = $1
  }
  
  # Check and update min for "core"
  if ($2 != "NA" && ($2 < min_core || min_core == "NA")) {
    min_core = $2
    min_core_gca = $1
  }
  
  # Check and update min for "softcore"
  if ($3 != "NA" && ($3 < min_softcore || min_softcore == "NA")) {
    min_softcore = $3
    min_softcore_gca = $1
  }
  
  # Check and update min for "dispensable"
  if ($4 != "NA" && ($4 < min_dispensable || min_dispensable == "NA")) {
    min_dispensable = $4
    min_dispensable_gca = $1
  }
  
  # Check and update min for "private"
  if ($5 != "NA" && ($5 < min_private || min_private == "NA")) {
    min_private = $5
    min_private_gca = $1
  }
}
END {
  # Print the results
  print "Min core: " min_core " (GCA: " min_core_gca ")"
  print "Min softcore: " min_softcore " (GCA: " min_softcore_gca ")"
  print "Min dispensable: " min_dispensable " (GCA: " min_dispensable_gca ")"
  print "Min private: " min_private " (GCA: " min_private_gca ")"
}
' "$collapsed_file"
