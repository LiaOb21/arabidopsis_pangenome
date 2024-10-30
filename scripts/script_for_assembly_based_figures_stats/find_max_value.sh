#!/bin/bash

# Input file with the collapsed data
collapsed_file="pav_collapsed_output.txt"

# Use awk to find the maximum values for each class and track the corresponding GCA
awk '
NR > 1 {  # Skip the header row
  # Check and update max for "core"
  if ($2 != "NA" && $2 > max_core) {
    max_core = $2
    max_core_gca = $1
  }
  
  # Check and update max for "softcore"
  if ($3 != "NA" && $3 > max_softcore) {
    max_softcore = $3
    max_softcore_gca = $1
  }
  
  # Check and update max for "dispensable"
  if ($4 != "NA" && $4 > max_dispensable) {
    max_dispensable = $4
    max_dispensable_gca = $1
  }
  
  # Check and update max for "private"
  if ($5 != "NA" && $5 > max_private) {
    max_private = $5
    max_private_gca = $1
  }
}
END {
  # Print the results
  print "Max core: " max_core " (GCA: " max_core_gca ")"
  print "Max softcore: " max_softcore " (GCA: " max_softcore_gca ")"
  print "Max dispensable: " max_dispensable " (GCA: " max_dispensable_gca ")"
  print "Max private: " max_private " (GCA: " max_private_gca ")"
}
' "$collapsed_file"
