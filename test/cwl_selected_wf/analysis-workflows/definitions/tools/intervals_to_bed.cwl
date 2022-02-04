#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
requirements:
- class: DockerRequirement
  dockerPull: "ubuntu:bionic"
- class: ResourceRequirement
  ramMin: 4000
- class: InitialWorkDirRequirement
  listing:
  - entryname: 'intervals_to_bed.pl'
    entry: |
      use feature qw(say);

      for my $line (<>) {
          chomp $line;

          next if substr($line,0,1) eq '@'; #skip header lines

          my ($chrom, $start, $stop) = split(/\t/, $line);
          say(join("\t", $chrom, $start-1, $stop));
      }
hints:
  NetworkAccess:
    networkAccess: true
  LoadListingRequirement:
    loadListing: deep_listing
inputs:
  interval_list:
    type: File
    inputBinding:
      position: 1
baseCommand: ['/usr/bin/perl', 'intervals_to_bed.pl']
stdout: "interval_list.bed"
outputs:
  interval_bed:
    type: stdout
