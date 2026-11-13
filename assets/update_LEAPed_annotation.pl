#!/usr/bin/env perl

=head1 NAME

  update_LEAPed_annotation.pl script

=head1 SYNOPSIS

  perl update_LEAPed_annotation.pl -in <FILE> -dataset <DSNAME> -out <OUTFILE> <OPTIONS>

=head1 DESCRIPTION

  The input for this script is the output of the LEAP pipeline, a gff3 file with the updated transcript coordinates.
  It extracts transcript starts and ends, and LEAP tags, and checks that the transcript coordinate changes are consistent
  with the LEAP tags.

  If the updated transcript is not a MANE transcript, the script will modify the transcript coordinates accordingly
  and assign it the corresponding LEAP remarks in the database. For MANE transcripts,
  it first looks for an Ensembl_canonical_extended model with matching protected exon
  boundaries and CDS. It reuses that model, keeping the farther boundary at each end,
  or makes a copy if no matching model exists. The original MANE model is preserved.
  Multiple matching existing models are an error requiring review.

  The output is a set of three files:
  - <OUTFILE>_input_check.txt
  - <OUTFILE>_edits.log
  - <OUTFILE>_postedit_check.txt

=head1 OPTIONS

  -dataset          Species name in Otter as configured in the server_config/species.dat file
  -in               Input GFF file
  -out              Output file
  -chr              Only process data on this chromosome
  -dry_run          Do sanity checks and simulate the annotation update
  -write            Do sanity checks and actually update the annotation in the database

=head1 CONTACT

  gencode-help@ebi.ac.uk

=cut


use strict;
use warnings;
use Getopt::Long;
use Bio::Otter::Lace::Defaults;
use Bio::Otter::Server::Config;
use Bio::Vega::Exon;
use Sys::Hostname;
use Try::Tiny;
use Bio::Otter::ServerAction::Script::Region;
use Bio::Otter::Server::Support::Local;
use List::Util qw[min max];
$| = 1;

my $infile;
my $dataset_name;
my $author_name;
my $only_chr;
my $outfile;
my $dry_run;
my $write;

&GetOptions(
            'in=s'            => \$infile,
            'dataset=s'       => \$dataset_name,
            'author=s'        => \$author_name,
            'chr=s'           => \$only_chr,
            'out=s'           => \$outfile,
            'dry_run!'        => \$dry_run,
            'write!'          => \$write,
            );

unless ($infile and $dataset_name){
  die "You must provide an input file and a dataset name!\n";
}

if (!$dry_run and !$write){
  print "\nNeither -dry_run nor -write flags have been set.\nThis will be a dry run.\n\n";
}
if ($dry_run and $write){
  print "\nUsing both -dry_run and -write flags at the same time.\nThis will be a dry run.\n\n";
  $write = 0;
}

#Default output file name
unless ($outfile){
  my ($s, $m, $h, $d, $mo, $y) = localtime();
  $outfile = join("_", $y+1900, map {sprintf("%02d", $_)} ($mo, $d, $h, $m, $s));
}

#Create attribute objects
my $leap_5_att = Bio::EnsEMBL::Attribute->new(-code => 'remark', -value => 'LEAP_5');
my $leap_3_att = Bio::EnsEMBL::Attribute->new(-code => 'remark', -value => 'LEAP_3');
my $canon_ext_att = Bio::EnsEMBL::Attribute->new(-code => 'remark', -value => 'Ensembl_canonical_extended');


#Read input file and store transcript coordinates, exon coordinates and LEAPed transcript ids in hashes
print "Reading input GFF3 file... ";

my %transcript_coords;
my %exon_coords;
my %leap_tags;
if ($infile =~ /\.gz$/){
  open (IN, "zcat $infile | ") or die "Can't open input file $infile: $!";
}
else{
  open (IN, "cat $infile | ") or die "Can't open input file $infile: $!";
}
while (<IN>){
  next if /^#/;
  chomp;
  my @cols = split(/\t/);
  my $chr = $cols[0];
  $chr =~ s/^chr//;
  if ($only_chr){
    next unless $chr eq $only_chr;
  }
  if ($chr eq "Y"){
    next if $_ =~ /_PAR_Y/; #skip Y PAR annotation in the input file as the havana database doesn't store this annotation
  }
  my @atts = split(/;/, $cols[8]);
  my ($start, $end, $strand) = @cols[3,4,6];
  my $transcript_id;
  my $exon_number;
  foreach my $att (@atts){
    my ($field, $value) = split(/=/, $att);
    if ($field eq "transcript_id"){
      $transcript_id = $value;
      $transcript_id =~ s/\.\d+//; #remove version
    }
    elsif ($field eq "leap_tag"){
      for (split(/,/, $value)){
        if ($_ =~ /^LEAP_(5|3)$/){
          $leap_tags{$_}{$transcript_id} = 1;
        }
      }
    }
    elsif ($field eq "exon_number"){
      $exon_number = $value;
    }
  }
  next if $transcript_id =~ /^OTT/; #skip OTT transcripts
  if ($cols[2] eq "transcript"){
    $transcript_coords{$transcript_id}{'chr'} = $chr;
    $transcript_coords{$transcript_id}{'start'} = $start;
    $transcript_coords{$transcript_id}{'end'} = $end;
    $transcript_coords{$transcript_id}{'strand'} = $strand;
  }
  elsif ($cols[2] eq "exon"){
    $exon_coords{$transcript_id}{$exon_number}{'start'} = $start;
    $exon_coords{$transcript_id}{$exon_number}{'end'} = $end;
  }
}
close (IN);

print "DONE\n";


#Connect to the havana database
my $dataset = Bio::Otter::Server::Config->SpeciesDat->dataset($dataset_name);
my $otter_dba = $dataset->otter_dba or die "can't get db adaptor\n";
my $ga = $otter_dba->get_GeneAdaptor();
my $ta = $otter_dba->get_TranscriptAdaptor();
my $csa = $otter_dba->get_CoordSystemAdaptor();
my $assembly_version = $csa->fetch_by_rank(1)->version;

#Check that the coordinates updates are consistent
print "Checking consistency of coordinate updates... ";

open (CHK, ">".$outfile."_input_check.txt") or die "Can't open file ".$outfile."_input_check.txt:$!";
print CHK join("\t",
                "transcript_id",
                "chr",
                "strand",
                "current_start",
                "updated_start",
                "current_end",
                "updated_end",
                "LEAP_5_tag",
                "LEAP_3_tag",
                "skipped_transcript",
                "warnings",
                )."\n";

#Loop over transcripts, get gene stable ids and check that the coordinates are consistent with the LEAP tags
my %transcripts_to_update;
my %genes_to_update;
my %new_gene_coords;
foreach my $transcript_id (keys %transcript_coords){
  my $skip_transcript = 0;
  my @warnings;
  my $transcript = $ta->fetch_by_stable_id($transcript_id);
  if ($transcript){
    my $gene = $transcript->get_Gene;
    $new_gene_coords{$gene->stable_id}{'start'} = $gene->start;
    $new_gene_coords{$gene->stable_id}{'end'} = $gene->end;

    #Sanity check: first/last exons of transcript in database must be present in the input gff3 file
    my @exons = @{$transcript->get_all_Exons};
    my $exon_count = scalar(@exons);
    if (!$exon_coords{$transcript_id}{1} or ($transcript->strand == 1 and $exon_coords{$transcript_id}{1}{'end'} != $exons[0]->end)
        or ($transcript->strand == -1 and $exon_coords{$transcript_id}{1}{'start'} != $exons[0]->start)){
      push(@warnings, "missing_or_wrong_first_exon");
      $skip_transcript = 1;
    }
    if (!$exon_coords{$transcript_id}{$exon_count} or ($transcript->strand == 1 and $exon_coords{$transcript_id}{$exon_count}{'start'} != $exons[$#exons]->start)
        or ($transcript->strand == -1 and $exon_coords{$transcript_id}{$exon_count}{'end'} != $exons[$#exons]->end)){
      push(@warnings, "missing_or_wrong_last_exon");
      $skip_transcript = 1;
    }

    #Sanity check: transcript start/end changes in the input gff3 file with respect to the database must come with the corresponding LEAP tag
    if ($transcript_coords{$transcript_id}{'start'} != $transcript->start){
      if ($transcript_coords{$transcript_id}{'start'} < $transcript->start){
        $transcripts_to_update{$transcript_id} = 1;
        $genes_to_update{$gene->stable_id} = 1;
        $new_gene_coords{$gene->stable_id}{'start'} = min($new_gene_coords{$gene->stable_id}{'start'}, $transcript->start);
        if ($transcript->strand == 1 and !($leap_tags{"LEAP_5"}{$transcript_id})){
          push(@warnings, "missing_LEAP_5_tag_start_change");
          $skip_transcript = 1;
        }
        elsif ($transcript->strand == -1 and !($leap_tags{"LEAP_3"}{$transcript_id})){
          push(@warnings, "missing_LEAP_3_tag_start_change");
          $skip_transcript = 1;
        }
      }
      else{
        push(@warnings, "unexpected_start_change");
        $skip_transcript = 1;
      }
    }
    if ($transcript_coords{$transcript_id}{'end'} != $transcript->end){
      if ($transcript_coords{$transcript_id}{'end'} > $transcript->end){
        $transcripts_to_update{$transcript_id} = 1;
        $genes_to_update{$gene->stable_id} = 1;
        $new_gene_coords{$gene->stable_id}{'end'} = max($new_gene_coords{$gene->stable_id}{'end'}, $transcript->end);
        if ($transcript->strand == 1 and !($leap_tags{"LEAP_3"}{$transcript_id})){
          push(@warnings, "missing_LEAP_3_tag_end_change");
          $skip_transcript = 1;
        }
        elsif ($transcript->strand == -1 and !($leap_tags{"LEAP_5"}{$transcript_id})){
          push(@warnings, "missing_LEAP_5_tag_end_change");
          $skip_transcript = 1;
        }
      }
      else{
        push(@warnings, "unexpected_end_change");
        $skip_transcript = 1;
      }
    }
  }
  #If transcript stable id from gff3 file is not current in the database, skip the transcript
  else{
    push(@warnings, "transcript_not_found_in_database");
    $skip_transcript = 1;
  }
  #Delete skipped transcripts from the hash before going on
  if ($skip_transcript){
    delete($transcripts_to_update{$transcript_id});
  }
  #Print results of the sanity check in the <OUTFILE>_check.txt file
  print CHK join("\t",
                $transcript_id,
                $transcript_coords{$transcript_id}{'chr'},
                $transcript_coords{$transcript_id}{'strand'},
                defined($transcript) ? $transcript->start : "na",
                $transcript_coords{$transcript_id}{'start'},
                defined($transcript) ? $transcript->end : "na",
                $transcript_coords{$transcript_id}{'end'},
                defined($leap_tags{"LEAP_5"}{$transcript_id}) ? "LEAP_5" : "na",
                defined($leap_tags{"LEAP_3"}{$transcript_id}) ? "LEAP_3" : "na",
                $skip_transcript ? "yes" : "no",
                join(", ", @warnings) || "na",
                )."\n";
}
print "DONE\n";


#Go on to edit the annotation
open (ED, ">".$outfile."_edits.log") or die "Can't open file ".$outfile."_edits.log:$!";
if ($write){
  print "Updating the gene annotation ... ";
  open (PEC, ">".$outfile."_postedit_check.txt") or die "Can't open file ".$outfile."_postedit_check.txt:$!";
}
else{
  #If no -write flag, the script still prepares the annotation update and may output some useful warning messages
  print "Simulating the gene annotation update...\n";
}

#Loop over list of genes to be altered
#Log file: <OUTFILE>_edits.log
my %gene_update_plans;
foreach my $gene_id (sort keys %genes_to_update){
  print ED "\n#######\n\nUpdating gene $gene_id ...\n";
  my $gene = $ga->fetch_by_stable_id($gene_id);

  #Process gene updates and store new gene annotation in the database
  my $msg = process_gene($gene);
  if ($msg =~ /lock ok,write ok,unlock ok/){
    print ED "\nRESULT: gene ".$gene->stable_id." was modified\n";
  }
  elsif ($msg =~ /unchanged/){
    print ED "\nRESULT: gene ".$gene->stable_id." was not modified\n";
  }
  else{
    print ED "\nRESULT: gene ".$gene->stable_id." could not be unlocked\n";
  }


  #If the -write flag was set, carry out some post-edits sanity checks
  #and print the output in the <OUTFILE>_postedit_check.txt file
  if ($write){
    my $updated_gene = $ga->fetch_by_stable_id($gene_id);
    # Compare against planned union endpoints, which may be farther than either
    # individual input row. Also verify that original MANE structures stayed fixed.
    foreach my $finding (leap_check_uploaded_gene($updated_gene, $gene_update_plans{$gene_id} || [])){
      print PEC join("\t", $gene_id, @$finding)."\n";
    }
  }
}
print "DONE\n";

if ($write){
  close (PEC);
}

#############


sub process_gene {
  my ($gene) = shift;
  my $g_msg = "unchanged";

  #Get Otter region for the gene
  #The region has to span the whole host gene to stop this from being "spliced out" as an incomplete gene.
  #This would prevent the new annotation from being added to the host gene.
  my $local_server = Bio::Otter::Server::Support::Local->new( otter_dba => $otter_dba );
  $local_server->set_params(
        dataset => $dataset_name,
        chr     => $gene->seq_region_name,
        start   => $new_gene_coords{$gene->stable_id}{'start'},
        end     => $new_gene_coords{$gene->stable_id}{'end'},
        cs      => "chromosome",
        csver   => $assembly_version,
  );
  my $region_action = Bio::Otter::ServerAction::Script::Region->new_with_slice($local_server);
  my $region = $region_action->get_region; #Bio::Vega::Region object

  #Update existing gene: host gene that provided its stable id to the gene being stored
  if (scalar(grep {$_->stable_id eq $gene->stable_id} $region->genes) == 1){
    my ($db_gene) = grep {$_->stable_id eq $gene->stable_id} $region->genes;

    # Plan the whole gene before changing objects. This prevents both duplicate
    # copies and a later input row undoing a reused model's farther endpoints.
    my $plans = leap_plan_updates($db_gene, \%transcripts_to_update, \%transcript_coords);
    $gene_update_plans{$gene->stable_id} = $plans;

    #Edit transcripts
    my $changed_transcript = 0; #Flag indicates if at least one transcript was modified
    foreach my $plan (@$plans){
      print ED "PLAN: ".$plan->{'action'}." target=".$plan->{'target_id'}.
        " sources=".join(',', @{$plan->{'source_ids'}}).
        " endpoints=".$plan->{'start'}."-".$plan->{'end'}."\n";
      next if $plan->{'action'} eq 'already_covered';
      my $db_tr = $plan->{'transcript'};
      if ($plan->{'action'} eq 'create_canonical_extended'){
        $db_tr = $db_tr->new_dissociated_copy;
        $db_gene->add_Transcript($db_tr);
      }
      print ED "\n#Transcript ".($db_tr->dbID ? $db_tr->stable_id : $db_tr->stable_id." (copy)")."\n";
      my $old_start = $db_tr->seq_region_start;
      my $old_end = $db_tr->seq_region_end;
      my $new_start = $plan->{'start'};
      my $new_end = $plan->{'end'};

      #Calculate the slice offset, i.e. start of the Otter region w.r.t. the chromosome start.
      #Update coordinates of the exon object taking into account this slice offset,
      #or the exon coordinates in the havana database will be wrong
      my $slice_offset = $region->slice->start - 1;

      #Store exon and translation coordinates to compare them with those after the transcript update
      my @current_tr_exons = @{$db_tr->get_all_Exons};
      print ED "BEFORE - exons: ".join('+', map {$_->vega_hashkey} @current_tr_exons)."\n";
      my $old_translation_start = $db_tr->coding_region_start + $slice_offset;
      my $old_translation_end = $db_tr->coding_region_end + $slice_offset;
      print ED "BEFORE - translation: $old_translation_start - $old_translation_end\n";
      #Store translation sequence for sanity check after the update (it must not change)
      my $old_translation = $db_tr->translation->seq;


      #Update transcript coordinates, which involves updating or replacing terminal exons and updating
      #the translation start/end coordinates in the exons.
      #There are several options when an exon has to be updated: update the current exon, create a new exon
      #or replace the exon with another existing exon.

      #It has been observed that if an exon is updated/created in any transcript to match an existing exon
      #in a MANE transcript, the new exon may replace the exon of the MANE transcript when the annotation
      #is saved in the database, also updating the last modification date of the MANE transcript.
      #This is not desirable, so the script now searches for an existing matching MANE exon and replaces
      #the exon to be updated with that exon, which prevents the update of the MANE transcript.

      if ($old_start ne $new_start){
        #First exon, positive strand
        if ($db_tr->strand == 1){
          print ED "Updating transcript start (+) from $old_start to $new_start\n";
          #By default, we want to update the current start exon
          my $new_start_exon = $db_tr->start_Exon;
          #But...
          #a) We may want to create a new exon to be updated if the current exon is shared with a MANE transcript,
          #   or when editing a canonical extended model (including a MANE copy).
          if ($plan->{'canonical'} or exon_in_mane_transcript($db_tr->start_Exon, $db_gene) or
              scalar grep {$_->value =~/^(MANE_select|MANE_plus_clinical)$/} @{$db_tr->get_all_Attributes('remark')}
          ){
            $new_start_exon = $db_tr->start_Exon->new_dissociated_copy;
            $new_start_exon->stable_id("");
            $new_start_exon->version("");
          }
          #b) If the updated exon coordinates would match an existing MANE exon, use this exon instead of a new one
          elsif (my $mane_exon = find_matching_exon_in_mane_transcript($new_start - $slice_offset, $db_tr->start_Exon->end, $db_tr->start_Exon->phase, $db_tr->start_Exon->end_phase, $db_gene)){
            $new_start_exon = $mane_exon;
          }

          #Set new exon start and swap the old and new exons
          $new_start_exon->start($new_start - $slice_offset);
          $db_tr->swap_exons($db_tr->start_Exon, $new_start_exon);
          #The "swap_exons" method updates the exons in the translation object,
          #but not the translation start position within the start exon.
          #We must do this now if the translation starts in the first exon.
          if ($db_tr->translation->start_Exon == $new_start_exon){
            my $old_seq_start = $db_tr->translation->start;
            my $new_seq_start = $db_tr->translation->start + $old_start - $new_start;
            print ED "Updating translation seq_start (+) from $old_seq_start to $new_seq_start\n";
            $db_tr->translation->start($new_seq_start);
            #If it is a single-exon CDS, edit the translation end position too
            if ($db_tr->translation->end_Exon == $new_start_exon){
              my $old_seq_end = $db_tr->translation->end;
              my $new_seq_end = $db_tr->translation->end + $old_start - $new_start;
              print ED "Updating translation seq_end (+) from $old_seq_end to $new_seq_end\n";
              $db_tr->translation->end($new_seq_end);
            }
          }
          #Add the LEAP_5 remark attribute
          leap_add_remark_once($db_tr, $leap_5_att);
        }

        #Last exon, negative strand
        elsif ($db_tr->strand == -1){
          print ED "Updating transcript start (-) from $old_start to $new_start\n";
          #By default, we want to update the current end exon
          my $new_end_exon = $db_tr->end_Exon;
          #But...
          #a) We may want to create a new exon to be updated if the current exon is shared with a MANE transcript,
          #   or when editing a canonical extended model (including a MANE copy).
          if ($plan->{'canonical'} or exon_in_mane_transcript($db_tr->end_Exon, $db_gene) or
              scalar grep {$_->value =~/^(MANE_select|MANE_plus_clinical)$/} @{$db_tr->get_all_Attributes('remark')}
          ){
            $new_end_exon = $db_tr->end_Exon->new_dissociated_copy;
            $new_end_exon->stable_id("");
            $new_end_exon->version("");
          }
          #b) If the updated exon coordinates would match an existing MANE exon, use this exon instead of a new one
          elsif (my $mane_exon = find_matching_exon_in_mane_transcript($new_start - $slice_offset, $db_tr->end_Exon->end, $db_tr->end_Exon->phase, $db_tr->end_Exon->end_phase, $db_gene)){
            $new_end_exon = $mane_exon;
          }

          #Set new exon start and swap the old and new exons
          $new_end_exon->start($new_start - $slice_offset);
          $db_tr->swap_exons($db_tr->end_Exon, $new_end_exon);
          #No need to update the translation end position within the end exon on the negative strand
          #Add the LEAP_3 remark attribute
          leap_add_remark_once($db_tr, $leap_3_att);
        }
        $changed_transcript = 1;
      }

      if ($old_end ne $new_end){
        #Last exon, positive strand
        if ($db_tr->strand == 1){
          print ED "Updating transcript end (+) from $old_end to $new_end\n";
          #By default, we want to update the current end exon
          my $new_end_exon = $db_tr->end_Exon;
          #But...
          #a) We may want to create a new exon to be updated if the current exon is shared with a MANE transcript,
          #   or when editing a canonical extended model (including a MANE copy).
          if ($plan->{'canonical'} or exon_in_mane_transcript($db_tr->end_Exon, $db_gene) or
              scalar grep {$_->value =~/^(MANE_select|MANE_plus_clinical)$/} @{$db_tr->get_all_Attributes('remark')}
          ){
            $new_end_exon = $db_tr->end_Exon->new_dissociated_copy;
            $new_end_exon->stable_id("");
            $new_end_exon->version("");
          }
          #b) If the updated exon coordinates would match an existing MANE exon, use this exon instead of a new one
          elsif (my $mane_exon = find_matching_exon_in_mane_transcript($db_tr->end_Exon->start, $new_end - $slice_offset, $db_tr->end_Exon->phase, $db_tr->end_Exon->end_phase, $db_gene)){
            $new_end_exon = $mane_exon;
          }

          #Set new exon start and swap the old and new exons
          $new_end_exon->end($new_end - $slice_offset);
          $db_tr->swap_exons($db_tr->end_Exon, $new_end_exon);
          #no need to update the translation end position within the end exon on the positive strand
          #Add LEAP_3 remark attribute
          leap_add_remark_once($db_tr, $leap_3_att);
        }

        #First exon, negative strand
        elsif ($db_tr->strand == -1){
          print ED "Updating transcript end (-) from $old_end to $new_end\n";
          #By default, we want to update the current start exon
          my $new_start_exon = $db_tr->start_Exon;
          #But...
          #a) We may want to create a new exon to be updated if the current exon is shared with a MANE transcript,
          #   or when editing a canonical extended model (including a MANE copy).
          if ($plan->{'canonical'} or exon_in_mane_transcript($db_tr->start_Exon, $db_gene) or
              scalar grep {$_->value =~/^(MANE_select|MANE_plus_clinical)$/} @{$db_tr->get_all_Attributes('remark')}
          ){
            $new_start_exon = $db_tr->start_Exon->new_dissociated_copy;
            $new_start_exon->stable_id("");
            $new_start_exon->version("");
          }
          #b) If the updated exon coordinates would match an existing MANE exon, use this exon instead of a new one
          elsif (my $mane_exon = find_matching_exon_in_mane_transcript($db_tr->start_Exon->start, $new_end - $slice_offset, $db_tr->start_Exon->phase, $db_tr->start_Exon->end_phase, $db_gene)){
            $new_start_exon = $mane_exon;
          }

          #Set new exon end and swap the old and new exons
          $new_start_exon->end($new_end - $slice_offset);
          $db_tr->swap_exons($db_tr->start_Exon, $new_start_exon);
          #The "swap_exons" method updates the exons in the translation object,
          #but not the translation start position within the start exon.
          #We must do this now if the translation starts in the first exon.
          if ($db_tr->translation->start_Exon == $new_start_exon){
            my $old_seq_start = $db_tr->translation->start;
            my $new_seq_start = $db_tr->translation->start + $new_end - $old_end;
            print ED "Updating translation seq_start (+) from $old_seq_start to $new_seq_start\n";
            $db_tr->translation->start($new_seq_start);
            #If single-exon CDS, edit the translation end position too
            if ($db_tr->translation->end_Exon == $new_start_exon){
              my $old_seq_end = $db_tr->translation->end;
              my $new_seq_end = $db_tr->translation->end + $new_end - $old_end;
              print ED "Updating translation seq_end (-) from $old_seq_end to $new_seq_end\n";
              $db_tr->translation->end($new_seq_end);
            }
          }
          #Add LEAP_5 remark attribute
          leap_add_remark_once($db_tr, $leap_5_att);
        }
        $changed_transcript = 1;
      }
      my @updated_tr_exons = @{$db_tr->get_all_Exons};
      print ED "AFTER - exons ".join('+', map {$_->vega_hashkey} @updated_tr_exons)."\n";
      if (scalar @current_tr_exons != scalar @updated_tr_exons){
        print ED "WARNING: number of exons has changed!\n";
      }
      for (my $i=0; $i++; $i <= $#updated_tr_exons){
        if (($i == 0 and $db_gene->strand == 1) or ($i == $#updated_tr_exons and $db_gene->strand == -1)){
          if ($updated_tr_exons[$i]->start > $current_tr_exons[$i]->end){
            print ED "WARNING: first exon start is greater than before!\n";
          }
          if ($updated_tr_exons[$i]->end != $current_tr_exons[$i]->end){
            print ED "WARNING: first exon end has changed!\n";
          }
        }
        elsif (($i == $#updated_tr_exons and $db_gene->strand == 1) or ($i == 0 and $db_gene->strand == -1)){
          if ($updated_tr_exons[$i]->end < $current_tr_exons[$i]->end){
            print ED "WARNING: last exon end is smaller than before!\n";
          }
          if ($updated_tr_exons[$i]->start != $current_tr_exons[$i]->start){
            print ED "WARNING: last exon start has changed!\n";
          }
        }
        else{
          if ($updated_tr_exons[$i]->start != $current_tr_exons[$i]->start or $updated_tr_exons[$i]->end != $current_tr_exons[$i]->end){
            print ED "WARNING: internal exon coordinates have changed!\n";
          }
        }
      }
      my $new_translation_start = $db_tr->coding_region_start + $slice_offset;
      my $new_translation_end = $db_tr->coding_region_end + $slice_offset;
      print ED "AFTER - translation: $new_translation_start - $new_translation_end\n";

      unless ($old_translation_start == $new_translation_start  and $old_translation_end == $new_translation_end){
        print ED "WARNING: translation coordinates have changed!\n";
      }

      if ($old_translation ne $db_tr->translation->seq){
        print ED "WARNING: translation sequence has changed!\n";
      }

      #If this is the copy of a MANE transcript, remove its stable id so that it is created as a new transcript in the database
      #Also remove the MANE attribute remark and add an 'Ensembl_canonical_extended' remark
      if (my ($att) = grep {$_->value =~/^(MANE_select|MANE_plus_clinical)$/} @{$db_tr->get_all_Attributes('remark')}){
        $db_tr->stable_id("");
        $db_tr->version("");
        $db_tr->translation->stable_id("");
        $db_tr->translation->version("");
        my $attribs = $db_tr->{'attributes'};
        @{$attribs} = grep {$_->value ne "MANE_select" and $_->value ne "MANE_plus_clinical"} @{$attribs};
        @{$attribs} = grep {!($_->value =~ /^NoMerge:/)} @{$attribs};
        $db_tr->add_Attributes($canon_ext_att);
        $db_tr->add_Attributes(Bio::EnsEMBL::Attribute->new(-code => 'hidden_remark', -value => "Copy of ".$att->value." transcript, extended by LEAP"));
        #Add a new transcript name since name duplicates will prevent Otter from opening a session
        my $new_transcript_name = get_new_transcript_name($db_gene, $otter_dba);
        @{$attribs} = grep {$_->code ne "name"} @{$attribs};
        $db_tr->add_Attributes(Bio::EnsEMBL::Attribute->new(-code => 'name', -value => $new_transcript_name));
      }
    }

    #Write region
    if ($changed_transcript){
      $local_server->authorized_user($gene->gene_author->name); # preserve authorship
      if ($write){
        $g_msg = write_gene_region($region_action, $region);
      }
      print ED $g_msg."\n";
    }
  }
  else{
    print ED "More than one gene with stable id ".$gene->stable_id."!!!\n";
  }

  return $g_msg;
}


#Given a region (slice of the genome with gene annotation modified by this script), lock the region in the database
#so that no one else can edit it while it is being processed, and update the annotation in the database
sub write_gene_region {
  my ($region_action, $region) = @_;

  my @msg;
  my $lock;
  try {
    $region_action->server->add_param( hostname => hostname );
    $lock = $region_action->lock_region;
    push @msg, 'lock ok';
  }
  catch {
    my ($err) = ($_ =~ m/^MSG: (Failed to lock.*)$/m);
    push @msg, "lock failed: '$err'";
    push @msg, "lock failed: '$_'";
  };

  if ($lock) {
    my $new_region;
    try {
      $region_action->server->set_params( data => $region, locknums => $lock->{locknums} );
      $new_region = $region_action->write_region;
      push @msg, 'write ok';
    }
    catch {
      my $err = $_;
      chomp $err;
      push @msg, "write failed: '$err'";
    };
    $region_action->server->set_params( locknums => $lock->{locknums} );
    $region_action->unlock_region;
    push @msg, 'unlock ok';
  }

  return join(',', @msg);
}


sub get_new_transcript_name {
  my ($gene, $dba) = @_;
  #Get a list of existing clone-based transcript names
  #in order to infer the clone-based gene name and the next transcript number.
  #The current gene name may have been replaced with a symbol, so use this as the default name
  #but replace it with a clone-based name if available in the transcript names
  my $gene_name = $gene->get_all_Attributes('name')->[0]->value;
  my $max = 0;
  foreach my $tr (@{$gene->get_all_Transcripts}){
    if (scalar @{$tr->get_all_Attributes('name')}){
      my ($root, $n) = $tr->get_all_Attributes('name')->[0]->value =~ /^(.+)-(\d{3,})$/;
      if ($n and $n > $max){
        $max = $n;
      }
      if ($root){
        $gene_name = $root;
      }
    }
  }

  my $new_tr_name;
  my $name_in_db;
  do{
    $new_tr_name = sprintf("%s-%03d", $gene_name, ++$max);
    #Sanity check - avoid duplication with past or present transcripts
    my $dbcon = $dba->dbc;
    my $sth = $dbcon->prepare("SELECT t.stable_id FROM transcript t, transcript_attrib ta
                                WHERE t.transcript_id=ta.transcript_id AND ta.value=\"$new_tr_name\"");
    $sth->execute();
    my @tsi = $sth->fetchrow_array;
    $sth->finish();
    $name_in_db = scalar @tsi;
  }
  while ($name_in_db);

  return $new_tr_name;
}


#Find if an exon is also attached to a MANE transcript
sub exon_in_mane_transcript {
  my ($query_exon, $gene) = @_;
  foreach my $tr (@{$gene->get_all_Transcripts}){
    if (grep {$_->value =~/^(MANE_select|MANE_plus_clinical)$/} @{$tr->get_all_Attributes('remark')}){
      foreach my $exon (@{$tr->get_all_Exons}){
        if ($exon->stable_id eq $query_exon->stable_id){
          return 1;
        }
      }
    }
  }
  return;
}


#Given a set of start, end, phase and end_phase values for an updated exon,
#find if an existing exon in a MANE transcript has the same features
sub find_matching_exon_in_mane_transcript {
  my ($start, $end, $phase, $end_phase, $gene) = @_;
  #print "START=$start, END=$end, PHASE=$phase, END_PHASE=$end_phase\n";
  foreach my $tr (@{$gene->get_all_Transcripts}){
    if (grep {$_->value =~/^(MANE_select|MANE_plus_clinical)$/} @{$tr->get_all_Attributes('remark')}){
      foreach my $exon (@{$tr->get_all_Exons}){
        if ($exon->start == $start and $exon->end == $end and $exon->phase == $phase and $exon->end_phase == $end_phase){
          return $exon;
        }
      }
    }
  }
  return;
}


# BEGIN CANONICAL_EXTENDED_HELPERS
# These helpers use only core Perl and transcript objects. They are also tested
# independently of the Otter installation, without connecting to a database.
use Scalar::Util qw(refaddr);

sub leap_mane_remark {
  my ($transcript) = @_;
  my ($attribute) = grep {$_->value =~ /^(MANE_select|MANE_plus_clinical)$/}
    @{$transcript->get_all_Attributes('remark')};
  return $attribute ? $attribute->value : '';
}

sub leap_is_canonical_extended {
  my ($transcript) = @_;
  return scalar grep {lc($_->value) eq 'ensembl_canonical_extended'}
    @{$transcript->get_all_Attributes('remark')};
}

sub leap_structure_key {
  my ($transcript) = @_;
  my @exons = sort {$a->start <=> $b->start || $a->end <=> $b->end}
    @{$transcript->get_all_Exons};
  die "Transcript has no exons: ".$transcript->stable_id."\n" unless @exons;
  # The region used during editing may have a different slice origin from a
  # subsequently fetched database model. Normalize protected positions to genome.
  my $offset = $transcript->seq_region_start - $transcript->start;
  my @key = ($transcript->seq_region_name, $transcript->strand, scalar @exons);
  foreach my $coordinate ($transcript->coding_region_start, $transcript->coding_region_end){
    push @key, defined($coordinate) ? $coordinate + $offset : 'no_CDS';
  }
  for my $index (0 .. $#exons){
    my $exon = $exons[$index];
    push @key, $index > 0 ? $exon->start + $offset : 'outer_start';
    push @key, $index < $#exons ? $exon->end + $offset : 'outer_end';
    push @key, map {defined($_) ? $_ : '.'} ($exon->phase, $exon->end_phase);
  }
  return join('|', @key);
}

sub leap_add_remark_once {
  my ($transcript, $attribute) = @_;
  return if grep {$_->value eq $attribute->value} @{$transcript->get_all_Attributes('remark')};
  $transcript->add_Attributes($attribute);
}

sub leap_plan_updates {
  my ($gene, $selected, $coordinates) = @_;
  my @original = @{$gene->get_all_Transcripts};
  my @canonical = grep {leap_is_canonical_extended($_)} @original;
  foreach my $transcript (@canonical){
    die "Canonical extended transcript is also marked MANE: ".$transcript->stable_id."\n"
      if leap_mane_remark($transcript);
  }
  my %plans;
  foreach my $source (sort {$a->stable_id cmp $b->stable_id} @original){
    my $source_id = $source->stable_id;
    next unless $selected->{$source_id};
    my $desired = $coordinates->{$source_id};
    die "Missing or invalid proposed endpoints for $source_id\n"
      unless $desired && defined($desired->{'start'}) && defined($desired->{'end'})
      && $desired->{'start'} =~ /^\d+$/ && $desired->{'end'} =~ /^\d+$/
      && $desired->{'start'} > 0 && $desired->{'end'} >= $desired->{'start'};
    die "Proposed endpoints would shorten source transcript $source_id\n"
      if $desired->{'start'} > $source->seq_region_start || $desired->{'end'} < $source->seq_region_end;
    my $mane = leap_mane_remark($source);
    my $target = $source;
    my $structure = leap_structure_key($source);
    my $new_copy = 0;
    if ($mane){
      my @matches = grep {leap_structure_key($_) eq $structure} @canonical;
      die "Multiple matching canonical extended transcripts for $source_id: ".
        join(',', sort map {$_->stable_id} @matches)."\n" if @matches > 1;
      if (@matches){
        $target = $matches[0];
      }
      else{
        $new_copy = 1;
      }
    }
    my $key = $new_copy ? 'copy:'.$structure : 'existing:'.refaddr($target);
    unless ($plans{$key}){
      $plans{$key} = {
        transcript => $target, target_id => $new_copy ? 'new_copy_of_'.$source_id : $target->stable_id,
        new_copy => $new_copy, structure => $structure,
        start => $target->seq_region_start, end => $target->seq_region_end,
        original_start => $target->seq_region_start, original_end => $target->seq_region_end,
        source_ids => [], mane_sources => [], canonical => $new_copy || leap_is_canonical_extended($target),
      };
    }
    my $plan = $plans{$key};
    $plan->{'start'} = $desired->{'start'} if $desired->{'start'} < $plan->{'start'};
    $plan->{'end'} = $desired->{'end'} if $desired->{'end'} > $plan->{'end'};
    push @{$plan->{'source_ids'}}, $source_id;
    if ($mane){
      push @{$plan->{'mane_sources'}}, {id => $source_id, remark => $mane, structure => $structure,
        start => $source->seq_region_start, end => $source->seq_region_end};
    }
  }
  my @result;
  foreach my $key (sort keys %plans){
    my $plan = $plans{$key};
    $plan->{'action'} = $plan->{'new_copy'} ? 'create_canonical_extended'
      : ($plan->{'start'} == $plan->{'original_start'} && $plan->{'end'} == $plan->{'original_end'})
      ? 'already_covered' : $plan->{'canonical'} ? 'update_canonical_extended' : 'update_transcript';
    push @result, $plan;
  }
  return \@result;
}

sub leap_check_uploaded_gene {
  my ($gene, $plans) = @_;
  my @transcripts = @{$gene->get_all_Transcripts};
  my %by_id = map {$_->stable_id => $_} @transcripts;
  my @findings;
  foreach my $plan (@$plans){
    foreach my $mane (@{$plan->{'mane_sources'}}){
      my $original = $by_id{$mane->{'id'}};
      my $unchanged = $original && $original->seq_region_start == $mane->{'start'}
        && $original->seq_region_end == $mane->{'end'}
        && leap_structure_key($original) eq $mane->{'structure'}
        && leap_mane_remark($original) eq $mane->{'remark'};
      push @findings, [$mane->{'id'}, $mane->{'remark'}.($unchanged ? '_not_updated' : '_wrongly_updated_or_missing')];
    }
    my @targets = $plan->{'new_copy'}
      ? grep {leap_is_canonical_extended($_) && !leap_mane_remark($_)
          && leap_structure_key($_) eq $plan->{'structure'}
          && $_->seq_region_start == $plan->{'start'} && $_->seq_region_end == $plan->{'end'}} @transcripts
      : ($by_id{$plan->{'target_id'}} ? ($by_id{$plan->{'target_id'}}) : ());
    my $ok = @targets == 1 && $targets[0]->seq_region_start == $plan->{'start'}
      && $targets[0]->seq_region_end == $plan->{'end'}
      && leap_structure_key($targets[0]) eq $plan->{'structure'}
      && (!$plan->{'canonical'} || (leap_is_canonical_extended($targets[0]) && !leap_mane_remark($targets[0])));
    push @findings, [@targets == 1 ? $targets[0]->stable_id : $plan->{'target_id'},
      ($ok ? 'planned_endpoints_and_structure_match' : 'ERROR_planned_upload_mismatch').
      ',action='.$plan->{'action'}.',sources='.join('+', @{$plan->{'source_ids'}})];
  }
  return @findings;
}
