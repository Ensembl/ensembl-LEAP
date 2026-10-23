#!/hps/software/users/ensembl/ensw/swenv/env/prod-perl-5-26/bin/perl -w

# Create a GTF or GFF file from an Ensembl core schema database

use strict;
use Getopt::Long;
use Gencode::Default;
use Gencode::Ensembl2GTF_new;
use Gencode::Ensembl2GFF3_new;
use Bio::EnsEMBL::SeqFeature;
use FileHandle;
$|=1;


my $dbhost    = undef;
my $dbport    = 3306;
my $dbuser    = undef;
my $dbpass    = undef;
my $dbname    = undef;
my $dnadbhost = undef;
my $dnadbport = 3306;
my $dnadbuser = undef;
my $dnadbpass = undef;
my $dnadbname = undef;
my $chromosome = undef;
my $start      = undef;
my $end        = undef;
my $outfile   = undef;
my $verbose    = 0;
my $source     = undef;
my $biotype    = undef;
my $analysis   = undef;
my $no_version = 0;
my $assembly_name;
my $vega_only;
my $allow_all_artifacts;
my $gene_list;
my $transcript_list;
my $transform_to_asv;
my $format;
my $add_chr;
my $mane;
my $include_non_ref;

&GetOptions(
            'host=s'          => \$dbhost,
            'port=s'          => \$dbport,
            'user=s'          => \$dbuser,
            'pass=s'          => \$dbpass,
            'dbname=s'        => \$dbname,
            'dnahost=s'       => \$dnadbhost,
            'dnaport=s'       => \$dnadbport,
            'dnauser=s'       => \$dnadbuser,
            'dnapass=s'       => \$dnadbpass,
            'dnadbname=s'     => \$dnadbname,            
            'file=s'          => \$outfile,
            'chr=s'           => \$chromosome,
            'start=s'         => \$start,
            'end=s'           => \$end,
            'source=s'        => \$source,
            'biotype=s'       => \$biotype,
            'analysis=s'      => \$analysis,
            'vega'            => \$vega_only,
            'artifacts'       => \$allow_all_artifacts,
            'no_version'      => \$no_version,
            'verbose!'        => \$verbose,
            'asv=s'           => \$assembly_name,
            'list=s'          => \$gene_list,
            'tlist=s'         => \$transcript_list,
            'new_asv=s'       => \$transform_to_asv,
            'format=s'        => \$format,
            'add_chr!'        => \$add_chr,
            'mane=s'          => \$mane,
            'non_ref!'        => \$include_non_ref,
           );


my $usage =  <<_USAGE_;

This script dumps annotation from an Ensembl-like database into a GTF file

perl db2gtf.pl -host ... -port ... -user ... -pass ... -dbname ... [-asv ...] [-chr...] [-start...] [-end...] [-source ...] [-biotype ...] [-analysis ...] [-vega] [-new_asv ...] [-no_version] [-add_chr] [-list] [-tlist] -format ... -file ...

 -host        database host [required]
 -port        database port (default: 3306) [optional]
 -user        database user [required]
 -pass        database password (default: no password) [optional]
 -dbname      database name [required]
 -dnahost     DNA database host [optional]
 -dnaport     DNA database port (default: 3306) [optional]
 -dnauser     DNA database user [optional]
 -dnapass     DNA database password (default: no password) [optional]
 -dnadbname   DNA database name [optional] 
 -asv         genome assembly version [optional]
 -chr         only dump annotation on this chromosome [optional]
 -start       only dump annotation from this genomic coordinate [optional]
 -end         only dump annotation up to this genomic coordinate [optional]
 -source      restrict to this gene source [optional]
 -biotype     restrict to this gene biotype [optional]
 -analysis    restrict to this gene analysis name [optional]
 -vega        exclude gene or transcripts with a "not for VEGA" attribute
 -artifacts   include all genes and transcripts with an "artifact" biotype
 -list        restrict to this gene stable id list [optional]
 -tlist       restrict to this transcript stable id list [optional]
 -new_asv     transform genes to this genome assembly version [optional]
 -no_version  do not append a version number to the stable ids [optional]
 -mane        add RefSeq dbxref attributes based on this MANE file [optional]
 -non_ref     include non-reference sequence regions (patches and haplotypes)
 -add_chr     prepend "chr" to the chromosome name
 -format      output file format, either gtf or gff3 [required]
 -file        output file name [required]

_USAGE_

die $usage unless ($dbhost and $dbuser and $dbname and $format and $outfile);


#connect to databases
my $dba = Gencode::Default->dbconnect($dbhost, $dbport, $dbname, $dbuser, $dbpass);
my $slice_adaptor = $dba->get_SliceAdaptor();
my $gene_adaptor = $dba->get_GeneAdaptor();
my $transcript_adaptor = $dba->get_TranscriptAdaptor();
my $sfa = $dba->get_SimpleFeatureAdaptor();

if ($dnadbhost and $dnadbuser and $dnadbname){
  my $dnadba = Gencode::Default->dbconnect($dnadbhost, $dnadbport, $dnadbname, $dnadbuser, $dnadbpass);
  $dba->dnadb($dnadba);
}

my %selected_genes;
if ($gene_list){
  open (LIST, $gene_list) or die "can't open $gene_list";
  while (<LIST>){
    chomp;
    $selected_genes{$_} = 1;
  }
  close (LIST);
}
my %printed_genes; #To avoid duplications in some particular cases

my %selected_transcripts;
if ($transcript_list){
  open (LIST, $transcript_list) or die "can't open $transcript_list";
  while (<LIST>){
    chomp;
    $selected_transcripts{$_} = 1;
  }
  close (LIST);
}

my %ens2refseq;
if ($mane){
  open (MANE, $mane) or die "can't open $mane";
  while (<MANE>){
    chomp;
    my @fs = split(/\t/);
    my $rs_tr = $fs[5];
    my $rs_prot = $fs[6];
    my $ens_tr = $fs[7];
    my $ens_prot = $fs[8];
    $ens2refseq{$ens_tr} = $rs_tr;
    $ens2refseq{$ens_prot} = $rs_prot; 
  }
  close (MANE);
}

my $gtffp = new FileHandle;
$gtffp->open(">$outfile") or die "Unable to open $outfile for write";
$gtffp->autoflush(1);

if ($format =~ "gff3"){
  print $gtffp "##gff-version 3\n"; #required as the GFF3 file first line
}

#fetch genes from merged database, print to GTF file
my @slicenames = ();
if($chromosome){
  #my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $onlychrom );
  push(@slicenames, $chromosome);
}
elsif ($dbname =~ /loutre/){ #seq_region_id order matches reference chromosome names
  @slicenames = map {$_->seq_region_name} grep {$_->coord_system->is_default} @{$slice_adaptor->fetch_all('chromosome', undef, 0, 1)};
}
else{
  if ($include_non_ref or scalar keys %selected_genes or scalar keys %selected_transcripts){
    @slicenames = map {$_->seq_region_name} grep {$_->coord_system->is_default} @{sort_slices_by_name($slice_adaptor->fetch_all('toplevel', undef, 1, 1))};
  }
  else{
    @slicenames = map {$_->seq_region_name} grep {$_->coord_system->is_default} @{sort_slices_by_name($slice_adaptor->fetch_all('toplevel', undef, 0, 1))};            
  }
}
my $ccdslist;
my $write_utr = 1;
my $coding_types = '^(protein_coding|nonsense_mediated_decay|non_stop_decay|polymorphic_pseudogene|protein_coding_LoF|IG_.+_gene|TR_.+_gene)$';

foreach my $slicename (@slicenames){
  next if $dbname =~ /loutre_(human|mouse)/ and !($slicename =~ /chr.+-38/);
  print STDERR "Checking chromosome $slicename\n";
  my $slice = $slice_adaptor->fetch_by_region("toplevel", $slicename, $start, $end, 1, $assembly_name) or die "can't get slice $slicename\n";
  print STDERR $slice->seq_region_name.":".$slice->start."-".$slice->end."\n"; 
  #Fetch slice genes
  #my $genes = $slice->get_all_Genes(undef, undef, 1); no lazy-loading
  my $genes;
  if (scalar keys %selected_genes and scalar keys %selected_genes < 2000){ #apparently, the default maximum stable id list size is 2000 in the fetch_all_by_stable_id_list method
    my @gsids = keys %selected_genes;
    $genes = $gene_adaptor->fetch_all_by_stable_id_list(\@gsids, $slice);
  }
  elsif (scalar keys %selected_transcripts and scalar keys %selected_transcripts < 2000){
    my @tsids = keys %selected_transcripts;
    my $trs = $transcript_adaptor->fetch_all_by_stable_id_list(\@tsids, $slice);
    my %seen;
    foreach my $tr (@$trs){
      next unless $tr->is_current;
      my $gene = $tr->get_Gene;
      unless ($seen{$gene->stable_id}){
        push(@{$genes}, $gene);
        $seen{$gene->stable_id} = 1;
      }
    }
  }  
  else{
    $genes = $slice->get_all_Genes();
  }
  print STDERR ($genes ? scalar(@$genes) : "0")." genes\n";
  my $localcoords = 0;
  my $trim_to_slice = 0;
  my $include_codons = 1;
  my $v_ta;
  #Get PAR regions (assembly exception features) associated to the slice
  my @par_regions = grep {$_->type eq "PAR"} @{$slice->get_all_AssemblyExceptionFeatures};
  
  my $seqname = $slice->seq_region_name;
  $seqname =~ s/-38//;
  if ($add_chr and !($seqname=~/^chr/) and !($seqname=~/^(HG|HS|KI|GL)/)){
    $seqname = "chr".$seqname;
  }

  if ($format eq "gff3"){
    print $gtffp "##sequence-region $seqname ".$slice->start." ".$slice->end."\n";
  }
  next unless $genes;
  
  GENE:
  foreach my $gene (@$genes){
    next unless $gene->is_current;
    if ($gene_list){
      next unless $selected_genes{$gene->stable_id};
    }
    if ($transcript_list){ #Skip gene unless one of its transcripts is on the list
      next unless scalar (grep {$selected_transcripts{$_->stable_id}} @{$gene->get_all_Transcripts});
    }
    next if $printed_genes{$gene->stable_id};
    unless ($allow_all_artifacts){
      next if $gene->biotype eq "artifact" and scalar @{$gene->get_all_Attributes('artef_dupl')} == 0;
    }

    if ($transform_to_asv){
      my $gene2 = $gene->transform($gene->coord_system_name, $transform_to_asv);
      unless ($gene2){
        print "\nCouldn't transform gene ".$gene->stable_id."\n\n";
        next GENE;
      }
      $gene = $gene2;
      #my $gene2 = $gene->transform($gene->coord_system_name, $transform_to_asv);
      #print "***".$gene2->stable_id."\t".$gene2->seq_region_name."\t".$gene2->start."\t".$gene2->end."***\n";
    }
    else{
      $gene = $gene->transform ("toplevel");
    }

    my $gene_name = get_name($gene);
    my $gene_id = $gene->stable_id;
    my $gene_version = $gene->version;
    $gene_id .= ".".$gene_version unless $no_version;
    my @extra_tags = ();
    #Add PAR tag to tag list if approppriate
    my $gene_on_par = 0;
    if (scalar(@par_regions) and is_on_PAR_region($gene, \@par_regions)){
      push(@extra_tags, "PAR");
      $gene_on_par = 1;
    }
    #Make PAR-specific gene id
    my $gene_id2 = $gene_id; #Keep a copy of the original id for GFF3
    if ($gene_on_par){
      $gene_id .= "_PAR_Y";
    }
    my $sourcename = $gene->source;
    if ($source){
      next unless $sourcename eq $source;
    }

    my $e = 0;
    my $gene_strand = _get_strand($gene);
    my $gene_biotype = $gene->biotype;
    if ($biotype){
      next unless $gene_biotype eq $biotype;
    }
    #my $gene_status = $gene->status;
    #$gene_status = "unknown" unless $gene_status;
    my $gene_status = "unknown";
    if ($analysis){
      next unless $gene->analysis->logic_name eq $analysis;
    }
    
    #my $elevel = _get_level($gene);

    #Create gene line
    my $gene_line = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
                            $seqname,
                            $sourcename,
                            'gene',
                            $gene->seq_region_start,
                            $gene->seq_region_end,
                            '.',
                            $gene_strand,
                            '.',
                           );

    #Add attributes
    if ($format eq "gtf"){
      $gene_line .= sprintf("gene_id \"%s\"; gene_type \"%s\"; gene_name \"%s\";",   
                            $gene_id,
                            $gene_biotype,
                            $gene_name,
                            #$elevel,
                           );
    }
    elsif ($format eq "gff3"){
      $gene_line .= sprintf("ID=%s;gene_id=%s;gene_type=%s;gene_name=%s;",     
                            $gene_id, #ID
                            $gene_id2,
                            $gene_biotype,
                            $gene_name,
                            #$elevel,
                           );
    }


    #normal gene line
    foreach my $attrib (@{$gene->get_all_Attributes('remark')}){
      if ($attrib->value eq "not for VEGA"){
        next GENE if $vega_only;
        push(@extra_tags, "tag \"not_for_VEGA\";");
      }
    }

    #Add gene tags
    if (scalar @extra_tags){
      if ($format eq "gtf"){
        foreach my $tag (@extra_tags){
          $gene_line .= " tag \"$tag\";";
        }
      }
      elsif ($format eq "gff3"){
        #Cluster all 'tag' attributes
        $gene_line .= "tag=".join(",", @extra_tags).";";
      }
    }

    #Remove trailing semicolon from GFF3 line if present
    if ($format eq "gff3" and $gene_line =~ /;$/){
      $gene_line =~ s/;$//;
    } 
    print $gtffp $gene_line."\n";



    TRANSCRIPT:
    foreach my $transcript (@{$gene->get_all_Transcripts}){
      if ($transcript_list){
        next unless $selected_transcripts{$transcript->stable_id};
      }
      my @ccdsids;
      my $extras;
      @extra_tags = ();
      $e = 0;
      my $ccdstagset = 0;
      my $stop_codon = undef;
      my $write_seleno = 0;
      my $cons_pseudos;
      my $elevel = 2;
      
      #Change source to transcript source if different
      if ($sourcename ne $transcript->source){
        $sourcename = $transcript->source;
      }
      
      my $transcript_name = _get_name($transcript);
      my $transcript_id = $transcript->stable_id;
      #Deal with NULL transcript ids (Ensembl RNAseq)
      if (!$transcript_id or $transcript_id eq "NULL"){
        $transcript_id = $gene_id;
        $transcript_id =~ s/ROUGHG/ROUGHT/;
      }
      unless ($allow_all_artifacts){
        next if $transcript->biotype eq "artifact" and scalar @{$gene->get_all_Attributes('artef_dupl')} == 0;
      }
      my $transcript_version = $transcript->version;
      $transcript_id .= ".".$transcript_version unless $no_version;

      my $transcript_biotype = $transcript->biotype;
      #my $transcript_status = $transcript->status;
      #$transcript_status = "unknown" unless $transcript_status;
      my $transcript_status = "unknown";
        
      #Filter out "not for VEGA" transcripts   
      foreach my $attrib (@{$transcript->get_all_Attributes('remark')}){
        if ($attrib->value eq "not for VEGA"){
          next TRANSCRIPT if $vega_only;
        }
      }
      
      #Add MANE tags
      if (scalar @{$transcript->get_all_Attributes('MANE_Select')}){
        push(@extra_tags, "MANE_Select");
      }
      elsif (scalar @{$transcript->get_all_Attributes('MANE_Plus_Clinical')}){
        push(@extra_tags, "MANE_Plus_Clinical");
      }
           
      #Add other annotation remarks
      my $remarks = transcript_remarks($transcript);
      push(@extra_tags, @$remarks);

      my $transcript_strand = _get_strand($transcript);

      if($transcript_strand ne $gene_strand){
        warn "\nSTRAND MISMATCH BETWEEN $transcript_id AND $gene_id!";
      }

      print STDERR "Transcript: $sourcename, $transcript_id\n" if($verbose);

      my $protein_id;
      if ($transcript->translation){
        $protein_id = $transcript->translation->stable_id;
        $protein_id .= ".".$transcript->translation->version unless $no_version;
      }

     #Make PAR-specific transcript id
      my $transcript_id2 = $transcript_id; #Keep a copy of the original id for GFF3
      if ($gene_on_par){
        $transcript_id .= "_PAR_Y";
      }

      #Create transcript line
      my $transcript_line = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
                               $seqname,
                               $sourcename,
                               'transcript',
                               $transcript->seq_region_start,
                               $transcript->seq_region_end,
                               '.',
                               $transcript_strand,
                               '.',
                              );

      #Add attributes
      if ($format eq "gtf"){
        $transcript_line .= sprintf("gene_id \"%s\"; transcript_id \"%s\"; ".
                               "gene_type \"%s\"; gene_name \"%s\"; ".
                               "transcript_type \"%s\"; transcript_name \"%s\";",
                               $gene_id,
                               $transcript_id,
                               $gene_biotype,
                               $gene_name,
                               $transcript_biotype,
                               $transcript_name,
                               #$elevel,
                              );
      }
      elsif ($format eq "gff3"){
        $transcript_line .= sprintf("ID=%s;Parent=%s;gene_id=%s;transcript_id=%s;".
                               "gene_type=%s;gene_name=%s;".
                               "transcript_type=%s;transcript_name=%s;",
                               $transcript_id, #ID
                               $gene_id, #Parent
                               $gene_id2,
                               $transcript_id2,
                               $gene_biotype,
                               $gene_name,
                               $transcript_biotype,
                               $transcript_name,
                               #$elevel,
                              );
      }

      my $transcript_line_extra = "";
      #Add transcript tags
      if (scalar @extra_tags){
        if ($format eq "gtf"){
          foreach my $tag (@extra_tags){
            $transcript_line_extra .= " tag \"$tag\";";
          }
        }
        elsif ($format eq "gff3"){
          $transcript_line_extra .= "tag=".join(",", @extra_tags).";";
        }
      }

      if ($protein_id){
        if ($format eq "gtf"){
          $transcript_line_extra .= " protein_id \"".$protein_id."\";";
        }
        elsif ($format eq "gff3"){
          $transcript_line_extra .= "protein_id=$protein_id;";
        }
      }
                 
      #For the MANE FTP files, add RefSeq transcript dbxref
      if ($mane){
        if ($format eq "gtf"){
          #db_xref "RefSeq:<accession.version>"
          $transcript_line_extra .= " db_xref \"RefSeq:".$ens2refseq{$transcript_id}."\";";
        }          
        elsif ($format eq "gff3"){
          #Dbxref=RefSeq:<accession.version>
          $transcript_line_extra .= "Dbxref=RefSeq:".$ens2refseq{$transcript_id}.";";
        }      
      }
      
      #Remove trailing ';' from the GFF3 file
      if ($format eq "gff3" and $transcript_line_extra =~ /;$/){
        $transcript_line_extra =~ s/;$//;
      }
      
      #Print transcript line
      $transcript_line .= $transcript_line_extra;
      print $gtffp $transcript_line."\n";
      
      #For the MANE FTP files, add RefSeq protein dbxref (CDS lines)
      if ($mane){
        if ($protein_id and $ens2refseq{$protein_id}){
          if ($format eq "gtf"){
            $transcript_line_extra .= " db_xref \"RefSeq:".$ens2refseq{$protein_id}."\";";
          }
          elsif ($format eq "gff3"){
            $transcript_line_extra .= ";Dbxref=RefSeq:".$ens2refseq{$protein_id}; 
            #provided two dbxref attributes, but only one will be printed depending on the feature type (hack in the Ensembl2GTF_new and Ensembl2GFF3_new modules)
          }
        }
      }

      #Pass transcript attributes to transcript child features
      if ($format eq "gtf"){
        $extras = $transcript_line_extra;
      }
      elsif ($format eq "gff3"){
        $extras = $transcript_line_extra;
      }

      if ($format eq "gtf"){
        write_transcript_gtf($gtffp, $slice, $gene, $transcript, $seqname, $sourcename,
                             $gene_id, $gene_biotype, $gene_status, $gene_version, $gene_name, 
                             $transcript_id, $transcript_strand, $transcript_biotype, $transcript_status, $transcript_version, $transcript_name, 
                             $extras, $assembly_name, $coding_types, $write_seleno);
      }
      elsif ($format eq "gff3"){
        write_transcript_gff3($gtffp, $slice, $gene, $transcript, $seqname, $sourcename,
                              $gene_id, $gene_biotype, $gene_status, $gene_version, $gene_name, 
                              $transcript_id, $transcript_strand, $transcript_biotype, $transcript_status, $transcript_version, $transcript_name, 
                              $extras, $assembly_name, $coding_types, $write_seleno);
      }
      $printed_genes{$gene->stable_id} = 1;
    }
  }
}

$gtffp->close() or die "Problem closing $outfile.\n";


########   ########   ########   ########   ########





#sort the chromosome names
sub _sortbychrnum {

  my @awords = split /_/,$a;
  my @bwords = split /_/,$b;

  my $anum = $awords[0];
  my $bnum = $bwords[0];

  $anum =~ s/chr//;
  $bnum =~ s/chr//;

  if ($anum !~ /^[0-9]*$/) {
    if ($bnum !~ /^[0-9]*$/) {
      return $anum cmp $bnum;
    } else {
      return 1;
    }
  }
  if ($bnum !~ /^[0-9]*$/) {
    return -1;
  }

  if ($anum <=> $bnum) {
    return $anum <=> $bnum;
  } else {
    if ($#awords == 0) {
      return -1;
    } elsif ($#bwords == 0) {
      return 1;
    } else {
	return $awords[1] cmp $bwords[1];
      }
  }
}

#get ordered array of chromosome names
sub _get_all_slice_names {
  my ($slice_adaptor, $coordsys, $version, $remove_NTs) = @_;

  my @slices = ();

  foreach my $chr (@{$slice_adaptor->fetch_all($coordsys, $version, undef, undef, undef)}){
    if($remove_NTs && ($chr->seq_region_name =~ /_NT_/)){ next; }
    push @slices, $chr->seq_region_name;
  }

  @slices = sort _sortbychrnum @slices;

  return \@slices;
}

#get gene/transcript name
sub get_name{
  my $object = shift;

  if (defined($object->display_xref) and $object->display_xref->display_id){
    return $object->display_xref->display_id;
  }
  elsif (defined($object->external_name)){
    return $object->external_name;
  }
  elsif ($object->get_all_Attributes('name')->[0]){
    return $object->get_all_Attributes('name')->[0]->value;
  }
  return _get_id($object);
}

