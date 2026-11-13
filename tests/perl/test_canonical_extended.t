use strict;
use warnings;
use Test::More;
use File::Temp qw(tempfile);
use Storable qw(dclone);
use Scalar::Util qw(refaddr);

my $script = shift or die 'Expected uploader script path';
open my $input, '<', $script or die $!;
my $text = do {local $/; <$input>};
my ($helpers) = $text =~ /# BEGIN CANONICAL_EXTENDED_HELPERS\n(.*)\z/s;
die 'Missing helper block' unless $helpers;
# Load the actual uploader helpers and editing routine, without executing its
# database connection or main program. All database-facing objects below are mocks.
eval $helpers; die $@ if $@;

{
  package MockAttribute;
  sub new { my ($class, $value, $code) = @_; bless {value=>$value, code=>$code || 'remark'}, $class }
  sub value { $_[0]{value} }
  sub code { $_[0]{code} }
  package Bio::EnsEMBL::Attribute;
  sub new { my ($class, %args) = @_; MockAttribute->new($args{-value}, $args{-code}) }
  package MockExon;
  sub new {my ($class, $start, $end, $id) = @_; bless {start=>$start,end=>$end,id=>$id,phase=>0,end_phase=>0,version=>1}, $class}
  sub start {$_[0]{start}=$_[1] if @_>1; $_[0]{start}}
  sub end {$_[0]{end}=$_[1] if @_>1; $_[0]{end}}
  sub phase {$_[0]{phase}}
  sub end_phase {$_[0]{end_phase}}
  sub stable_id {$_[0]{id}=$_[1] if @_>1; $_[0]{id}}
  sub version {$_[0]{version}=$_[1] if @_>1; $_[0]{version}}
  sub new_dissociated_copy {Storable::dclone($_[0])}
  sub vega_hashkey {join(':', $_[0]{id}, $_[0]{start}, $_[0]{end})}
  package MockTranslation;
  sub start_Exon {$_[0]{start_exon}}
  sub end_Exon {$_[0]{end_exon}}
  sub start {$_[0]{start}=$_[1] if @_>1; $_[0]{start}}
  sub end {$_[0]{end}=$_[1] if @_>1; $_[0]{end}}
  sub seq {'PEPTIDE'}
  sub stable_id {$_[0]{id}=$_[1] if @_>1; $_[0]{id}}
  sub version {$_[0]{version}=$_[1] if @_>1; $_[0]{version}}
  package MockTranscript;
  sub new {
    my ($class, $id, $start, $end, $strand, $remark) = @_;
    my @exons = (MockExon->new($start, 150, $id.'e1'), MockExon->new(200, 250, $id.'e2'), MockExon->new(300, $end, $id.'e3'));
    @exons = reverse @exons if $strand == -1;
    my $translation = bless {start_exon=>$exons[0], end_exon=>$exons[-1],
      start=>$strand == 1 ? 121-$start : $end-320+1,
      end=>$strand == 1 ? 21 : 150-120+1, id=>'protein_'.$id, version=>1}, 'MockTranslation';
    bless {id=>$id,exons=>\@exons,strand=>$strand,chr=>'12',offset=>0,
      attributes=>[$remark ? MockAttribute->new($remark) : (), MockAttribute->new($id, 'name')],
      translation=>$translation, dbID=>1, version=>1}, $class;
  }
  sub stable_id {$_[0]{id}=$_[1] if @_>1; $_[0]{id}}
  sub version {$_[0]{version}=$_[1] if @_>1; $_[0]{version}}
  sub dbID {$_[0]{dbID}}
  sub start {my @x = sort {$a->start <=> $b->start} @{$_[0]{exons}}; $x[0]->start}
  sub end {my @x = sort {$a->end <=> $b->end} @{$_[0]{exons}}; $x[-1]->end}
  sub seq_region_start {$_[0]->start + $_[0]{offset}}
  sub seq_region_end {$_[0]->end + $_[0]{offset}}
  sub seq_region_name {$_[0]{chr}}
  sub strand {$_[0]{strand}}
  sub get_all_Exons {$_[0]{exons}}
  sub start_Exon {$_[0]{exons}[0]}
  sub end_Exon {$_[0]{exons}[-1]}
  sub translation {$_[0]{translation}}
  sub coding_region_start {
    my $t = $_[0]{translation};
    $_[0]{strand} == 1 ? $t->start_Exon->start + $t->start - 1 : $t->end_Exon->end - $t->end + 1;
  }
  sub coding_region_end {
    my $t = $_[0]{translation};
    $_[0]{strand} == 1 ? $t->end_Exon->start + $t->end - 1 : $t->start_Exon->end - $t->start + 1;
  }
  sub get_all_Attributes {my ($self, $code)=@_; [grep {$_->code eq $code} @{$self->{attributes}}]}
  sub add_Attributes {my $self=shift; push @{$self->{attributes}}, @_}
  sub new_dissociated_copy {my $new=Storable::dclone($_[0]); $new->{dbID}=undef; $new}
  sub swap_exons {
    my ($self, $old, $new) = @_;
    @{$self->{exons}} = map {Scalar::Util::refaddr($_) == Scalar::Util::refaddr($old) ? $new : $_} @{$self->{exons}};
    for my $key ('start_exon','end_exon') {
      $self->{translation}{$key} = $new if Scalar::Util::refaddr($self->{translation}{$key}) == Scalar::Util::refaddr($old);
    }
  }
  package MockGene;
  sub new {my $class=shift; bless {transcripts=>[@_]}, $class}
  sub stable_id {'gene1'}
  sub strand {$_[0]{transcripts}[0]->strand}
  sub seq_region_name {'12'}
  sub get_all_Transcripts {$_[0]{transcripts}}
  sub add_Transcript {push @{$_[0]{transcripts}}, $_[1]}
  sub gene_author {bless {}, 'MockAuthor'}
  sub get_all_Attributes {[MockAttribute->new('GENE', 'name')]}
  package MockAuthor;
  sub name {'original_author'}
  package MockSlice;
  sub start {$_[0]{start}}
  package MockRegion;
  sub genes {$_[0]{gene}}
  sub slice {bless {start=>($_[0]{gene}{offset} || 0)+1}, 'MockSlice'}
  package Bio::Otter::Server::Support::Local;
  sub new {bless {}, shift}
  sub set_params {}
  sub authorized_user {}
  package Bio::Otter::ServerAction::Script::Region;
  sub new_with_slice {bless {}, shift}
  sub get_region {$main::test_region}
}

sub model {MockTranscript->new(@_)}
sub plan {
  my ($models, $desired) = @_;
  my %selected = map {$_=>1} keys %$desired;
  leap_plan_updates(MockGene->new(@$models), \%selected, $desired);
}
my $m = model('mane',100,350,1,'MANE_select');
my $c = model('canon',90,370,1,'Ensembl_canonical_extended');
subtest 'same or already longer canonical: no new model' => sub {
  for my $end (350,370){
    my $plans = plan([$m,$c], {mane=>{start=>90,end=>$end}});
    is(scalar @$plans,1);
    is($plans->[0]{action},'already_covered');
    is($plans->[0]{target_id},'canon');
    is_deeply([$plans->[0]{start},$plans->[0]{end}],[90,370]);
  }
};
subtest 'new candidate extends existing model; farther endpoints retained' => sub {
  my $plans=plan([$m,$c], {mane=>{start=>80,end=>360}});
  is($plans->[0]{action},'update_canonical_extended');
  is_deeply([$plans->[0]{start},$plans->[0]{end}],[80,370]);
};
subtest 'canonical own input and two MANE sources coalesce, independent of order' => sub {
  my $plus=model('plus',100,350,1,'MANE_plus_clinical');
  for my $models ([$m,$c,$plus],[$plus,$c,$m]){
    my $plans=plan($models,{mane=>{start=>80,end=>350},plus=>{start=>100,end=>390},canon=>{start=>90,end=>375}});
    is(scalar @$plans,1);
    is_deeply([$plans->[0]{start},$plans->[0]{end}],[80,390]);
    is_deeply([sort @{$plans->[0]{source_ids}}],[qw(canon mane plus)]);
  }
};
subtest 'no matching canonical: prepare one copy for equivalent MANE models' => sub {
  my $plus=model('plus',100,350,1,'MANE_plus_clinical');
  my $plans=plan([$m,$plus],{mane=>{start=>80,end=>350},plus=>{start=>100,end=>390}});
  is(scalar @$plans,1);
  is($plans->[0]{action},'create_canonical_extended');
  is_deeply([$plans->[0]{start},$plans->[0]{end}],[80,390]);
};
subtest 'protected structure includes internal exons, every splice-facing boundary and CDS' => sub {
  for my $mutation (
    sub {$_[0]{exons}[0]->end(149)}, sub {$_[0]{exons}[1]->start(201)},
    sub {$_[0]{exons}[1]->end(249)}, sub {$_[0]{exons}[2]->start(301)},
    sub {$_[0]{translation}->start($_[0]{translation}->start+1)},
    sub {$_[0]{translation}->end($_[0]{translation}->end-1)},
    sub {$_[0]{exons}[1]{phase}=1}, sub {$_[0]{exons}[1]{end_phase}=1},
    sub {$_[0]{chr}='1'}, sub {$_[0]{strand}=-1}, sub {splice @{$_[0]{exons}},1,1}
  ){
    my $other=Storable::dclone($c); $mutation->($other);
    my $plans=plan([$m,$other],{mane=>{start=>80,end=>380}});
    is($plans->[0]{action},'create_canonical_extended');
  }
};
subtest 'multiple matching existing models stop planning before mutations' => sub {
  my $other=model('other',80,380,1,'Ensembl_canonical_extended');
  eval {plan([$m,$c,$other],{mane=>{start=>80,end=>390}})};
  like($@,qr/Multiple matching canonical extended.*canon,other/);
  is($c->seq_region_start,90);
};
subtest 'lowercase canonical tag accepted and repeated remarks avoided' => sub {
  my $lower=model('lower',90,370,1,'ensembl_canonical_extended');
  is(plan([$m,$lower],{mane=>{start=>90,end=>370}})->[0]{action},'already_covered');
  leap_add_remark_once($lower,MockAttribute->new('LEAP_5')) for 1..3;
  is(scalar(grep {$_->value eq 'LEAP_5'} @{$lower->get_all_Attributes('remark')}),1);
};
subtest 'reject source shortening and invalid proposed coordinates' => sub {
  for my $coords ({start=>101,end=>380},{start=>80,end=>349},{start=>0,end=>380},{start=>'x',end=>380}){
    eval {plan([$m,$c],{mane=>$coords})};
    like($@,qr/Proposed endpoints would shorten|Missing or invalid/);
  }
};
subtest 'genomic structure independent of slice offset' => sub {
  my $other=Storable::dclone($c);
  $_->start($_->start-50), $_->end($_->end-50) for @{$other->{exons}};
  $other->{offset}=50;
  is(leap_structure_key($other),leap_structure_key($c));
};

# Bind the real process_gene routine to local mock services and edit objects.
my (%transcripts_to_update,%transcript_coords,%new_gene_coords,%gene_update_plans);
my ($otter_dba,$dataset_name,$assembly_version,$write)=(undef,'test','GRCh38',1);
my $leap_5_att=MockAttribute->new('LEAP_5');
my $leap_3_att=MockAttribute->new('LEAP_3');
my $canon_ext_att=MockAttribute->new('Ensembl_canonical_extended');
our ($test_region,$write_count);
my ($edits,$edits_path)=tempfile();
open ED,'>', $edits_path or die $!;
my ($process)=$text =~ /(sub process_gene \{.*?)(?=\n#Given a region)/s;
eval $process; die $@ if $@;
my ($exon_helpers)=$text =~ /(sub exon_in_mane_transcript \{.*?)(?=\n# BEGIN CANONICAL_EXTENDED_HELPERS)/s;
eval $exon_helpers; die $@ if $@;
sub get_new_transcript_name {'GENE-new-copy'}
sub write_gene_region {
  $write_count++;
  # Simulate ID assignment only; never use a real database.
  my $n=0;
  for my $tr (@{$test_region->{gene}->get_all_Transcripts}){
    $tr->stable_id('assigned_'.++$n) unless $tr->stable_id;
  }
  return 'lock ok,write ok,unlock ok';
}
sub edit {
  my ($gene,$desired)=@_;
  %transcript_coords=%$desired;
  %transcripts_to_update=map {$_=>1} keys %$desired;
  %new_gene_coords=(gene1=>{start=>1,end=>1000});
  %gene_update_plans=();
  $test_region=bless {gene=>$gene},'MockRegion';
  $write_count=0;
  process_gene($gene);
  my @findings=leap_check_uploaded_gene($gene,$gene_update_plans{gene1});
  ok(!grep({$_->[1] =~ /ERROR|wrongly/} @findings),'post-upload checks agree');
}
subtest 'actual editor: reuse both strands at either or both ends, preserve MANE and CDS' => sub {
  for my $strand (1,-1){
    for my $new ([80,370],[90,390],[80,390]){
      my $mane=model('mane',100,350,$strand,'MANE_select');
      my $canonical=model('canon',90,370,$strand,'Ensembl_canonical_extended');
      my $gene=MockGene->new($mane,$canonical);
      my $original_key=leap_structure_key($canonical);
      edit($gene,{mane=>{start=>$new->[0],end=>$new->[1]}});
      is(scalar @{$gene->get_all_Transcripts},2,'no new transcript');
      is($canonical->stable_id,'canon','canonical ID retained');
      is($canonical->translation->stable_id,'protein_canon','protein ID retained');
      is_deeply([$canonical->seq_region_start,$canonical->seq_region_end],$new);
      is_deeply([$mane->seq_region_start,$mane->seq_region_end],[100,350]);
      is(leap_structure_key($canonical),$original_key,'CDS and protected exons preserved');
      is($write_count,1,'one gene write');
      # The next upload of the same MANE candidate must be a no-op.
      edit($gene,{mane=>{start=>$new->[0],end=>$new->[1]}});
      is($write_count,0,'repeated upload skips existing matching model');
      is(scalar @{$gene->get_all_Transcripts},2);
    }
  }
};
subtest 'actual editor: retain farther endpoint at opposite end' => sub {
  for my $strand (1,-1){
    my $mane=model('mane',100,350,$strand,'MANE_select');
    my $canonical=model('canon',90,390,$strand,'Ensembl_canonical_extended');
    my $gene=MockGene->new($mane,$canonical);
    edit($gene,{mane=>{start=>80,end=>370}});
    is_deeply([$canonical->seq_region_start,$canonical->seq_region_end],[80,390]);
  }
};
subtest 'actual editor: no match creates exactly one canonical copy, retaining original MANE' => sub {
  for my $strand (1,-1){
    my $mane=model('mane',100,350,$strand,'MANE_plus_clinical');
    my $gene=MockGene->new($mane);
    edit($gene,{mane=>{start=>80,end=>390}});
    is(scalar @{$gene->get_all_Transcripts},2);
    my ($copy)=grep {$_ != $mane} @{$gene->get_all_Transcripts};
    ok(leap_is_canonical_extended($copy));
    is(leap_mane_remark($copy),'');
    is(leap_mane_remark($mane),'MANE_plus_clinical');
    is($mane->stable_id,'mane');
    is($copy->translation->stable_id,'','new copy does not reuse protein ID');
    is_deeply([$copy->seq_region_start,$copy->seq_region_end],[80,390]);
    is_deeply([$mane->seq_region_start,$mane->seq_region_end],[100,350]);
    edit($gene,{mane=>{start=>80,end=>390}});
    is(scalar @{$gene->get_all_Transcripts},2,'newly created copy reused on repeat upload');
    is($write_count,0);
  }
};
subtest 'actual editor: canonical input cannot undo merged MANE extension' => sub {
  my $mane=model('mane',100,350,1,'MANE_select');
  my $canonical=model('canon',90,370,1,'Ensembl_canonical_extended');
  my $gene=MockGene->new($mane,$canonical);
  edit($gene,{canon=>{start=>90,end=>375},mane=>{start=>80,end=>390}});
  is_deeply([$canonical->seq_region_start,$canonical->seq_region_end],[80,390]);
  is($write_count,1);
};
subtest 'actual editor: unchanged covered models do not write, ordinary models update in place' => sub {
  my $mane=model('mane',100,350,1,'MANE_select');
  my $canonical=model('canon',90,370,1,'Ensembl_canonical_extended');
  edit(MockGene->new($mane,$canonical),{mane=>{start=>90,end=>370}});
  is($write_count,0);
  my $ordinary=model('ordinary',100,350,1,'unrelated_remark');
  my $gene=MockGene->new($ordinary);
  edit($gene,{ordinary=>{start=>80,end=>390}});
  is(scalar @{$gene->get_all_Transcripts},1);
  is($ordinary->stable_id,'ordinary');
  is_deeply([$ordinary->seq_region_start,$ordinary->seq_region_end],[80,390]);
};
subtest 'actual editor: single exon CDS remains fixed on both strands' => sub {
  for my $strand (1,-1){
    my @models=(model('mane',100,350,$strand,'MANE_select'),
                model('canon',90,370,$strand,'Ensembl_canonical_extended'));
    for my $tr (@models){
      my $exon=MockExon->new($tr->start,$tr->end,$tr->stable_id.'single');
      $tr->{exons}=[$exon];
      $tr->{translation}{start_exon}=$exon;
      $tr->{translation}{end_exon}=$exon;
      $tr->{translation}{start}=$strand == 1 ? 120-$exon->start+1 : $exon->end-320+1;
      $tr->{translation}{end}=$strand == 1 ? 320-$exon->start+1 : $exon->end-120+1;
    }
    my $gene=MockGene->new(@models);
    edit($gene,{mane=>{start=>80,end=>390}});
    is_deeply([$models[1]->coding_region_start,$models[1]->coding_region_end],[120,320]);
    is_deeply([$models[1]->seq_region_start,$models[1]->seq_region_end],[80,390]);
    is_deeply([$models[0]->seq_region_start,$models[0]->seq_region_end],[100,350]);
  }
};
subtest 'actual editor: slice offsets and shared exons do not alter original models' => sub {
  for my $strand (1,-1){
    my $mane=model('mane',100,350,$strand,'MANE_select');
    my $canonical=model('canon',100,350,$strand,'Ensembl_canonical_extended');
    # Share actual exon objects with MANE; the reused canonical must detach ends.
    $canonical->{exons}=[@{$mane->{exons}}];
    $canonical->{translation}{start_exon}=$canonical->{exons}[0];
    $canonical->{translation}{end_exon}=$canonical->{exons}[-1];
    my $gene=MockGene->new($mane,$canonical);
    for my $exon (@{$mane->{exons}}){
      $exon->start($exon->start-50); $exon->end($exon->end-50);
    }
    $_->{offset}=50 for ($mane,$canonical);
    $gene->{offset}=50;
    edit($gene,{mane=>{start=>80,end=>390}});
    is_deeply([$canonical->seq_region_start,$canonical->seq_region_end],[80,390]);
    is_deeply([$mane->seq_region_start,$mane->seq_region_end],[100,350]);
    is_deeply([$canonical->coding_region_start+50,$canonical->coding_region_end+50],[120,320]);
  }
};
subtest 'post-upload check detects changed or missing canonical and altered MANE' => sub {
  my $mane=model('mane',100,350,1,'MANE_select');
  my $canonical=model('canon',90,370,1,'Ensembl_canonical_extended');
  my $gene=MockGene->new($mane,$canonical);
  my $plans=leap_plan_updates($gene,{mane=>1},{mane=>{start=>90,end=>370}});
  $canonical->{exons}[0]->start(91);
  my @bad=leap_check_uploaded_gene($gene,$plans);
  ok(grep({$_->[1] =~ /ERROR_planned_upload_mismatch/} @bad));
  $gene->{transcripts}=[$mane];
  $mane->{exons}[0]->end(151);
  @bad=leap_check_uploaded_gene($gene,$plans);
  ok(grep({$_->[1] =~ /wrongly_updated_or_missing/} @bad));
  ok(grep({$_->[1] =~ /ERROR_planned_upload_mismatch/} @bad));
};
close ED;
unlink $edits_path;
done_testing();
