package NMDFrameshift;

use strict;
use warnings;
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use Bio::Seq;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub version { return '3.1'; }
sub feature_types { return ['Transcript']; }
sub variant_feature_types { return ['VariationFeature']; }

sub get_header_info {
    return {
        WildtypeCDS       => "Wildtype CDS (from start to native stop codon)",
        WildtypeProtein   => "Wildtype protein sequence (stop excluded)",
        MutantCDS         => "Mutant CDS (to first new in-frame stop, inclusive)",
        MutantProtein     => "Mutant protein sequence (truncated at new stop)",
        NMDEscape         => "NMD escape prediction: NMDER1 (start-proximal <150nt), NMDER2 (long exon >400nt from splice site), NMDER3 (distal: last exon or last 50nt of penultimate exon), NMD (no escape), NoPTC (no premature termination codon found)",
    };
}

sub run {
    my ($self, $tva) = @_;
    my @ocs = @{ $tva->get_all_OverlapConsequences || [] };

    return {} unless grep { $_->SO_term eq 'frameshift_variant' } @ocs;
    return {} if grep { $_->SO_term =~ /splice/ } @ocs;

    my $tv = $tva->transcript_variation;
    my $tr = $tv->transcript;

    my $codon_table = _get_codon_table($tr);

    my $wt_cds = defined($tr->{_variation_effect_feature_cache})
               ? $tr->{_variation_effect_feature_cache}->{translateable_seq}
               : $tr->translateable_seq;

    return {} unless defined $wt_cds && length($wt_cds) > 0;
    $wt_cds =~ s/\-//g;
    return {} if $wt_cds eq '';

    my $wt_seq_obj = Bio::Seq->new(-seq => $wt_cds, -moltype => 'dna');
    my $wt_protein = $wt_seq_obj->translate(undef, undef, 0, $codon_table)->seq();
    $wt_protein =~ s/\*.*//;

    my $cds_start = $tv->cds_start // 0;
    my $cds_end   = $tv->cds_end   // 0;

    my $cds_len = length($wt_cds);
    $cds_start = $cds_start < 1 ? 1 : $cds_start;
    $cds_end   = $cds_end > $cds_len ? $cds_len : $cds_end;

    my $head = ($cds_start > 1) ? substr($wt_cds, 0, $cds_start - 1) : '';
    my $tail = ($cds_end < $cds_len) ? substr($wt_cds, $cds_end) : '';

    my $mutation = $tva->feature_seq // '';
    $mutation =~ s/\-//g;

    my $mut_cds_core = $head . $mutation . $tail;
    return {} if !defined($mut_cds_core) || $mut_cds_core eq '';

    my $utr3_seq = "";
    if (my $utr3 = $tr->three_prime_utr) {
        $utr3_seq = ($utr3->seq() // '');
        $utr3_seq =~ s/\-//g;
    }

    my $full_mutant_dna = $mut_cds_core . $utr3_seq;

    my $mut_seq_obj = Bio::Seq->new(-seq => $full_mutant_dna, -moltype => 'dna');
    my $mut_protein_full = $mut_seq_obj->translate(undef, undef, 0, $codon_table)->seq();

    my $has_ptc = 0;
    my $ptc_aa_index = undef; # Store amino acid index instead
    my $mut_cds = $mut_cds_core;
    my $mut_protein = $mut_protein_full;
    $mut_protein =~ s/\*$//;

    if ($mut_protein_full =~ /\*/) {
        $ptc_aa_index = index($mut_protein_full, '*');
        $has_ptc = 1;
        $mut_cds = substr($full_mutant_dna, 0, $ptc_aa_index * 3 + 3);
        $mut_protein = substr($mut_protein_full, 0, $ptc_aa_index);
    }

    my $nmd_escape = "NoPTC";
    if ($has_ptc && defined $ptc_aa_index) {
        # Pass the full mutant DNA and AA index to the new function
        $nmd_escape = $self->_determine_nmd_escape_v2($tr, $full_mutant_dna, $ptc_aa_index, $codon_table);
    }

    return {
        WildtypeCDS       => $wt_cds,
        WildtypeProtein   => $wt_protein,
        MutantCDS         => $mut_cds,
        MutantProtein     => $mut_protein,
        NMDEscape         => $nmd_escape,
    };
}


sub _determine_nmd_escape_v2 {
    my ($self, $tr, $full_mutant_dna, $ptc_aa_index, $codon_table) = @_;

    my $cds_genomic_start;
    my $cds_genomic_end;
    my $strand = $tr->strand;
    
    if ($strand == 1) {
        $cds_genomic_start = $tr->coding_region_start;
        $cds_genomic_end = $tr->coding_region_end;
    } else {
        $cds_genomic_start = $tr->coding_region_end;
        $cds_genomic_end = $tr->coding_region_start;
    }

    my $target_nt_pos_in_transcript = $ptc_aa_index * 3 + 1; # 1-based position in the artificial transcript sequence

    my @exons = sort { $strand == 1 ? $a->start <=> $b->start : $b->start <=> $a->start } 
                @{ $tr->get_all_Exons };

    my $current_transcript_pos = 1;
    my $ptc_genomic_pos = undef;
    my $ptc_exon = undef;
    my $ptc_exon_index = -1;

    foreach my $i (0..$#exons) {
        my $exon = $exons[$i];
        my $exon_length = $exon->length;

        my $exon_start_in_transcript = $current_transcript_pos;
        my $exon_end_in_transcript = $current_transcript_pos + $exon_length - 1;

        # Check if the PTC falls within this exon
        if ($target_nt_pos_in_transcript >= $exon_start_in_transcript && 
            $target_nt_pos_in_transcript <= $exon_end_in_transcript) {
            
            $ptc_exon = $exon;
            $ptc_exon_index = $i;
            
            # Calculate the exact genomic position of the PTC
            my $offset_within_exon = $target_nt_pos_in_transcript - $exon_start_in_transcript;
            if ($strand == 1) {
                $ptc_genomic_pos = $exon->start + $offset_within_exon;
            } else {
                $ptc_genomic_pos = $exon->end - $offset_within_exon;
            }
            last;
        }

        $current_transcript_pos += $exon_length;
    }

    return "NMD" unless defined $ptc_genomic_pos;

    my $num_exons = scalar @exons;

    # Rule 1: Start-proximal PTC (<150 nt from CDS start)
    if ($ptc_aa_index < 50) {
        return "NMDER1";
    }

    # Rule 3: Distal region
    # (a) PTC in last exon → NMD escape
    if ($ptc_exon_index == $num_exons - 1) {
        return "NMDER3";
    }
    
    # (b) PTC in penultimate exon AND within last 50 nt of that exon
    if ($num_exons >= 2 && $ptc_exon_index == $num_exons - 2) {
        my $dist_from_end = $ptc_exon->length - ($ptc_genomic_pos - $ptc_exon->start) - 1;
        if ($strand == -1) {
            $dist_from_end = $ptc_genomic_pos - $ptc_exon->start;
        }
        if ($dist_from_end <= 50) {
            return "NMDER3";
        }
    }

    # Rule 2: Long exon with PTC >400 nt from splice donor
    if ($ptc_exon->length > 400) {
        my $dist_to_splice_donor = $ptc_exon->length - ($ptc_genomic_pos - $ptc_exon->start) - 1;
        if ($strand == -1) {
            $dist_to_splice_donor = $ptc_genomic_pos - $ptc_exon->start;
        }
        if ($dist_to_splice_donor > 400) {
            return "NMDER2";
        }
    }

    return "NMD";
}

sub _get_codon_table {
    my $tr = shift;
    my $codon_table = 1;
    if (defined($tr->{_variation_effect_feature_cache})) {
        $codon_table = $tr->{_variation_effect_feature_cache}->{codon_table} || 1;
    } else {
        my ($attrib) = @{$tr->slice->get_all_Attributes('codon_table')};
        $codon_table = $attrib ? $attrib->value || 1 : 1;
    }
    return $codon_table;
}

1;