
# c_mps = speed of light (m/s)
set c_mps 299792458.
# c_mmps = speed of light (mm/s)
set c_mmps [expr $c_mps * 1000.]
# lamda = wavelength
# f = frequency
# lambda = c / f
# f = c / lambda
set band_list [list vhf1 vhf2 uhf1 uhf2]
set band_low(vhf1) 54000000.
set band_low(vhf2) 174000000.
set band_low(uhf1) 470000000.
set band_low(uhf2) 698000000.
set band_high(vhf1) 88000000.
set band_high(vhf2) 216000000.
set band_high(uhf1) 692000000.
# channel 51 is deprecated in US. 692 to 698mhz
set band_high(uhf2) 890000000.

# set up wavelength ranges
set ref_first 0
set ref $ref_first
foreach band $band_list {
    set power_of_two 1
    set power_of_2_arr($ref) 0
    # propose
    set lambda_low [expr { $c_mmps / ( $band_low($band) * ( $power_of_two * 1.) ) } ]
    set lambda_high [expr { $c_mmps / ( $band_high($band) * ( $power_of_two * 1.) ) } ]
    while { $lambda_low > 1 && $lambda_high > 1 } {
        set lambda_lowf($ref) $lambda_low
        set lambda_highf($ref) $lambda_high
        set band_arr($ref) $band
        set power_of_2_arr($ref) $power_of_two
        # create new proposal
        incr power_of_two
        incr ref
        set lambda_low [expr { $c_mmps / ( $band_low($band) * ( $power_of_two * 1.) ) } ]
        set lambda_high [expr { $c_mmps / ( $band_high($band) * ( $power_of_two * 1.) ) } ]
        #	puts "band $band ref $ref 2^${power_of_two} lambda_low $lambda_lowf($ref) lambda_high $lambda_highf($ref)"
        #	puts -nonewline "."
    }
}
# ref_last is one more than last reference, for faster looping checks
set ref_last $ref
# code historical reference B
# check band ranges against an interval of theoretical antenna lengths (ITAL)

set ital_mm_list [list 310.0 320.0 130.0 110.0 100.0 90.0 50.0 6.0 25.0 65.0 11.0 28.5 2.44 1.9 0.51 0.64 1.63 2.05 2.375]
#set ital_chc_list [list 9 9 9 9 9 9 9 9 9 L-bar L-bar hangar tube 24awg 22awg 14awg 12awg balunloop1]
# set ITAL range:
set ital_mm_first 0.1
set ital_mm_last 980.
set ital_mm_interval 0.1
puts "check band ranges against ITAL from ${ital_mm_first} to ${ital_mm_last}mm"
# convert to natural number referencing for arrays
set ital_first 0
set ital_last [expr { int( ( $ital_mm_last - $ital_mm_first ) / $ital_mm_interval ) } ]
incr ital_last
set ital_interval 1
set report_lists [list ]

for {set ital_i $ital_first} {$ital_i < $ital_last} {incr ital_i $ital_interval } {
    # convert to ital length (mm)
    set ital_mm_i [expr { ( $ital_i * 1. ) * $ital_mm_interval + $ital_mm_first } ]

    #puts "ital_i $ital_i ital_mm_i $ital_mm_i"
    set ital_hits 0
    foreach band_i $band_list {
        set ital_hits_arr($band_i) 0
    }

    set ital_band_stats_lists [list ]
    set ital_q_low_list [list ]
    set ital_q_high_list [list ]
    set ital_pot_list [list ]
    set ital_band_list [list ]

    # bands and fractional bands are in range of ref_first to ref_last
    # check bands against antenna length ital_mm_i
    for {set band_i $ref_first} { $band_i < $ref_last} {incr band_i} {

        # Is ital_mm_i in band range?
        if { $ital_mm_i > $lambda_highf($band_i) && $ital_mm_i < $lambda_lowf($band_i) } {
            # record hit
            incr ital_hits
            incr ital_hits_arr($band_arr($band_i))
            # quantify harmonic vs. boundary conditions 
            # collect %diff of band high/low vs ital for statistical analysis
            # ital_low is percent diff between band_low and antenna length
            set ital_low [expr { ( $lambda_lowf($band_i) - $ital_mm_i ) / $ital_mm_i } ]
            lappend ital_q_low_list $ital_low
            # ital_high is percent diff between band_high and antenna length
            set ital_high [expr { ( $ital_mm_i - $lambda_highf($band_i)  ) / $ital_mm_i } ]
            lappend ital_q_high_list $ital_high
            # record power of 2 for geometric averaging
            lappend ital_pot_list $power_of_2_arr($band_i)
            # record band
            lappend ital_band_list $band_arr($band_i)
        }
        # output something to show process status
        #puts "ital_i $ital_i band $band_arr($band_i) power_of_two $power_of_2_arr($band_i) hit: $ital_hits "
    }

    if { $ital_hits > 0 } {
        # calculate statistics
        set ital_q_low 0.
        set ital_q_high 0.
        # shift factors into reverse to maximize significance with lower power of two factor
        set ital_pot_max [lindex [lsort -integer $ital_pot_list ] end]
        set ital_pot_sum 0
        # for geometric averaging based on power of two factors (pot)
        set ital_sum_factors_k 0.
        foreach ital_pot_i $ital_pot_list {
            set ital_sum_factors_k [expr { $ital_sum_factors_k + $ital_pot_i } ]
        }
        for { set i 0 } {$i < $ital_hits} {incr i} {
            set ital_pot_i [lindex $ital_pot_list $i]
            set ital_q_low [lindex $ital_q_low_list $i]
            set ital_q_high [lindex $ital_q_high_list $i]
            set ital_band_i [lindex $ital_band_list $i]
            set ital_q_low_adjusted [expr { ( $ital_pot_max - $ital_pot_i) * $ital_q_low / $ital_sum_factors_k } ]
            set ital_q_high_adjusted [expr { ( $ital_pot_max - $ital_pot_i) * $ital_q_high / $ital_sum_factors_k } ]
            # record statistics for this band
            set band_stats_list [list $ital_band_i $ital_pot_i [format "%1.3f" $ital_q_low_adjusted] [format "%1.3f" $ital_q_high_adjusted] $ital_hits_arr(${ital_band_i})]
            lappend ital_band_stats_lists $band_stats_list
        }

        # analyze statistics for this band via sorting
        # lambda_HF is smaller wavelength
        # lambda_LF is larger wavelength
        # ital_mm test is positive if lambda_hf < ital_mm < lambda_lf
        
        # In general, antenna length should be less than wavelength,
        # such as 1/2 or 1/4 are common, with reason: resonance
        
        # suggesting that ital s/b < lambda_LF, so lets use this for 1st ranking criteria
        # where %diff lambda_lf-ital should approach 0..
        # how might this be expressed in terms of lamda_hf wavelength?
        # %diff lamda_hf-ital should approach some large value.. to show ital is bias toward longer waves.

        # another approach might be to rank best as cases where power_of_two for each band is close to same..
        set ital_band_stats_sorted_2_lists [lsort -index 2 -real -decreasing $ital_band_stats_lists ]
        # lets try by band for any insights
        set p2_uhf1 ""
        set lf_uhf1 ""
        set hf_uhf1 ""
        set p2_uhf2 ""
        set lf_uhf2 ""
        set hf_uhf2 ""
        set p2_vhf1 ""
        set lf_vhf1 ""
        set hf_vhf1 ""
        set p2_vhf2 ""
        set lf_vhf2 ""
        set hf_vhf2 ""
        set row_i 0
        foreach band_ii $band_list {
            set q_hits_${band_ii} 0
        }
        set rows_count [llength $ital_band_stats_lists]
        while { ( $lf_uhf1 eq "" || $hf_uhf1 eq "" || $lf_uhf2 eq "" || $hf_uhf2 eq "" || $lf_vhf1 eq "" || $hf_vhf1 eq "" || $lf_vhf2 eq "" || $hf_vhf2 eq "" ) && $row_i < $rows_count } {
            set row_list [lindex $ital_band_stats_sorted_2_lists $row_i]
            set band [lindex $row_list 0]
            set q_hits_${band} [lindex $row_list 4]
            switch -- $band {
                uhf1 { 
                    if { $lf_uhf1 eq "" && $hf_uhf1 eq "" } {
                        set p2_uhf1 [lindex $row_list 1]
                        set lf_uhf1 [lindex $row_list 2]
                        set hf_uhf1 [lindex $row_list 3]
                    }
                }
                uhf2 { 
                    if { $lf_uhf2 eq "" && $hf_uhf2 eq "" } {
                        set p2_uhf2 [lindex $row_list 1]
                        set lf_uhf2 [lindex $row_list 2]
                        set hf_uhf2 [lindex $row_list 3]
                    }
                }
                vhf1 { 
                    if { $lf_vhf1 eq "" && $hf_vhf1 eq "" } {
                        set p2_vhf1 [lindex $row_list 1]
                        set lf_vhf1 [lindex $row_list 2]
                        set hf_vhf1 [lindex $row_list 3]
                    }
                }
                vhf2 {
                    if { $lf_vhf2 eq "" && $hf_vhf2 eq "" } {
                        set p2_vhf2 [lindex $row_list 1]
                        set lf_vhf2 [lindex $row_list 2]
                        set hf_vhf2 [lindex $row_list 3]
                    }
                }
                default { 
                    puts "Error (206): Band '$band' does not exist. This shouldn't happen."
                }
            }
            incr row_i
        }
        # Code moved to end for historical reference (A)
        # create report table, one row per ital
        if { [lsearch -real -exact $ital_mm_list $ital_mm_i] > -1 } {
            set tested_p "* "
        } else {
            set tested_p "  "
        }
        set $ital_mm_i [format "%7.1f" $ital_mm_i]

        # Let's add a couple of numbers to determine relative band responsiveness
        # using principles from linear regression and least squares
        
        # using 3 instead of 4, because uhf2 is deprecated
        set hits_avg_3pt [expr { $ital_hits / 3. } ]
        set sigma_d2 [expr { sqrt( pow( $q_hits_vhf1 - $hits_avg_3pt , 2.) +  pow( $q_hits_vhf2 - $hits_avg_3pt , 2.) + pow( $q_hits_uhf1 - $hits_avg_3pt , 2.) ) / $hits_avg_3pt } ]
        set sigma_d2 [format "%7.5f" $sigma_d2]

        set report_row_list [list $tested_p $ital_i $ital_mm_i $ital_hits $sigma_d2 $q_hits_vhf1 $p2_vhf1 $lf_vhf1 $hf_vhf1 $q_hits_vhf2 $p2_vhf2 $lf_vhf2 $hf_vhf2 $q_hits_uhf1 $p2_uhf1 $lf_uhf1 $hf_uhf1 $q_hits_uhf2 $p2_uhf2 $lf_uhf2 $hf_uhf2 ]
        #puts $report_row_list
        lappend report_lists $report_row_list
    }
}

# report statistics
# tested_p means antenna has test results already in system to check against.
set report_titles_list [list *_p ital_i ital_mm q_hits sigma_d2 qh_vhf1 p2_vhf1 lf_vhf1 hf_vhf1 qh_vhf2 p2_vhf2 lf_vhf2 hf_vhf2 qh_uhf1 p2_uhf1 lf_uhf1 hf_uhf1 qh_uhf2 p2_uhf2 lf_uhf2 hf_uhf2 ]
set report_lists [linsert $report_lists 0 $report_titles_list]

set formatStr2 { % 3s %  7s  % 8s   % 6s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s  % 7s  % 7s}
set formatStr3 { % 3s %  7d  % 8.3f   % 6d   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s   % 7s  % 7s  % 7s}
set row_i 0

# write a copy to file
set fid [open "antenna-calc-results-${ital_mm_first}-${ital_mm_last}.txt" w]

foreach row_list $report_lists {
    set tested_p [lindex $row_list 0]
    set ital_i [lindex $row_list 1]
    set ital_mm [lindex $row_list 2]
    set q_hits [lindex $row_list 3]
    set sigma_d2 [lindex $row_list 4]
    set qh_vhf1 [lindex $row_list 5]
    set p2_vhf1 [lindex $row_list 6]
    set lf_vhf1 [lindex $row_list 7]
    set hf_vhf1 [lindex $row_list 8]
    set qh_vhf2 [lindex $row_list 9]
    set p2_vhf2 [lindex $row_list 10]
    set lf_vhf2 [lindex $row_list 11]
    set hf_vhf2 [lindex $row_list 12]
    set qh_uhf1 [lindex $row_list 13]
    set p2_uhf1 [lindex $row_list 14]
    set lf_uhf1 [lindex $row_list 15]
    set hf_uhf1 [lindex $row_list 16]
    set qh_uhf2 [lindex $row_list 17]
    set p2_uhf2 [lindex $row_list 18]
    set lf_uhf2 [lindex $row_list 19]
    set hf_uhf2 [lindex $row_list 20]

    if { $row_i eq 0 } {
        puts [format $formatStr2 $tested_p $ital_i $ital_mm $q_hits $sigma_d2 $qh_vhf1 $p2_vhf1 $lf_vhf1 $hf_vhf1 $qh_vhf2 $p2_vhf2 $lf_vhf2 $hf_vhf2 $qh_uhf1 $p2_uhf1 $lf_uhf1 $hf_uhf1 $qh_uhf2 $p2_uhf2 $lf_uhf2 $hf_uhf2 ]
        puts $fid [format $formatStr2 $tested_p $ital_i $ital_mm $q_hits $sigma_d2 $qh_vhf1 $p2_vhf1 $lf_vhf1 $hf_vhf1 $qh_vhf2 $p2_vhf2 $lf_vhf2 $hf_vhf2 $qh_uhf1 $p2_uhf1 $lf_uhf1 $hf_uhf1 $qh_uhf2 $p2_uhf2 $lf_uhf2 $hf_uhf2 ]
    } else {
        puts [format $formatStr3 $tested_p $ital_i $ital_mm $q_hits $sigma_d2 $qh_vhf1 $p2_vhf1 $lf_vhf1 $hf_vhf1 $qh_vhf2 $p2_vhf2 $lf_vhf2 $hf_vhf2 $qh_uhf1 $p2_uhf1 $lf_uhf1 $hf_uhf1 $qh_uhf2 $p2_uhf2 $lf_uhf2 $hf_uhf2 ]
        puts $fid [format $formatStr3 $tested_p $ital_i $ital_mm $q_hits $sigma_d2 $qh_vhf1 $p2_vhf1 $lf_vhf1 $hf_vhf1 $qh_vhf2 $p2_vhf2 $lf_vhf2 $hf_vhf2 $qh_uhf1 $p2_uhf1 $lf_uhf1 $hf_uhf1 $qh_uhf2 $p2_uhf2 $lf_uhf2 $hf_uhf2 ]
    }
    incr row_i
}
close $fid

# reference A historical code
if { 0 } {
    if { [llength $ital_band_stats_lists] > 1 } {
        # add data from next best ranking
        set ital_q_high_3_list [lindex $ital_band_stats_sorted_3_lists end-1]
        append ital_q_high_adj "([lindex $ital_q_high_3_list 3])"
        append ital_q_high_adj_band "([lindex $ital_q_high_3_list 0]:[lindex $ital_q_high_3_list 1])"
        set ital_q_low_2_list [lindex $ital_band_stats_sorted_2_lists end-1]
        append ital_q_low_adj "([lindex $ital_q_low_2_list 2])"
        append ital_q_low_adj_band "([lindex $ital_q_low_2_list 0]:[lindex $ital_q_low_2_list 1])"
        if { [llength $ital_band_stats_lists] > 2 } {
            # add data from next best ranking
            set ital_q_high_3_list [lindex $ital_band_stats_sorted_3_lists end-2]
            append ital_q_high_adj "([lindex $ital_q_high_3_list 3])"
            append ital_q_high_adj_band "([lindex $ital_q_high_3_list 0]:[lindex $ital_q_high_3_list 1])"
            set ital_q_low_2_list [lindex $ital_band_stats_sorted_2_lists end-2]
            append ital_q_low_adj "([lindex $ital_q_low_2_list 2])"
            append ital_q_low_adj_band "([lindex $ital_q_low_2_list 0]:[lindex $ital_q_low_2_list 1])"
            
        }
    }
}
# code historical reference B 
if { 0 } {
    # check ranges against length of known, optimized antennas by experiment (loae) in mm and channel count (chc).
    puts "check ranges against LOAE"
    set loae_mm_list [list 310 320 130 110 100 90 50 6 25 65 11 28.542]
    set loae_chc_list [list 9 9 9 9 9 9 9 9 9 L-bar L-bar hangar]
    set ant_count [llength $loae_mm_list]
    set ant_num 0
    set ant_i 0
    foreach loae_i $loae_mm_list {
        puts "loae_i $loae_i"
        set loae_hits_arr($loae_i) 0
        set loae_last_ref($loae_i) 0
        puts "ref_first $ref_first ref_last $ref_last"
        set last_band_ref ""
        for {set i $ref_first} { $i < $ref_last} {incr i} {
            if { $loae_i > $lambda_highf($i) && $loae_i < $lambda_lowf($i) } {
                if { $loae_hits_arr($loae_i) < 1 || $last_band_ref ne $band_arr($i)} {
                    puts "loae_i $loae_i loae_chc [lindex $loae_chc_list $ant_i] hit: $loae_hits_arr($loae_i) last_ref $loae_last_ref($loae_i) band $band_arr($i) power_of_two $power_of_2_arr($i) lambda_low $lambda_lowf($i) lambda_high $lambda_highf($i)"
                    set last_band_ref $band_arr($i)
                }
                set loae_last_ref($loae_i) $i
                incr loae_hits_arr($loae_i)
            }
        }
        puts "loae_i $loae_i loae_chc [lindex $loae_chc_list $ant_i] hit: $loae_hits_arr($loae_i) last_ref $loae_last_ref($loae_i) band $band_arr($loae_last_ref($loae_i)) power_of_two $power_of_2_arr($loae_last_ref($loae_i)) lambda_low $lambda_lowf($loae_last_ref($loae_i)) lambda_high $lambda_highf($loae_last_ref($loae_i))"
        puts ""
        incr ant_i
    }
}
