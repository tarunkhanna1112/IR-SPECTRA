proc gaussian {} {
	set f [open "template.log" "r"]
	set data [read $f]
	close $f

	set n_atoms [input1]
	set n_modes [expr { (3 * $n_atoms) - 6 }]

	 #  Condition to sort the file
	set nmode1 5.1700

	set g [open "vecs_gaussian" "w"]
	puts $g "Atom	AN		X			Y			Z"
	set k -1

	puts "[llength $data]"
	while { $k < [llength $data] } {
		set k1 $k
		if { [lindex $data $k] == $nmode1 } {
			for {set nm 1} {$nm <= $n_modes} {incr nm 3} {
				puts $g "frequency = [lindex $data $k]"
				while { [lindex $data $k1] != 1 } {
					incr k1
				}

				if { [lindex $data $k1] == 1 } {
					for {set i $k1} {$i < [expr { $k1 + (11 * $n_atoms) } ]} {incr i 11} {	
						puts $g "[lindex $data $i]	[lindex $data [expr { $i + 1 }]]		[lindex $data [expr { $i + 2 }]]		[lindex $data [expr { $i + 3 }]]		[lindex $data [expr { $i + 4 }]]"
					}
				}

				puts $g ""
				puts $g "frequency = [lindex $data [expr { $k + 1 } ]]"
				for {set i $k1} {$i < [expr { $k1 + (11 * $n_atoms) } ]} {incr i 11} {	
					puts $g "[lindex $data $i]	[lindex $data [expr { $i + 1 }]]		[lindex $data [expr { $i + 5 }]]		[lindex $data [expr { $i + 6 }]]		[lindex $data [expr { $i + 7 }]]"
				}

				puts $g ""
				puts $g "frequency = [lindex $data [expr { $k + 2 } ]]"
				for {set i $k1} {$i < [expr { $k1 + (11 * $n_atoms) } ]} {incr i 11} {	
					puts $g "[lindex $data $i]	[lindex $data [expr { $i + 1 }]]		[lindex $data [expr { $i + 8 }]]		[lindex $data [expr { $i + 9 }]]		[lindex $data [expr { $i + 10 }]]"
				}
			puts $g ""
			set k [expr { $i + 8 }]  
			set k1 $k
			}
			set k [llength $data]
		
		}
		incr k
	}
	close $g	
}

proc amber {} {
	set n_atoms [input1]
	set n_modes [expr { (3 * $n_atoms) - 6 }]
	set f [open "[input2]" "r"]
	set data [read $f]
	close $f
	set g [open "vecs_amber" "w"]
	puts $g "ATOM NUMBER		X			Y			Z"
	set k 0
	set tmo 0
	while { $k != [llength $data] } {
		if { [lindex $data $k] == "****" } {
			incr tmo
			if { $tmo > 6 } {
				puts $g "freq = [lindex $data [expr { $k + 2 }]]"
				set j 1
				for {set i [expr { $k + 3 }]} {$i <= [expr { (3 *$n_atoms) + ($k + 2) }]} {incr i 3} {
					puts $g "$j		[lindex $data $i]		[lindex $data [expr { $i + 1 }]]		[lindex $data [expr { $i + 2 }]]"
					incr j
				}
			puts $g ""
			}
		}
	incr k
	}
close $g
}

proc max_fluc_amber {} {
	set n_atoms [input1]
	set n_modes [expr { (3 * $n_atoms) - 6 }]
	set f [open "vecs_amber" "r"]
	set data [read $f]
	close $f
	set g [open "max_fluc_amber" "w"]
	puts $g "freq	atom_number	max_fluc"
	set fg [open "amber_fluc" "w"]
	set k 0
	while { $k != [llength $data] } {
		set max_fluc 0.0
		if { [lindex $data $k] == "freq" } {
			puts $fg " # frequency = [lindex $data [expr { $k + 2 }]]"	
			for {set i [expr { $k + 4 }]} {$i < [expr { (4 *$n_atoms) + ($k + 4) }]} {incr i 4} {
				# calculation of maximum fluctuation and the fluctuating atom
				set x [lindex $data $i] 
				set y [lindex $data [expr { $i + 1 }]]
				set z [lindex $data [expr { $i + 2 }]]
				set fluc [expr { sqrt(($x * $x) + ($y * $y) + ($z * $z)) }]
				puts $fg "[lindex $data [expr { $i -1 } ]] $fluc"
				if { $fluc > $max_fluc } {
					set max_fluc $fluc
					set atno [lindex $data [expr { $i - 1 }]]
				}
			}
		puts $g "[lindex $data [expr { $k + 2 }]] $atno $max_fluc"
		puts $fg ""
		}
	incr k
	}
close $g
close $fg
}

proc max_fluc_gaussian {} {
	set n_atoms [input1]
	set n_modes [expr { (3 * $n_atoms) - 6 }]
	set f [open "vecs_gaussian" "r"]
	set data [read $f]
	close $f
	set g [open "max_fluc_gaussian" "w"]
	set fg [open "gaussian_fluc" "w"]
	puts $g "freq	atom_number	max_fluc"
	set k 0
	while { $k != [llength $data] } {
		set max_fluc 0.0
		if { [lindex $data $k] == "frequency" } {	
			puts $fg " # frequency = [lindex $data [expr { $k + 2 }]]"
			for {set i [expr { $k + 5 }]} {$i < [expr { (5 *$n_atoms) + ($k + 5) }]} {incr i 5} {
				# calculation of maximum fluctuation and the fluctuating atom
				set x [lindex $data $i] 
				set y [lindex $data [expr { $i + 1 }]]
				set z [lindex $data [expr { $i + 2 }]]
				set fluc [expr { sqrt(($x * $x) + ($y * $y) + ($z * $z)) }]
				puts $fg "[lindex $data [expr { $i - 2 } ]] $fluc"
				if { $fluc > $max_fluc } {
					set max_fluc $fluc
					set atno [lindex $data [expr { $i - 2 }]]
				}
			}
		puts $fg ""
		puts $g "[lindex $data [expr { $k + 2 }]] $atno $max_fluc"
		}
	incr k
	}
close $g
close $fg
}

proc amber_visualisation {} {
	set f [open "template.log" "r"]
	set data [read $f]
	close $f

	set h [open "vecs_amber" "r"]
	set data1 [read $h]
	close $h

	set p [open "[input3]" "r"]
	set prm [read $p]
	close $p

	set pd [open "[input4]" "r"]
	set pdb [read $pd]
	close $pd

	set parid 0

	set pr 0
	while { [lindex $prm $pr] != "MASS" } {
		incr pr
	}
	incr pr 2

	set n_atoms [input1]
	set n_modes [expr { (3 * $n_atoms) - 6 }]

	#  Condition to sort the file
	set nmode1 5.1700



																			# ***********CREATING FILE 1************

	puts " step 1 of 7 "
  set g1 [open "1" "w"]
	set stng0 [string first "Symbolic Z-matrix:"  $data]
	set data2 [string replace $data $stng0 end ""]
	puts $g1 "$data2"
	close $g1

																			# ***********CREATING FILE 2************

	puts " step 2 of 7 "
  set g2 [open "2" "w"]
	puts $g2 " Symbolic Z-matrix:"
 	puts $g2 " Charge =  0 Multiplicity = 1"
	for {set i 1} {$i <= $n_atoms} {incr i} {
		while {[lindex $pdb $parid] != $i || [lindex $pdb [expr { $parid + 1 }]] == "File" } {
			incr parid
		}
		puts $g2 " [lindex $pdb [expr { $parid + 9 } ]] 	[lindex $pdb [expr { $parid + 4 }]] 	[lindex $pdb [expr { $parid + 5 }]] 	[lindex $pdb [expr { $parid + 6 }]]"
	}
	puts $g2 ""
	puts $g2 " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad"
	close $g2
	set g2 [open "2" "r"]
	set data3 [read $g2]
	close $g2
		
		
																			# ***********CREATING FILE 3************

	puts " step 3 of 7 "
  set g3 [open "1" "w"]
	set stng0 [string first "Standard orientation:"  $data]
	set data4_1 [string replace $data $stng0 end ""]
	set stng0 [string first "Initialization pass."  $data4_1]
	set data4 [string replace $data4_1 0 $stng0 ""]
	puts $g3 "$data4"
	close $g3

																			# ***********CREATING FILE 4************

	# space variables

	set sp2_1(4) ""
	set sp2_1(3) " "
	set sp2_1(2) "  "
	set sp2_1(1) "   "
	set sp2_2(2) ""
	set sp2_2(1) " "
	set sp2_3(11) ""
	set sp2_3(10) " "
	set sp2_3(9) "  "
	set sp2_3(8) "   " 

	puts " step 4 of 7 " 

	set g4 [open "4" "w"]

	puts $g4 " ---------------------------------------------------------------------"
	
	puts $g4 " Center     Atomic      Atomic             Coordinates (Angstroms)"
 	puts $g4 " Number     Number       Type             X           Y           Z"
	puts $g4 " ---------------------------------------------------------------------"

	set pr1 $pr
	set parid 0
	for {set i 1} {$i <= $n_atoms} {incr i} {
		while { [lindex $pdb $parid] != $i || [lindex $pdb [expr { $parid + 1 }]] == "File" } {
			incr parid
		}
		set coordx [lindex $pdb [expr { $parid + 4 }]]
		set coordx [format "%.6f" $coordx]
		set coordy [lindex $pdb [expr { $parid + 5 }]]
		set coordy [format "%.6f" $coordy]
		set coordz [lindex $pdb [expr { $parid + 6 }]]
		set coordz [format "%.6f" $coordz]

		set mas [lindex $prm $pr1]
		if { $mas != 1.0 } {
			set mas [expr { $mas / 2.0 }]
		} 
		set mas [format "%.0f" $mas]
		set strl1 [string length $i]
		set strl2 [string length $mas]
		set strl3x [string length $coordx]
		set strl3y [string length $coordy]
		set strl3z [string length $coordz]
		puts $g4 "   $sp2_1($strl1)$i         $sp2_2($strl2)$mas           0     $sp2_3($strl3x)$coordx $sp2_3($strl3y)$coordy $sp2_3($strl3z)$coordz"	
		incr pr1
	}
	puts $g4 " ---------------------------------------------------------------------"
	close $g4
	set g4 [open "4" "r"]
	set data5 [read $g4]
	close $g4
																			# ***********CREATING FILE 5************

	puts " step 5 of 7 " 

	set g5 [open "5" "w"]
	set stng0 [string first "Rotational constants (GHZ):"  $data]
	set data6_1 [string replace $data 0 $stng0 ""]
	set stng0 [string first "                      1                      2                      3"  $data6_1]
	set data6 [string replace $data6_1 $stng0 end ""]
	puts $g5 "$data6"
	close $g5

																			# ***********CREATING FILE 6************


	puts " step 6 of 7 " 

	set g [open "6" "w"]

	set k -1
	set kk 0
	set dokk1 0
	set okk1 0
	# SPACE VARIABLES

	# for string 1

	set sp1(1) "	"
	set sp1(2) " "
	set sp1(3) ""

	# for string 2

	set sp21(6) "    "
	set sp21(7) "   "
	set sp21(8) "  "
	set sp21(9) " "
	set sp22(6) "  "
	set sp22(7) " "
	set sp23(6) "  "
	set sp23(7) " "
	set sp24(6) "   "
	set sp24(7) "  "
	set sp24(8) " "
	set sp24(9) ""
	# for string 4

	set sp4(1) " "
	set sp4(2) ""
	set sp4n(1) "  "
	set sp4n(2) " "
	set sp4n(3) ""
	set sprc4(4) " "
	set sprc4(5) ""

	while { $k < [llength $data] } {
		set k1 $k
		if { [lindex $data $k] == $nmode1 } {
			for {set nm 1} {$nm <= $n_modes} {incr nm 3} {
				puts "			**** Analyzing mode $nm [expr { $nm + 1 }] and [expr { $nm + 2 }] ****"

				# FORMING STRING 1

				set stl1 [string length $nm]
				set stl2 [string length [expr { $nm + 1 }]]
				set stl3 [string length [expr { $nm + 2 }]]
				puts $g "\t\t\t\t\t\t\t\t\t\t$sp1($stl1)$nm\t\t\t\t\t\t\t\t\t\t$sp1($stl2)[expr { $nm + 1 }]\t\t\t\t\t\t\t\t\t\t$sp1($stl3)[expr { $nm + 2 }]"
				puts $g "\t\t\t\t\t\t\t\t\t\t\tA\t\t\t\t\t\t\t\t\t\t\t A\t\t\t\t\t\t\t\t\t\t\tA"

				# FORMING STRING 2
				
				set dk $k
				for {set t 0} {$t < 12} {incr t 3} {
					if { $t == 0 } {
						for {set ff 0} {$ff < 3} {incr ff} {
							while { [lindex $data1 $kk] != "freq" } {
								incr kk
							}
							incr kk
							set newf($ff) [expr { double(round(10000*[lindex $data1 [expr { $kk + 1 }]])) /10000 }]	
							if { [string length $newf($ff)] < 6 } {
								set newf($ff) [format "%.4f" $newf($ff)]
							}						
						}
						set stlf1 [string length $newf(0)]
						set stlf2 [string length $newf(1)]
						set stlf3 [string length $newf(2)]
						puts $g " Frequencies --  $sp21($stlf1)$newf(0)  \t\t\t\t\t  $sp21($stlf2)$newf(1)  \t\t\t\t\t $sp21($stlf3)$newf(2)"
					} 
					if { $t == 3 } {
						set dk [expr { $dk + $t + 1 }]
						set rm1 0.0001
						set rm2 0.0001
						set rm3 0.0001
						set stlm1 [string length $rm1]
						set stlm2 [string length $rm2]
						set stlm3 [string length $rm3]
						puts $g " Red. masses --    $sp22($stlm1)$rm1  \t\t\t\t\t\t  $sp22($stlm2)$rm2  \t\t\t\t\t\t $sp22($stlm3)$rm3"	
					}
					if { $t == 6 } {
						set dk [expr { $dk + $t }]
						set fc1 0.0001
						set fc2 0.0001
						set fc3 0.0001
						set stlfc1 [string length $fc1]
						set stlfc2 [string length $fc2]
						set stlfc3 [string length $fc3]
						puts $g " Frc consts  --    $sp23($stlfc1)$fc1  \t\t\t\t\t\t  $sp23($stlfc2)$fc2  \t\t\t\t\t\t $sp23($stlfc3)$fc3"						
					}
					if { $t == 6 } {
						set dk [expr { $dk + $t }]
						set ir1 0.0001
						set ir2 0.0001
						set ir3 0.0001
						set stlir1 [string length $ir1]
						set stlir2 [string length $ir2]
						set stlir3 [string length $ir3]
						puts $g " IR Inten    --   $sp24($stlir1)$ir1  \t\t\t\t\t\t $sp24($stlir2)$ir2  \t\t\t\t\t  $sp24($stlir3)$ir3"						
					}
				}

				# FORMING STRING 3

				puts $g "  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z"

				# FORMING STRING 4		

				set pr1 $pr
				set kk1 $dokk1
				for {set i $k1} {$i < [expr { $k1 + (11 * $n_atoms) } ]} {incr i 11} {
					set achan [expr { ($i - $k1) / 11 }]
					set an [expr { $achan + 1 }]
					set stlan [string length $an]
					set am [lindex $prm $pr1]
					if { $am != 1.0 } {
						set am [expr { $am / 2.0 }]
					} 
					set am [format "%.0f" $am]
					set stlam [string length $am]
					incr pr1
					for {set ff 0} {$ff < 3} {incr ff} {
						while { [lindex $data1 $kk1] != "freq" } {
							incr kk1
						}
						incr kk1
						set xr($ff) [expr { double(round(100*[lindex $data1 [expr { $kk1 + 3 + (4*$achan) }]])) /100 }]
						set stlch [string length $xr($ff)]
						if { $stlch < 4 } {
							set xr($ff) [format "%.2f" $xr($ff)]
						}

						set yr($ff) [expr { double(round(100*[lindex $data1 [expr { $kk1 + 4 + (4*$achan) }]])) /100 }]
						set stlch [string length $yr($ff)]
						if { $stlch < 4 } {
							set yr($ff) [format "%.2f" $yr($ff)]
						}

						set zr($ff) [expr { double(round(100*[lindex $data1 [expr { $kk1 + 5 + (4*$achan) }]])) /100 }]
						set stlch [string length $zr($ff)]
						if { $stlch < 4 } {
							set zr($ff) [format "%.2f" $zr($ff)]
						}
					}
						set stlxr0 [string length $xr(0)]
						set stlyr0 [string length $yr(0)]
						set stlzr0 [string length $zr(0)] 
						set stlxr1 [string length $xr(1)]
						set stlyr1 [string length $yr(1)]
						set stlzr1 [string length $zr(1)]
						set stlxr2 [string length $xr(2)]
						set stlyr2 [string length $yr(2)]
						set stlzr2 [string length $zr(2)]
						puts $g "   $sp4n($stlan)$an  $sp4($stlam)$am    $sprc4($stlxr0)$xr(0)  $sprc4($stlyr0)$yr(0)  $sprc4($stlzr0)$zr(0)    $sprc4($stlxr1)$xr(1)  $sprc4($stlyr1)$yr(1)  $sprc4($stlzr1)$zr(1)    $sprc4($stlxr2)$xr(2)  $sprc4($stlyr2)$yr(2)  $sprc4($stlzr2)$zr(2)"
						set dokk1 $kk1
						set kk1 $okk1
				}
			set okk1 $dokk1	
			set k [expr { $i + 8 }]  
			set k1 $k
			}
			set k [llength $data]		
		}
		incr k
	}
	close $g	
	set g [open "6" "r"]
	set data7 [read $g]
	close $g

# 																		******************* CREATING FILE 7 ************************

	puts " Step 7 of 7 "
	set g7 [open "7" "w"]
	set stng0 [string first " - Thermochemistry -" $data]
	set data8_1 [string replace $data 0 $stng0 ""]
	set stng0 [string first "GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad" $data8_1]
	set data8 [string replace $data8_1 $stng0 end ""]
	puts $g7 "$data8"
	close $g7

#																			***************** CONCATINATION **********************

	set data9 [concat $data2 $data3 $data4 $data5 $data6 $data7 $data8]
	set g8 [open "amber_freq.log" "w"]
	puts $g8 "$data9"
	close $g8

}

proc plot {} {
	set p [open "plot.gp" "w"]
	puts $p "plot 'max_fluc_amber' u 1:2 w lp"
	puts $p "pause 5"
	close $p
	exec gnuplot plot.gp
}

proc input1 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f
	return [lindex $data 0]
}

proc input2 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f
	return [lindex $data 1]
}

proc input3 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f
	return [lindex $data 2]
	
}

proc input4 {} {
	set f [open "input" "r"]
	set data [read $f]
	close $f

	set file1 [lindex $data 2]
	set file2 [lindex $data 3]

	set gf [open "input_cpptraj" "w"]
 
	puts $gf "trajin $file2"
	puts $gf "trajout mod_pdb.pdb"
	puts $gf "go"

	close $gf

	exec cpptraj -p $file1 -i input_cpptraj

	file delete input_cpptraj
	return mod_pdb.pdb
}

proc delete_file {} {
	for {set i 1} {$i < 8} {incr i} {
		file delete $i
	}
	file delete mod_pdb.pdb
}

# gaussian
amber
max_fluc_amber
# max_fluc_gaussian
amber_visualisation
delete_file
#plot

				
				 

























