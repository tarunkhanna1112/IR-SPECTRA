proc input_files {} {

	puts " **** PROCEDURE :: CREATING INPUT FILES **** "

	set f [open "[input1]" "r"]
	set data [read $f]
	close $f
	
	set f2 [open "[input3]" "r"]
	set dat [read $f2]
	close $f2

	set g [open "atomic_mass" "w"]
	set h [open "cartesian_coord" "w"]
	set m [open "cartesian_force_constant" "w"]
	set n [open "cartesian_fc_matrix" "w"]

	set k 10

	# ATOMIC WEIGHTS

	puts " writing the atomic weights "

	
	while { [lindex $data $k] != "weights" || [lindex $data [expr { $k - 1 }]] != "atomic"  || [lindex $data [expr { $k - 2 }]] != "Real" } {
		incr k
	}	
	set num_elem [lindex $data [expr { $k + 3 }]]
	set num_coord $num_elem
	incr k 4

	for {set i $k} {$i < [expr { $k + $num_elem }]} {incr i} {
		puts $g "{ [lindex $data $i] }"
	}		

	 # CARTESIAN COORDINATES

	puts " writing the cartesian coordinates :: reading from the mol2 file "
	set kk 0
	while {[lindex $dat $kk] != "@<TRIPOS>ATOM"} {
		incr kk
	}
	set atnum 1
	while { [lindex $dat $kk] != "@<TRIPOS>BOND"} {
		if {[lindex $dat $kk] == $atnum} {
			puts $h "{ [lindex $dat [expr { $kk + 2 }]] [lindex $dat [expr {$kk + 3}]] [lindex $dat [expr {$kk + 4}]] }"
			incr atnum
		}
		incr kk
	}

	# CARTESIAN FORCE CONSTANTS LOWER TRAIANGULAR MATRIX

	puts " writing the cartesian force constants "

	while {[lindex $data $k] != "Constants" || [lindex $data [expr { $k - 1 }]] != "Force" || [lindex $data [expr { $k - 2 }]] != "Cartesian"} {
		incr k
	}	
	set num_elem [lindex $data [expr { $k + 3 }]]
	incr k 4
	
	set l 0
	set rl $k
	set j 1
	puts $m "*"
	for {set i $k} {$i < [expr { $k + $num_elem }] } {incr i} {
		if { $i >= [expr { $rl + $l }] } {
			puts $m "[lindex $data $i]"
			incr j
			puts $m "**"
			incr l
			set rl [expr { $i + 1 }]
			if { $rl < [expr { $k  + $num_elem }] } {
				puts $m "*"
			}
		} else {
			puts $m "[lindex $data $i]"
		}
	}

	close $g 
	close $h
	close $m

	set m [open "cartesian_force_constant" "r"]
	set data [read $m]
	close $m

	set k 0

	while { $k < [llength $data] } {
		set dum [open "dummy" "w"]
		if { [lindex $data $k] == "*" } {
			incr k 
	 		while { [lindex $data $k] != "**" && $k < [llength $data] } {
				puts $dum "[lindex $data $k]" 
		 		incr k 
			}
		}
		close $dum
		set dum [open "dummy" "r"]
		set data1 [read $dum]
		close $dum
		puts $n " { $data1 } "
		incr k
	}
	close $n
}

proc hessian_full {} {
	puts " **** PROCEDURE :: TRANSFORMING THE UPPERTRAINGULAR HESIAN TO FULL HESSIAN **** "

	package require math::linearalgebra

	# IN THIS SECTION WE WILL TRANSFORM THE UPPER TRAINGULAR HESSIAN READ FROM CHECKPOINT FILE OF GAUSSIAN TO FULL MATRIX

	set n_atoms [input2]
	set n_rows [expr { $n_atoms * 3 }]

	set f [open "cartesian_fc_matrix" "r"]
	set data [read $f]
	close  $f

	set g [open "cart_fc_matrix_full" "w"]

	for {set i 0} {$i < $n_rows} {incr i} {
		set m1 [lindex $data $i]
		set m2 $m1
		set col $i
		while { $col < [expr { $n_rows - 1 }] } {
			set m2 [linsert $m2 end [lindex $data [expr { $col + 1 }] $i]]
			incr col
		}
		puts $g "{ $m2 }"
	}
	close $g
}

proc bonding_info {a} {

	#puts " **** PROCEDURE :: READING THE BONDING INFORMATION FROM THE MOL2 FILE OF THE COMPOUND **** "

	# IN THIS SECTION WE WILL STORE THE BONDING INFORMATION OF THE RESIDUE IN AN ARRAY 

	set elem $a

	set n_atoms [input2]

	set f [open "[input3]" "r"]
	set data [read $f]
	close $f
	set k 0
	while { [lindex $data $k] != $n_atoms } { 
		incr k
	}
	set n_bonds [lindex $data [expr { $k + 1}]]
	
	while { [lindex $data $k] != "@<TRIPOS>BOND" } {
		incr k
	}
	incr k

	# storing the information for the bond pairs in the format atom(atom_number,bonded_atom_number)

	for {set j 1} {$j <= $n_atoms} {incr j} {
		set l 1
		for {set i $k} {$i < [expr { $k + (4*$n_bonds) }]} {incr i 4} {
			if { [lindex $data [expr { $i + 1 }]] == $j } {
				set atom($j,$l) [lindex $data [expr {$i + 2}]]
				incr l
			}
			if { [lindex $data [expr { $i + 2 }]] == $j } {
				set atom($j,$l) [lindex $data [expr {$i + 1}]]
				incr l
			}
		}
		set atom($j,$l) "end"
	}
	set l 1
	set dat ""
	while { $atom($elem,$l) != "end" } {
		set dat [linsert $dat end $atom($elem,$l)]
		incr l
	} 
	return [list $dat]
}

proc bond_length_AB_2 {A B} {

	# THIS SECTION WILL CALCULATE THE VECTOR BETWEEN THE ATOM A AND ATOM B AND GIVES OUT THE LIST

	package require math::linearalgebra

	set f [open "cartesian_coord" "r"]
	set data [read $f]
	close $f

	set elem1 [expr { $A - 1 }]
	set elem2 [expr { $B - 1 }]
 
	set vec1x [lindex $data $elem1 0]
	set vec1y [lindex $data $elem1 1]
	set vec1z [lindex $data $elem1 2]

	set vec2x [lindex $data $elem2 0]
	set vec2y [lindex $data $elem2 1]
	set vec2z [lindex $data $elem2 2]

	set vecx [expr { $vec1x - $vec2x }] 
	set vecy [expr { $vec1y - $vec2y }]
	set vecz [expr { $vec1z - $vec2z }]
		
	set bl [expr { ($vecx*$vecx) + ($vecy*$vecy) + ($vecz*$vecz) }]
	
	return $bl

}

proc AB_vectors {A B} {

	# THIS SECTION WILL CALCULATE THE VECTOR BETWEEN THE ATOM A AND ATOM B AND GIVES OUT THE LIST

	package require math::linearalgebra

	set f [open "cartesian_coord" "r"]
	set data [read $f]
	close $f

	set elem1 [expr { $A - 1 }]
	set elem2 [expr { $B - 1 }]
 
	set vec1x [lindex $data $elem1 0]
	set vec1y [lindex $data $elem1 1]
	set vec1z [lindex $data $elem1 2]

	set vec2x [lindex $data $elem2 0]
	set vec2y [lindex $data $elem2 1]
	set vec2z [lindex $data $elem2 2]

	set vecx [expr { $vec1x - $vec2x }] 
	set vecy [expr { $vec1y - $vec2y }]
	set vecz [expr { $vec1z - $vec2z }]

	set vec [list $vecx $vecy $vecz]

	set norm [::math::linearalgebra::unitLengthVector $vec]

	return $norm

}

proc plane_ABC { A B C D } {

	# IN THIS SECTION WE CALCULATE THE PLANE OF THE ATOMS A B AND C FOR BOND ANGLE FORCE CONSTANT CALCULATION

	package require math::linearalgebra

	if { $D == 0 } {
		set AB_vec [AB_vectors $A $B]
		set CB_vec [AB_vectors $C $B]
		
		# NORMAL TO PLANE ABC

		set normal [::math::linearalgebra::crossproduct $CB_vec $AB_vec]
		set unit_normal [::math::linearalgebra::unitLengthVector $normal]
		return [list $unit_normal]
	} else { 
		set AB_vec [AB_vectors $A $B]
		set BC_vec [AB_vectors $B $C]
		set DC_vec [AB_vectors $D $C]
		set CB_vec [AB_vectors $C $B]

		# NORMAL TO PLANE ABC

		set normal [::math::linearalgebra::crossproduct $CB_vec $AB_vec]
		set unit_normal_ABC [::math::linearalgebra::unitLengthVector $normal]

		# NORMAL TO PLANE BCD

		set normal [::math::linearalgebra::crossproduct $DC_vec $BC_vec]
		set unit_normal_BCD [::math::linearalgebra::unitLengthVector $normal]

		return [list $unit_normal_ABC $unit_normal_BCD]
	}
}

proc eigenvalue_AB {} {

	# IN THIS SECTION WE WILL DETERMINE THE EIGENVALUE AND EIGENVECTOR OF K_AB MATRIX MAKING USE OF THE LARGEST EIGENVALUE CALCULATED THROUGH TCL LIB

	package require math::linearalgebra
	package require math::polynomials
		
	set f [open "dummy" "r"]
	set data [read $f]
	close $f

	set a11 [lindex $data 0 0]
	set a12 [lindex $data 0 1]
	set a13 [lindex $data 0 2]
	set a21 [lindex $data 1 0]
	set a22 [lindex $data 1 1]
	set a23 [lindex $data 1 2]
	set a31 [lindex $data 2 0]
	set a32 [lindex $data 2 1]
	set a33 [lindex $data 2 2]
	set lambda1 [lindex $data 3]

	set evec1 [lindex $data 4] 

	# CHARACTERSTIC EQUATION OF MATRIX A :: lambda^3 + a*lambda^2 + b*lambda + c = 0

	set a [expr { ($a11 + $a22 + $a33) * -1.0 }]
	set b [expr { ($a11*$a22) + ($a33*$a11) + ($a22*$a33) - ($a23*$a32) - ($a21*$a12) - ($a31*$a13) }]
	set c [expr { (($a11*$a22*$a33) - ($a32*$a23*$a11) - ($a12*$a21*$a33) + ($a23*$a31*$a12) + ($a13*$a21*$a32) - ($a22*$a31*$a13)) * -1.0 }]
	
	set ch_poly [::math::polynomials::polynomial [list $c $b $a 1.0]]

	set leval_poly [::math::polynomials::polynomial [list [expr { -1*$lambda1 }] 1.0]]	

	set rem_poly [::math::polynomials::divPolyn $ch_poly $leval_poly]
	
	set g [open "dummy1" "w"]
	puts $g "$rem_poly"
	close $g

	# SOLVING THE SECOND ORDER POLYNOMIAL FOR THE OTHER TW0 EIGENVECTORS a*lambda^2 + b*lambda + c = 0

	set g [open "dummy1" "r"]
	set data1 [read $g]
	close $g

	set a [lindex $data1 1 0]
	set b [lindex $data1 1 1]
	set c [lindex $data1 1 2]

	set check [expr { ($b*$b) - (4*$a*$c) }]
	if { $check < 0 } {
		
		#puts "***WARNING****EIGENVALUE FOUND TO BE COMPLEX CHECK FOR THE DIRECTION OF LARGEST EIGENVECTOR"
		set eigenvalue [list $lambda1]
		set eigenvector [list $evec1]		
		return [list $eigenvalue $eigenvector]		
	} else {
		set root1 [expr { ((-1*$b) + (sqrt($check))) / 2.0 }]
		set root2 [expr { ((-1*$b) - (sqrt($check))) / 2.0 }]
	}	
	if { $root1 > 0 && $root2 > 0 } {
		set evec2 [eigenvector_AB $root1]
		set evec3 [eigenvector_AB $root2]
		set eigenvalue [list $lambda1 $root1 $root2]
		set eigenvector [list $evec1 $evec2 $evec3]		

		return [list $eigenvalue $eigenvector]		
	} 
	if { $root1 > 0 && $root2 < 0 } {
		set evec2 [eigenvector_AB $root1]
		set eigenvalue [list $lambda1 $root1 $root2]
		set eigenvector [list $evec1 $evec2]		

		return [list $eigenvalue $eigenvector]		
	} 
	if { $root1 < 0 && $root2 > 0 } {
		set evec3 [eigenvector_AB $root2]
		set eigenvalue [list $lambda1 $root1 $root2]
		set eigenvector [list $evec1 $evec3]		

		return [list $eigenvalue $eigenvector]		
	} 
	if { $root1 < 0 && $root2 < 0 } {
		set eigenvalue [list $lambda1 $root1 $root2]
		set eigenvector [list $evec1]		

		return [list $eigenvalue $eigenvector]		
	} 

}	

proc eigenvector_AB {root} {
	# IN THIS SECTION WE WILL DETERMINE EIGENVECTOR OF K_AB MATRIX 

	package require math::linearalgebra
	
	set f [open "dummy" "r"]
	set data [read $f]
	close $f

	set row1 [lindex $data 0]
	set row2 [lindex $data 1]
	set row3 [lindex $data 2]
	
	set A_matrix [list $row1 $row2 $row3]

	set lrow1 [list $root 0.0 0.0]
	set lrow2 [list 0.0 $root 0.0]
	set lrow3 [list 0.0 0.0 $root]

	set lambda_matrix [list $lrow1 $lrow2 $lrow3]

	set char_matrix [::math::linearalgebra::sub $A_matrix $lambda_matrix]

	set det [::math::linearalgebra::det $char_matrix]

	set svd [::math::linearalgebra::determineSVD $char_matrix]

	set g [open "dummy1" "w"]
	puts $g "$svd"
	close $g

	set g [open "dummy1" "r"]
	set data1 [read $g]
	close $g

	set v_mat [lindex $data1 2]

	set vec [::math::linearalgebra::transpose $v_mat]

	set g [open "dummy1" "w"]
	puts $g "$vec"
	close $g

	set g [open "dummy1" "r"]
	set data1 [read $g]
	close $g

	set evec [lindex $data1 2]	

	return $evec
}

proc force_constantAB { A B } {

	# CALCULATION OF BOND FORCE CONSTANT BETWEEN ATOM A AND ATOM B :: FORCE ON ATOM A DUE TO DISPLACEMENT OF ATOM B

	set f [open "cart_fc_matrix_full" "r"]
	set data [read $f]
	close $f

	set h [open "dummy" "w"]
					
	set row [expr { ($A-1) * 3 }]
	set col [expr { ($B-1) * 3 }]
		
	# READING THE K MATRIX FOR PAIR A B

	for {set k 0} {$k < 3} {incr k} {
		puts $h "*"
		for {set k1 0} {$k1 < 3} {incr k1} {
			puts $h "[lindex $data [expr { $k + $row }] [expr {$k1 + $col}]]"		
		}
		puts $h "**"
	}
					
	close $h
		
	set m [open "dummy" "r"]
	set data3 [read $m]
	close $m

	set n [open "dummy_matrix" "w"]

	set k 0

	while { $k < [llength $data3] } {
		set dum [open "dummy1" "w"]
		if { [lindex $data3 $k] == "*" } {
			incr k 
				while { [lindex $data3 $k] != "**" && $k < [llength $data3] } {
				puts $dum "[lindex $data3 $k]" 
					incr k 
			}
		}
		close $dum
		set dum [open "dummy1" "r"]
		set data4 [read $dum]
		close $dum
		puts $n " { $data4 } "
		incr k
	}
	close $n		
					
	set h [open "dummy_matrix" "r"]
	set data3 [read $h]
	close $h

	set kab [::math::linearalgebra::scale -2016.783 $data3]
			
	set kab1 [::math::linearalgebra::largesteigen $kab]
		
	set dum [open "dummy" "w"]
	puts $dum "$kab $kab1"
	close $dum

	set kab [eigenvalue_AB]
				
	set dum [open "dummy" "w"]
	puts $dum "$kab"
	close $dum

	set dum [open "dummy" "r"]
	set data3 [read $dum]
	close $dum

	set eigvec [open "eigenvector" "w"]
	set eigval [open "eigenvalue" "w"]

	puts $eigvec [lindex $data3 1]
	puts $eigval [lindex $data3 0]

	close $eigvec
	close $eigval

	set eigvec [open "eigenvector" "r"]
	set dat_ev [read $eigvec]
	close $eigvec
	set eigval [open "eigenvalue" "r"]
	set dat_eval [read $eigval]
	close $eigval

	set vectorAB [AB_vectors $A $B]			
	set kab 0.0
	set listsize [llength $dat_ev]
	for {set ve 0} {$ve < $listsize} {incr ve} {
		set vec1 [::math::linearalgebra::dotproduct $vectorAB [lindex $dat_ev $ve]]
		set vec1 [expr { abs($vec1) } ]
		set kab [expr { $kab + ($vec1*[lindex $dat_eval $ve]) }]
	}
	return [list $A $B $kab]	
}
	
proc seminario {} {

	puts " **** PROCEDURE :: IN THIS PROCEDURE WE WILL USE SEMINARIO EXPRESSION FOR FORCE CONSTANT TO DETERMINE THE FORCE FIELD PARAMETERS **** "

	package require math::linearalgebra

	# IN THIS SECTION WE WILL USE SEMINARIO METHOD TO CALCULATE FORCE FIELD PARAMETERS : A-B-C-D

	set f [open "cart_fc_matrix_full" "r"]
	set data [read $f]
	close $f

	set n_atoms [input2]
	set n_rows [expr { $n_atoms * 3 }]
	
	# WRINTING THE FORCE CONSTANTS 

	set bf [open "results.bond_fc" "w"]
	set af [open "results.angle_fc" "w"]
	set df [open "results.dihedral_fc" "w"]

	for {set i 1} {$i <= $n_atoms} {incr i} {

		puts " 		**** ANALYZING ATOM $i OF $n_atoms **** 	"

		branch $i

		# CALCULATING THE FORCE FIELD PARAMETERS FOR ALL CONNECTIONS TO ATOM i ACCORDING TO SEMINARIO EXPRESSION OF FORCE CONSTANTS

		set g [open "branch" "r"]
		set data2 [read $g]
		close $g

		set k 0
		while { [lindex $data2 $k] != "A" } {
			incr k
		}
		incr k 

		set A_index $k
		set A_length [llength [lindex $data2 $k]]

		while { [lindex $data2 $k] != "B" } {
			incr k
		}
		incr k 

		set B_index $k
		set B_length [llength [lindex $data2 $k]]

		while { [lindex $data2 $k] != "C" } {
			incr k
		}
		incr k 

		set C_index $k
		set C_length 0
		for {set j $k} {$j < [expr { $k + $B_length }]} {incr j} {
			set C_length [expr { $C_length + [llength [lindex $data2 $k]] }]
		}

		while { [lindex $data2 $k] != "D" } {
			incr k
		}
		incr k 

		set D_index $k
		set D_length 0
		for {set j $k} {$j < [expr { $k + $C_length }]} {incr j} {
			set D_length [expr { $D_length + [llength [lindex $data2 $k]] }]
		}

		set A $i

		# BONDED TERM

		for {set nb 0} {$nb < $B_length} {incr nb} {

			set B [lindex $data2 $B_index $nb]
			set k_AB [force_constantAB $A $B]
			puts $bf "{ $k_AB }"

			# ANGLE TERMS

			set BC_pairs [llength [lindex $data2 [expr { $nb + $C_index }]]]

			for {set nc 0} {$nc < $BC_pairs} {incr nc} {

				set C [lindex $data2 [expr { $C_index + $nb }] $nc]
				
				set normal [plane_ABC $A $B $C 0]

				set dum [open "dummy" "w"]
				puts $dum "$normal"
				close $dum

				set dum [open "dummy" "r"]
				set data3 [read $dum]
				close $dum
			
				set vec_AB [AB_vectors $A $B]
				set vec_CB [AB_vectors $C $B]
	
				set PA [::math::linearalgebra::crossproduct [lindex $data3 0] $vec_AB]
				set PC [::math::linearalgebra::crossproduct $vec_CB [lindex $data3 0]]

				set R_AB2 [bond_length_AB_2 $A $B]
				set R_CB2 [bond_length_AB_2 $C $B]

				set k_AB [force_constantAB $A $B]

				set eigenv_AB [open "eigenvalue" "r"]
				set e_AB [read $eigenv_AB]
				close $eigenv_AB

				set eigenve_AB [open "eigenvector" "r"]
				set ev_AB [read $eigenve_AB]
				close $eigenve_AB

				set k_CB [force_constantAB $C $B]

				set eigenv_CB [open "eigenvalue" "r"]
				set e_CB [read $eigenv_CB]
				close $eigenv_CB

				set eigenve_CB [open "eigenvector" "r"]
				set ev_CB [read $eigenve_CB]
				close $eigenve_CB

				set term1 0.0
				set term2 0.0

				set num_val [llength $ev_AB]
				for {set val 0} {$val < $num_val} {incr val} {
					set vec_term1 [lindex $ev_AB $val]
					set dot1 [::math::linearalgebra::dotproduct $PA $vec_term1]
					set dot1 [expr { abs($dot1) }]
					set term1 [expr { $term1 + ([lindex $e_AB $val] * $dot1) }]
				}

				set num_val [llength $ev_CB]
				for {set val 0} {$val < $num_val} {incr val} {
					set vec_term2 [lindex $ev_CB $val]
					set dot2 [::math::linearalgebra::dotproduct $PC $vec_term2]
					set dot2 [expr { abs($dot2) }]
					set term2 [expr { $term2 + ([lindex $e_CB $val] * $dot2) }]
				}
				
				set term1 [expr { $term1 * $R_AB2 }]
				set term2 [expr { $term2 * $R_CB2 }]

				set term1 [expr { 1.0 / $term1 }]
				set term2 [expr { 1.0 / $term2 }]
				set k_theta_ABC [expr { $term1 + $term2 }]
				set k_theta_ABC [expr { 1.0 / $k_theta_ABC }]
			
				puts $af "{ $A $B $C $k_theta_ABC }"

				# DIHEDRALS AND IMPROPER TERMS

				set BCD_pairs [llength [lindex $data2 [expr { $nc + $D_index }]]]

				for {set nd 0} {$nd < $BCD_pairs} {incr nd} {

					# DIHEDRAL TERMS
					set D [lindex $data2 [expr { $D_index + $nc }] $nd]
	
					set normal [plane_ABC $A $B $C $D]

					set dum [open "dummy" "w"]
					puts $dum "$normal"
					close $dum

					set dum [open "dummy" "r"]
					set data3 [read $dum]
					close $dum
			
					set vec_AB [AB_vectors $A $B]
					set vec_BC [AB_vectors $B $C]
					set vec_CD [AB_vectors $C $D]
						

					set cp1 [::math::linearalgebra::crossproduct $vec_AB $vec_BC]
					set cp1 [::math::linearalgebra::dotproduct $cp1 $cp1]
					set cp2 [::math::linearalgebra::crossproduct $vec_BC $vec_CD]
					set cp2 [::math::linearalgebra::dotproduct $cp2 $cp2]

					set R_BA2 [bond_length_AB_2 $A $B]
					set R_CD2 [bond_length_AB_2 $C $D]

					set k_DC [force_constantAB $D $C]

					set eigenv_DC [open "eigenvalue" "r"]
					set e_DC [read $eigenv_DC]
					close $eigenv_DC
	
					set eigenve_DC [open "eigenvector" "r"]
					set ev_DC [read $eigenve_DC]
					close $eigenve_DC

					set term1 0.0
					set term2 0.0

					set num_val [llength $ev_AB]
					for {set val 0} {$val < $num_val} {incr val} {
						set vec_term1 [lindex $ev_AB $val]
						set dot1 [::math::linearalgebra::dotproduct [lindex $data3 0] $vec_term1]
						set dot1 [expr { abs($dot1) }]
						set term1 [expr { $term1 + ([lindex $e_AB $val] * $dot1) }]
					}

					set num_val [llength $ev_DC]
					for {set val 0} {$val < $num_val} {incr val} {
						set vec_term2 [lindex $ev_DC $val]
						set dot2 [::math::linearalgebra::dotproduct [lindex $data3 1] $vec_term2]
						set dot2 [expr { abs($dot2) }]
						set term2 [expr { $term2 + ([lindex $e_DC $val] * $dot2) }]
					}
				
					set term1 [expr { $term1 * $R_AB2 * $cp1 }]
					set term2 [expr { $term2 * $R_CD2 * $cp2 }]

					set term1 [expr { 1.0 / $term1 }]
					set term2 [expr { 1.0 / $term2 }]
					set k_phi_ABCD [expr { $term1 + $term2 }]
					set k_phi_ABCD [expr { 1.0 / $k_phi_ABCD }]
			
					puts $df "{ $A $B $C $D $k_phi_ABCD }"
				}			 
			}
		incr D_index $nc
		}

	# GETTING THE INFORMATION OF THE ELEMENTS WHICH ORGINATES AT BRANCH B AND INCLUDE ELEMENT A AS THE B BRANCH OF THAT TREE

		set B_A_length $B_length

		for {set b_branch 0} {$b_branch < $B_A_length} {incr b_branch} {
			
			set A1($b_branch) [lindex $data2 $B_index $b_branch]

			branch $A1($b_branch)

			# CALCULATING THE FORCE FIELD PARAMETERS FOR ALL CONNECTIONS TO ATOM i ACCORDING TO AMBER FORCE FIELD

			set g [open "branch" "r"]
			set data2b [read $g]
			close $g

			set k 0
			while { [lindex $data2b $k] != "A" } {
				incr k
			}
			incr k 

			set A_index $k
			set A_length [llength [lindex $data2b $k]]

			while { [lindex $data2b $k] != "B" } {
				incr k
			}
			incr k 

			set B_index $k
			set B_length [llength [lindex $data2b $k]]

			while { [lindex $data2b $k] != "C" } {
				incr k
			}
			incr k 

			set C_index $k
			set C_length 0
			for {set j $k} {$j < [expr { $k + $B_length }]} {incr j} {
				set C_length [expr { $C_length + [llength [lindex $data2b $k]] }]
			}

			while { [lindex $data2b $k] != "D" } {
				incr k
			}
			incr k 

			set D_index $k
			set D_length 0
			for {set j $k} {$j < [expr { $k + $C_length }]} {incr j} {
				set D_length [expr { $D_length + [llength [lindex $data2b $k]] }]
			}

				# BONDED TERM

			for {set nb 0} {$nb < $B_length} {incr nb} {

				set B [lindex $data2 $B_index $nb]

				if { $B == $A } {

					# ANGLE TERMS

					set BC_pairs [llength [lindex $data2 [expr { $nb + $C_index }]]]

					for {set nc 0} {$nc < $BC_pairs} {incr nc} {

						set C [lindex $data2 [expr { $C_index + $nb }] $nc]

						set count 0
						for {set chk 0} {$chk < $b_branch} {incr chk} {
							if { $A1($chk) == $C } {
								incr count
							}
						}

						if { $count == 0 } {
				
							set normal [plane_ABC $A $B $C 0]

							set dum [open "dummy" "w"]
							puts $dum "$normal"
							close $dum

							set dum [open "dummy" "r"]
							set data3 [read $dum]
							close $dum
			
							set vec_AB [AB_vectors $A $B]
							set vec_CB [AB_vectors $C $B]
	
							set PA [::math::linearalgebra::crossproduct [lindex $data3 0] $vec_AB]
							set PC [::math::linearalgebra::crossproduct $vec_CB [lindex $data3 0]]

							set R_AB2 [bond_length_AB_2 $A $B]
							set R_CB2 [bond_length_AB_2 $C $B]

							set k_AB [force_constantAB $A $B]

							set eigenv_AB [open "eigenvalue" "r"]
							set e_AB [read $eigenv_AB]
							close $eigenv_AB

							set eigenve_AB [open "eigenvector" "r"]
							set ev_AB [read $eigenve_AB]
							close $eigenve_AB

							set k_CB [force_constantAB $C $B]

							set eigenv_CB [open "eigenvalue" "r"]
							set e_CB [read $eigenv_CB]
							close $eigenv_CB

							set eigenve_CB [open "eigenvector" "r"]
							set ev_CB [read $eigenve_CB]
							close $eigenve_CB

							set term1 0.0
							set term2 0.0

							set num_val [llength $ev_AB]
							for {set val 0} {$val < $num_val} {incr val} {
								set vec_term1 [lindex $ev_AB $val]
								set dot1 [::math::linearalgebra::dotproduct $PA $vec_term1]
								set dot1 [expr { abs($dot1) }]
								set term1 [expr { $term1 + ([lindex $e_AB $val] * $dot1) }]
							}

							set num_val [llength $ev_CB]
							for {set val 0} {$val < $num_val} {incr val} {
								set vec_term2 [lindex $ev_CB $val]
								set dot2 [::math::linearalgebra::dotproduct $PC $vec_term2]
								set dot2 [expr { abs($dot2) }]
								set term2 [expr { $term2 + ([lindex $e_CB $val] * $dot2) }]
							}
				
							set term1 [expr { $term1 * $R_AB2 }]
							set term2 [expr { $term2 * $R_CB2 }]

							set term1 [expr { 1.0 / $term1 }]
							set term2 [expr { 1.0 / $term2 }]
							set k_theta_ABC [expr { $term1 + $term2 }]
							set k_theta_ABC [expr { 1.0 / $k_theta_ABC }]
			
							puts $af "{ $A $B $C $k_theta_ABC }"
						}

						# DIHEDRALS AND IMPROPER TERMS

						set BCD_pairs [llength [lindex $data2 [expr { $nc + $D_index }]]]

						for {set nd 0} {$nd < $BCD_pairs} {incr nd} {

							# DIHEDRAL TERMS
							set D [lindex $data2 [expr { $D_index + $nc }] $nd]
	
							set normal [plane_ABC $A $B $C $D]

							set dum [open "dummy" "w"]
							puts $dum "$normal"
							close $dum

							set dum [open "dummy" "r"]
							set data3 [read $dum]
							close $dum
			
							set vec_AB [AB_vectors $A $B]
							set vec_BC [AB_vectors $B $C]
							set vec_CD [AB_vectors $C $D]
						

							set cp1 [::math::linearalgebra::crossproduct $vec_AB $vec_BC]
							set cp1 [::math::linearalgebra::dotproduct $cp1 $cp1]
							set cp2 [::math::linearalgebra::crossproduct $vec_BC $vec_CD]
							set cp2 [::math::linearalgebra::dotproduct $cp2 $cp2]

							set R_BA2 [bond_length_AB_2 $A $B]
							set R_CD2 [bond_length_AB_2 $C $D]

							set k_DC [force_constantAB $D $C]

							set eigenv_DC [open "eigenvalue" "r"]
							set e_DC [read $eigenv_DC]
							close $eigenv_DC
	
							set eigenve_DC [open "eigenvector" "r"]
							set ev_DC [read $eigenve_DC]
							close $eigenve_DC

							set term1 0.0
							set term2 0.0

							set num_val [llength $ev_AB]
							for {set val 0} {$val < $num_val} {incr val} {
								set vec_term1 [lindex $ev_AB $val]
								set dot1 [::math::linearalgebra::dotproduct [lindex $data3 0] $vec_term1]
								set dot1 [expr { abs($dot1) }]
								set term1 [expr { $term1 + ([lindex $e_AB $val] * $dot1) }]
							}

							set num_val [llength $ev_DC]
							for {set val 0} {$val < $num_val} {incr val} {
								set vec_term2 [lindex $ev_DC $val]
								set dot2 [::math::linearalgebra::dotproduct [lindex $data3 1] $vec_term2]
								set dot2 [expr { abs($dot2) }]
								set term2 [expr { $term2 + ([lindex $e_DC $val] * $dot2) }]
							}
				
							set term1 [expr { $term1 * $R_AB2 * $cp1 }]
							set term2 [expr { $term2 * $R_CD2 * $cp2 }]

							set term1 [expr { 1.0 / $term1 }]
							set term2 [expr { 1.0 / $term2 }]
							set k_phi_ABCD [expr { $term1 + $term2 }]
							set k_phi_ABCD [expr { 1.0 / $k_phi_ABCD }]
			
							puts $df "{ $A $B $C $D $k_phi_ABCD }"
						}			 
					}
				incr D_index $nc
				}
			}
		}
	}
	close $af
	close $df
	close $bf	
}

proc branch {A} {

	# STORING PAIRING INFORMATION OF 'A' ATOMS CALLED  'B' ATOMS

	set pair [bonding_info $A]

	set g1 [open "branch" "w"]
	set g [open "dummy" "w"]
	puts $g "$pair"
	puts $g1 "A"
	puts $g1 "{$A}"
	puts $g1 "B"
	puts $g1 "$pair"
	close $g

	set g [open "dummy" "r"]
	set data2 [read $g]
	close $g
		
	set num_pairs_A [llength [lindex $data2 0]]

	# STORING PAIRING INFORMATION OF 'A' ATOM NEIGHBOURS 'B'
	
	for {set j 0} {$j < $num_pairs_A} {incr j} {
		set elem_B($j) [lindex $data2 0 $j]
	}
			
	set g [open "dummy" "w"]

	# STORING PAIR INFORMATION OF 'B' ATOMS CALLED 'C' ATOMS

	puts $g1 "C"
	for {set j 0} {$j < $num_pairs_A} {incr j} {
		set pair_B($j) [bonding_info $elem_B($j)]
		puts $g "$pair_B($j)"
	}
	close $g

	set g [open "dummy" "r"]
	set data2 [read $g]
	close $g

	# STORING PAIRING INFORMATION OF 'B' ATOM NEIGHBOURS 'C'

	for {set j 0} {$j < $num_pairs_A} {incr j} {
		set num_pairs_B($j) [llength [lindex $data2 $j]]
	}
	for {set j 0} {$j < $num_pairs_A} {incr j} {
		for {set d 0} {$d < $num_pairs_B($j)} {incr d} {
			set elem_C($j,$d) [lindex $data2 $j $d]
		}
	}	

	# PREVENTING THE BRANCH FROM GOING BACKWARD

	for {set j 0} {$j < $num_pairs_A} {incr j} {
		set latom($j) ""
		for {set d 0} {$d < $num_pairs_B($j)} {incr d} {
			if { [lindex $data2 $j $d] != $A } {
				set latom($j) [linsert $latom($j) end [lindex $data2 $j $d]]
			}
		}
		puts $g1 "{$latom($j)}"
	} 

	set g [open "dummy" "w"]
	for {set j 0} {$j < $num_pairs_A} {incr j} {
		puts $g "{$latom($j)}"
	}
	close $g

	set g [open "dummy" "r"]
	set data2 [read $g]
	close $g

	# STORING AGAIN THE PAIRING INFORMATION OF 'B' ATOM NEIGHBOURS 'C' WITH ONLY FORWARD MOVEMENT OF THE BRANCH

	for {set j 0} {$j < $num_pairs_A} {incr j} {
		set num_pairs_B($j) [llength [lindex $data2 $j]]
	}
	for {set j 0} {$j < $num_pairs_A} {incr j} {
		for {set d 0} {$d < $num_pairs_B($j)} {incr d} {
			set elem_C($j,$d) [lindex $data2 $j $d]
		}
	}	

	set g [open "dummy" "w"]	
		

	# STORING PAIR INFORMATION OF 'C' ATOMS CALLED 'D' ATOMS

	puts $g1 "D"

	for {set j 0} {$j < $num_pairs_A} {incr j} {
		for {set d 0} {$d < $num_pairs_B($j)} {incr d} {
			set pair_C($j,$d) [bonding_info $elem_C($j,$d)]
			puts $g "$pair_C($j,$d)"
		}
	}
	close $g
		
	set g [open "dummy" "r"]
	set data2 [read $g]
	close $g

	# STORING PAIRING INFORMATION OF 'C' ATOM NEIGHBOURS 'D'

	set dstart 0
	for {set j 0} {$j < $num_pairs_A} {incr j} {
		set d $dstart
		for {set d $dstart} {$d < [expr { $dstart + $num_pairs_B($j) }]} {incr d} {
			set num_pairs_C($j,$d) [llength [lindex $data2 $d]]
		}
		set dstart [expr { $dstart + $num_pairs_B($j) }]
	}

	# PREVENTING THE BRANCH FROM GOING BACKWARD :: REMOVING ANY OVERLAP WITH NODE B

	set dstart 0
	for {set j 0} {$j < $num_pairs_A} {incr j} {
		set d $dstart
		for {set d $dstart} {$d < [expr { $dstart + $num_pairs_B($j) }]} {incr d} {
			set latom1 ""
			for {set k1 0} {$k1 < $num_pairs_C($j,$d)} {incr k1} {
				set elem_D($j,$d,$k1) [lindex $data2 $d $k1]
				set count 0
				for {set k2 0} {$k2 < $num_pairs_A} {incr k2} {
					if { $elem_B($k2) == $elem_D($j,$d,$k1) } {
						incr count
					}
				}
				if { $count == 0} {
					set latom1 [linsert $latom1 end [lindex $data2 $d $k1]]
				}
			}
			puts $g1 "{$latom1}"
		}
		set dstart [expr { $dstart + $num_pairs_B($j) }]
	} 
	close $g1	
}

	
proc unique_comb {} {

	# WRITING THE FRCMOD FILE COMPATABLE WITH AMBER

	set f [open "bonded" "r"]
	set data [read $f]
	close $f
	set f1 [open "fc.bonded" "w"]
	puts $f1 "{ [lindex $data 0 0]-[lindex $data 0 1] [expr { [lindex $data 0 2] / 2.0 }] [lindex $data 0 3] }"

	set g [open "angle" "r"]
	set data1 [read $g]
	close $g
	set g1 [open "fc.angle" "w"]
	puts $g1 "{ [lindex $data1 0 0]-[lindex $data1 0 1]-[lindex $data1 0 2] [expr { [lindex $data1 0 3] / 2.0 }] [lindex $data1 0 4] } "

	set h [open "dihedrals" "r"]
	set data2 [read $h]
	close $h
	set h1 [open "fc.dihedrals" "w"]
	puts $h1 "{ [lindex $data2 0 0]-[lindex $data2 0 1]-[lindex $data2 0 2]-[lindex $data2 0 3] 0.0 0.0 0.0 1.0 }"

	# removing the repitation in bond terms

	set nterms [llength $data]
	for {set i 1} {$i < $nterms} {incr i} {
		set term1 [lindex $data $i 0]
		set term2 [lindex $data $i 1]
		set count 0
		for {set j [expr { $i - 1 }]} {$j >= 0} {set j [expr { $j - 1 }]} {
			set cterm1 [lindex $data $j 0] 
			set cterm2 [lindex $data $j 1]
			if { $cterm1 == $term1 && $cterm2 == $term2 } {
				incr count
			} elseif { $cterm1 == $term2 && $cterm2 == $term1 } {
				incr count
			} 
		}
		if { $count == 0 } {
			puts $f1 "{ [lindex $data $i 0]-[lindex $data $i 1] [expr { [lindex $data $i 2] / 2.0 }] [lindex $data $i 3] }"
		}
	}
	close $f1

	# removing the repitation in angle terms

	set nterms [llength $data1]
	for {set i 1} {$i < $nterms} {incr i} {
		set term1 [lindex $data1 $i 0]
		set term2 [lindex $data1 $i 1]
		set term3 [lindex $data1 $i 2]
		set count 0
		for {set j [expr { $i - 1 }]} {$j >= 0} {set j [expr { $j - 1 }]} {
			set cterm1 [lindex $data1 $j 0] 
			set cterm2 [lindex $data1 $j 1]
			set cterm3 [lindex $data1 $j 2]
			if { $cterm1 == $term1 && $cterm2 == $term2 && $cterm3 == $term3 } {
				incr count
			} elseif { $cterm1 == $term1 && $cterm2 == $term3 && $cterm3 == $term2 } {
				incr count
			} elseif { $cterm1 == $term2 && $cterm2 == $term1 && $cterm3 == $term3 } {
				incr count
			} elseif { $cterm1 == $term2 && $cterm2 == $term3 && $cterm3 == $term1 } {
				incr count
			} elseif { $cterm1 == $term3 && $cterm2 == $term1 && $cterm3 == $term2 } {
				incr count
			} elseif { $cterm1 == $term3 && $cterm2 == $term2 && $cterm3 == $term1 } {
				incr count
			} 
		}
		if { $count == 0 } {
			puts $g1 "{ [lindex $data1 $i 0]-[lindex $data1 $i 1]-[lindex $data1 $i 2] [expr { [lindex $data1 $i 3] / 2.0 }] [lindex $data1 $i 4] } "
		}
	}
	close $g1

	# removing the repititions in dihedral terms
	
	set nterms [llength $data2]
	for {set i 1} {$i < $nterms} {incr i} {
		set term1 [lindex $data2 $i 0]
		set term2 [lindex $data2 $i 1]
		set term3 [lindex $data2 $i 2]
		set term4 [lindex $data2 $i 3]
		set count 0
		for {set j [expr { $i - 1 }]} {$j >= 0} {set j [expr { $j - 1 }]} {
			set cterm1 [lindex $data2 $j 0] 
			set cterm2 [lindex $data2 $j 1]
			set cterm3 [lindex $data2 $j 2]
			set cterm4 [lindex $data2 $j 3]
			if { $cterm1 == $term1 && $cterm2 == $term2 && $cterm3 == $term3 && $cterm4 == $term4} {
				incr count
			} elseif { $cterm1 == $term1 && $cterm2 == $term2 && $cterm3 == $term4 && $term4 == $cterm3 } {
				incr count
			} elseif { $cterm1 == $term1 && $cterm2 == $term3 && $cterm3 == $term2 && $term4 == $cterm4 } {
				incr count
			} elseif { $cterm1 == $term1 && $cterm2 == $term3 && $cterm3 == $term4 && $term4 == $cterm2 } {
				incr count
			} elseif { $cterm1 == $term1 && $cterm2 == $term4 && $cterm3 == $term2 && $term4 == $cterm3 } {
				incr count
			} elseif { $cterm1 == $term1 && $cterm2 == $term4 && $cterm3 == $term3 && $term4 == $cterm2 } {
				incr count
			} elseif { $cterm1 == $term2 && $cterm2 == $term1 && $cterm3 == $term3 && $term4 == $cterm4 } {
				incr count
			} elseif { $cterm1 == $term2 && $cterm2 == $term1 && $cterm3 == $term4 && $term4 == $cterm3 } {
				incr count
			} elseif { $cterm1 == $term2 && $cterm2 == $term3 && $cterm3 == $term1 && $term4 == $cterm4 } {
				incr count
			} elseif { $cterm1 == $term2 && $cterm2 == $term3 && $cterm3 == $term4 && $term4 == $cterm1 } {
				incr count
			} elseif { $cterm1 == $term2 && $cterm2 == $term4 && $cterm3 == $term1 && $term4 == $cterm3 } {
				incr count
			} elseif { $cterm1 == $term2 && $cterm2 == $term4 && $cterm3 == $term3 && $term4 == $cterm1 } {
				incr count
			} elseif { $cterm1 == $term3 && $cterm2 == $term1 && $cterm3 == $term2 && $term4 == $cterm4 } {
				incr count
			} elseif { $cterm1 == $term3 && $cterm2 == $term1 && $cterm3 == $term4 && $term4 == $cterm4 } {
				incr count
			} elseif { $cterm1 == $term3 && $cterm2 == $term2 && $cterm3 == $term4 && $term4 == $cterm4 } {
				incr count
			} elseif { $cterm1 == $term3 && $cterm2 == $term2 && $cterm3 == $term1 && $term4 == $cterm4 } {
				incr count
			} elseif { $cterm1 == $term3 && $cterm2 == $term4 && $cterm3 == $term1 && $term4 == $cterm4 } {
				incr count
			} elseif { $cterm1 == $term3 && $cterm2 == $term4 && $cterm3 == $term2 && $term4 == $cterm4 } {
				incr count
			} elseif { $cterm1 == $term4 && $cterm2 == $term1 && $cterm3 == $term2 && $term4 == $cterm3 } {
				incr count
			} elseif { $cterm1 == $term4 && $cterm2 == $term1 && $cterm3 == $term3 && $term4 == $cterm2 } {
				incr count
			} elseif { $cterm1 == $term4 && $cterm2 == $term2 && $cterm3 == $term1 && $term4 == $cterm3 } {
				incr count
			} elseif { $cterm1 == $term4 && $cterm2 == $term2 && $cterm3 == $term3 && $term4 == $cterm1 } {
				incr count
			} elseif { $cterm1 == $term4 && $cterm2 == $term3 && $cterm3 == $term2 && $term4 == $cterm1 } {
				incr count
			} elseif { $cterm1 == $term4 && $cterm2 == $term3 && $cterm3 == $term1 && $term4 == $cterm2 } {
				incr count
			} 
		}
		if { $count == 0 } {
			puts $h1 "{ [lindex $data2 $i 0]-[lindex $data2 $i 1]-[lindex $data2 $i 2]-[lindex $data2 $i 3] 0.0 0.0 0.0 1.0 } "
		}
	}
	close $h1
}

proc amber {} {

	# IN THIS PROCEDURE WE WILL WRITE THE PARAMETERS IN AMBER'S FORMAT

	package require math::linearalgebra

	set pi 3.14

	amber_atom_names

	set f1 [open "results.bond_fc" "r"]
	set data1 [read $f1]
	close $f1
	set ff1 [open "bonded" "w"]

	set f2 [open "results.angle_fc" "r"]
	set data2 [read $f2]
	close $f2
	set ff2 [open "angle" "w"]
		
	set f3 [open "results.dihedral_fc" "r"]
	set data3 [read $f3]
	close $f3
	set ff3 [open "dihedrals" "w"]

	set g [open "amber_atom_types" "r"]
	set data [read $g]
	close $g

	# CHANGING THE BONDED TERMS

	set nterms [llength $data1]

	for {set i 0} {$i < $nterms} {incr i} {
		set term1 [lindex $data1 $i 0]	
		set term2 [lindex $data1 $i 1]
		set ab_len [bond_length_AB_2 $term1 $term2]
		set ab_len [expr { sqrt($ab_len) }]
		set newt1 [lindex $data 0 [expr { $term1 - 1 }]]
		set newt2 [lindex $data 0 [expr { $term2 - 1 }]]
		
		puts $ff1 "{ $newt1 $newt2 [lindex $data1 $i 2] $ab_len }"
	}

	# CHANGING THE ANGLE TERMS

	set nterms [llength $data2]

	for {set i 0} {$i < $nterms} {incr i} {
		set term1 [lindex $data2 $i 0]	
		set term2 [lindex $data2 $i 1]
		set term3 [lindex $data2 $i 2]
		set eq_ang1 [AB_vectors $term1 $term2]
		set eq_ang2 [AB_vectors $term3 $term2]
		set eq_ang [::math::linearalgebra::dotproduct $eq_ang1 $eq_ang2]
		set eq_ang [expr { ((acos($eq_ang)) * 180.0) / $pi }]
		set newt1 [lindex $data 0 [expr { $term1 - 1 }]]
		set newt2 [lindex $data 0 [expr { $term2 - 1 }]]
		set newt3 [lindex $data 0 [expr { $term3 - 1 }]]
		
		puts $ff2 "{ $newt1 $newt2 $newt3 [lindex $data2 $i 3] $eq_ang }"
	}

	# CHANGING THE DIHEDRAL TERMS

	set nterms [llength $data3]

	set di [open "harmonic_dih" "w"]

	for {set i 0} {$i < $nterms} {incr i} {
		set term1 [lindex $data3 $i 0]	
		set term2 [lindex $data3 $i 1]
		set term3 [lindex $data3 $i 2]
		set term4 [lindex $data3 $i 3]
		set dihe [plane_ABC $term1 $term2 $term3 $term4]

		set dum [open "dummy" "w"]
		puts $dum "$dihe"
		close $dum

		set dum [open "dummy" "r"]
		set datdum [read $dum]
		close $dum
 
		set eq_dh1 [lindex $datdum 0]
		set eq_dh2 [lindex $datdum 1]

		set eq_dh [::math::linearalgebra::dotproduct $eq_dh1 $eq_dh2]
		set eq_dh [expr { ((acos($eq_dh)) * 180.0) / $pi }]

		set newt1 [lindex $data 0 [expr { $term1 - 1 }]]
		set newt2 [lindex $data 0 [expr { $term2 - 1 }]]
		set newt3 [lindex $data 0 [expr { $term3 - 1 }]]
		set newt4 [lindex $data 0 [expr { $term4 - 1 }]]
		
		puts $ff3 "{ $newt1 $newt2 $newt3 $newt4 [lindex $data3 $i 4] $eq_dh }"
		puts $di "{ [lindex $data3 $i 4] $eq_dh }"
	}	
	close $ff1
	close $ff2
	close $ff3
	close $di
}

proc amber_frcmod {} {
	
	# CREATING THE FRCMOD FILE FOR AMBER INPUT

	puts "		**** ENTERING PROCEDURE amber_frcmod ****		"

	set f1 [open "fc.bonded" "r"]
	set data1 [read $f1]
	close $f1

	set f2 [open "fc.angle" "r"]
	set data2 [read $f2]
	close $f2

	set f3 [open "fc.dihedrals" "r"]
	set data3 [read $f3]
	close $f3

	set g [open "frcmod_vf" "w"]

	puts $g "Remark line goes here"

	# PUTTING MASS IN

	puts $g "MASS"

	puts $g ""

	# PUTTING BONDED PARAMETER IN

	# SPACE VARIABLES

	set sp(5) ""
	set sp(4) " "
	set sp(3) "  "
	
	puts $g "BOND"
	
	set num_term [llength $data1]

	for {set i 0} {$i < $num_term} {incr i} {
		set stl [string length [lindex $data1 $i 0]]
		puts $g "[lindex $data1 $i 0]$sp($stl)			[format "%.2f" [lindex $data1 $i 1]]			[format "%.2f" [lindex $data1 $i 2]]"
	}
	puts $g ""

	# PUTTING ANGLE PARAMETER IN

	# SPACE VARIABLES

	set sp(8) ""
	set sp(7) " "
	set sp(6) "  "
	set sp(5) "   "
	
	puts $g "ANGLE"
	
	set num_term [llength $data2]

	for {set i 0} {$i < $num_term} {incr i} {
		set stl [string length [lindex $data2 $i 0]]
		puts $g "[lindex $data2 $i 0]$sp($stl)			[format "%.2f" [lindex $data2 $i 1]]			[format "%.2f" [lindex $data2 $i 2]]"
	}
	puts $g ""

	# PUTTING DIHEDRAL PARAMETER IN

	# SPACE VARIABLES

	set sp(11) ""
	set sp(10) " "
	set sp(9) "  "
	set sp(8) "   "
	set sp(7) "    "
	
	puts $g "DIHE"
	
	set num_term [llength $data3]

	for {set i 0} {$i < $num_term} {incr i} {
		set stl [string length [lindex $data3 $i 0]]
		puts $g "[lindex $data3 $i 0]$sp($stl)			[format "%.2f" [lindex $data3 $i 1]]			[format "%.2f" [lindex $data3 $i 2]]			[format "%.2f" [lindex $data3 $i 3]]			[format "%.2f" [lindex $data3 $i 4]]"
	}
	puts $g ""

	# PUTTING IMPROPER TERMS IN

	puts $g "IMPROPER"
	puts $g ""

	# PUTTING NON BONDED PARAMETERS IN

	puts $g "NONBOND"
	puts $g ""

	close $g
}

proc amber_atom_names {} {

	# IN THIS SECTION WE WILL STORE THE AMBER NAMES OF THE ATOMS IN A FILE NAMED AMBER NAMES

	set f [open "[input3]" "r"]
	set data [read $f]
	close $f
	set k 0

	while { [lindex $data $k] != "@<TRIPOS>ATOM" } {
		incr k
	}
	incr k

	set at_list ""
	set j 1
	while { [lindex $data $k] != "@<TRIPOS>BOND" } {
		if { [lindex $data $k] == "$j" } {
			set at_list [linsert $at_list end [lindex $data [expr {$k + 5}]]]
			incr j
		}
		incr k
	}
	set g [open "amber_atom_types" "w"]
	puts $g "[list $at_list]"
	close $g
}

proc input1 {} {
	set f  [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 0]
}

proc input2 {} {
	set f  [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 1]
}

proc input3 {} {
	set f  [open "input" "r"]
	set data [read $f]
	close $f

	return [lindex $data 2]
}

# IN THE NEXT 2 PROCEDURES WE WILL FIT THE HARMONIC POTENTIAL TO THE COSINE POTENTIAL IMPLEMENTED IN AMBER AND OTHER SIMULATION PACKAGES AS THE POTENTIAL TO MODEL DIHEDRALS 

proc data_generation {} {
	
	# IN THIS PROCEDURE WE WILL FIT THE HARMONIC POTENTIAL DETERMINED ABOVE TO THE CORRESPONDING COSINE POTENTIAL :: V = K(1+cos(nx-xeq)) for x E (0,pi)

	puts "	***** CHANGING THE DIHEDRALS TO THE AMBER READIBLE FORMAT *****	"

	package require math::linearalgebra

	set pi 3.14

	set f [open "harmonic_dih" "r"]
	set data [read $f]
	close $f

	set h [open "cosine_dh" "w"]


	# GENERATING THE DATA POINTS FOR EACH DIHEDRAL USING HARMONIC POTENTIAL :: 1/2*k*(x-xeq)^2

	set num_dh [llength $data]

	for {set i 0} {$i < $num_dh} {incr i} {
		set g [open "fitting_data" "w"]
		for {set j 0} {$j < 360} {incr j 2} {
			set x $j
			set k [lindex $data $i 0]
			set xeq [lindex $data $i 1]
			set x_xeq [expr { $x - $xeq }]
			set x_xeq2 [expr { $x_xeq * $x_xeq }]
			set E [expr { ($k * $x_xeq2) / (2.0 * 57.29 * 57.29) }]
			puts $g "$x,$E" 
		}
		close $g
		set parm [fitting]
		puts $h "$parm"
	}
	close $h		
}
	 
proc fitting {} {
	
	# FITTING THE DATA FROM HARMONIC TO COSINE DIHEDRAL USING Levenberg–Marquardt algorithm

	set f [open "fitting_data" "r"]
	set data [read $f]
	close $f

	# INITIAL GUESS FOR K AND n

	set K 10.0
	set n 2.0
	
}
	
proc delete_files {} {
	
	# THIS PROCEDURE WILL DELETE ALL THE FILES REQUIRED DURING CALCULATIONS ONLY

	file delete branch
	file delete atomic_mass
	file delete cartesian_coord
	file delete cartesian_fc_matrix
	file delete dummy
	file delete dummy1
	file delete eigenvalue
	file delete dummy_matrix
	file delete eigenvector
	file delete cartesian_force_constant
	file delete cart_fc_matrix_full
	file delete results.bond_fc
	file delete results.angle_fc
	file delete results.dihedral_fc
	file delete amber_atom_types
	file delete bonded
	file delete angle
	file delete dihedrals
	file delete harmonic_dih
	file delete fc.bonded
	file delete fc.angle
	file delete fc.dihedrals
}

input_files
hessian_full
seminario
amber
unique_comb
amber_frcmod
delete_files
#data_generation
