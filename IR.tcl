puts "DO YOU WANT TO VISUALISE NORMAL MODES? Y/N"
set ans1 [gets stdin]

if { $ans1 == "y" || $ans1 == "Y" } {

	puts "		************************************************"
	puts "		AMBER NORMAL MODE VISUALISATION USING GAUSS VIEW"
	puts "		************************************************"

	puts ""

	puts "		BY TARUN KHANNA"
	puts "		IMPERIAL COLLEGE LONDON,U.K."
	puts ""

	puts "		*** MAKE SURE THAT YOU HAVE RUN NMODE ANALYSIS IN AMBER BEFORE USING THIS CODE **"
	puts ""
	puts ""

	puts "		## ENTER THE PRMTOP OF THE MOLECULE"
	set opt1 [gets stdin]
	exec ls $opt1

	puts ""

	puts "		## ENTER THE PDB OF THE MOLECULE"
	set opt2 [gets stdin]
	exec ls $opt2

	puts ""

	set f [open "input" "w"]

	set g [open "$opt1" "r"]
	set data [read $g]
	close $g

	set k 0

	while { [lindex $data $k] != "%FORMAT(10I8)" } {
		incr k
	}
	incr k
	set natom [lindex $data $k]

	exec ls vecs

	puts $f "$natom vecs $opt1 $opt2"

	close $f

	exec tclsh nmode_with_visualisation_v2.tcl | tee log1

	file delete input

	puts "		### DONE: 'amber_freq.log' can now be visualised using GAUSS VIEW (IGNORE ANY MASSAGE WHICH POPS UP ON OPENING THIS FILE WOTH GAUSS VIEW) ###"
	
}


puts "DO YOU WANT TO DERIVE METAL BONDED AND ANGLE PARAMETERS? y/n"
set ans2 [gets stdin]

if { $ans2 == "Y" || $ans2 == "y" } {
	
	puts "		**********************************************************************************"
	puts "		DERIVING METAL PARAMTERS USING THE HESSIAN MATRIX OF GAUSSIAN FREUENCY CALCULATION"
	puts "		\[USING THE SEMINARIO METHOD\]"
	puts "		**********************************************************************************"

	puts ""

	puts "		BY TARUN KHANNA"
	puts "		IMPERIAL COLLEGE LONDON,U.K."
	puts ""

	puts "		*** MAKE SURE THAT YOU HAVE RUN QM FREQUENCY CALCULATION USING GAUSSIAN SOFTWARE PACKAGE BEFORE USING THIS CODE **"
	puts ""
	puts ""

	puts "		## ENTER THE FORMATTED CHECKPOINT FILE (USE formchk COMMAND TO CONVERT THE GAUSSIAN CHECKPOINT FILE INTO FORMATTED CHECKPOINT) "
	set opt1 [gets stdin]
	exec ls $opt1

	puts ""

	puts "		## ENTER THE MOL2 FILE OF THE MOLECULE"
	set opt2 [gets stdin]
	exec ls $opt2

	puts ""

	set f [open "input" "w"]

	set g [open "$opt2" "r"]
	set data [read $g]
	close $g

	set natom [lindex $data 2]

	puts $f "$opt1 $natom $opt2"

	close $f

	exec tclsh force_constants_freq.tcl | tee log2

	file delete input

	puts ""
	puts ""
	puts "		### DONE: 'frcmod_vf' file in the folder contains the metal bond and angle parameters ###"
}






















