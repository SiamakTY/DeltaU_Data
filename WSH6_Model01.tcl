
# This model simulates the nonlinear behavior of WSH6 experimental wall tested by 
# Dazio et al. (2009) under quasi-static cyclic loads.
   
# Elements used in this model: 
# a) Zero length element for the base shear spring
# b) Displacement-based beam-column elements for the unsupported section of the wall
# c) Rigid element (elasticBeamColumn with a large elasticity modulus) for the wall-loading beam intersection

# Number of integration points per each beam-column element: 2 
# The utilized integration technique: Gauss-Legendre

# Material models: 
# Reinforcing Steel Material + MinMax wrapper material for the longitudinal bars (the associated effect of buckling & low-cycle fatigue is taken into account)
# Concrete04 (Popovics Concrete) material model for plain & confined concrete
# 
# Buckling length has been approximated using the approach by (Dhakal & Maekawa, 2002).
# Coffin-Manson fatigue model parameters (Cf & alpha) were evaluated employing a two-step approach: 
# i) Material constants for the Koh-Stephens (1991) fatigue model were estimated using the empirical method developed by (Tripathi et al., 2018).
# ii) The Coffin-Manson material constants were subsequently approximated using the empirical expression by (Brown & Kunnath, 2000).
#
# Number of generated nodes : 7
# Nodes 1 and 2 are located at coordinates (0,0), and are used for generating the shear spring at the wall base.
# Nodes 6 and 7 are used to generate the top rigid element.
# Number of generated elements (including the zero-length element) : 6

# Local y and z axes:
#	   .----------------------------------------------.
#	   |                                              |
# West |          y <----------.                      | East
#(Left)|                       |                      |(Right)
#	   '-----------------------|----------------------'
#                              |
#                              z
# For more information, see sectionCoordinateGuide.pdf

# The simulated hysteretic behavior may be compared with the experimental values using the lateral load-displacement data provided in the WSH6_measured.csv file.

# Code developer: Siamak TAHAEI YAGHOUBI 
# Contact e-mail address: tahaeiyaghoubi@itu.edu.tr
# Last Modification Date: 24.11.2025

wipe
set resultsDir "C:/.../" ; # Path to the directory of output files (must be modified by the user)
file mkdir $resultsDir
set pi 3.1415926535897931


logFile $resultsDir/logFile.out

# Analysis parameters
set numIncrs 150 ; # Number of displacement steps to achieve the target displacement by the loading protocol
set tol 1.0e-6 ; # Tolerance value for the convergence test
set factorLP 2.0 ; # A coefficient multiplied by the plastic hinge length to result in the length of the first (i.e., base) displacement-based element
set intType "Legendre" ; # The utilized numerical integration technique 
set numIntgrPts 2 ; # Number of integration points along each displacement-based element
# disps: Target displacement values (in meters) for the cyclic pushover analysis, specified by the loading protocol:
set disps [list 0.0094934  -0.0085847    0.010402  -0.0098859     0.02578   -0.032152    0.025413   -0.026388    0.037778   -0.039238    0.044824   -0.039611     0.05121   -0.051527    0.051833   -0.050955    0.064699   -0.064479    0.065321    -0.06381    0.077615   -0.078141    0.077579   -0.077794    0.094703   -0.090815]

#Geometry and section parameters
set wallWidth 150.0e-3 ; # Wall thickness (m)
set wallDepth 2000.0e-3 ; # Wall length (m)
set Hcl 3.990 ; # Clear height: the unsupported height of the wall (from top of the base block to the bottom of the loading beam) in meters 
set Heff 4.52 ; # Effective height: the distance from base block top to the approximate impact point of the lateral load in meters
set covery 21.0e-3 ; # Concrete cover in the local y (in-plane) direction (m)
set coverz 21.0e-3 ; # Concrete cover in the local z (out-of-plane) direction (m)
set hoopLength 341.0e-3 ; # Edge confinement lateral reinforcement length in the local y (i.e., in-plane) direction (m)
set Lsp 0.135 ; # Strain penetration length in meters, per (Priestley et al, 2007)
set LP 0.6706 ; # Length of the plastic hinge in meters, per(Kazaz, 2013)
set HcLsp [expr $Hcl + $Lsp] ; # Clear/unsupported wall height + length of strain penetration

# Define gravity loads
set P 1476.0e3 ; # The applied axial load value in Newtons (derived from (Dazio et al., 2009))
set PW [expr 1.0*($wallWidth*$wallDepth*$Hcl + 0.681)*2500.0*9.81] ; # Weight of the wall + the loading beam in Newtons, where 0.681 (m3) is the approximate volume of the loading beam.

# Bar buckling parameters
set bLength1 [expr 2.0*50.0e-3] ; # Buckling length (m) of edge phi 12 bars as per the method by (Dhakal & Maekawa, 2002). Note that edge lateral reinforcement spacing is 50 mm.
set bLength2 [expr 15.0*8.0e-3] ; # Assumed buckling length (m) of web phi 8 bars
set bLength3 [expr 50.0e-3] ; # Buckling length (m) of edge phi 8 bars as per the method by (Dhakal & Maekawa, 2002). Note that edge lateral reinforcement spacing is 50 mm.

# Low cycle fatigue parameters
set Cf1 0.25658 ; # Fatigue ductility coefficient (edge phi 12 bars)
set alpha1 0.51509 ; # Fatigue ductility exponent (edge phi 12 bars)
set Cd1 [expr 1.25*$Cf1] ; # Cyclic strength reduction constant approximated as 1.25*fatigue ductility coefficient (edge phi 12 bars) 
set Cf2 0.16982 ; # Fatigue ductility coefficient (web phi 8 bars)
set alpha2 0.59501 ; # Fatigue ductility exponent (web phi 8 bars)
set Cd2 [expr 1.25*$Cf2]; # Cyclic strength reduction constant approximated as 1.25*fatigue ductility coefficient (web phi 8 bars)
set Cf3 0.27143 ; # Fatigue ductility coefficient (edge phi 8 bars)
set alpha3 0.50734 ; # Fatigue ductility exponent (edge phi 8 bars)
set Cd3 [expr 1.25*$Cf3]; # Cyclic strength reduction constant approximated as 1.25*fatigue ductility coefficient (edge phi 8 bars)

# Recommended values are used for the following parameters based on (https://opensees.berkeley.edu/wiki/index.php/Reinforcing_Steel_Material): 		
set beta1 1.0 ; # Amplification factor for the buckled stress strain curve (edge phi 12 bars)
set r1    0.4 ; # Buckling reduction factor (edge phi 12 bars)
set gamma1 0.5 ; # Buckling constant (edge phi 12 bars)
set beta2 $beta1 ; # Amplification factor for the buckled stress strain curve (web phi 8 bars)
set r2    $r1 ;    # Buckling reduction factor (web phi 8 bars)     
set gamma2 $gamma1 ; # Buckling constant (web phi 8 bars)
set beta3 $beta1 ; # Amplification factor for the buckled stress strain curve (edge phi 8 bars)
set r3    $r1 ;    # Buckling reduction factor (edge phi 8 bars)     
set gamma3 $gamma1 ; # Buckling constant (edge phi 8 bars)

# Other steel properties
set dbar1 12.0e-3 ; # Diameter of edge phi 12 bars (m)
set dbar2 8.0e-3 ; # Diameter of phi 8 bars (m)
set fy1 576.0e6 ;  # Yield stress of phi 12 bars (Pascals)
set fy2 583.7e6 ;  # Yield stress of phi 8 bars (Pascals)
set E1  235.0e9 ;  # Young's modulus of phi 12 bars (Pascals); derived from (Dazio et al., 1999, sec. 2.5.4)
set E2  201.35e9 ;  # Young's modulus of phi 8 bars (Pascals); derived from (Dazio et al., 1999, sec. 2.5.4)
set eRup1 0.075 ; # Strain corresponding to tensile strength of phi 12 bars (mm/mm); derived from (Dazio et al., 1999, sec. 2.5.4) 
set eRup2 0.080 ; # Strain corresponding to tensile strength of phi 8 bars (mm/mm); derived from (Dazio et al., 1999, sec. 2.5.4) 
set fu1 674.9e6 ;  # Tensile strength of phi 12 bars (Pascals)  
set fu2 714.4e6 ;  # Tensile strength of phi 8 bars (Pascals) 
set Esh1 5763.4033e6 ; # Tangent at initial strain hardening of phi 12 bars, estimated following the method by (Yun & Gardner, 2017)
set Esh2 6308.5151e6 ; # Tangent at initial strain hardening of phi 8 bars, estimated following the method by (Yun & Gardner, 2017)
set esh1 0.03 ;   # Strain at initiation of strain hardening curve of phi 12 bars, estimated following the method by (Yun & Gardner, 2017)
set esh2 0.0267 ; # Strain at initiation of strain hardening curve of phi 8 bars, estimated following the method by (Yun & Gardner, 2017)

# Concrete material properties
set fc 45.6e6 ; # Plain concrete compressive strength (Pascals)
set fcc 54.2369e6 ; # Confined concrete compressive strength in Pascals per (Yassin, 1994)
set ep0TH -0.002 ; # Theoretical peak strain of plain concrete per (Yassin, 1994)
set epc0TH -0.0023788 ; # Theoretical peak strain of confined concrete per (Yassin, 1994)
set epuTH -0.0034255 ; #  Theoretical ultimate strain of plain concrete per (Yassin, 1994)
set epcuTH -0.032935 ; # Theoretical ultimate strain of confined concrete per (Yassin, 1994)
set ft 2.2429e6 ; # Concrete tensile strength in Pascals per (Wong & Vecchio, 2006)
set Ec 31960.7358e6 ; # Initial Young's modulus of concrete in Pascals per (ACI 318-19)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Model calibration process
set epc0M -0.0029850 ; # Modified peak strain of the CONFINED concrete. Can be tuned for better displacement capacity calibration. 
set ep0M [expr ($epc0M/$epc0TH)*$ep0TH] ; # Evaluating the modified peak strain of PLAIN concrete using the modified peak strain of confined concrete.
set epuM [expr $epuTH + ($ep0M - $ep0TH)] ; #  Evaluating the modified ultimate strain of PLAIN concrete.
set epcuM [expr $epcuTH + ($epc0M - $epc0TH)] ; # Evaluating the modified ultimate strain of CONFINED concrete.
#  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Continuing the analysis using the calibrated concrete parameters (i.e., ep0M, epuM, epc0M, epcuM)
puts "Weight of the wall is estimated as [expr 0.001*ceil(1000.0*0.001*$PW)] kN "
set ALR [expr 0.01*ceil(100.0*100.0*($P + $PW)/($wallWidth*$wallDepth*$fc))] ; # Axial load ratio
puts "Axial load ratio \[(applied vertical load + wall weight)/gross area\] is: $ALR%"

set clHDiv [expr ceil(($HcLsp - $factorLP*$LP)/($factorLP*$LP))]; # Determining the number of elements that the clear wall height must be divided into.
puts "The unsupported (i.e., clear) wall height has been divided into [expr $clHDiv + 1] elements" 
set intersectPoint [expr int($clHDiv + 3)] ; # The point RC wall intersects with the loading beam. 
#puts "The RC wall intersects with the loading beam @ node: $intersectPoint" ; # Activate this if you wish!
set impactNode [expr int($intersectPoint + 1)] ; # The point lateral load is applied to the wall, for example, using hydraulic jacks.

# Printing the calibrated concrete properties
puts " " ; # Empty line for the sake of neatness
puts "Calibrated peak strain of plain concrete is: [expr 1.0e-6*ceil(1.0e6*$ep0M)]"
puts "Calibrated ultimate strain value of plain concrete is [expr 1.0e-6*ceil(1.0e6*$epuM)]"
puts "Calibrated peak strain of confined concrete is: [expr 1.0e-6*ceil(1.0e6*$epc0M)]"
puts "Calibrated ultimate strain value of confined concrete is [expr 1.0e-6*ceil(1.0e6*$epcuM)]" 
puts " " ; # Empty line for the sake of neatness
puts [string repeat "=" 35] ; # For the sake of neatness 

# Create ModelBuilder (with two-dimensions and 3 DOF/node)
model basic -ndm 2 -ndf 3
# Create nodes
#    tag        X       Y
node  1       0.0     0.0
node  2		  0.0	  0.0
puts "A shear spring will be defined using Nodes 1 & 2 @ (0.0,0.0) coordinates" 
set node3Height [expr $factorLP*$LP]
puts "Height of node 3 is: [expr 0.001*ceil(1000.0*$node3Height)] (m)"
node  3       0.0    $node3Height 
for {set ii 4} {$ii <= [expr $intersectPoint]} {incr ii} {
	set nodeHeight [expr $node3Height + ($ii - 3)*(($HcLsp - $node3Height)/$clHDiv)]
	puts "Height of node $ii is: [expr 0.001*ceil(1000.0*$nodeHeight)] (m)"
	node $ii  0.0 $nodeHeight
}
node $ii 0.0 [expr $Lsp + $Heff] ;
puts "Height of the lateral loading impact point + strain penetration length (node $ii) is [expr 0.001*ceil(1000.0*($Lsp + $Heff))] (m)"
puts [string repeat "-" 25]

# Fix support at the base of the wall
#    tag   DX   DY   RZ
fix   1     1    1    1

# Define materials for the nonlinear wall model
# ---------------------------------------------
# CONCRETE
# uniaxialMaterial Concrete04    $matTag   $fc    $ec     $ecu   $Ec   <$fct $et>   <$beta>
uniaxialMaterial   Concrete04		1     -$fcc  $epc0M  $epcuM  $Ec    $ft  0.00045  0.1 ; # Confined concrete
uniaxialMaterial   Concrete04		2     -$fc   $ep0M   $epuM   $Ec    $ft  0.00045  0.1 ; # Plain concrete

set lsr1 [expr $bLength1/$dbar1] ; # Bar slenderness (edge phi 12 bars)
set lsr2 [expr $bLength2/$dbar2] ; # Bar slenderness (web phi 8 bars)
set lsr3 [expr $bLength3/$dbar2] ; # Bar slenderness (edge phi 8 bars)

# STEEL
# uniaxialMaterial ReinforcingSteel $matTag $fy  $fu  $Es $Esh  $esh  $eult < -GABuck $lsr $beta $r $gama > < -DMBuck $lsr < $alpha >> < -CMFatigue $Cf $alpha $Cd > < -IsoHard <$a1 <$limit> > > 
uniaxialMaterial ReinforcingSteel   3     $fy1 $fu1 $E1 $Esh1 $esh1 $eRup1 -GABuck $lsr1 $beta1 $r1 $gamma1 -CMFatigue $Cf1 $alpha1 $Cd1 ; # edge phi 12 bars
uniaxialMaterial ReinforcingSteel   4     $fy2 $fu2 $E2 $Esh2 $esh2 $eRup2 -GABuck $lsr2 $beta2 $r2 $gamma2 -CMFatigue $Cf2 $alpha2 $Cd2 ; # web phi 8 bars
uniaxialMaterial ReinforcingSteel   5     $fy2 $fu2 $E2 $Esh2 $esh2 $eRup2 -GABuck $lsr3 $beta3 $r3 $gamma3 -CMFatigue $Cf3 $alpha3 $Cd3 ; # edge phi 8 bars

uniaxialMaterial MinMax  6       3    -min  -1.0e16   -max   $eRup1 ; # edge phi 12 bars
uniaxialMaterial MinMax  7       4    -min  -1.0e16   -max   $eRup2 ; # web phi 8 bars
uniaxialMaterial MinMax  8       5    -min  -1.0e16   -max   $eRup2 ; # edge phi 8 bars

# Shear spring
set shearStiffness [expr 0.02*$Ec*$wallWidth*$wallDepth] ; # The effective shear stiffness per PEER/ATC 72-1
# Materials 9 & 10 are used for defining the zero length shear spring that works only in the lateral load direction (i.e., global X direction)
# uniaxialMaterial Elastic   $matTag      $E             <$eta> <$Eneg> 
uniaxialMaterial   Elastic      9         $shearStiffness
uniaxialMaterial   Elastic      10        1.0e16 ; # To constrain the spring in vertical and out-of-plane directions (i.e., global Y & Z directions)

# Define cross-section for the nonlinear wall
# ------------------------------------------
# Bar area values derived from the introduced bar diameters
set As1   [expr $pi*($dbar1**2)/4.0];     # Area of phi 12 bars
set As2   [expr $pi*($dbar2**2)/4.0];     # Area of phi 8 bars

# Defining some auxiliary points for generating the concrete fibers
# puts " "
set i1y [expr $wallDepth/2.0 - $covery - $hoopLength] ;
set i1z [expr -1.0*($wallWidth/2.0 - $coverz)] ;
# puts "i1 point coordinates: ($i1y,$i1z)" ; # Activate to see the coordinates on screen

set j1y [expr $wallDepth/2.0 - $covery] ;
set j1z [expr -1.0*($wallWidth/2.0 - $coverz)] ;
# puts "j1 point coordinates: ($j1y,$j1z)" ; # Activate to see the coordinates on screen

set k1y [expr $wallDepth/2.0 - $covery] ;
set k1z [expr $wallWidth/2.0 - $coverz] ;
# puts "k1 point coordinates: ($k1y,$k1z)" ; # Activate to see the coordinates on screen

set l1y [expr $wallDepth/2.0 - $covery - $hoopLength] ;
set l1z [expr $wallWidth/2.0 - $coverz] ;
# puts "l1 point coordinates: ($l1y,$l1z)" ; # Activate to see the coordinates on screen

# puts " "
set i2y [expr -1.0*$j1y] ;
set i2z  $i1z 
# puts "i2 point coordinates: ($i2y,$i2z)" ; # Activate to see the coordinates on screen

set j2y [expr -1.0*$i1y] ;
set j2z  $j1z 
# puts "j2 point coordinates: ($j2y,$j2z)" ; # Activate to see the coordinates on screen

set k2y [expr -1.0*$l1y] ;
set k2z  $k1z;
# puts "k2 point coordinates: ($k2y,$k2z)" ; # Activate to see the coordinates on screen

set l2y [expr -1.0*$k1y] ;
set l2z  $l1z ; 
# puts "l2 point coordinates: ($l2y,$l2z)" ; # Activate to see the coordinates on screen

# puts " "
set i3y [expr -1.0*$wallDepth/2.0] ;
set i3z [expr -1.0*$wallWidth/2.0] ;
# puts "i3 point coordinates: ($i3y,$i3z)" ; # Activate to see the coordinates on screen

set j3y [expr -1.0*$i3y] ;
set j3z  $i3z ;
# puts "j3 point coordinates: ($j3y,$j3z)" ; # Activate to see the coordinates on screen

set k3y [expr -1.0*$i3y] ;
set k3z [expr -1.0*($wallWidth/2.0 - $coverz)]
# puts "k3 point coordinates: ($k3y,$k3z)" ; # Activate to see the coordinates on screen

set l3y $i3y ;
set l3z $k3z
# puts "l3 point coordinates: ($l3y,$l3z)" ; # Activate to see the coordinates on screen

# puts " "
set i4y  $i3y ;
set i4z  [expr -1.0*$l3z] ;
# puts "i4 point coordinates: ($i4y,$i4z)" ; # Activate to see the coordinates on screen

set j4y  $j3y ;
set j4z  [expr -1.0*$k3z] ;
# puts "j4 point coordinates: ($j4y,$j4z)" ; # Activate to see the coordinates on screen

set k4y  $k3y ;
set k4z  [expr -1.0*$j3z] ;
# puts "k4 point coordinates: ($k4y,$k4z)" ; # Activate to see the coordinates on screen

set l4y $l3y ;
set l4z [expr -1.0*$i3z] ;
# puts "l4 point coordinates: ($l4y,$l4z)" ; # Activate to see the coordinates on screen

set webDepth [expr $wallDepth - 2.0*$covery - 2.0*$hoopLength]  
# puts "web width is: $webDepth (m)"; puts "" ; # Activate to see the estimated web depth (plain concrete material model will be used in this region) 

# Creating the fiber section
section Fiber 1 {

		# Create the concrete core fibers
	#   patch quad $matTag $numSubdivIJ $numSubdivJK   $yI $zI     $yJ $zJ         $yK  $zK         $yL $zL
		patch quad    1     32            1          $i1y $i1z     $j1y $j1z       $k1y $k1z        $l1y $l1z
		patch quad    1     32            1          $i2y $i2z     $j2y $j2z       $k2y $k2z        $l2y $l2z

		# Create the concrete cover fibers (top, bottom, left, right, web)
		patch quad    2  [expr int(ceil($wallDepth/($hoopLength/32.0)))]  1      $i3y $i3z     $j3y $j3z       $k3y $k3z        $l3y $l3z
		patch quad    2  [expr int(ceil($wallDepth/($hoopLength/32.0)))]  1      $i4y $i4z     $j4y $j4z       $k4y $k4z        $l4y $l4z
		patch quad    2  [expr int(ceil($covery/($hoopLength/32.0)))]     1      $j1y $j1z     $k3y $k3z       $j4y $j4z        $k1y $k1z
		patch quad    2  [expr int(ceil($covery/($hoopLength/32.0)))]     1      $l3y $l3z     $i2y $i2z       $l2y $l2z        $i4y $i4z
		patch quad    2  [expr int(ceil($webDepth/($hoopLength/32.0)))]   1      $j2y $j2z     $i1y $i1z       $l1y $l1z        $k2y $k2z

		# Create the reinforcing fibres 
		#(Boundary elements)
	#   layer straight $matTag $numFiber $areaFiber $yStart     $zStart    $yEnd      $zEnd
		layer straight 6  3         $As1      770.0e-3   45.0e-3    970.0e-3    45.0e-3  ; # Left edge lower bars (3*phi 12)
		layer straight 6  3         $As1     -970.0e-3   45.0e-3   -770.0e-3    45.0e-3  ; # Right edge lower bars (3*phi 12)
		layer straight 6  3         $As1      770.0e-3  -45.0e-3    970.0e-3   -45.0e-3  ; # Left edge upper bars (3*phi 12)
		layer straight 6  3         $As1	 -970.0e-3  -45.0e-3   -770.0e-3   -45.0e-3  ; # Right edge upper bars (3*phi 12)
		
		layer straight 8  1         $As2	  645.0e-3  47.0e-3    645.0e-3  47.0e-3  ; # Left edge lower bar (1*phi 8)
		layer straight 8  1         $As2	 -645.0e-3  47.0e-3   -645.0e-3  47.0e-3  ; # Right edge lower bar (1*phi 8)
		layer straight 8  1         $As2	  645.0e-3 -47.0e-3    645.0e-3 -47.0e-3  ; # Left edge upper bar (1*phi 8)
		layer straight 8  1         $As2	 -645.0e-3 -47.0e-3   -645.0e-3 -47.0e-3  ; # Right edge upper bar (1*phi 8)
		
		#(Web bars)
		layer straight 7  4         $As2	 145.0e-3   47.0e-3    645.0e-3   47.0e-3  ; # Left lower web bars (4*phi 8)
		layer straight 7  4         $As2	-645.0e-3   47.0e-3   -145.0e-3   47.0e-3  ; # Right lower web bars (4*phi 8)
		layer straight 7  4         $As2	 145.0e-3  -47.0e-3    645.0e-3  -47.0e-3  ; # Left upper web bars (4*phi 8)
		layer straight 7  4         $As2	 -645.0e-3  -47.0e-3   -145.0e-3  -47.0e-3 ; # Right upper web bars (4*phi 8)
		layer straight 7  2         $As2      0.00      -47.0e-3    0.00       47.0e-3 ; # Middle web bars (2*phi 8)
}

# Define the wall element
#                tag
geomTransf Linear 1

#element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2 ...<-doRayleigh $rFlag> <-orient $x1 $x2 $x3 $yp1 $yp2 $yp3>
element zeroLength  1000    1      2      -mat 9 10 10    -dir 1 2 3 ; # Base shear spring that accounts for shear deformations
puts "A zero length element was created between nodes 1 & 2"

# A stack of line displacement-based beam column elements for the HcLsp height:
for {set ii 1} {$ii <= [expr $intersectPoint - 2]} {incr ii} {
	#element dispBeamColumn $eleTag  $iNode         $jNode          $numIntgrPts $secTag $transfTag <-mass $massDens> <-cMass> <-integration $intType>
	 element  dispBeamColumn $ii    [expr $ii + 1]  [expr $ii + 2]  $numIntgrPts   1       1         -integration $intType
	puts "Ele$ii generated between node[expr $ii + 1] & node[expr $ii + 2]"
} 
# The subsequently-defined rigid element (an element with elasticity defined as 25 times the uncracked concrete Young's modulus) 
# simulates the behaviour of RC wall-loading beam connection
#element elasticBeamColumn $eleTag      $iNode                $jNode       $A                            $E               $Iz                                    $transfTag <-mass $massDens> <-cMass>
element elasticBeamColumn  $ii          [expr $impactNode - 1]  $impactNode    [expr $wallWidth*$wallDepth]  [expr 25.0*$Ec]  [expr $wallWidth*($wallDepth**3)/12.0]     1
puts "Rigid Ele$ii generated between Node [expr $impactNode - 1] & Node $impactNode" 

puts [join [lrepeat 2 [string repeat "-" 25]] "\n"]; 

# Define gravity loads
# Create a Plain load pattern with a Linear TimeSeries
timeSeries Linear 1
pattern Plain 1 1 {

        # Create nodal load at node2
	#     nd            FX      FY       MZ
	load  $impactNode   0.0  [expr -$P]  0.0 ; # The applied axial load 
	load  $impactNode   0.0  [expr -$PW] 0.0 ; # Weight of the wall & loading beam
}
# Create the system of equation, a sparse solver with partial pivoting
system BandGeneral

# Create the constraint handler, the transformation method
constraints Transformation

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer RCM

# Create the convergence test, the norm of the residual with a tolerance of
# 1e-6 and a max number of iterations of 10
test NormDispIncr $tol  10  3

# Create the solution algorithm, a Newton-Raphson algorithm
algorithm Newton

# Create the integration scheme, the LoadControl scheme using steps of 0.1
integrator LoadControl 0.1

# Create the analysis object
analysis Static

# Perform the analysis
analyze 10
#====================================================================
# Set the gravity loads to be constant & reset the time in the domain
loadConst -time 0.0

# Difine some recorders
recorder Node -file $resultsDir/TopDisp.out -node $impactNode -dof 1 disp ; # Top displacement
recorder Node -file $resultsDir/RBase.out  -time -node 1 -dof 1 2 reaction ; # Base reactions

# Create recorders to monitor stresses and strains in the wall
recorder Element -file $resultsDir/StWest.out -ele 1 section 1 fiber   970.0e-3 45.0e-3 6 stressStrain ; # Longitudinal steel @ the extreme west
recorder Element -file $resultsDir/StEast.out -ele 1 section 1 fiber  -970.0e-3 45.0e-3 6 stressStrain ; # Longitudinal steel @ the extreme east
recorder Element -file $resultsDir/ConCoreWest.out -ele 1 section 1 fiber [expr $wallDepth/2.0 - $covery]        0.0e-3 1 stressStrain ; # Core/confined concrete @ the extreme west
recorder Element -file $resultsDir/ConCoreEast.out -ele 1 section 1 fiber [expr -1.0*($wallDepth/2.0 - $covery)] 0.0e-3 1 stressStrain ; # Core/confined concrete @ the extreme east
recorder Element -file $resultsDir/ConCoverWest.out -ele 1 section 1 fiber [expr $wallDepth/2.0]                 0.0e-3 2 stressStrain ; # Cover/plain concrete @ the far west
recorder Element -file $resultsDir/ConCoverEast.out -ele 1 section 1 fiber [expr -1.0*$wallDepth/2.0]            0.0e-3 2 stressStrain ; # Cover/plain concrete @ the far east
print -file $resultsDir/modelFile.out

# Cyclic pushover analysis
set patTag 1 ; # pattern tag
set targetDispNo [llength $disps] ; # Number of target displacements (i.e., displacement amplitudes) of the loading protocol

for {set ii 0} {$ii < $targetDispNo} {incr ii} { 
        
		
		set maxU [lindex $disps $ii] ; # Target displacement
		puts "Target delta value = [expr 0.001*ceil(1000.0*1000.0*$maxU)] mm; target drift ratio = [expr 0.001*ceil(1000.0*100.0*$maxU/$Heff)]%"
		set patTag [expr int($patTag + 1)] ;

		set dU [expr $maxU/$numIncrs]
		
		pattern Plain $patTag "Linear" {
					load $impactNode 1.0 0.0 0.0
				}
				
		# integrator DisplacementControl $node           $dof        $incr <$numIter $Î”Umin $Î”Umax>
		integrator DisplacementControl   $impactNode      1           $dU    

		constraints Transformation
		numberer RCM
		system BandGeneral
		test NormDispIncr $tol  1000
		algorithm NewtonLineSearch

		set currentDisp 0.0;
		set ok 0
		
		if {$maxU >= 0.0} {
			while {$ok == 0 && $currentDisp < $maxU} {

				set ok [analyze 1]

				if {$ok != 0} {
					puts "Newton Line Search failed ... trying ModifiedNewton -initial"
					test NormDispIncr $tol  1000 
					algorithm ModifiedNewton -initial
					set ok [analyze 1]
					if {$ok == 0} {
					puts "ModifiedNewton -initial worked ... back to Newton Line Search"
					test NormDispIncr $tol  1000 
					algorithm NewtonLineSearch
					}  
				}
				if {$ok != 0} {
					puts "ModifiedNewton -initial failed also... trying KrylovNewton"
					test NormDispIncr $tol  1000 
					algorithm KrylovNewton
					set ok [analyze 1]
					if {$ok == 0} {
						puts "KrylovNewton worked ... back to Newton Line Search"
						test NormDispIncr $tol  1000 
						algorithm NewtonLineSearch
					}
				}
				if {$ok != 0} {
					puts "KrylovNewton failed also... trying Broyden"
					test NormDispIncr $tol  1000 
					algorithm Broyden 500
					set ok [analyze 1]
					if {$ok == 0} {
						puts "Broyden worked ... back to Newton Line Search"
						test NormDispIncr 1.0e-12  1000 
						algorithm NewtonLineSearch
					}
				
				}
				if {$ok != 0} {
					puts "Broyden failed also... trying BFGS"
					test NormDispIncr [expr 100.0*$tol]  1000 
					algorithm BFGS
					set ok [analyze 1]
					if {$ok == 0} {
						puts "BFGS worked ... back to Newton Line Search"
						test NormDispIncr 1.0e-12  1000 
						algorithm NewtonLineSearch 
					} elseif {$ok != 0} {
					puts "All algorithms FAILED"
					}
				
				}

				set currentDisp [nodeDisp $impactNode  1]
			}
		}
		
		if {$maxU < 0.0} {
			while {$ok == 0 && $currentDisp > $maxU} {

				set ok [analyze 1]

				if {$ok != 0} {
					puts "Newton Line Search failed ... trying ModifiedNewton -initial"
					test NormDispIncr $tol  1000 
					algorithm ModifiedNewton -initial
					set ok [analyze 1]
					if {$ok == 0} {
					puts "ModifiedNewton -initial worked ... back to Newton Line Search"
					test NormDispIncr $tol  1000 
					algorithm NewtonLineSearch
					}  
				}
				if {$ok != 0} {
					puts "ModifiedNewton -initial failed also... trying KrylovNewton"
					test NormDispIncr $tol  1000 
					algorithm KrylovNewton
					set ok [analyze 1]
					if {$ok == 0} {
						puts "KrylovNewton worked ... back to Newton Line Search"
						test NormDispIncr $tol  1000 
						algorithm NewtonLineSearch
					}
				}
				if {$ok != 0} {
					puts "KrylovNewton failed also... trying Broyden"
					test NormDispIncr $tol  1000 
					algorithm Broyden 500
					set ok [analyze 1]
					if {$ok == 0} {
						puts "Broyden worked ... back to Newton Line Search"
						test NormDispIncr $tol  1000 
						algorithm NewtonLineSearch
					}
				}
				if {$ok != 0} {
					puts "Broyden failed also... trying BFGS"
					test NormDispIncr [expr 100.0*$tol]  1000 
					algorithm BFGS
					set ok [analyze 1]
					if {$ok == 0} {
						puts "BFGS worked ... back to Newton Line Search"
						test NormDispIncr 1.0e-12  1000 
						algorithm NewtonLineSearch 
					} elseif {$ok != 0} {
					puts "All algorithms FAILED"
					}
				
				}
				set currentDisp [nodeDisp $impactNode  1]
			}
		}
		

		if {$ok == 0} {
		  puts "Cyclic analysis completed SUCCESSFULLY";
		} else {
		  puts "Cyclic analysis FAILED";    
		}
		print node $impactNode
}

# References:
# ACI Committee 318. (2019). Building Code Requirements for Structural Concrete (ACI 318-19) and Commentary. Farmington Hills, MI, USA: American Concrete Institute.
# Brown, J., & Kunnath, S. K. (2000). Low cycle fatigue behavior of longitudinal reinforcement in reinforced concrete bridge columns. Orlando, Florida, USA: Multidisciplinary Center for Earthquake Engineering Research (MCEER).
# Coffin Jr, L. F. (1954). A study of the effects of cyclic thermal stresses on a ductile metal. Transactions of the American Society of Mechanical engineers, 76(6), 931-949.
# Dazio, A., Beyer, K., & Bachmann, H. (2009). Quasi-static cyclic tests and plastic hinge analysis of rc structural walls. Engineering Structures, 31(7), 1556–1571.
# Dazio, A., Wenk, T., & Bachmann, H. (1999). Versuche an Stahlbetontragwänden unter zyklisch-statischer Einwirkung (Vol. 239). ETH Zurich.
# Dhakal, R. P., & Maekawa, K. (2002b). Reinforcement stability and fracture of cover concrete in reinforced concrete members. Journal of Structural Engineering, 128(10), 1253-1262.
# Kazaz, İ. (2013). Analytical study on plastic hinge length of structural walls. Journal of Structural Engineering, 139(11), 1938-1950.
# Koh, S. K., & Stephens, R. I. (1991). Mean stress effects on low cycle fatigue for a high strength steel. Fatigue Fracture of Engineering Materials and Structures, 14(4), 413-428.
# Manson, S. S. (1953). Behavior of materials under conditions of thermal stress (Vol. 2933). National Advisory Committee for Aeronautics.
# PEER/ATC 72-1. (2010). Modeling and Acceptance Criteria for Seismic Design and Analysis of Tall Buildings. Redwood City, CA, USA: Pacific Earthquake Engineering Research Center, Applied Technology Council.
# Popovics, S. (1973). A numerical approach to the complete stress strain curve for concrete. Cement and Concrete Research, 3(5), 583-599.
# Priestley, M. J., Calvi, G. M., & Kowalsky, M. J. (2007). Displacement Based Seismic Design of Structures (1 ed.). Pavia, Italy: IUSS Press.
# Tripathi, M., Dhakal, R. P., Dashti, F., & Massone, L. M. (2018). Low-cycle fatigue behaviour of reinforcing bars including the effect of inelastic buckling. Construction and Building Materials, 190, 1226–1235.
# Wong and Vecchio, VecTor2 and FormWorks Users Manual. University of Toronto, 2006
# Yassin, M. H. M. (1994). Nonlinear analysis of prestressed concrete structures under monotonic and cyclic loads. University of California, Berkeley.
# Yun, X., & Gardner, L. (2017). Stress-strain curves for hot-rolled steels. Journal of Constructional Steel Research, 133, 36-46