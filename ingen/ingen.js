// default values of options
var gle_kind="optimal";

function mkOption(curval,value,label,style)
{
   style = typeof style !== "undefined" ? style : ""
   document.write("<option value='"+value+"' "+(curval==value?"selected='selected'":"")+" style='"+style+"'>"+label+"</option>");
}

//generate a drop down for units
function mkUnit(curval,value,mode)
{
  document.write("<select name='"+value+"' onchange='this.form.submit()'> ");
  switch (mode) {
  case "time":
    mkOption(curval,"ps",   "picosec.");
    mkOption(curval,"fs",   "femtosec.");
    mkOption(curval,"s",    "second");
    mkOption(curval,"aut",   "at.time");
    break;
  case "freq":
    mkOption(curval,"ps1",   "ps^-1");
    mkOption(curval,"thz",   "THz");
    mkOption(curval,"thzrad",   "THz rad");
    mkOption(curval,"cm1",     "cm^-1");
    mkOption(curval,"auw",    "at.freq");
    break;
  case "nrg":
    mkOption(curval,"ev",   "eV");
    mkOption(curval,"k",    "K");
    mkOption(curval,"aue",    "at.energy");
    break;
  }
  document.write("</select>");
}

function mkField(curval,value,label,uval,units)
{
  document.write("<div class='label'>"+label+"</div>&nbsp;<input class='number' type='text' name='"+value+"' value='"+curval+"' onchange='this.form.submit()'></input>");
  mkUnit(uval,"u"+value,units);
  document.write("<br/>\n");
}

function getOption(varname)
{
   var obox=document.getElementById(varname);
   return obox.options[obox.selectedIndex].value;
}


function gleUpdate()
{   
   // Sets the input form data
   gle_kind=getOption("gle_kind");
   gle_in=document.getElementById("gle_input");
   switch (gle_kind)
   { 
       case "optimal":
          httpGet("ingen/gle_optimal.html","gle_input",false);
          break;
       case "quantum":
          gle_in.innerHTML="QUANTUM";
          break;
          
   }
}


function unitname(unit)
{
  switch (unit) {
  case "ps": return "picoseconds";
  case "fs": return "femtoseconds";
  case "s": return "seconds";
  case "aut": return "atomic time units";
  case "ps1": 
  case "thz": return "THz";
  case "thzrad": return "THz*rad";
  case "cm1": return "cm^-1";
  case "auw": return "atomic frequency units";
  case "k": return "K";
  case "ev": return "eV";
  case "aue": return "atomic energy units"; 
  default: return "Unknown unit $unit";
}
}

function conv_u2i(val,unit,mode)
{
  cc=1.;
  switch (mode) {
  case "time":
    switch(unit) {
      case "ps": cc=1.;             break;
      case "fs": cc=1e-3;             break;
      case "s":  cc=1.e12;         break; 
      case "aut": cc=2.4188843e-05;  break;
      default: document.write( " ERROR! Wrong unit for mode time <br/>" );
    }
  break;
  case "freq":  
    switch(unit) {
      case "ps1": cc=1.;             break;
      case "thz": cc=1.;             break;
      case "thzrad": cc=1./6.2831853;             break;
      case "cm1":  cc=0.029979246;         break; 
      case "auw": cc=6579.6839;  break;
      default: document.write( " ERROR! Wrong unit for mode freq <br/>");
    }
  break;
  case "nrg":  
    switch(unit) {
      case "k": cc=1.;             break;
      case "ev": cc=11604.506;             break;
      case "kcalmol": cc=503.55573;             break;
      case "kjmol":  cc=120.27222;         break; 
      case "aue": cc=315774.67;  break;
      default: document.write( " ERROR! Wrong unit for mode nrg <br/>");
    }
  break;
  }
  return cc*val;   
}

function conv_i2u(val,unit,mode)
{
  cc=conv_u2i(1.0,unit,mode);
  return val/cc;
}

/*
<?php
include("ingen/units.php");
//Parses POST and prepares handy variables to be used later
$gle_kind=$_POST['kind'];     if($gle_kind=="") $gle_kind="optimal";
$out_mode=$_POST['outmode'];  if($out_mode=="") $out_mode="raw";
$multi_mat=0;

$OUTS=$DESC="";
$DESC="# Shall generate a description of the parameters, sometime\n";
?>

<div class="columns">
<div class="lcolumn">
<b>GLE type: &nbsp;&nbsp; </b>
<select name="kind" onchange='this.form.submit()'> 
<?php
option($gle_kind,"optimal",  "Optimal sampling"  );
option($gle_kind,"quantum",  "Quantum thermostat"  );
option($gle_kind,"piglet",   "PIGLET"  );
option($gle_kind,"pimd",     "PI+GLE"  );
option($gle_kind,"respa",    "GLE-RESPA"  );
option($gle_kind,"delta",    "Delta thermostat"  );
option($gle_kind,"exp2",     "Exponential 2x2"        );
option($gle_kind,"kdelta",   "Delta-like K(&omega;)"  );
option($gle_kind,"smart",    "Smart sampling", "color:red" );
?>
</select> <br/>
<b>GLE options:</b><br/>
<?php
switch($gle_kind) {
case "exp2":
  include("ingen/gle_exp2.php");
  break;
case "kdelta":
  include("ingen/gle_kdelta.php");
  break;
case "respa":
  include("ingen/gle_respa.php");
  break;
case "optimal":
  include("ingen/gle_optimal.php");
  break;
case "quantum":
  include("ingen/gle_quantum.php");
  break;
case "delta":
  include("ingen/gle_delta.php");
  break;
case "pimd":
  include("ingen/gle_pimd.php");
  break;
case "piglet":
  include("ingen/gle_piglet.php");
  break;
case "smart":
  include("ingen/gle_smart.php");
  break;  
}
?>
</div>
<div class="rcolumn"> <p>
<?php
switch($gle_kind) {
case "exp2":
echo <<<EOT
     Exponentially-damped noise, corresponding to a memory kernel &gamma; exp(-t/&tau;)/&tau;.
     &gamma; sets the strength of the friction, and &tau; a cutoff on the power spectrum of the
     noise, which will reduced the contribution from frequencies larger than 1/&tau;.
     Enforces fluctuation-dissipation theorem, so canonical sampling will result.
EOT;
break; 
case "kdelta":
echo <<<EOT
     &delta;-like memory kernel: the Fourier transform of the memory kernel is peaked at 
     frequency &omega;0, has width &Delta;&omega; and strength &gamma;. Modes with 
     frequency close to &omega;0 will be most affected by the GLE.
     Enforces fluctuation-dissipation theorem, so canonical sampling will result: not to be
     confused with the &delta; thermostat, which selectively excites a few normal modes.
EOT;
break; 
case "respa":
echo <<<EOT
     Memory kernel suitable to stabilize RESPA multiple time step dynamics
     with a reduced effect on diffusion and low-frequency dynamics.  
     Frequencies above &omega;F will be damped with an effective friction 
     &gamma;&infin;, while the effective friction becomes negligible at
     much lower frequencies. The default values have been tested for biological systems
     with flexible water and an outer time step of 12fs. 
     Enforces fluctuation-dissipation theorem, so canonical sampling will result.
EOT;
break; 
case "optimal":
echo <<<EOT
     Optimal-sampling GLE: tries to minimize the correlation time for sampling
     potential or total energy for normal modes over a broad frequency range, 
     geometrically centered at &omega;0. 
     Enforces fluctuation-dissipation theorem, so canonical sampling will result.
EOT;
break; 
case "smart":
echo <<<EOT
     Smart-sampling GLE: similar to optimal sampling, but tries to improve the sampling
     efficiency for slower collective motion -- with the idea that fast modes will be 
     sampled many times anyway. One needs to specify a rough estimate of the maximum 
     frequency present in the system being studied, and estimate for how long the 
     simulation will be run. The chosen parameters will yield a good compromise 
     between sampling as efficiently as possible the slowest motion within the reach of the
     simulation, and sampling efficiently the faster modes. 
     Enforces fluctuation-dissipation theorem, so canonical sampling will result.
EOT;
break; 
case "quantum":
echo <<<EOT
     "Quantum thermostat": models nuclear quantum effects by pinning the 
     stationary distribution for a harmonic oscillator to the quantum expectation
     value. This is a non-equilibrium GLE, breaking fluctuation-dissipation theorem.
     Quantum fluctuations are fitted up to a given cutoff in frequency
     (which depends on temperature, due to scaling laws). 
     Both configurational properties and the momentum distribution make physical 
     sense, but bear in mind that this is an approximate model, and the accuracy
     may vary depending on anharmonicity. Nevertheless, it was proven to give 
     a very significant improvement over a purely classical description in many 
     systems, and at no cost.
     Compare with path integral results on a smaller or simpler model if 
     quantitative accuracy is very important! Note in particular that weakly-coupled
     parameters suffer a lot from anharmonicity, and shouldn't in general be used.
EOT;
break; 
case "pimd":
echo <<<EOT
     "PI+GLE": to be used in conjunction with imaginary time path integral MD 
     to accelerate convergence to the quantum mechanical results. 
     It is an extension of the quantum thermostat to allow systematic improvement
     of the accuracy of the modelling of nuclear quantum effects.
     This is a non-equilibrium GLE, breaking fluctuation-dissipation theorem.
     Quantum fluctuations are fitted up to a given cutoff in frequency
     (which depends on temperature, due to scaling laws). 
     Only configurational properties make sense. Also sampling efficiency has been
     optimized, and should be comparable to the outcome of an "optimal sampling" 
     classical simulation. 
EOT;
break; 
case "piglet":
echo <<<EOT
     "PIGLET": to be used in conjunction with imaginary time path integral MD 
     to accelerate convergence to the quantum mechanical results. Also converges
     efficiently the quantum kinetic energy - so it is in general to be preferred to
     PI+GLE.  <br/>
EOT;
break; 
case "delta":
echo <<<EOT
     &delta; thermostat: only normal modes within a narrow frequency window are set to 
     a large temperature, while all the others are kept frozen at a much lower T. 
     This can be used to selectively excite normal modes. Note that an unconventional
     ensemble will be sampled, and care must be taken to extract quantitative information.
EOT;
break; 
}
?>
</p>
</div>
</div>
<div class="columns"> </div>
<div class="columns">
<div class="lcolumn">
<b>Output format: &nbsp;&nbsp; </b>
<select name="outmode" onchange='this.form.submit()'> 
<?php
option($out_mode,"raw",      "Raw matrices"        );
option($out_mode,"cp2k",     "CP2K input section"  );
option($out_mode,"cpmd",     "CPMD input files"  );
option($out_mode,"aims",     "FHI-AIMS input section"  );
option($out_mode,"ipi",      "i-PI input section"  );
?>
</select><br/>
<?php
switch($out_mode) {
case "raw":
  include("ingen/out_raw.php");
  break;
case "cp2k":
  include("ingen/out_cp2k.php");
  break;
case "cpmd":
  include("ingen/out_cpmd.php");
  break;
case "aims":
  include("ingen/out_aims.php");
  break;  
case "ipi":
  include("ingen/out_ipi.php");
  break;
}
?>
</div>
<div class="rcolumn">
<b>Generated parameters:</b><br/>
<pre style="overflow:auto; width:100%; height:300px; border:1px solid black;">
<?php
echo $OUTS;
?>
</pre>
</div>
</div>
<div class="columns"> </div>
*/
