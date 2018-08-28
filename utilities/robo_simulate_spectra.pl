#! /usr/bin/env perl

use warnings;
use strict;

my $pi = 4 * atan2(1,1);

use Getopt::Long;


my @control = ();
my %opt = ( #output
	    'wave_min' => 0,
	    'wave_max' => 10000,
	    'wave_delta' => 0.25,
	    'sigma'    => 0.01,
	    'control'  => \@control,
    );

GetOptions(\%opt,
	   "wave_min|w=s",
	   "wave_max|W=s",
	   "wave_delta|d=s",
	   "sigma|s=s",
	   "control|C=s@",
	   "help|h"
    );

if ((exists($opt{help}))||($#ARGV == -1)) {
    print "robo_simulate_spectra.pl [options] <input> <output_base>\n";
    print "     -w --wave_min          Minimum Wavelength    ($opt{wave_min})\n";
    print "     -W --wave_max          Maximum Wavelength    ($opt{wave_max})\n";
    print "     -d --wave_delta        Delta Wavelength      ($opt{wave_delta})\n";
    print "     -s --sigma             Gaussian Noise sigma  ($opt{sigma})\n";
    exit();
}
my $input = shift(@ARGV);
open(I,"$input") || die "Cannot open input file: $input: $!\n";
my $output = shift(@ARGV);
if (defined($output)) {
    open(OS,">${output}.sim.spect") || die "Cannot open output file: $output: $!\n";
    open(OC,">${output}.sim.cont")  || die "Cannot open output file: $output: $!\n";
    open(OT,">${output}.sim.slines")|| die "Cannot open output file: $output: $!\n";
    open(OL,">${output}.sim.lines") || die "Cannot open output file: $output: $!\n";
}
else {
    *OS = *STDOUT;
    *OL = *STDOUT;
}

my @M;
my @S;
my @F;
my @O;
my @L;
while (<I>) {
    chomp;
    unless ($_ =~ /^#/) {
	$_ =~ s/^\s+//;
	my ($mean,$sigma,$EQW,$leftovers) = split /\s+/;
	# unless (defined($offset)) {
	#     $offset = 0.0;
	# }
	push @M,$mean;
	push @S,$sigma;
	push @F,($EQW / (-1000.0));
	push @L,$leftovers;
    }
}

srand();
## Create continuum control
push @control, 1.0;
unshift @control, 1.0;
# 1.0 control points 1.0
my @control_x = ();
my @control_m = ();
for (my $i = 0; $i <= $#control; $i++) {
    if ($i == 0) { 
	$control_x[$i] = $opt{wave_min}; 
    }
    elsif ($i == $#control) { $control_x[$i] = $opt{wave_max}; }
    else {
	$control_x[$i] = $i * ($opt{wave_max} - $opt{wave_min}) / $#control + $opt{wave_min};
    }
#    print STDERR "$i $control_x[$i] $control[$i]\n";
}
for (my $i = 0; $i <= $#control; $i++) {
    if ($i == 0) { 
	$control_m[$i] = ($control[$i+1] - $control[$i]) / ( $control_x[$i+1] - $control_x[$i]);
    }
    if ($i == $#control) {
	$control_m[$i] = ($control[$i] - $control[$i-1]) / ( $control_x[$i] - $control_x[$i-1]);
    }
    else {
	$control_m[$i] = (0.5 * ($control[$i] - $control[$i-1]) / ( $control_x[$i] - $control_x[$i-1]) +
			  0.5 * ($control[$i+1] - $control[$i]) / ( $control_x[$i+1] - $control_x[$i]));
    }
}

my $t_i = 0;
for (my $x = $opt{wave_min}; $x <= $opt{wave_max}; $x += $opt{wave_delta}) {
    my $y = 1.0 + rand_gaussian() * ($opt{sigma});

    if ($x > $control_x[$t_i + 1]) { $t_i ++; }

    my $t = ($x - $control_x[$t_i]) / ($control_x[$t_i + 1] - $control_x[$t_i]);
    
    my $O = ((2 * $t**3 - 3 * $t**2 + 1) * $control[$t_i] +
	     ($t**3 - 2 * $t**2 + $t) * $control_m[$t_i] +
	     (3 * $t**2 - 2 * $t**3) * $control[$t_i+1] +
	     ($t**3 - $t**2) * $control_m[$t_i+1]);

    my $L = 0.0;
    for (my $i = 0; $i <= $#M; $i++) {
	if (abs(($x - $M[$i])/$S[$i]) < 8.0) {
	    my $l = gaussian($x,$M[$i],$S[$i],$F[$i]);
	    if ($l / $F[$i] >= 0) {
		$L += $l;
	    }
	}
    }
    $y += $L;
    $y *= $O;
    print OS "$x $y\n";
    print OT "$x $L\n";
    print OC "$x $O\n";
}
for (my $i = 0; $i <= $#M; $i++) {
    print OL "$M[$i] $S[$i] $F[$i] $L[$i]\n";
}

sub rand_gaussian {
    my ($v1,$v2,$s);
    do {
	$v1 = 2.0 * (rand()) - 1.0;
	$v2 = 2.0 * (rand()) - 1.0;
	$s = $v1 * $v1 + $v2 * $v2;
    } while ($s >= 1.0);
    if ($s == 0.0) {
	return(0.0);
    }
    else {
	return($v1 * sqrt(-2.0 * log($s) / $s));
    }
}
    

sub gaussian {
    my $x = shift;
    my $m = shift;
    my $s = shift;
    my $F = shift;

    return($F / ($s * sqrt(2 * $pi)) * exp(-0.5 * (($x - $m)/$s)**2.0));
}
