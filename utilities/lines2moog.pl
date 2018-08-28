#! /usr/bin/perl -w

use Getopt::Long;


## These are the set of possible command line options, and they're set
## (=>) to default values.
%opt = ( #name
	 'temperature' => 4660,
	 'logg' => 0.75,
	 'feh' => -2.4,
	 'vt'  => 2.05,
	 'model' => 'KUR',
	 'autosynthesis' => 'cno',
	 'other_parameter' => 0,
	 'fwhm' => 10,
	 'break_low' => 10,
	 'break_high' => 9000,
	 'parm1' => 0.02,
	 'parm2' => 1.05,
	 'parm3' => 1.05,
	 'parm4' => 0.10,
	 'iterations' => 15,
	 'flag' => 'l',
	 'flux_calibrated' => 0,
	 'pos_eqw_only' => 0

    );
## This call reads the command line, and stores any of the options it
## finds into the correct place.  The | says "this is a short alias
## for the option" and the =i says "integer!" =s says "string!".  If
## there's no =, then it sets it as a boolean

GetOptions(\%opt,
	   'name|N=s',
	   'temperature|T=i',
	   'logg|G=i',
	   'feh|F=i',
	   'vt|v=i',
	   'model|M=s',
	   'autosynthesis|A=s',
	   'other_parameter=i',
	   'fwhm=i',
	   'break_low=i',
	   'break_high=i',
	   'parm1|x=i',
	   'parm2|y=i',
	   'parm3|z=i',
	   'parm4|q=i',
	   'iterations|I=i',
	   'flag=s',
	   'flux_calibrated',
	   'pos_eqw_only',
	   'help|h'
    );


## We now read the filename from the command line, and set it to the
## "name" option.
my $input = shift(@ARGV);
unless(defined($opt{name})) {
    unless (defined($input)) {
	$input = '';
    }
    $opt{name} = $input;
## Then we strip off the .lines extension so we just have teh base name.
    $opt{name} =~ s%^.*/%%;
    $opt{name} =~ s%.lines$%%;
}

## List a help entry if the help option was found above.
if (defined($opt{help})) {
    print "Usage: lines2moog.pl [options] <robospect.lines>\n";
    print "  -h, --help                This help            (found)\n";
    print "               STAR PARAMETERS\n";
    print "  -N, --name <star name>                         ($opt{name})\n";
    print "  -T, --temperature <T>                          ($opt{temperature})\n";
    print "  -G, --logg <logg>                              ($opt{logg})\n";
    print "  -F, --feh  <Fe/H>                              ($opt{feh})\n";
    print "  -v, --vt   <v_t>                               ($opt{vt})\n";
    print "               MOOG FITTING PARAMETERS\n";
    print "  -M, --model <model>                            ($opt{model})\n";
    print "  -A, --autosynthesis <this thing>               ($opt{autosynthesis})\n";
    print "      --other_parameter <no clue>                ($opt{other_parameter})\n";
    print "      --fwhm <maybe?>                            ($opt{fwhm})\n";
    print "      --break_low <x>                            ($opt{break_low})\n";
    print "      --break_high <x>                           ($opt{break_high})\n";
    print "  -x, --parm1 <monkey>                           ($opt{parm1})\n";
    print "  -y, --parm2 <chicken>                          ($opt{parm2})\n";
    print "  -z, --parm3 <veal>                             ($opt{parm3})\n";
    print "  -q, --parm4 <cheese>                           ($opt{parm4})\n";
    print "  -I, --iterations <N_iter>                      ($opt{iterations})\n";
    print "      --flag <flag for lines>                    ($opt{flag})\n";
    print "      --pos_eqw_only                             ($opt{pos_eqw_only})\n";
    print "      --flux_calibrated                          ($opt{flux_calibrated}\n";
    exit(10);  ## And then exit without doing anything.
}


## Construct the output file by converting .lines to .moog.input
my $output = $input;
$output =~ s/.robolines$/.moog.input/;
## Open the input file for reading or stop and complain.
open(L,$input) || die "Could not open input lines file: $input\n";  
## Open the output file for writing (>) or stop and complain.
open(O,">$output") || die "Could not open output MOOG input file: $output\n";

## Print header information based on what MOOG expects the first line to be.
print O "$opt{name} $opt{temperature} $opt{logg} $opt{feh} $opt{vt} ";
print O "$opt{model} $opt{autosynthesis} $opt{other_parameter} $opt{fwhm} ";
print O "$opt{break_low} $opt{break_high} ";
print O "$opt{parm1} $opt{parm2} $opt{parm3} $opt{parm4} $opt{iterations}\n";

## Read through the input file
while(<L>) {
    chomp; ## Cut off any newlines
    unless ($_ =~ /^#/) {  ## If this line is NOT a comment, then:
	## Name all of our variables we'll use in this unless-block
	my ($x0,$mean,$sigma,$flux,$eta,$dmean,$dsigma,$dflux,$deta,
	    $fwhm_mean,$fwhm_sigma,$fwhm_flux,
	    $eqw,$deqw,$chi,$flags,$group,undef,
	    $comment_x0,$comment_species,$comment_xp,$comment_loggf,
	    @comment_remainder);
	## Then use the split command to separate the line we've read
	## from the input list into all of the variables.
	($x0,$mean,$sigma,$flux,$eta,$dmean,$dsigma,$dflux,$deta,
	 $fwhm_mean,$fwhm_sigma,$fwhm_flux,
	 $eqw,$deqw,$chi,$flags,$group,undef,
	 $comment_x0,$comment_species,$comment_xp,$comment_loggf,
	 @comment_remainder) = split /\s+/;
	## Convert the text flags into a numeric value for perl
	$flags = hex($flags);
	unless (($comment_x0 eq 'Auto-found')||  ## Unless this line was auto-found,
		($flags & 0x00000f)) {           ## or it has a bad flag value:
	    my $this_line;  ## Define the output line variable
	    ## Format a line of MOOG output(input?) based on the
	    ## values we found above in the ROBOSPECT output lines.
	    if (($opt{pos_eqw_only})&&($eqw < 0)) {
		next;
	    }
	    if ($opt{flux_calibrated}) {
		$this_line = sprintf("% 10.2f     % 3.1g     % 5.2g     % 5.2g     % 25.1g  %s\n",
				     $comment_x0,$comment_species,$comment_xp,$comment_loggf,
				     $eqw,$opt{flag});
	    }
	    else {
		$this_line = sprintf("% 10.2f     % 3.1f     % 5.2f     % 5.2f     % 25.1f  %s\n",
				     $comment_x0,$comment_species,$comment_xp,$comment_loggf,
				     $eqw,$opt{flag});
	    }
	    ## Add this line to an array of hashes of hashes.  Yikes.
	    ## Think of it like this: We know that we have a bunch of
	    ## species we want to sort by, and each species can have a
	    ## line at a different place.  So we declare that we want
	    ## to store something by species, and we call it
	    ## $output_data{THIS_SPECIES}.  Next, we say for each
	    ## species we have a wavelength, which is just:
	    ## $output_data{THIS_SPECIES}{THIS_WAVELENGTH} (in perl
	    ## shorthand).  Now, Maybe we have lots of lines at this
	    ## combination, so we add $this_line (that we defined
	    ## above) to an array that is now named: @{
	    ## $output_data{THIS_SPECIES}{THIS_WAVELENGTH} }.  This is
	    ## the perl long-hand way of doing it, but it's clearer
	    ## for me to read.  I'm pushing something onto an array,
	    ## but the name of that array is held by that
	    ## hash-of-hashes.  
	    push @{ $output_data{$comment_species}{$comment_x0} }, $this_line;
	} ## End unless-this-line-was-good
    } ## End unless-this-line-was-comment
} ## End read this line. Loop!
close(L); ## Close input ROBOSPECT .lines file

foreach my $s (sort {$a <=> $b} (keys %output_data)) {  ## foreach numerically sorted key of %output_data (these keys are the THIS_SPECIES above)
    foreach my $x (sort {$a <=> $b} (keys %{ $output_data{$s} })) { ## foreach numerically sorted key of %{ $output_data{THIS_SPECIES} } (these keys are the THIS_WAVELENGTH above)
	foreach my $l (@{ $output_data{$s}{$x} }) { ## foreach element of the array we constructed at each THIS_SPECIES/THIS_WAVELENGTH combo
	    print O "$l"; ## print out the line we constructed above
	}
    }
}
close(O);  ## Close output MOOG .input file


    
